//! Main HAWK-512 verification loop.
//!
//! Port of `Zh(verify_inner)` (hawk_vrfy.c:1585-2154).

use crate::error::HawkError;
use crate::keygen::mp31::{mp_add, mp_montymul};
use crate::verify::consts::{
    BITS_LIM01, BITS_LIMS0, BITS_LIMS1, HAWK_512_MAX_TNORM, P1, P1_0I, P1_M16, P1_R3, P2, P2_0I,
    P2_M16, P2_R3,
};
use crate::verify::fx32::{fx32_fft, fx32_ifft, fx32_of, fx32_rint};
use crate::verify::helpers::make_t1;
use crate::verify::mp::{mp_div, mp_poly_to_ntt, mp_poly_to_ntt_autoadj, GM_P1, GM_P2};
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;

/// tbmask: returns 0xFFFF_FFFF if x has its high bit set (i.e., x as i32 < 0),
/// else 0. Port of `tbmask` (ng_inner.h).
#[inline]
fn tbmask(x: u32) -> u32 {
    (x as i32 >> 31) as u32
}

/// Inner verification function.
///
/// Inputs:
///   - `msg`: the message being verified
///   - `pub_bytes`: encoded public key (1024 bytes for HAWK-512)
///   - `q00`: decoded public-key polynomial (first n/2 coefficients are the
///     auto-adjoint half; q00[n/2..n] should be zero from `decode_public`)
///   - `q01`: decoded public-key polynomial (n coefficients)
///   - `salt`: 24-byte signature salt
///   - `s1`: n-element i16 signature polynomial
///
/// Returns `Ok(())` if the signature is cryptographically valid, or
/// `Err(HawkError::InvalidSignature)` if verification fails at any check.
///
/// Port of `Zh(verify_inner)` (hawk_vrfy.c:1585-2154).
pub fn verify_inner(
    msg: &[u8],
    pub_bytes: &[u8],
    q00: &[i16],
    q01: &[i16],
    salt: &[u8],
    s1: &[i16],
) -> Result<(), HawkError> {
    let logn: u32 = 9;
    let n = 512usize;
    let hn = 256usize;
    let max_tnorm = HAWK_512_MAX_TNORM;

    // -------------------------------------------------------------------------
    // Step 1: Compute hashes.
    //
    //   hpub = SHAKE256(pub_bytes)[..32]   (hpub_len = 2^(logn-4) = 32)
    //   hm   = SHAKE256(msg || hpub)[..64]
    //   h    = SHAKE256(hm || salt)[..n/4] = h0 || h1
    //
    // Port of hawk_vrfy.c:1665-1707.
    // -------------------------------------------------------------------------
    let hpub_len: usize = 1 << (logn - 4); // 32
    let mut hpub = [0u8; 32];
    {
        let mut sh = Shake256::default();
        sh.update(pub_bytes);
        sh.finalize_xof().read(&mut hpub[..hpub_len]);
    }

    let mut hm = [0u8; 64];
    {
        let mut sh = Shake256::default();
        sh.update(msg);
        sh.update(&hpub[..hpub_len]);
        sh.finalize_xof().read(&mut hm);
    }

    // n/4 = 128 bytes; split into h0 (first n/8=64 bytes) and h1 (next 64).
    let mut h_bytes = [0u8; 128];
    {
        let mut sh = Shake256::default();
        sh.update(&hm);
        sh.update(salt);
        sh.finalize_xof().read(&mut h_bytes[..n / 4]);
    }
    let h0 = &h_bytes[..n / 8]; // 64 bytes
    let h1 = &h_bytes[n / 8..n / 4]; // 64 bytes

    // -------------------------------------------------------------------------
    // Step 2: Compute t1 = h1 - 2*s1 in fx32 + FFT.
    //
    // Also check the sym-break condition on t1.
    // Port of hawk_vrfy.c:1724-1746.
    // -------------------------------------------------------------------------
    let sh_t1: i32 = 29 - (1 + BITS_LIMS1[logn as usize] as i32); // 29 - 10 = 19
    let sh_q00: i32 = 29 - 9; // bits_lim00[9] = 9 → sh_q00 = 20
    let sh_q01: i32 = 29 - BITS_LIM01[logn as usize] as i32; // 29 - 12 = 17

    let mut ft1 = vec![0u32; n];
    let mut csb: u32 = 0xFFFF_FFFF;
    for u in (0..n).step_by(8) {
        let mut hb = h1[u >> 3] as u32;
        for v in 0..8 {
            let w_s1 = s1[u + v] as i32 as u32;
            let w = (hb & 1).wrapping_sub(w_s1 << 1);
            ft1[u + v] = fx32_of(w as i32, sh_t1 as u32);
            // Sym-break: if first nonzero value has high bit set, reject.
            if (csb & w) >> 31 != 0 {
                return Err(HawkError::InvalidSignature);
            }
            csb &= !tbmask(w.wrapping_neg());
            hb >>= 1;
        }
    }
    if csb != 0 {
        // t1 is entirely zero — not valid.
        return Err(HawkError::InvalidSignature);
    }
    fx32_fft(logn, &mut ft1);

    // -------------------------------------------------------------------------
    // Step 3: Convert q00 (auto-adjoint) to fx32 + FFT.
    //
    // Port of hawk_vrfy.c:1785-1806.
    // -------------------------------------------------------------------------
    if q00[0] < 0 {
        return Err(HawkError::InvalidSignature);
    }
    let cst_q00 = q00[0] as i32;

    // fq00 has only n/2 real slots: fq00[0..hn] = real parts.
    // fq00[hn..n] would be imaginary parts (zero for auto-adjoint), but we
    // allocate n words so the FFT can work in-place at full size.
    let mut fq00 = vec![0u32; n];
    fq00[0] = 0; // constant term cleared (injected back later as cstup)
    fq00[hn] = 0;
    for u in 1..hn {
        let z = fx32_of(q00[u] as i32, sh_q00 as u32);
        fq00[u] = z;
        fq00[n - u] = z.wrapping_neg();
    }
    fx32_fft(logn, &mut fq00);

    // -------------------------------------------------------------------------
    // Step 4: Convert q01 to fx32 + FFT.
    //
    // Port of hawk_vrfy.c:1841-1844.
    // -------------------------------------------------------------------------
    let mut fq01 = vec![0u32; n];
    for u in 0..n {
        fq01[u] = fx32_of(q01[u] as i32, sh_q01 as u32);
    }
    fx32_fft(logn, &mut fq01);

    // -------------------------------------------------------------------------
    // Step 5: fq01 <- (q01 * t1) / q00 in FFT representation.
    //
    // This is the most delicate part: fixed-point complex division.
    // Port of hawk_vrfy.c:1896-1932.
    //
    // The constant term of q00 was cleared before FFT. We now inject it back
    // as an additive correction:
    //   cstup = cst_q00 << (sh_q00 - (logn - 1))
    //         = cst_q00 << (20 - 8) = cst_q00 << 12
    // (This is because the FFT-domain constant corresponds to the sum of all
    //  coefficients, and the scaling factor cancels the downscaling of FFT.)
    // -------------------------------------------------------------------------
    // cstup = (uint32_t)cstq00 << (sh_q00 - (logn - 1))
    let cstup: u32 = (cst_q00 as u32) << (sh_q00 as u32 - (logn - 1));

    for u in 0..hn {
        // M(c, d) = (c as i32 as i64 * d as i32 as i64) as u64
        let q01_re = fq01[u];
        let q01_im = fq01[u + hn];
        let t1_re = ft1[u];
        let t1_im = ft1[u + hn];

        // x = q01 * t1 (complex product, i64 arithmetic)
        let x_re: u64 = ((q01_re as i32 as i64) * (t1_re as i32 as i64)
            - (q01_im as i32 as i64) * (t1_im as i32 as i64)) as u64;
        let x_im: u64 = ((q01_re as i32 as i64) * (t1_im as i32 as i64)
            + (q01_im as i32 as i64) * (t1_re as i32 as i64)) as u64;

        // Extract signs, then work with absolute values.
        let sx_re: u64 = (x_re as i64 >> 63) as u64;
        let sx_im: u64 = (x_im as i64 >> 63) as u64;
        let x_re_abs: u64 = (x_re ^ sx_re).wrapping_sub(sx_re);
        let x_im_abs: u64 = (x_im ^ sx_im).wrapping_sub(sx_im);

        let x_re_hi = (x_re_abs >> 32) as u32;
        let x_re_lo = x_re_abs as u32;
        let x_im_hi = (x_im_abs >> 32) as u32;
        let x_im_lo = x_im_abs as u32;

        // Divisor = fq00[u] + cstup (must be in (0, 2^30)).
        let w00: u32 = cstup.wrapping_add(fq00[u]);
        // Reject on three conditions:
        //   1. `w00 == 0` (would divide by zero) or `w00 >= 2^30`
        //      (out of the valid range for a valid pub key); the unsigned
        //      comparison `(w00 - 1) >= 0x3FFF_FFFF` catches both since
        //      `w00 = 0` wraps to `0xFFFFFFFF >= 0x3FFFFFFF`.
        //   2. `x_re_hi >= w00`: the 64/32 quotient would overflow u32.
        //      This matters for API correctness (u64/u32 division doesn't
        //      panic in Rust — the quotient simply exceeds u32 and gets
        //      truncated later, producing garbage).
        //   3. Same for `x_im_hi`.
        //
        // Port of hawk_vrfy.c:1893-1897. All three conditions indicate
        // malformed / malicious public key or signature and must trigger
        // rejection.
        if (w00.wrapping_sub(1)) >= 0x3FFF_FFFF || x_re_hi >= w00 || x_im_hi >= w00 {
            return Err(HawkError::InvalidSignature);
        }

        // Unsigned 64-by-32 division: (hi:lo) / w00.
        let y_re = (((x_re_hi as u64) << 32) | x_re_lo as u64) / w00 as u64;
        let y_im = (((x_im_hi as u64) << 32) | x_im_lo as u64) / w00 as u64;

        // Re-apply signs.
        let y_re = y_re as u32;
        let y_im = y_im as u32;
        let sx_re32 = sx_re as u32;
        let sx_im32 = sx_im as u32;
        fq01[u] = (y_re ^ sx_re32).wrapping_sub(sx_re32);
        fq01[u + hn] = (y_im ^ sx_im32).wrapping_sub(sx_im32);
    }
    fx32_ifft(logn, &mut fq01);

    // After iFFT, fq01 holds (q01*t1)/q00 in coefficient form.

    // -------------------------------------------------------------------------
    // Step 6: s0 = round(h0/2 + (q01*t1)/q00) and t0 = h0 - 2*s0.
    //
    // Port of hawk_vrfy.c:1948-1965.
    // -------------------------------------------------------------------------
    // sh_s0 = sh_t1 + sh_q01 - sh_q00 - (logn - 1)
    let sh_s0: i32 = sh_t1 + sh_q01 - sh_q00 - (logn as i32 - 1);
    let lim_s0 = 1i32 << BITS_LIMS0[logn as usize];
    let mut t0 = vec![0i16; n];
    for u in (0..n).step_by(8) {
        let mut h0b = h0[u >> 3] as u32;
        for v in 0..8 {
            let bit = h0b & 1;
            // w = h0_bit/2 (scaled) + fq01[u+v]  (scaled)
            let w = fx32_of(bit as i32, sh_s0 as u32).wrapping_add(fq01[u + v]);
            // Round to nearest integer at precision sh_s0+1.
            let z = fx32_rint(w, (sh_s0 + 1) as u32);
            if z < -lim_s0 || z >= lim_s0 {
                return Err(HawkError::InvalidSignature);
            }
            // t0[u+v] = h0_bit - 2*s0   (both as i16)
            let w2 = bit.wrapping_sub((z as u32) << 1);
            t0[u + v] = w2 as i32 as i16;
            h0b >>= 1;
        }
    }

    // -------------------------------------------------------------------------
    // Step 7: Dual-prime NTT rounds to compute n*sqnorm_Q(t).
    //
    // Uses t0, t1 (recomputed from h1+s1), q00, q01 in NTT form.
    // Accumulates n*sqnorm_Q(t)/2 modulo P1 and P2; they must agree.
    // Port of hawk_vrfy.c:2008-2134.
    // -------------------------------------------------------------------------
    let mut c1 = vec![0u32; n];
    let mut c2 = vec![0u32; n];
    let mut tnorm: u32 = 0;

    for i in 0..2usize {
        let (p, p0i, r3, m16, gm): (u32, u32, u32, u32, &[u32]) = if i == 0 {
            (P1, P1_0I, P1_R3, P1_M16, GM_P1.as_slice())
        } else {
            (P2, P2_0I, P2_R3, P2_M16, GM_P2.as_slice())
        };

        // c2 <- t1 = h1 - 2*s1 (recomputed from s1 + h1).
        let mut c2hi = s1.to_vec(); // starts as s1
        make_t1(logn, &mut c2hi, h1);
        mp_poly_to_ntt(logn, &mut c2, &c2hi, p, p0i, gm);

        // c1 <- q00 (auto-adjoint half — only first hn coefficients used).
        mp_poly_to_ntt_autoadj(logn, &mut c1, q00, p, p0i, gm);

        // c1 <- 1/c1   (Montgomery's trick to invert hn values with 1 true inversion).
        // After this, c1[0..hn] = element-wise inverse of q00 in NTT form.
        //
        // Algorithm (hawk_vrfy.c:2045-2057):
        //   bx = c1[0]
        //   c1[hn] = bx
        //   for u in 1..hn: bx = montymul(bx, c1[u]); c1[u+hn] = bx
        //   bx = mp_div(1, bx, ...)
        //   for u in (1..hn).rev():
        //     ix = montymul(bx, c1[u+hn-1], ...); bx = montymul(bx, c1[u], ...); c1[u] = ix
        //   c1[0] = bx
        let mut bx = c1[0];
        c1[hn] = bx;
        for u in 1..hn {
            bx = mp_montymul(bx, c1[u], p, p0i);
            c1[u + hn] = bx;
        }
        bx = mp_div(1, bx, p, p0i, m16);
        for u in (1..hn).rev() {
            let ix = mp_montymul(bx, c1[u + hn - 1], p, p0i);
            bx = mp_montymul(bx, c1[u], p, p0i);
            c1[u] = ix;
        }
        c1[0] = bx;

        // c2 <- c2 * c1 = t1/q00 (element-wise)
        // nnacc = Tr_p(t1 * adj(t1) / q00)  (in double-anti-Montgomery repr)
        //
        // Port of hawk_vrfy.c:2062-2071.
        let mut nnacc: u32 = 0;
        for u in 0..hn {
            let x1 = c2[u];
            let x2 = c2[(n - 1) - u];
            let qx = c1[u];
            let x1 = mp_montymul(x1, qx, p, p0i);
            nnacc = mp_add(nnacc, mp_montymul(x1, x2, p, p0i), p);
            let x2 = mp_montymul(x2, qx, p, p0i);
            c2[u] = x1;
            c2[(n - 1) - u] = x2;
        }
        // nnacc is now in double-anti-Montgomery representation (divided by R^2).

        // c1 <- q01 in NTT form.
        // Port of hawk_vrfy.c:2080.
        mp_poly_to_ntt(logn, &mut c1, q01, p, p0i, gm);

        // c1 <- c1 * c2 * R3 = q01 * t1/q00, adjusted to normal Montgomery repr.
        // Port of hawk_vrfy.c:2083-2086.
        for u in 0..n {
            c1[u] = mp_montymul(mp_montymul(c1[u], c2[u], p, p0i), r3, p, p0i);
        }

        // c2 <- t0 in NTT form.
        // Port of hawk_vrfy.c:2089.
        mp_poly_to_ntt(logn, &mut c2, &t0, p, p0i, gm);

        // c2 <- c2 + c1 = t0 + q01*t1/q00 = e.
        // Port of hawk_vrfy.c:2092-2094.
        for u in 0..n {
            c2[u] = mp_add(c2[u], c1[u], p);
        }

        // c1 <- q00 (auto-adjoint) again.
        // Port of hawk_vrfy.c:2101.
        mp_poly_to_ntt_autoadj(logn, &mut c1, q00, p, p0i, gm);

        // nnacc += Tr_p(q00 * e * adj(e))
        // Port of hawk_vrfy.c:2104-2110.
        for u in 0..hn {
            let x1 = c2[u];
            let x2 = c2[(n - 1) - u];
            let qx = c1[u];
            nnacc = mp_add(
                nnacc,
                mp_montymul(qx, mp_montymul(x1, x2, p, p0i), p, p0i),
                p,
            );
        }

        // Convert from double-anti-Montgomery to normal by multiplying with R3.
        // Port of hawk_vrfy.c:2114.
        nnacc = mp_montymul(nnacc, r3, p, p0i);

        if i == 0 {
            tnorm = nnacc;
        } else if tnorm != nnacc {
            // Values disagree mod P1 and P2 → tnorm > P2 → too large.
            return Err(HawkError::InvalidSignature);
        }
    }

    // -------------------------------------------------------------------------
    // Step 8: Final bound check.
    //
    // tnorm must be divisible by 2^(logn-1) = hn, and
    // (tnorm >> (logn-1)) <= max_tnorm.
    //
    // Port of hawk_vrfy.c:2153.
    // -------------------------------------------------------------------------
    if (tnorm & ((hn as u32) - 1)) != 0 {
        return Err(HawkError::InvalidSignature);
    }
    if (tnorm >> (logn - 1)) > max_tnorm {
        return Err(HawkError::InvalidSignature);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::keygen::HawkKeypair;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn verify_rejects_mutated_signature() {
        let mut rng = ChaCha20Rng::from_seed([7u8; 32]);
        let kp = HawkKeypair::generate(&mut rng);
        let sig = kp.secret.sign(b"hello", &mut rng).unwrap();
        let mut bytes = sig.to_bytes().unwrap();
        bytes[0] ^= 0xAA;
        let mutated = crate::sign::HawkSignature::from_bytes(&bytes);
        if let Ok(mutated) = mutated {
            let r = kp.public.verify(b"hello", &mutated);
            assert!(r.is_err(), "Should reject mutated sig");
        }
        // If from_bytes fails, decode-level rejection is also acceptable.
    }

    #[test]
    fn verify_rejects_wrong_message() {
        let mut rng = ChaCha20Rng::from_seed([11u8; 32]);
        let kp = HawkKeypair::generate(&mut rng);
        let sig = kp.secret.sign(b"hello", &mut rng).unwrap();
        let r = kp.public.verify(b"goodbye", &sig);
        assert!(r.is_err(), "Should reject wrong message");
    }

    #[test]
    fn verify_accepts_valid_signature() {
        let mut rng = ChaCha20Rng::from_seed([42u8; 32]);
        let kp = HawkKeypair::generate(&mut rng);
        let mut srng = ChaCha20Rng::from_seed([99u8; 32]);
        let sig = kp.secret.sign(b"test message", &mut srng).unwrap();
        assert!(
            kp.public.verify(b"test message", &sig).is_ok(),
            "Should accept valid signature"
        );
    }
}
