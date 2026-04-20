//! HAWK keypair generation top-level orchestrator.
//!
//! Port of c-reference/hawk-512/ng_hawk.c:425-686 (`Hawk_keygen`).

use crate::keygen::fxp::{vect_fft, vect_ifft};
use crate::keygen::fxr::{Fxr, FXR_ZERO};
use crate::keygen::mp31::{
    mp_add, mp_intt, mp_mkgm, mp_mkgmigm, mp_montymul, mp_norm, mp_ntt, mp_set,
};
use crate::keygen::ntru::solve_ntru;
use crate::keygen::ntru_profile::SOLVE_HAWK_512;
use crate::keygen::poly::poly_sqnorm;
use crate::keygen::primes::PRIMES;
use crate::keygen::q_derive::{make_q001, BITS_LIM00, BITS_LIM01, BITS_LIM11};
use crate::keygen::regen_fg::{hawk_regen_fg, HAWK_REGEN_FG_SEED_LEN};

/// Maximum norm floor check for (f,g). Value 2080 is from ng_hawk.c:462
/// for HAWK-512 (logn=9).
const HAWK_512_L2LOW: u32 = 2080;

/// Maximum constant-term threshold for 1/(f*adj(f) + g*adj(g)). Value
/// 4294967 in fxr_of_scaled32 form is ~1/1000. From ng_hawk.c:463.
const HAWK_512_D0HIGH: Fxr = Fxr(4294967);

/// Output of a successful HAWK-512 keygen.
pub struct HawkKeygenOutput {
    /// Secret polynomial f.
    pub f: Vec<i8>,
    /// Secret polynomial g.
    pub g: Vec<i8>,
    /// Secret polynomial F.
    pub f_cap: Vec<i8>,
    /// Secret polynomial G.
    pub g_cap: Vec<i8>,
    /// Public Q-polynomial q00 (auto-adjoint, int16).
    pub q00: Vec<i16>,
    /// Public Q-polynomial q01 (int16).
    pub q01: Vec<i16>,
    /// Public Q-polynomial q11 (auto-adjoint, int32).
    pub q11: Vec<i32>,
    /// Seed (24 bytes) that was used to derive (f, g).
    pub seed: Vec<u8>,
}

/// Generate a HAWK-512 keypair from an RNG that supplies 24-byte seeds.
///
/// `rng_bytes` is a callback or closure that fills a byte buffer with random
/// data (simulating the C's `rng(rng_context, buf, len)`). The function may
/// call it multiple times if the initial seed attempt fails sanity checks.
///
/// Returns the keypair output on success, or Err(()) if tmp is too small
/// or logn is unsupported.
///
/// Port of `Hawk_keygen` (ng_hawk.c:425-686) for HAWK-512 (logn=9) only.
pub fn hawk_keygen_512<F>(mut rng_bytes: F) -> HawkKeygenOutput
where
    F: FnMut(&mut [u8]),
{
    let logn: u32 = 9;
    let n = 1usize << logn;
    let hn = n >> 1;
    let seed_len = HAWK_REGEN_FG_SEED_LEN; // 24 for HAWK-512
    let l2low = HAWK_512_L2LOW;
    let d0high = HAWK_512_D0HIGH;
    let lim00 = 1i32 << BITS_LIM00[logn as usize];
    let lim01 = 1i32 << BITS_LIM01[logn as usize];
    let lim11 = 1i32 << BITS_LIM11[logn as usize];

    let mut seed = vec![0u8; seed_len];
    let mut f = vec![0i8; n];
    let mut g = vec![0i8; n];

    let p0 = PRIMES[0].p;
    let p0i = PRIMES[0].p0i;
    let r2_0 = PRIMES[0].r2;
    let p1 = PRIMES[1].p;
    let p1_0i = PRIMES[1].p0i;

    loop {
        // Step 1: Sample (f, g) from a fresh seed.
        rng_bytes(&mut seed);
        hawk_regen_fg(&mut f, &mut g, &seed);

        // Step 2: Parity check.
        if parity(&f) != 1 || parity(&g) != 1 {
            continue;
        }

        // Step 3: Squared-norm check.
        let norm2 = poly_sqnorm(&f).wrapping_add(poly_sqnorm(&g));
        if norm2 < l2low {
            continue;
        }

        // Step 4: Invertibility mod (X^n+1, PRIMES[0].p) + build q00 in tmp.
        let mut tmp = vec![0u32; 20 * n]; // Generous buffer.
        let mut invertible = true;
        {
            let (t1, rest) = tmp.split_at_mut(n);
            let (t2, rest2) = rest.split_at_mut(n);
            let (t3, rest3) = rest2.split_at_mut(n);
            let (t4, _) = rest3.split_at_mut(n);
            mp_mkgmigm(logn, t1, t2, PRIMES[0].g, PRIMES[0].ig, p0, p0i);
            for u in 0..n {
                t3[u] = mp_set(f[u] as i32, p0);
                t4[u] = mp_set(g[u] as i32, p0);
            }
            mp_ntt(logn, t3, t1, p0, p0i);
            mp_ntt(logn, t4, t1, p0, p0i);
            for u in 0..n {
                let prod_f = mp_montymul(t3[u], t3[(n - 1) - u], p0, p0i);
                let prod_g = mp_montymul(t4[u], t4[(n - 1) - u], p0, p0i);
                let x = mp_add(prod_f, prod_g, p0);
                if x == 0 {
                    invertible = false;
                    break;
                }
                t1[u] = mp_montymul(r2_0, x, p0, p0i);
            }
        }
        if !invertible {
            continue;
        }
        // iNTT to get plain t1 = f*adj(f) + g*adj(g).
        {
            let (t1, rest) = tmp.split_at_mut(n);
            let (t2, _) = rest.split_at_mut(n);
            mp_intt(logn, t1, t2, p0, p0i);
            for u in 0..n {
                t1[u] = mp_norm(t1[u], p0) as u32;
            }
        }

        // Step 5: Invertibility mod PRIMES[1].p.
        let mut invertible2 = true;
        {
            let (t1, rest) = tmp.split_at_mut(n);
            let (t2, rest2) = rest.split_at_mut(n);
            let (t3, _) = rest2.split_at_mut(n);
            for u in 0..n {
                t2[u] = mp_set(t1[u] as i32, p1);
            }
            mp_mkgm(logn, t3, PRIMES[1].g, p1, p1_0i);
            mp_ntt(logn, t2, t3, p1, p1_0i);
            for u in 0..n {
                if t2[u] == 0 {
                    invertible2 = false;
                    break;
                }
            }
        }
        if !invertible2 {
            continue;
        }

        // Step 6: Constant-term check on 1/(f*adj(f) + g*adj(g)).
        // Convert t1 (plain signed) to Fxr, FFT, invert first half, zero second,
        // iFFT, check rt1[0] against d0high.
        let rt1_check_passes = {
            let mut rt1: Vec<Fxr> = (0..n).map(|u| Fxr::of(tmp[u] as i32)).collect();
            vect_fft(logn, &mut rt1);
            for u in 0..hn {
                rt1[u] = rt1[u].inv();
            }
            for u in hn..n {
                rt1[u] = FXR_ZERO;
            }
            vect_ifft(logn, &mut rt1);
            // fxr_lt(d0high, rt1[0]) = d0high < rt1[0] → reject.
            !d0high.lt(rt1[0])
        };
        if !rt1_check_passes {
            continue;
        }

        // Step 7: Solve the NTRU equation.
        let err = solve_ntru(&SOLVE_HAWK_512, logn, &f, &g, &mut tmp);
        if err != 0 {
            continue;
        }

        // Step 8: Extract F and G from tmp (packed as 2n int8 at tmp[0..n/2] u32).
        // The C sets `int8_t *tF = (int8_t *)tt32; int8_t *tG = tF + n;`.
        // So tmp[0..n*bytes] is F, then tmp[n*bytes..2n*bytes] is G.
        // n bytes = n/4 u32 words (for HAWK-512, n=512 → 128 u32 words per F/G).
        //
        // Unpack: for each i8, extract from tmp[i/4] via byte-shift.
        let mut f_cap = vec![0i8; n];
        let mut g_cap = vec![0i8; n];
        for u in 0..n {
            let byte_idx = u;
            let word_idx = byte_idx / 4;
            let byte_in_word = byte_idx % 4;
            let b = (tmp[word_idx] >> (byte_in_word * 8)) as u8;
            f_cap[u] = b as i8;
        }
        // G starts at byte offset n.
        for u in 0..n {
            let byte_idx = n + u;
            let word_idx = byte_idx / 4;
            let byte_in_word = byte_idx % 4;
            let b = (tmp[word_idx] >> (byte_in_word * 8)) as u8;
            g_cap[u] = b as i8;
        }

        // Step 9: Derive Q-polynomials via make_q001.
        // The C passes (uint32_t *)(tG + n) as tmp to make_q001, which is
        // byte offset 2*n = n/2 u32 words into tt32.
        // We allocate a fresh buffer here for simplicity.
        let mut q_tmp = vec![0u32; 20 * n];
        match make_q001(
            logn, lim00, lim01, lim11, &f, &g, &f_cap, &g_cap, &mut q_tmp,
        ) {
            Ok(q_out) => {
                // Step 10: Success. Return.
                return HawkKeygenOutput {
                    f: f.clone(),
                    g: g.clone(),
                    f_cap,
                    g_cap,
                    q00: q_out.q00,
                    q01: q_out.q01,
                    q11: q_out.q11,
                    seed: seed.clone(),
                };
            }
            Err(()) => {
                // Retry.
                continue;
            }
        }
    }
}

/// Helper: sum of coefficients mod 2.
fn parity(f: &[i8]) -> u32 {
    let mut pp: u32 = 0;
    for &c in f.iter() {
        pp = pp.wrapping_add(c as u8 as u32);
    }
    pp & 1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hawk_keygen_512_terminates_with_chacha_rng() {
        // Use ChaCha20 for determinism.
        use rand::{RngCore, SeedableRng};
        use rand_chacha::ChaCha20Rng;
        let mut rng = ChaCha20Rng::from_seed([42u8; 32]);
        let kp = hawk_keygen_512(|buf| rng.fill_bytes(buf));
        assert_eq!(kp.f.len(), 512);
        assert_eq!(kp.g.len(), 512);
        assert_eq!(kp.f_cap.len(), 512);
        assert_eq!(kp.g_cap.len(), 512);
        assert_eq!(kp.q00.len(), 512);
        assert_eq!(kp.q01.len(), 512);
        assert_eq!(kp.q11.len(), 512);
        assert_eq!(kp.seed.len(), 24);
    }

    #[test]
    fn parity_basic() {
        // All zeros → parity 0.
        let zeros = vec![0i8; 8];
        assert_eq!(parity(&zeros), 0);
        // Single odd element → parity 1.
        let mut one = vec![0i8; 8];
        one[0] = 1;
        assert_eq!(parity(&one), 1);
        // Two odd elements → parity 0.
        let mut two = vec![0i8; 8];
        two[0] = 1;
        two[1] = 1;
        assert_eq!(parity(&two), 0);
    }

    #[test]
    fn poly_sqnorm_basic() {
        use crate::keygen::poly::poly_sqnorm;
        let f = vec![1i8, -1, 2, 0];
        // 1 + 1 + 4 + 0 = 6
        assert_eq!(poly_sqnorm(&f), 6);
    }
}
