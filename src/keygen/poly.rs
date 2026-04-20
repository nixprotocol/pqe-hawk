//! Small-polynomial helpers for the keygen path.
//!
//! Port of the first few functions in c-reference/hawk-512/ng_poly.c.

use crate::keygen::mp31::mp_set;

/// Convert a polynomial of int8 coefficients to the [0, p) representation.
///
/// Port of `poly_mp_set_small` (ng_poly.c:3-12).
pub fn poly_mp_set_small(logn: u32, d: &mut [u32], f: &[i8], p: u32) {
    let n = 1usize << logn;
    debug_assert!(d.len() >= n);
    debug_assert!(f.len() >= n);
    for u in 0..n {
        d[u] = mp_set(f[u] as i32, p);
    }
}

// ============================================================
// Sub-task 15.19a: divrem31, poly_max_bitlength, poly_big_to_fixed,
//                  poly_sub_scaled (general case)
// ============================================================

/// Divide `x` by 31 using a magic-multiply reciprocal.
/// Returns `(q, r)` where `q = x / 31` and `r = x % 31`.
///
/// Port of the `DIVREM31` macro (ng_inner.h:1146-1152).
#[inline]
pub fn divrem31(x: u32) -> (u32, u32) {
    let q = x.wrapping_mul(67651) >> 21;
    let r = x - 31u32.wrapping_mul(q);
    (q, r)
}

/// Compute the maximum bit-length across all coefficients of a big-int polynomial.
///
/// `f` contains `n * flen` limbs laid out as `f[coeff_index + limb_index * n]`
/// (stride-`n` limb order). Returns 0 if `flen == 0`.
///
/// Port of `poly_max_bitlength` (ng_poly.c:55-91).
pub fn poly_max_bitlength(logn: u32, f: &[u32], flen: usize) -> u32 {
    use crate::keygen::mp31::tbmask;
    use crate::keygen::zint31::lzcnt;

    if flen == 0 {
        return 0;
    }
    let n = 1usize << logn;
    let mut t: u32 = 0;
    let mut tk: u32 = 0;
    for u in 0..n {
        // Sign bit of the top limb (limb index flen-1) for coefficient u.
        // In the C, `f` advances by 1 each outer iteration so `f[(flen-1) << logn]`
        // accesses coefficient u of the top limb.
        let top_limb = f[(flen - 1) * n + u];
        let m: u32 = (top_limb >> 30).wrapping_neg() & 0x7FFF_FFFF;

        // Find the top non-zero sign-adjusted word and its index.
        let mut c: u32 = 0;
        let mut ck: u32 = 0;
        for v in 0..flen {
            let w = f[v * n + u] ^ m;
            // nz = all-ones if w != 0, else 0
            let nz = ((w.wrapping_sub(1)) >> 31).wrapping_sub(1);
            c ^= nz & (c ^ w);
            ck ^= nz & (ck ^ (v as u32));
        }

        // Replace (t, tk) with (c, ck) if (ck, c) > (tk, t) lexicographically.
        let rr = tbmask((tk.wrapping_sub(ck)) | ((tk ^ ck).wrapping_sub(1) & t.wrapping_sub(c)));
        t ^= rr & (t ^ c);
        tk ^= rr & (tk ^ ck);
    }

    31u32
        .wrapping_mul(tk)
        .wrapping_add(32u32.wrapping_sub(lzcnt(t)))
}

/// Convert a big-int polynomial to a fixed-point (Q32.32) polynomial.
///
/// Each coefficient of `f` (stored with stride `n`, `len` limbs per coeff) is
/// converted to a `Fxr` value scaled by `2^(-sc)`. If `len == 0`, `d` is zeroed.
///
/// Port of `poly_big_to_fixed` (ng_poly.c:94-170).
pub fn poly_big_to_fixed(
    logn: u32,
    d: &mut [crate::keygen::fxr::Fxr],
    f: &[u32],
    len: usize,
    sc: u32,
) {
    use crate::keygen::fxr::Fxr;

    let n = 1usize << logn;

    if len == 0 {
        for x in d[..n].iter_mut() {
            *x = Fxr(0);
        }
        return;
    }

    // Split sc = 31*sch + scl with scl in 1..31.
    let (mut sch, mut scl) = divrem31(sc);
    // If scl == 0, adjust so scl becomes 31 and sch decreases by 1.
    let z = scl.wrapping_sub(1) >> 31; // 1 if scl == 0, else 0
    sch = sch.wrapping_sub(z);
    scl |= 31 & z.wrapping_neg();

    // Masked indices for the three limbs we need (sch-1, sch, sch+1).
    let t0 = sch.wrapping_sub(1) & 0xFFFFFF;
    let t1 = sch & 0xFFFFFF;
    let t2 = sch.wrapping_add(1) & 0xFFFFFF;

    for u in 0..n {
        let mut w0: u32 = 0;
        let mut w1: u32 = 0;
        let mut w2: u32 = 0;

        for v in 0..len {
            let w = f[v * n + u];
            let vv = (v as u32) & 0xFFFFFF;
            w0 |= w & ((vv ^ t0).wrapping_sub(1) >> 31).wrapping_neg();
            w1 |= w & ((vv ^ t1).wrapping_sub(1) >> 31).wrapping_neg();
            w2 |= w & ((vv ^ t2).wrapping_sub(1) >> 31).wrapping_neg();
        }

        // Sign-extend for limbs beyond the array (use the sign word).
        let ws: u32 = (f[(len - 1) * n + u] >> 30).wrapping_neg() >> 1;
        w0 |= ws & ((len as u32).wrapping_sub(sch) >> 31).wrapping_neg();
        w1 |= ws & ((len as u32).wrapping_sub(sch).wrapping_sub(1) >> 31).wrapping_neg();
        w2 |= ws & ((len as u32).wrapping_sub(sch).wrapping_sub(2) >> 31).wrapping_neg();

        // Sign-extend w2 to bit 31.
        w2 |= (w2 & 0x4000_0000) << 1;

        // Assemble the 64-bit Q32.32 value.
        let xl = (w0 >> (scl - 1)) | (w1 << (32 - scl));
        let xh = (w1 >> scl) | (w2 << (31 - scl));
        d[u] = Fxr::of_scaled32((xl as u64) | ((xh as u64) << 32));
    }
}

/// Subtract `k * f * 2^sc` from `F` (big-int polynomial).
///
/// General case only (logn >= 4). Both `F` and `f` use stride-n limb layout.
/// `sc` is the total bit-shift count; it is split into 31*sch + scl internally.
///
/// Port of `poly_sub_scaled` general case (ng_poly.c:484-498) plus prologue
/// (ng_poly.c:178-188). The logn=1/2/3 optimised switch arms are not ported
/// because HAWK-512 always uses logn=9.
pub fn poly_sub_scaled(
    logn: u32,
    f_big: &mut [u32],
    flen_f: usize,
    f: &[u32],
    flen: usize,
    k: &[i32],
    sc: u32,
) {
    use crate::keygen::zint31::zint_add_scaled_mul_small;

    if flen == 0 {
        return;
    }
    let (sch, scl) = divrem31(sc);
    let sch = sch as usize;
    if sch >= flen_f {
        return;
    }

    let n = 1usize << logn;
    let shifted_len = flen_f - sch;

    // Equivalent to `F += sch << logn` in C: skip the first sch*n limbs.
    let f_big = &mut f_big[sch * n..];

    // General case (ng_poly.c:484-498):
    for u in 0..n {
        let mut kf: i32 = -(k[u]);
        // x starts at f_big[u]; the C increments x by 1 each inner step,
        // which maps to x_off advancing by 1.
        let mut x_off = u;
        for v in 0..n {
            // Both x and y are addressed with stride n.
            zint_add_scaled_mul_small(
                &mut f_big[x_off..],
                shifted_len,
                &f[v..],
                flen,
                n,
                kf,
                0,
                scl,
            );
            if u + v == n - 1 {
                x_off = 0;
                kf = -kf;
            } else {
                x_off += 1;
            }
        }
    }
}

// ============================================================
// Sub-task 15.19b: poly_mp_set, poly_mp_norm, poly_big_to_small,
//                  poly_sub_scaled_ntt, poly_sub_kf_scaled_depth1
// ============================================================

/// Convert a polynomial from signed 31-bit two's complement to [0, p).
/// Port of `poly_mp_set` (ng_poly.c:14-24).
pub fn poly_mp_set(logn: u32, f: &mut [u32], p: u32) {
    use crate::keygen::mp31::mp_set;
    let n = 1usize << logn;
    for u in 0..n {
        let mut x = f[u];
        x |= (x & 0x4000_0000) << 1;
        f[u] = mp_set(x as i32, p);
    }
}

/// Convert a polynomial from [0, p) to signed 31-bit two's complement.
/// Port of `poly_mp_norm` (ng_poly.c:26-34).
pub fn poly_mp_norm(logn: u32, f: &mut [u32], p: u32) {
    use crate::keygen::mp31::mp_norm;
    let n = 1usize << logn;
    for u in 0..n {
        f[u] = (mp_norm(f[u], p) as u32) & 0x7FFF_FFFF;
    }
}

/// Convert a signed 31-bit polynomial to int8, checking that every
/// coefficient is within `[-lim, +lim]`. Returns `true` on success.
/// Port of `poly_big_to_small` (ng_poly.c:37-52).
pub fn poly_big_to_small(logn: u32, d: &mut [i8], s: &[u32], lim: i32) -> bool {
    let n = 1usize << logn;
    for u in 0..n {
        let mut x = s[u];
        x |= (x & 0x4000_0000) << 1;
        let z = x as i32;
        if z < -lim || z > lim {
            return false;
        }
        d[u] = z as i8;
    }
    true
}

/// NTT-based variant of `poly_sub_scaled`, faster for large logn. Uses
/// `tmp` as scratch (needs ~(tlen+2)*n u32 words).
///
/// Port of `poly_sub_scaled_ntt` (ng_poly.c:502-551).
pub fn poly_sub_scaled_ntt(
    logn: u32,
    f_big: &mut [u32],
    flen_f_big: usize,
    f: &[u32],
    flen: usize,
    k: &[i32],
    sc: u32,
    tmp: &mut [u32],
) {
    use crate::keygen::mp31::{mp_intt, mp_mkgmigm, mp_montymul, mp_ntt, mp_set};
    use crate::keygen::primes::PRIMES;
    use crate::keygen::zint31::{zint_rebuild_crt, zint_sub_scaled};

    let n = 1usize << logn;
    let tlen = flen + 1;
    let (sch, scl) = divrem31(sc);

    // Layout in tmp: [gm(n) | igm(n) | fk(tlen*n) | t1(n)].
    let (gm, rest1) = tmp.split_at_mut(n);
    let (igm, rest2) = rest1.split_at_mut(n);
    let (fk, t1) = rest2.split_at_mut(tlen * n);

    for u in 0..tlen {
        let p = PRIMES[u].p;
        let p0i = PRIMES[u].p0i;
        let r2 = PRIMES[u].r2;
        mp_mkgmigm(logn, gm, igm, PRIMES[u].g, PRIMES[u].ig, p, p0i);
        for v in 0..n {
            t1[v] = mp_set(k[v], p);
        }
        mp_ntt(logn, &mut t1[..n], gm, p, p0i);

        let fs_base = u << logn;
        for v in 0..n {
            fk[fs_base + v] = mp_montymul(mp_montymul(t1[v], f[fs_base + v], p, p0i), r2, p, p0i);
        }
        let ff_slice = &mut fk[fs_base..fs_base + n];
        mp_intt(logn, ff_slice, igm, p, p0i);
    }

    zint_rebuild_crt(fk, tlen, n, 1, true, t1);

    for u in 0..n {
        zint_sub_scaled(&mut f_big[u..], flen_f_big, &fk[u..], tlen, n, sch, scl);
    }
}

/// Specialized `poly_sub_scaled` for depth 1 with f as int8. Used by
/// `solve_NTRU_intermediate` when depth == 1 for faster operation with
/// better precision.
///
/// Port of `poly_sub_kf_scaled_depth1` (ng_poly.c:554-735).
pub fn poly_sub_kf_scaled_depth1(
    logn_top: u32,
    f_big: &mut [u32],
    fglen: usize,
    k: &mut [u32],
    sc: u32,
    f: &[i8],
    tmp: &mut [u32],
) {
    use crate::keygen::mp31::{
        mp_add, mp_half, mp_intt, mp_mkgm, mp_mkigm, mp_montymul, mp_norm, mp_ntt, mp_set, mp_sub,
        tbmask,
    };
    use crate::keygen::primes::PRIMES;

    let logn = logn_top - 1;
    let n = 1usize << logn;
    let hn = n >> 1;

    // Layout in tmp: [gm(n) | t1(n) | t2(n)].
    let (gm, rest) = tmp.split_at_mut(n);
    let (t1, t2) = rest.split_at_mut(n);

    // Step 1: convert F to RNS.
    if fglen == 1 {
        let p = PRIMES[0].p;
        for u in 0..n {
            let mut xf = f_big[u];
            xf |= (xf & 0x4000_0000) << 1;
            f_big[u] = mp_set(xf as i32, p);
        }
    } else {
        // fglen == 2
        let p0 = PRIMES[0].p;
        let p0_0i = PRIMES[0].p0i;
        let z0 = mp_half(PRIMES[0].r2, p0);
        let p1 = PRIMES[1].p;
        let p1_0i = PRIMES[1].p0i;
        let z1 = mp_half(PRIMES[1].r2, p1);
        for u in 0..n {
            let xl = f_big[u];
            let xh_raw = f_big[u + n];
            let xh = xh_raw | ((xh_raw & 0x4000_0000) << 1);
            let yl0 = xl.wrapping_sub(p0 & !tbmask(xl.wrapping_sub(p0)));
            let yh0 = mp_set(xh as i32, p0);
            let r0 = mp_add(yl0, mp_montymul(yh0, z0, p0, p0_0i), p0);
            let yl1 = xl.wrapping_sub(p1 & !tbmask(xl.wrapping_sub(p1)));
            let yh1 = mp_set(xh as i32, p1);
            let r1 = mp_add(yl1, mp_montymul(yh1, z1, p1, p1_0i), p1);
            f_big[u] = r0;
            f_big[u + n] = r1;
        }
    }

    // Step 2: for each prime slot, subtract (2^sc)*k*f from F.
    for u in 0..fglen {
        let p = PRIMES[u].p;
        let p0i = PRIMES[u].p0i;
        let r2 = PRIMES[u].r2;
        let r3 = mp_montymul(r2, r2, p, p0i);
        mp_mkgm(logn, gm, PRIMES[u].g, p, p0i);

        // Scale k by 2^sc and transform to NTT.
        let mut scv = mp_montymul(1u32 << (sc & 31), r2, p, p0i);
        let mut mm = sc >> 5;
        while mm > 0 {
            scv = mp_montymul(scv, r2, p, p0i);
            mm -= 1;
        }
        for v in 0..n {
            let x = mp_set(k[v] as i32, p);
            k[v] = mp_montymul(scv, x, p, p0i);
        }
        mp_ntt(logn, k, gm, p, p0i);

        // Transform F[u*n..] to NTT.
        let fu_base = u << logn;
        let fu = &mut f_big[fu_base..fu_base + n];
        mp_ntt(logn, fu, gm, p, p0i);

        // Load f (interleaved even/odd from logn_top) into t1, t2 and NTT.
        for v in 0..n {
            t1[v] = mp_set(f[(v << 1)] as i32, p);
            t2[v] = mp_set(f[(v << 1) + 1] as i32, p);
        }
        mp_ntt(logn, t1, gm, p, p0i);
        mp_ntt(logn, t2, gm, p, p0i);

        // Inner loop: compute f(x)*f(x) at each NTT point using
        // f(x) = fe(x^2) + x*fo(x^2), then subtract k*f*f from F.
        for v in 0..hn {
            let xe0 = t1[(v << 1)];
            let xe1 = t1[(v << 1) + 1];
            let xo0 = t2[(v << 1)];
            let xo1 = t2[(v << 1) + 1];
            let xv0 = gm[hn + v];
            let xv1 = p.wrapping_sub(xv0);
            let xe0 = mp_montymul(xe0, xe0, p, p0i);
            let xe1 = mp_montymul(xe1, xe1, p, p0i);
            let xo0 = mp_montymul(xo0, xo0, p, p0i);
            let xo1 = mp_montymul(xo1, xo1, p, p0i);
            let xf0 = mp_sub(xe0, mp_montymul(xo0, xv0, p, p0i), p);
            let xf1 = mp_sub(xe1, mp_montymul(xo1, xv1, p, p0i), p);
            let xkf0 = mp_montymul(mp_montymul(xf0, k[(v << 1)], p, p0i), r3, p, p0i);
            let xkf1 = mp_montymul(mp_montymul(xf1, k[(v << 1) + 1], p, p0i), r3, p, p0i);
            fu[(v << 1)] = mp_sub(fu[(v << 1)], xkf0, p);
            fu[(v << 1) + 1] = mp_sub(fu[(v << 1) + 1], xkf1, p);
        }

        // Inverse NTT F back to RNS.
        mp_mkigm(logn, t1, PRIMES[u].ig, p, p0i);
        mp_intt(logn, fu, t1, p, p0i);

        // If there is a next iteration, restore k to its pre-scaled form.
        if (u + 1) < fglen {
            mp_intt(logn, k, t1, p, p0i);
            // Divide k by 2^sc: multiply by 2^(-sc) = 2^(32 - (sc & 31)) * R^(-(sc>>5)).
            let mut scv = 1u32 << (sc.wrapping_neg() & 31);
            let mut mm = sc >> 5;
            while mm > 0 {
                scv = mp_montymul(scv, 1, p, p0i);
                mm -= 1;
            }
            for v in 0..n {
                // Note: C stores the raw int32->uint32 cast (no & 0x7FFF_FFFF) so that
                // the next iteration's mp_set(*(int32_t*)&k[v]) can recover the negative
                // value. We do the same: cast to u32 without masking.
                k[v] = mp_norm(mp_montymul(scv, k[v], p, p0i), p) as u32;
            }
        }
    }

    // Step 3: convert RNS back to plain signed 31-bit integers.
    if fglen == 1 {
        let p = PRIMES[0].p;
        for u in 0..n {
            f_big[u] = (mp_norm(f_big[u], p) as u32) & 0x7FFF_FFFF;
        }
    } else {
        // fglen == 2: CRT reconstruction with two primes.
        let p0 = PRIMES[0].p;
        let p1 = PRIMES[1].p;
        let p1_0i = PRIMES[1].p0i;
        let s = PRIMES[1].s;
        let pp = (p0 as u64).wrapping_mul(p1 as u64);
        let hpp = pp >> 1;
        for u in 0..n {
            let x0 = f_big[u];
            let x1 = f_big[u + n];
            let x0m1 = x0.wrapping_sub(p1 & !tbmask(x0.wrapping_sub(p1)));
            let y = mp_montymul(mp_sub(x1, x0m1, p1), s, p1, p1_0i);
            let mut z = (x0 as u64).wrapping_add((p0 as u64).wrapping_mul(y as u64));
            // Normalize to [-pp/2, pp/2].
            z = z.wrapping_sub(pp & (hpp.wrapping_sub(z) >> 63).wrapping_neg());
            f_big[u] = (z as u32) & 0x7FFF_FFFF;
            f_big[u + n] = ((z >> 31) as u32) & 0x7FFF_FFFF;
        }
    }
}

// ============================================================
// Sub-task 15.24: poly_sqnorm
// ============================================================

/// Compute the squared norm of an int8 polynomial.
/// Port of `poly_sqnorm` (ng_poly.c:813-823).
pub fn poly_sqnorm(f: &[i8]) -> u32 {
    let mut s: u32 = 0;
    for &c in f.iter() {
        let x = c as i32;
        s = s.wrapping_add((x.wrapping_mul(x)) as u32);
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::keygen::primes::PRIMES;

    #[test]
    fn poly_mp_set_small_basic() {
        // Small polynomial [1, -1, 2, 0] at logn=2.
        let coeffs: [i8; 4] = [1, -1, 2, 0];
        let mut d = [0u32; 4];
        let p = PRIMES[0].p;
        poly_mp_set_small(2, &mut d, &coeffs, p);
        // 1 mod p = 1
        assert_eq!(d[0], 1);
        // -1 mod p = p - 1
        assert_eq!(d[1], p - 1);
        // 2 mod p = 2
        assert_eq!(d[2], 2);
        // 0 mod p = 0
        assert_eq!(d[3], 0);
    }

    // === 15.19a tests ===

    #[test]
    fn divrem31_basic() {
        assert_eq!(divrem31(0), (0, 0));
        assert_eq!(divrem31(31), (1, 0));
        assert_eq!(divrem31(30), (0, 30));
        assert_eq!(divrem31(62), (2, 0));
        assert_eq!(divrem31(100), (3, 7));
        assert_eq!(divrem31(1000), (32, 8));
    }

    #[test]
    fn divrem31_invariant() {
        // q * 31 + r == x and 0 <= r < 31.
        // The magic multiply is valid for x < 63506 (x * 67651 < 2^32).
        for x in [0u32, 1, 30, 31, 32, 62, 99, 100, 255, 999, 1000, 63000] {
            let (q, r) = divrem31(x);
            assert_eq!(q * 31 + r, x, "divrem31({}) failed", x);
            assert!(r < 31, "remainder {} >= 31 for input {}", r, x);
        }
    }

    #[test]
    fn poly_max_bitlength_empty() {
        let f = [0u32; 512];
        assert_eq!(poly_max_bitlength(9, &f, 0), 0);
    }

    #[test]
    fn poly_max_bitlength_all_zero() {
        let f = [0u32; 512];
        assert_eq!(poly_max_bitlength(9, &f, 1), 0);
    }

    #[test]
    fn poly_big_to_fixed_zero_len() {
        let f = [0u32; 512];
        let mut d = [crate::keygen::fxr::Fxr(0); 512];
        poly_big_to_fixed(9, &mut d, &f, 0, 0);
        for fx in d.iter() {
            assert_eq!(fx.0, 0);
        }
    }

    // === 15.19b tests ===

    #[test]
    fn poly_mp_set_norm_roundtrip() {
        let p = PRIMES[0].p;
        // Use logn=3 (n=8). Construct original values as signed 31-bit words.
        // The encoding: non-negative values are stored as-is (bit 30 = 0).
        // Negative values: bit 30 is the sign bit (two's complement within 31 bits).
        //   -1 as 31-bit two's comp = 0x7FFFFFFF (all ones in 31 bits).
        //   -2 as 31-bit two's comp = 0x7FFFFFFE.
        // Only values whose magnitude < p/2 (≈1073736704) roundtrip cleanly.
        // 0x7FFFFFFF encodes -1, 0x7FFFFFFE encodes -2 (bit-30 sign bit).
        // All magnitudes here are small enough to survive mp_norm.
        let original = [0u32, 1, 2, 127, 0x7FFFFFFF, 0x7FFFFFFE, 0x7FFFFF80, 10];
        let mut f8 = original;
        poly_mp_set(3, &mut f8, p);
        // All coefficients must be in [0, p).
        for i in 0..8 {
            assert!(f8[i] < p, "index {} out of range: {}", i, f8[i]);
        }
        poly_mp_norm(3, &mut f8, p);
        // Round-trip must recover the original 31-bit encoding.
        for i in 0..8 {
            assert_eq!(f8[i], original[i], "roundtrip failed at index {}", i);
        }
    }

    #[test]
    fn poly_big_to_small_accepts_inside_limit() {
        // Encode: 5=5, 10=10, -2=0x7FFFFFFE, 3=3 as 31-bit two's complement.
        let s = [5u32, 10, 0x7FFFFFFE, 3];
        let mut d = [0i8; 4];
        let ok = poly_big_to_small(2, &mut d, &s, 10);
        assert!(ok);
        assert_eq!(d, [5, 10, -2, 3]);
    }

    #[test]
    fn poly_big_to_small_rejects_out_of_range() {
        // 20 > lim=10, so the function must return false.
        let s = [5u32, 20, 3, 3];
        let mut d = [0i8; 4];
        let ok = poly_big_to_small(2, &mut d, &s, 10);
        assert!(!ok);
    }
}
