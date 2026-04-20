//! HAWK public-key polynomial (Q-poly) derivation from the secret key.
//!
//! Given (f, g, F, G), compute:
//!   - q00 = f*adj(f) + g*adj(g)      (auto-adjoint, i16 coefficients)
//!   - q01 = F*adj(f) + G*adj(g)      (i16 coefficients)
//!   - q11 = F*adj(F) + G*adj(G)      (auto-adjoint, i32 coefficients)
//!
//! Port of c-reference/hawk-512/ng_hawk.c:274-423 (`make_q001`).

use crate::keygen::mp31::{mp_add, mp_intt, mp_mkgm, mp_mkigm, mp_montymul, mp_norm, mp_ntt};
use crate::keygen::poly::poly_mp_set_small;
use crate::keygen::primes::PRIMES;

/// Maximum bit-size of q00 coefficients (excluding q00[0]).
/// Port of `bits_lim00` (ng_hawk.c:242-244).
pub const BITS_LIM00: [i8; 11] = [0, 0, 0, 0, 0, 0, 0, 0, 9, 9, 10];

/// Maximum bit-size of q01 coefficients.
/// Port of `bits_lim01` (ng_hawk.c:245-247).
pub const BITS_LIM01: [i8; 11] = [0, 0, 0, 0, 0, 0, 0, 0, 11, 12, 14];

/// Maximum bit-size of q11 coefficients (excluding q11[0]).
/// Port of `bits_lim11` (ng_hawk.c:248-250).
pub const BITS_LIM11: [i8; 11] = [0, 0, 0, 0, 0, 0, 0, 0, 13, 15, 17];

/// Output of [`make_q001`]: three Q-polynomials plus success flag.
///
/// The C packs q00 (i16), q01 (i16), q11 (i32) into a u32 buffer via pointer
/// casts. This Rust struct exposes them as clean typed arrays.
pub struct MakeQ001Output {
    pub q00: Vec<i16>,
    pub q01: Vec<i16>,
    pub q11: Vec<i32>,
}

/// Compute Q-polynomials from (f, g, F, G).
///
/// On success, returns `Ok(MakeQ001Output)` with q00, q01, q11 all of length n.
/// On failure (coefficient out of range or q00 not invertible), returns `Err(())`.
///
/// `lim00`, `lim01`, `lim11` are `1 << bits_lim{00,01,11}[logn]`.
///
/// Port of `make_q001` (ng_hawk.c:274-423).
///
/// Returns `Result<_, ()>` because the reference C returns a bare int flag
/// with no error detail; the sole caller (`hawk_keygen`) interprets `Err`
/// as "retry with a fresh (f, g)" and carries no payload that a richer
/// error type would expose.
#[allow(clippy::result_unit_err)]
pub fn make_q001(
    logn: u32,
    lim00: i32,
    lim01: i32,
    lim11: i32,
    f: &[i8],
    g: &[i8],
    f_cap: &[i8],
    g_cap: &[i8],
    tmp: &mut [u32],
) -> Result<MakeQ001Output, ()> {
    let n = 1usize << logn;
    let hn = n >> 1;
    let p = PRIMES[0].p;
    let p0i = PRIMES[0].p0i;
    let r2 = PRIMES[0].r2;

    // Layout: tmp = [t1 | t2 | t3 | t4 | t5] where t1..t4 are n u32, t5 is hn u32.
    assert!(
        tmp.len() >= 4 * n + hn,
        "tmp too small: {} < {}",
        tmp.len(),
        4 * n + hn
    );

    // Step 1: mkgm into t1. Then set t2=f, t3=g, t4=F and NTT all three.
    {
        let (t1, rest) = tmp.split_at_mut(n);
        let (t2, rest2) = rest.split_at_mut(n);
        let (t3, rest3) = rest2.split_at_mut(n);
        let (t4, _) = rest3.split_at_mut(n);
        mp_mkgm(logn, t1, PRIMES[0].g, p, p0i);
        poly_mp_set_small(logn, t2, f, p);
        poly_mp_set_small(logn, t3, g, p);
        poly_mp_set_small(logn, t4, f_cap, p);
        mp_ntt(logn, t2, t1, p, p0i);
        mp_ntt(logn, t3, t1, p, p0i);
        mp_ntt(logn, t4, t1, p, p0i);
    }

    // Step 2: Compute t2 = F*adj(f) and t5 = q00 = f*adj(f) + g*adj(g)
    // (auto-adjoint, half-size). For u in 0..hn.
    {
        let (_t1, rest) = tmp.split_at_mut(n);
        let (t2, rest2) = rest.split_at_mut(n);
        let (t3, rest3) = rest2.split_at_mut(n);
        let (t4, t5) = rest3.split_at_mut(n);
        for u in 0..hn {
            let xf = t2[u];
            let xfa = t2[(n - 1) - u];
            let xg = t3[u];
            let xga = t3[(n - 1) - u];
            let x_f_cap = t4[u];
            let x_f_cap_a = t4[(n - 1) - u];

            let new_t2_u = mp_montymul(r2, mp_montymul(x_f_cap, xfa, p, p0i), p, p0i);
            let new_t2_mirror = mp_montymul(r2, mp_montymul(x_f_cap_a, xf, p, p0i), p, p0i);
            let xq00 = mp_montymul(
                r2,
                mp_add(
                    mp_montymul(xf, xfa, p, p0i),
                    mp_montymul(xg, xga, p, p0i),
                    p,
                ),
                p,
                p0i,
            );
            if xq00 == 0 {
                // q00 is not invertible mod X^n+1 mod p.
                return Err(());
            }
            t2[u] = new_t2_u;
            t2[n - 1 - u] = new_t2_mirror;
            t5[u] = xq00;
        }
    }

    // Step 3: t4 = G (set + NTT). Then update t2 = F*adj(f) + G*adj(g).
    {
        let (_t1, rest) = tmp.split_at_mut(n);
        let (_t2, rest2) = rest.split_at_mut(n);
        let (_t3, rest3) = rest2.split_at_mut(n);
        let (t4, _t5) = rest3.split_at_mut(n);
        poly_mp_set_small(logn, t4, g_cap, p);
    }
    {
        let (t1, rest) = tmp.split_at_mut(n);
        let (_t2, rest2) = rest.split_at_mut(n);
        let (_t3, rest3) = rest2.split_at_mut(n);
        let (t4, _t5) = rest3.split_at_mut(n);
        mp_ntt(logn, t4, t1, p, p0i);
    }
    {
        let (_t1, rest) = tmp.split_at_mut(n);
        let (t2, rest2) = rest.split_at_mut(n);
        let (t3, rest3) = rest2.split_at_mut(n);
        let (t4, _t5) = rest3.split_at_mut(n);
        for u in 0..hn {
            let xg = t3[u];
            let xga = t3[(n - 1) - u];
            let x_g_cap = t4[u];
            let x_g_cap_a = t4[(n - 1) - u];
            t2[u] = mp_add(
                mp_montymul(r2, mp_montymul(x_g_cap, xga, p, p0i), p, p0i),
                t2[u],
                p,
            );
            t2[n - 1 - u] = mp_add(
                mp_montymul(r2, mp_montymul(x_g_cap_a, xg, p, p0i), p, p0i),
                t2[n - 1 - u],
                p,
            );
        }
    }

    // Step 4: Reload F into t3 and NTT. Compute q11 = F*adj(F) + G*adj(G) in t3.
    {
        let (_t1, rest) = tmp.split_at_mut(n);
        let (_t2, rest2) = rest.split_at_mut(n);
        let (t3, _) = rest2.split_at_mut(n);
        poly_mp_set_small(logn, t3, f_cap, p);
    }
    {
        let (t1, rest) = tmp.split_at_mut(n);
        let (_t2, rest2) = rest.split_at_mut(n);
        let (t3, _) = rest2.split_at_mut(n);
        mp_ntt(logn, t3, t1, p, p0i);
    }
    {
        let (_t1, rest) = tmp.split_at_mut(n);
        let (_t2, rest2) = rest.split_at_mut(n);
        let (t3, rest3) = rest2.split_at_mut(n);
        let (t4, _t5) = rest3.split_at_mut(n);
        for u in 0..hn {
            let x_f_cap = t3[u];
            let x_f_cap_a = t3[(n - 1) - u];
            let x_g_cap = t4[u];
            let x_g_cap_a = t4[(n - 1) - u];
            let xq11 = mp_montymul(
                r2,
                mp_add(
                    mp_montymul(x_f_cap, x_f_cap_a, p, p0i),
                    mp_montymul(x_g_cap, x_g_cap_a, p, p0i),
                    p,
                ),
                p,
                p0i,
            );
            t3[u] = xq11;
            t3[(n - 1) - u] = xq11;
        }
    }

    // Step 5: Expand t5 (half-size) into t4 (full-size), then mkigm into t1,
    // iNTT t2, t3, t4.
    {
        let (t1, rest) = tmp.split_at_mut(n);
        let (_t2, rest2) = rest.split_at_mut(n);
        let (_t3, rest3) = rest2.split_at_mut(n);
        let (t4, t5) = rest3.split_at_mut(n);
        for u in 0..hn {
            t4[u] = t5[u];
            t4[n - 1 - u] = t5[u];
        }
        mp_mkigm(logn, t1, PRIMES[0].ig, p, p0i);
    }
    {
        let (t1, rest) = tmp.split_at_mut(n);
        let (t2, rest2) = rest.split_at_mut(n);
        let (t3, rest3) = rest2.split_at_mut(n);
        let (t4, _t5) = rest3.split_at_mut(n);
        mp_intt(logn, t2, t1, p, p0i);
        mp_intt(logn, t3, t1, p, p0i);
        mp_intt(logn, t4, t1, p, p0i);
    }

    // Step 6: Extract + range-check q00 (from t4), q01 (from t2), q11 (from t3).
    let mut q00 = vec![0i16; n];
    let mut q01 = vec![0i16; n];
    let mut q11 = vec![0i32; n];
    {
        let (_t1, rest) = tmp.split_at_mut(n);
        let (t2, rest2) = rest.split_at_mut(n);
        let (t3, rest3) = rest2.split_at_mut(n);
        let (t4, _t5) = rest3.split_at_mut(n);
        for u in 0..n {
            let xq00 = mp_norm(t4[u], p);
            let xq01 = mp_norm(t2[u], p);
            let xq11 = mp_norm(t3[u], p);
            if u == 0 {
                if !(-32768..=32767).contains(&xq00) {
                    return Err(());
                }
            } else {
                if xq00 <= -lim00 || xq00 >= lim00 {
                    return Err(());
                }
                if xq11 <= -lim11 || xq11 >= lim11 {
                    return Err(());
                }
            }
            if xq01 <= -lim01 || xq01 >= lim01 {
                return Err(());
            }
            q00[u] = xq00 as i16;
            q01[u] = xq01 as i16;
            q11[u] = xq11;
        }
    }

    Ok(MakeQ001Output { q00, q01, q11 })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bits_lim_tables_match_c() {
        // Cross-reference ng_hawk.c:242-250.
        assert_eq!(BITS_LIM00, [0, 0, 0, 0, 0, 0, 0, 0, 9, 9, 10]);
        assert_eq!(BITS_LIM01, [0, 0, 0, 0, 0, 0, 0, 0, 11, 12, 14]);
        assert_eq!(BITS_LIM11, [0, 0, 0, 0, 0, 0, 0, 0, 13, 15, 17]);
    }

    #[test]
    fn bits_lim_tables_len() {
        assert_eq!(BITS_LIM00.len(), 11);
        assert_eq!(BITS_LIM01.len(), 11);
        assert_eq!(BITS_LIM11.len(), 11);
    }

    #[test]
    fn hawk512_lim_values() {
        // logn=9: lim00 = 1<<9 = 512, lim01 = 1<<12 = 4096, lim11 = 1<<15 = 32768.
        let logn = 9usize;
        assert_eq!(1i32 << BITS_LIM00[logn], 512);
        assert_eq!(1i32 << BITS_LIM01[logn], 4096);
        assert_eq!(1i32 << BITS_LIM11[logn], 32768);
    }
}
