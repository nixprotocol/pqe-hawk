//! Multi-precision big-integer arithmetic for HAWK keygen.
//!
//! Big integers are stored as arrays of `u32` limbs, where each limb holds
//! only the low 31 bits (top bit cleared). Signed big-ints use two's
//! complement over the whole array (bit 30 of the last limb is the sign
//! bit). "Stride" is an addressing parameter used by the reference C to
//! interleave multiple big-ints in a single contiguous buffer.
//!
//! Port of the reference C in c-reference/hawk-512/ng_zint31.c. This file
//! contains the primitive operations (mul_small, mod_small_unsigned,
//! add_mul_small, norm_zero). Higher-level operations (CRT rebuild,
//! Bezout, scaled add/sub) come in later sub-tasks.

// Functions are pub to allow FFI cross-check tests (tests/cross_check.rs) to
// reach them; no stability guarantees beyond this crate's test suite.

/// Multiply the big-int `m` (stride 1) by a small u32 value `x` (assumed `< 2^31`).
/// Returns the carry.
/// Port of `zint_mul_small` (ng_zint31.c:22-32).
pub fn zint_mul_small(m: &mut [u32], x: u32) -> u32 {
    let mut cc: u32 = 0;
    for limb in m.iter_mut() {
        let z = (*limb as u64) * (x as u64) + cc as u64;
        *limb = (z as u32) & 0x7FFF_FFFF;
        cc = (z >> 31) as u32;
    }
    cc
}

/// Reduce the big-int stored at `d[0], d[stride], d[2*stride], ..., d[(len-1)*stride]`
/// modulo the small prime `p`.
///
/// Port of `zint_mod_small_unsigned` (ng_zint31.c:35-56). `d` is viewed as a
/// sequence of `len` limbs with a row stride.
pub fn zint_mod_small_unsigned(
    d: &[u32],
    len: usize,
    stride: usize,
    p: u32,
    p0i: u32,
    r2: u32,
) -> u32 {
    use crate::keygen::mp31::{mp_add, mp_half, mp_montymul, tbmask};
    let mut x: u32 = 0;
    let z = mp_half(r2, p);
    // Process limbs from highest-index to lowest (high limb first).
    let mut u = len;
    while u > 0 {
        u -= 1;
        let limb = d[u * stride];
        // w = limb - p; if limb < p (so the subtraction underflowed), add p back.
        let w0 = limb.wrapping_sub(p);
        let w = w0.wrapping_add(p & tbmask(w0));
        x = mp_montymul(x, z, p, p0i);
        x = mp_add(x, w, p);
    }
    x
}

/// Reduce a **signed** big-int modulo a small prime `p`.
/// `Rx = 2^(31*len) mod p` (precomputed by the caller via `mp_rx31`).
///
/// Port of `zint_mod_small_signed` (ng_inner.h:606-616).
pub fn zint_mod_small_signed(
    d: &[u32],
    len: usize,
    stride: usize,
    p: u32,
    p0i: u32,
    r2: u32,
    rx: u32,
) -> u32 {
    use crate::keygen::mp31::mp_sub;
    if len == 0 {
        return 0;
    }
    let z = zint_mod_small_unsigned(d, len, stride, p, p0i, r2);
    // Sign-adjust: subtract Rx if the high bit of the top limb is set.
    let sign_bit = d[(len - 1) * stride] >> 30;
    mp_sub(z, rx & sign_bit.wrapping_neg(), p)
}

/// Add `s * y` to the big-int `x` (stride-based addressing).
///
/// After this call, `x` has `len + 1` limbs: the low `len` limbs are updated,
/// and the final carry is written at `x[len * xstride]`. `y` has stride 1 and
/// length `len`. Caller must ensure `x.len() >= (len + 1) * xstride` (or
/// equivalently, `x[len * xstride]` is writable).
///
/// Port of `zint_add_mul_small` (ng_zint31.c:58-76).
pub fn zint_add_mul_small(x: &mut [u32], len: usize, xstride: usize, y: &[u32], s: u32) {
    let mut cc: u32 = 0;
    for u in 0..len {
        let xw = x[u * xstride];
        let yw = y[u];
        let z = (yw as u64) * (s as u64) + (xw as u64) + (cc as u64);
        x[u * xstride] = (z as u32) & 0x7FFF_FFFF;
        cc = (z >> 31) as u32;
    }
    x[len * xstride] = cc;
}

/// Conditionally negate the big-int `a` (stride 1) when `ctl == 1`; leave it
/// unchanged when `ctl == 0`. Constant-time.
///
/// Port of `zint_negate` (ng_zint31.c:200-219).
pub fn zint_negate(a: &mut [u32], ctl: u32) {
    let mut cc: u32 = ctl;
    let m: u32 = ctl.wrapping_neg() >> 1;
    for limb in a.iter_mut() {
        let aw = (*limb ^ m).wrapping_add(cc);
        *limb = aw & 0x7FFF_FFFF;
        cc = aw >> 31;
    }
}

/// Rebuild big-ints from their RNS representation across `PRIMES[0..xlen]`.
///
/// `xx` contains `num_sets` sets; each set holds `n` interleaved big-ints of
/// `xlen` limbs each. On input, `xx[u * n + v]` for set `k` at limb index `u`
/// contains the residue of big-int `v` modulo `PRIMES[u].p`. On output,
/// `xx` contains the full big-ints in standard (stride-`n`) limb representation.
///
/// If `normalize_signed` is true, outputs are normalized to `(-m/2, m/2]` via
/// `zint_norm_zero`; otherwise outputs are in `[0, m)` where `m` is the
/// product of the first `xlen` primes.
///
/// `tmp` must have length at least `xlen` (it holds the running product of
/// primes).
///
/// Port of `zint_rebuild_CRT` (ng_zint31.c:133-198).
pub fn zint_rebuild_crt(
    xx: &mut [u32],
    xlen: usize,
    n: usize,
    num_sets: usize,
    normalize_signed: bool,
    tmp: &mut [u32],
) {
    use crate::keygen::mp31::{mp_montymul, mp_sub};
    use crate::keygen::primes::PRIMES;

    assert!(tmp.len() >= xlen);
    assert!(xx.len() >= num_sets * xlen * n);

    let mut uu: usize = 0;
    tmp[0] = PRIMES[0].p;
    for u in 1..xlen {
        let p = PRIMES[u].p;
        let p0i = PRIMES[u].p0i;
        let r2 = PRIMES[u].r2;
        let s = PRIMES[u].s;
        uu += n;
        let mut kk: usize = 0;
        for _k in 0..num_sets {
            for v in 0..n {
                let xp = xx[kk + v + uu];
                let xq = zint_mod_small_unsigned(&xx[kk + v..], u, n, p, p0i, r2);
                let xr = mp_montymul(s, mp_sub(xp, xq, p), p, p0i);
                zint_add_mul_small(&mut xx[kk + v..], u, n, &tmp[..u], xr);
            }
            kk += n * xlen;
        }
        // Update tmp: multiply tmp[0..u] in place and store carry in tmp[u].
        // split_at_mut avoids aliasing tmp[0..u] with tmp[u].
        let (head, tail) = tmp.split_at_mut(u);
        tail[0] = zint_mul_small(head, p);
    }

    if normalize_signed {
        // Build a Vec copy of tmp[0..xlen] so we can borrow xx mutably at the
        // same time. (tmp and xx are separate slices so this is only needed if
        // the borrow checker requires it; an explicit copy is the cleanest way.)
        let p_vec: Vec<u32> = tmp[..xlen].to_vec();
        let mut kk: usize = 0;
        for _k in 0..num_sets {
            for v in 0..n {
                zint_norm_zero(&mut xx[kk + v..], xlen, n, &p_vec);
            }
            kk += n * xlen;
        }
    }
}

/// Normalize the big-int `x` (stride `xstride`) around 0: if `x > p/2`,
/// replace `x` with `x - p`. Operates on `len` limbs. `p` has stride 1.
///
/// Port of `zint_norm_zero` (ng_zint31.c:78-131). The function is constant-
/// time: it computes the comparison via a branchless limb-wise compare and
/// a conditional-subtract loop with a mask.
pub fn zint_norm_zero(x: &mut [u32], len: usize, xstride: usize, p: &[u32]) {
    use crate::keygen::mp31::tbmask;

    // Compare x with (p-1)/2, limb by limb from high to low.
    let mut r: u32 = 0;
    let mut bb: u32 = 0;
    // Walk from highest limb to lowest.
    let mut u = len;
    while u > 0 {
        u -= 1;
        let wx = x[u * xstride];
        // wp = (p[u] >> 1) | (bb << 30): shift p right by 1, OR in bb from above.
        let wp = (p[u] >> 1) | (bb << 30);
        bb = p[u] & 1;
        // cc = -1, 0, or 1 depending on wp vs wx.
        let cc = wp.wrapping_sub(wx);
        let cc = (cc.wrapping_neg() >> 31) | ((cc >> 31).wrapping_neg());
        // If r != 0 keep it; else replace with cc.
        r |= cc & (r & 1).wrapping_sub(1);
    }

    // If r == -1 (i.e. top bit set), do the conditional subtract.
    let m = tbmask(r);
    let mut cc: u32 = 0;
    for j in 0..len {
        let xw = x[j * xstride];
        let w = xw.wrapping_sub(p[j]).wrapping_sub(cc);
        cc = w >> 31;
        let new_xw = xw ^ (((w & 0x7FFF_FFFF) ^ xw) & m);
        x[j * xstride] = new_xw;
    }
}

// ============================================================
// Sub-task 15.10: zint_bezout + helpers + add_scaled / sub_scaled
// ============================================================

/// Given odd `x`, compute `-1/x mod 2^31` via Newton–Hensel lifting.
/// Port of `mp_ninv31` (ng_zint31.c:363-372).
#[inline]
fn mp_ninv31(x: u32) -> u32 {
    let mut y = 2u32.wrapping_sub(x);
    y = y.wrapping_mul(2u32.wrapping_sub(x.wrapping_mul(y)));
    y = y.wrapping_mul(2u32.wrapping_sub(x.wrapping_mul(y)));
    y = y.wrapping_mul(2u32.wrapping_sub(x.wrapping_mul(y)));
    y = y.wrapping_mul(2u32.wrapping_sub(x.wrapping_mul(y)));
    y.wrapping_neg() & 0x7FFF_FFFF
}

/// Count leading zeros in a 32-bit value (constant-time bit-twiddling).
/// Port of `lzcnt` / `lzcnt_nonzero` (ng_inner.h:666-698).
#[inline]
pub(crate) fn lzcnt(x: u32) -> u32 {
    use crate::keygen::mp31::tbmask;
    let mut x = x;
    let mut m = tbmask((x >> 16).wrapping_sub(1));
    let mut s = m & 16;
    x = (x >> 16) ^ (m & (x ^ (x >> 16)));
    m = tbmask((x >> 8).wrapping_sub(1));
    s |= m & 8;
    x = (x >> 8) ^ (m & (x ^ (x >> 8)));
    m = tbmask((x >> 4).wrapping_sub(1));
    s |= m & 4;
    x = (x >> 4) ^ (m & (x ^ (x >> 4)));
    m = tbmask((x >> 2).wrapping_sub(1));
    s |= m & 2;
    x = (x >> 2) ^ (m & (x ^ (x >> 2)));
    s.wrapping_add(2u32.wrapping_sub(x) & tbmask(x.wrapping_sub(3)))
}

/// Replace `a` with `(a*xa + b*xb) >> 31` and `b` with `(a*ya + b*yb) >> 31`.
/// Dropped bits must be zero (caller's responsibility). If a result is
/// negative it is negated in place. Returns a 2-bit mask: bit 0 = `a` was
/// negated, bit 1 = `b` was negated.
/// Port of `zint_co_reduce` (ng_zint31.c:236-265).
fn zint_co_reduce(a: &mut [u32], b: &mut [u32], xa: i64, xb: i64, ya: i64, yb: i64) -> u32 {
    let len = a.len();
    debug_assert_eq!(b.len(), len);

    let mut cca: i64 = 0;
    let mut ccb: i64 = 0;
    for u in 0..len {
        let wa = a[u];
        let wb = b[u];
        // Cast i64 to u64 is a bit-reinterpret (two's complement), then multiply.
        let za = (wa as u64)
            .wrapping_mul(xa as u64)
            .wrapping_add((wb as u64).wrapping_mul(xb as u64))
            .wrapping_add(cca as u64);
        let zb = (wa as u64)
            .wrapping_mul(ya as u64)
            .wrapping_add((wb as u64).wrapping_mul(yb as u64))
            .wrapping_add(ccb as u64);
        if u > 0 {
            a[u - 1] = (za as u32) & 0x7FFF_FFFF;
            b[u - 1] = (zb as u32) & 0x7FFF_FFFF;
        }
        // Signed right-shift: reinterpret bits as i64, then arithmetic shift.
        cca = (za as i64) >> 31;
        ccb = (zb as i64) >> 31;
    }
    a[len - 1] = (cca as u32) & 0x7FFF_FFFF;
    b[len - 1] = (ccb as u32) & 0x7FFF_FFFF;

    let nega = ((cca as u64) >> 63) as u32;
    let negb = ((ccb as u64) >> 63) as u32;
    zint_negate(a, nega);
    zint_negate(b, negb);
    nega | (negb << 1)
}

/// Finish modular reduction after a co-reduce step.
/// Port of `zint_finish_mod` (ng_zint31.c:278-310).
fn zint_finish_mod(a: &mut [u32], m: &[u32], neg: u32) {
    let len = a.len();
    debug_assert_eq!(m.len(), len);

    // Compare a (assumed non-negative) with m.
    let mut cc: u32 = 0;
    for u in 0..len {
        cc = (a[u].wrapping_sub(m[u]).wrapping_sub(cc)) >> 31;
    }

    // xm = 0x3FFFFFFF if neg = 0, or 0xFFFFFFFF (all-ones) mask if neg = 1.
    // But the C is: `uint32_t xm = -neg >> 1;`  which is 0x7FFFFFFF when neg=1.
    // ym masks out the whole operation when neg=0 and cc=1.
    let xm: u32 = neg.wrapping_neg() >> 1; // 0x7FFFFFFF if neg=1, 0 if neg=0
    let ym: u32 = (neg | (1u32.wrapping_sub(cc))).wrapping_neg(); // all-ones or 0
    cc = neg;
    for u in 0..len {
        let mw = (m[u] ^ xm) & ym;
        let aw = a[u].wrapping_sub(mw).wrapping_sub(cc);
        a[u] = aw & 0x7FFF_FFFF;
        cc = aw >> 31;
    }
}

/// Replace `a` with `(a*xa + b*xb) / 2^31 mod m` and `b` with
/// `(a*ya + b*yb) / 2^31 mod m` (Montgomery-style four-way reduction).
/// `m0i = -1/m[0] mod 2^31`. `m`, `a`, `b` all have the same length.
/// Port of `zint_co_reduce_mod` (ng_zint31.c:317-358).
fn zint_co_reduce_mod(
    a: &mut [u32],
    b: &mut [u32],
    m: &[u32],
    m0i: u32,
    xa: i64,
    xb: i64,
    ya: i64,
    yb: i64,
) {
    let len = a.len();
    debug_assert_eq!(b.len(), len);
    debug_assert_eq!(m.len(), len);

    let mut cca: i64 = 0;
    let mut ccb: i64 = 0;
    // fa, fb: low 31-bit Montgomery reduction factors.
    // (uint32_t)xa truncates int64 to low 32 bits — same as `xa as u32` in Rust.
    let fa: u32 = ((a[0].wrapping_mul(xa as u32))
        .wrapping_add(b[0].wrapping_mul(xb as u32))
        .wrapping_mul(m0i))
        & 0x7FFF_FFFF;
    let fb: u32 = ((a[0].wrapping_mul(ya as u32))
        .wrapping_add(b[0].wrapping_mul(yb as u32))
        .wrapping_mul(m0i))
        & 0x7FFF_FFFF;

    for u in 0..len {
        let wa = a[u];
        let wb = b[u];
        let za = (wa as u64)
            .wrapping_mul(xa as u64)
            .wrapping_add((wb as u64).wrapping_mul(xb as u64))
            .wrapping_add((m[u] as u64).wrapping_mul(fa as u64))
            .wrapping_add(cca as u64);
        let zb = (wa as u64)
            .wrapping_mul(ya as u64)
            .wrapping_add((wb as u64).wrapping_mul(yb as u64))
            .wrapping_add((m[u] as u64).wrapping_mul(fb as u64))
            .wrapping_add(ccb as u64);
        if u > 0 {
            a[u - 1] = (za as u32) & 0x7FFF_FFFF;
            b[u - 1] = (zb as u32) & 0x7FFF_FFFF;
        }
        cca = (za as i64) >> 31;
        ccb = (zb as i64) >> 31;
    }
    a[len - 1] = cca as u32;
    b[len - 1] = ccb as u32;

    zint_finish_mod(a, m, ((cca as u64) >> 63) as u32);
    zint_finish_mod(b, m, ((ccb as u64) >> 63) as u32);
}

/// Compute the Bezout relation `x*u - y*v = 1` for two odd positive big-ints
/// `x` and `y` of `len` limbs. Returns 1 if `gcd(x,y) = 1` (outputs `u`,`v`
/// are valid); returns 0 otherwise.
///
/// `tmp` must have length at least `4*len`.
///
/// Port of `zint_bezout` (ng_zint31.c:374-571).
pub fn zint_bezout(
    u: &mut [u32],
    v: &mut [u32],
    x: &[u32],
    y: &[u32],
    len: usize,
    tmp: &mut [u32],
) -> u32 {
    if len == 0 {
        return 0;
    }
    assert!(tmp.len() >= 4 * len);
    assert!(u.len() >= len);
    assert!(v.len() >= len);
    assert!(x.len() >= len);
    assert!(y.len() >= len);

    // Montgomery reduction constants.
    let x0i = mp_ninv31(x[0]);
    let y0i = mp_ninv31(y[0]);

    // Buffer layout: tmp = [u0 | v0 | a | b], each of length `len`.
    // u1 = u (output), v1 = v (output).
    //
    // Initial values:
    //   a = x,  u0 = 1,    v0 = 0
    //   b = y,  u1 = y,    v1 = x - 1
    {
        // Split tmp into four disjoint slices.
        let (u0v0ab, _) = tmp.split_at_mut(4 * len);
        let (u0v0, ab) = u0v0ab.split_at_mut(2 * len);
        let (u0, v0) = u0v0.split_at_mut(len);
        let (a, b) = ab.split_at_mut(len);

        a.copy_from_slice(&x[..len]);
        b.copy_from_slice(&y[..len]);
        u0[0] = 1;
        for i in 1..len {
            u0[i] = 0;
        }
        for i in 0..len {
            v0[i] = 0;
        }
        u[..len].copy_from_slice(&y[..len]); // u1 = y
        v[..len].copy_from_slice(&x[..len]); // v1 = x
        v[0] -= 1; // v1 = x - 1
    }

    // Main reduction loop: reduce by 31 bits per iteration.
    let mut num: u32 = 62 * (len as u32) + 31;
    while num >= 30 {
        // ---- Extract 63-bit approximations of top words of a and b ----
        let mut c0: u32 = 0xFFFF_FFFF;
        let mut c1: u32 = 0xFFFF_FFFF;
        let mut cp: u32 = 0xFFFF_FFFF;
        let mut a0: u32 = 0;
        let mut a1: u32 = 0;
        let mut b0: u32 = 0;
        let mut b1: u32 = 0;

        {
            let (_, ab) = tmp.split_at(2 * len);
            let a = &ab[..len];
            let b = &ab[len..2 * len];

            let mut j = len;
            while j > 0 {
                j -= 1;
                let aw = a[j];
                let bw = b[j];
                a1 ^= c1 & (a1 ^ aw);
                a0 ^= c0 & (a0 ^ aw);
                b1 ^= c1 & (b1 ^ bw);
                b0 ^= c0 & (b0 ^ bw);
                cp = c0;
                c0 = c1;
                c1 &= (((aw | bw).wrapping_add(0x7FFF_FFFF)) >> 31).wrapping_sub(1);
            }
        }

        // Align top 32 bits. s is a shift count (uses lzcnt on nonzero input).
        let s = lzcnt(a1 | b1 | ((cp & c0) >> 1));
        let mut ha = (a1 << s) | (a0 >> (31 - s));
        let mut hb = (b1 << s) | (b0 >> (31 - s));

        // If j <= 62, use non-aligned bits.
        ha ^= cp & (ha ^ a1);
        hb ^= cp & (hb ^ b1);

        // If j <= 31, clear upper bits entirely.
        ha &= !c0;
        hb &= !c0;

        // Assemble 63-bit approximations.
        let xa_approx: u64 = ((ha as u64) << 31) | (tmp[2 * len] as u64);
        let xb_approx: u64 = ((hb as u64) << 31) | (tmp[3 * len] as u64);

        // ---- 31-step inner binary GCD ----
        let mut xa = xa_approx;
        let mut xb = xb_approx;
        let mut fg0: u64 = 1;
        let mut fg1: u64 = 1u64 << 32;
        for _ in 0..31 {
            let a_odd: u64 = (xa & 1).wrapping_neg();
            let dx = xa.wrapping_sub(xb);
            let dx_sign: u64 = (dx as i64 >> 63) as u64;
            let swap = a_odd & dx_sign;
            let t1 = swap & (xa ^ xb);
            xa ^= t1;
            xb ^= t1;
            let t2 = swap & (fg0 ^ fg1);
            fg0 ^= t2;
            fg1 ^= t2;
            xa = xa.wrapping_sub(a_odd & xb);
            fg0 = fg0.wrapping_sub(a_odd & fg1);
            xa >>= 1;
            fg1 <<= 1;
        }

        // ---- Split update factors ----
        fg0 = fg0.wrapping_add(0x7FFF_FFFF_7FFF_FFFFu64);
        fg1 = fg1.wrapping_add(0x7FFF_FFFF_7FFF_FFFFu64);
        let f0: i64 = (fg0 & 0xFFFF_FFFF) as i64 - 0x7FFF_FFFF;
        let g0: i64 = (fg0 >> 32) as i64 - 0x7FFF_FFFF;
        let f1: i64 = (fg1 & 0xFFFF_FFFF) as i64 - 0x7FFF_FFFF;
        let g1: i64 = (fg1 >> 32) as i64 - 0x7FFF_FFFF;

        // ---- Apply update factors to a, b (co_reduce) ----
        let negab = {
            let (_, ab) = tmp.split_at_mut(2 * len);
            let (a, b) = ab.split_at_mut(len);
            zint_co_reduce(a, b, f0, g0, f1, g1)
        };

        // Adjust signs of factors for co_reduce_mod.
        let mut f0m = f0;
        let mut g0m = g0;
        let mut f1m = f1;
        let mut g1m = g1;
        f0m -= (f0m + f0m) & -((negab & 1) as i64);
        g0m -= (g0m + g0m) & -((negab & 1) as i64);
        f1m -= (f1m + f1m) & -(((negab >> 1) & 1) as i64);
        g1m -= (g1m + g1m) & -(((negab >> 1) & 1) as i64);

        // ---- Apply to u0/u1 mod y ----
        {
            let (u0, rest) = tmp.split_at_mut(len);
            let (v0, _) = rest.split_at_mut(len);
            // u0 and u[..len] are disjoint (u0 is in tmp, u is external).
            // y is a shared ref — safe.
            zint_co_reduce_mod(u0, &mut u[..len], &y[..len], y0i, f0m, g0m, f1m, g1m);
            // v0 and v[..len] are disjoint.
            zint_co_reduce_mod(v0, &mut v[..len], &x[..len], x0i, f0m, g0m, f1m, g1m);
        }

        num -= 31;
    }

    // b should contain gcd(x,y).  Accept only if gcd == 1 and both inputs odd.
    let mut r: u32 = tmp[3 * len] ^ 1; // b[0] ^ 1
    for j in 1..len {
        r |= tmp[3 * len + j]; // b[j]
    }
    r |= (x[0] & y[0] & 1) ^ 1;
    1u32.wrapping_sub((r | r.wrapping_neg()) >> 31)
}

/// Compute `x += k * (2^(31*sch + scl)) * y`.
///
/// `x` has `xlen` limbs, `y` has `ylen` limbs, both with stride `stride`.
/// Port of `zint_add_scaled_mul_small` (ng_zint31.c:574-616).
pub fn zint_add_scaled_mul_small(
    x: &mut [u32],
    xlen: usize,
    y: &[u32],
    ylen: usize,
    stride: usize,
    k: i32,
    sch: u32,
    scl: u32,
) {
    if ylen == 0 {
        return;
    }

    // Sign extension word for y.
    let ysign: u32 = (y[stride * (ylen - 1)] >> 30).wrapping_neg() >> 1;
    let mut tw: u32 = 0;
    let mut cc: i32 = 0;
    let mut yi: usize = 0; // index into y (pointer arithmetic emulation)
    let mut yrem = ylen;

    for u in (sch as usize)..xlen {
        // Get next word of (2^scl)*y.
        let wy: u32 = if yrem > 0 {
            let w = y[yi];
            yi += stride;
            yrem -= 1;
            w
        } else {
            ysign
        };
        let wys = ((wy << scl) & 0x7FFF_FFFF) | tw;
        tw = wy >> (31 - scl);

        let z: u64 = (wys as i64)
            .wrapping_mul(k as i64)
            .wrapping_add(x[u * stride] as i64)
            .wrapping_add(cc as i64) as u64;
        x[u * stride] = (z as u32) & 0x7FFF_FFFF;

        // Carry: signed right-shift of z by 31 bits.
        let ccu = (z >> 31) as u32;
        cc = ccu as i32;
    }
}

/// Compute `x -= y * 2^(31*sch + scl)`.
///
/// `x` has `xlen` limbs, `y` has `ylen` limbs, both with stride `stride`.
/// Port of `zint_sub_scaled` (ng_zint31.c:619-652).
pub fn zint_sub_scaled(
    x: &mut [u32],
    xlen: usize,
    y: &[u32],
    ylen: usize,
    stride: usize,
    sch: u32,
    scl: u32,
) {
    if ylen == 0 {
        return;
    }

    let ysign: u32 = (y[stride * (ylen - 1)] >> 30).wrapping_neg() >> 1;
    let mut tw: u32 = 0;
    let mut cc: u32 = 0;
    let mut yi: usize = 0;
    let mut yrem = ylen;

    for u in (sch as usize)..xlen {
        let wy: u32 = if yrem > 0 {
            let w = y[yi];
            yi += stride;
            yrem -= 1;
            w
        } else {
            ysign
        };
        let wys = ((wy << scl) & 0x7FFF_FFFF) | tw;
        tw = wy >> (31 - scl);

        let w = x[u * stride].wrapping_sub(wys).wrapping_sub(cc);
        x[u * stride] = w & 0x7FFF_FFFF;
        cc = w >> 31;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::keygen::primes::PRIMES;

    #[test]
    fn mul_small_by_zero() {
        let mut m = [5u32, 10, 15];
        let cc = zint_mul_small(&mut m, 0);
        assert_eq!(m, [0, 0, 0]);
        assert_eq!(cc, 0);
    }

    #[test]
    fn mul_small_by_one_is_identity() {
        let mut m = [5u32, 10, 15];
        let cc = zint_mul_small(&mut m, 1);
        assert_eq!(m, [5, 10, 15]);
        assert_eq!(cc, 0);
    }

    #[test]
    fn mul_small_produces_carry() {
        // m = [0x7FFFFFFF, 0, 0] (value = 2^31 - 1)
        // multiply by 2 -> value = 2^32 - 2 = 0x7FFFFFFE * 2^1 with carry.
        // Limb 0: (2^31 - 1) * 2 = 2^32 - 2 = 0xFFFFFFFE. Low 31 bits = 0x7FFFFFFE. Carry = 1.
        // Limb 1: 0 * 2 + 1 = 1. Low 31 = 1. Carry = 0.
        let mut m = [0x7FFFFFFFu32, 0, 0];
        let cc = zint_mul_small(&mut m, 2);
        assert_eq!(m, [0x7FFFFFFE, 1, 0]);
        assert_eq!(cc, 0);
    }

    #[test]
    fn mod_small_unsigned_basic() {
        // Encode the value 100 as a 1-limb big-int; reduce mod PRIMES[0].p.
        // Result should be 100 (since 100 < p).
        let d = [100u32];
        let p = PRIMES[0].p;
        let p0i = PRIMES[0].p0i;
        let r2 = PRIMES[0].r2;
        let result = zint_mod_small_unsigned(&d, 1, 1, p, p0i, r2);
        assert_eq!(result, 100);
    }

    #[test]
    fn mod_small_unsigned_reduces_multi_limb() {
        // Build a 2-limb big-int representing p (the prime itself). Its reduction
        // mod p must be 0.
        // p fits in 31 bits so as a 1-limb big-int its value is p (already less than 2^31).
        // Use a 2-limb representation: [p, 0] = p. Reduce mod p = 0.
        let p = PRIMES[0].p;
        let p0i = PRIMES[0].p0i;
        let r2 = PRIMES[0].r2;
        let d = [p, 0u32];
        let result = zint_mod_small_unsigned(&d, 2, 1, p, p0i, r2);
        assert_eq!(result, 0);
    }

    #[test]
    fn add_mul_small_basic() {
        // x = [0, 0], y = [7], s = 3. After call: x = [21, 0, 0] (one extra limb for carry).
        let mut x = [0u32, 0, 0];
        let y = [7u32];
        zint_add_mul_small(&mut x, 1, 1, &y, 3);
        assert_eq!(x, [21, 0, 0]);
    }

    #[test]
    fn add_mul_small_accumulates() {
        // x = [5, 0, 0], y = [7], s = 3. After: x[0] = 5 + 21 = 26. Carry = 0.
        let mut x = [5u32, 0, 0];
        let y = [7u32];
        zint_add_mul_small(&mut x, 1, 1, &y, 3);
        assert_eq!(x, [26, 0, 0]);
    }

    #[test]
    fn norm_zero_does_nothing_when_x_small() {
        // x = 1, p = 1000 (use a small-ish p). Since 1 < 500, no change.
        // But zint_norm_zero requires p to be real, so use the full p = PRIMES[0].p.
        // x = [1] (1 word). (p-1)/2 >> 0, so 1 < (p-1)/2, no change expected.
        let p_limb = PRIMES[0].p;
        let p = [p_limb];
        let mut x = [1u32];
        zint_norm_zero(&mut x, 1, 1, &p);
        assert_eq!(x, [1]);
    }

    #[test]
    fn norm_zero_subtracts_when_x_large() {
        // x = p-1, should become (p-1) - p = -1 mod 2^31 = 0x7FFFFFFF (top bit per two's complement).
        // Actually: after subtraction, x = -1 encoded as 0x7FFFFFFF (low 31 bits) with sign
        // bit set in bit 30 of that limb. Verify the low 31 bits.
        let p_limb = PRIMES[0].p;
        let p = [p_limb];
        let mut x = [p_limb - 1];
        zint_norm_zero(&mut x, 1, 1, &p);
        // (p-1) - p = -1. In 31-bit two's complement, -1 is 0x7FFFFFFF.
        assert_eq!(x, [0x7FFFFFFF]);
    }

    // --- zint_negate tests ---

    #[test]
    fn negate_identity_when_ctl_zero() {
        let mut a = [1u32, 2, 3, 4];
        zint_negate(&mut a, 0);
        assert_eq!(a, [1, 2, 3, 4]);
    }

    #[test]
    fn negate_flips_when_ctl_one() {
        // Negate [5, 0, 0, 0].
        // m = 0x7FFFFFFF. cc = 1.
        // limb 0: (5 ^ 0x7FFFFFFF) + 1 = 0x7FFFFFFA + 1 = 0x7FFFFFFB. low31 = 0x7FFFFFFB. cc = 0.
        // limbs 1-3: (0 ^ 0x7FFFFFFF) + 0 = 0x7FFFFFFF. low31 = 0x7FFFFFFF. cc = 0.
        let mut a = [5u32, 0, 0, 0];
        zint_negate(&mut a, 1);
        assert_eq!(a, [0x7FFFFFFB, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF]);
    }

    #[test]
    fn negate_twice_is_identity() {
        let mut a = [0x12345u32, 0x6789A, 0xBCDE, 0];
        let orig = a;
        zint_negate(&mut a, 1);
        zint_negate(&mut a, 1);
        assert_eq!(a, orig);
    }

    // --- zint_rebuild_crt tests ---

    #[test]
    fn rebuild_crt_single_limb_noop() {
        // xlen=1: the for loop (1..1) never runs, xx is untouched.
        let mut xx = [42u32; 8];
        let mut tmp = [0u32; 1];
        zint_rebuild_crt(&mut xx, 1, 4, 2, false, &mut tmp);
        assert_eq!(xx, [42u32; 8]);
    }

    #[test]
    fn rebuild_crt_two_limbs_small_value() {
        // RNS representation of 100 across PRIMES[0] and PRIMES[1].
        // num_sets=1, n=1: xx = [100 % p0, 100 % p1]. After rebuild: xx = [100, 0].
        let mut xx = [100u32 % PRIMES[0].p, 100u32 % PRIMES[1].p];
        let mut tmp = [0u32; 2];
        zint_rebuild_crt(&mut xx, 2, 1, 1, false, &mut tmp);
        assert_eq!(xx, [100, 0]);
    }

    // --- zint_bezout tests ---

    #[test]
    fn bezout_of_coprime_small_values() {
        // x = 3, y = 5. gcd = 1.
        // Bezout: 3*u - 5*v = 1 with 0 <= u < 5 and 0 <= v < 3.
        // u = 2, v = 1: 3*2 - 5*1 = 1. ✓
        let x = [3u32];
        let y = [5u32];
        let mut u = [0u32];
        let mut v = [0u32];
        let mut tmp = vec![0u32; 4];
        let r = zint_bezout(&mut u, &mut v, &x, &y, 1, &mut tmp);
        assert_eq!(r, 1);
        assert_eq!(u[0], 2);
        assert_eq!(v[0], 1);
    }

    #[test]
    fn bezout_of_non_coprime_returns_zero() {
        // x = 9, y = 15. gcd = 3. Return 0.
        let x = [9u32];
        let y = [15u32];
        let mut u = [0u32];
        let mut v = [0u32];
        let mut tmp = vec![0u32; 4];
        let r = zint_bezout(&mut u, &mut v, &x, &y, 1, &mut tmp);
        assert_eq!(r, 0);
    }

    #[test]
    fn bezout_rejects_even_inputs() {
        // x = 4, y = 5. x is even. Should return 0.
        let x = [4u32];
        let y = [5u32];
        let mut u = [0u32];
        let mut v = [0u32];
        let mut tmp = vec![0u32; 4];
        let r = zint_bezout(&mut u, &mut v, &x, &y, 1, &mut tmp);
        assert_eq!(r, 0);
    }

    #[test]
    fn add_scaled_mul_small_zero_y_is_noop() {
        let mut x = [5u32, 10, 15, 20];
        let y: [u32; 0] = [];
        zint_add_scaled_mul_small(&mut x, 4, &y, 0, 1, 3, 0, 0);
        assert_eq!(x, [5, 10, 15, 20]);
    }

    #[test]
    fn sub_scaled_zero_y_is_noop() {
        let mut x = [5u32, 10, 15, 20];
        let y: [u32; 0] = [];
        zint_sub_scaled(&mut x, 4, &y, 0, 1, 0, 0);
        assert_eq!(x, [5, 10, 15, 20]);
    }
}
