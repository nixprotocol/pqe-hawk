//! `GF(2)[X]/(X^n+1)` polynomial arithmetic for HAWK-512 signing.
//!
//! Port of the `bp_*` helpers from c-reference/hawk-512/hawk_sign.c:25-295.
//! These implement carryless (GF(2)) polynomial multiplication via Karatsuba
//! recursion, used by `basis_m2_mul` to compute (t0, t1) = B*h mod 2.
//!
//! Scratch-buffer requirements (see C macro expansions):
//!   bp_muladd_64  : 0 bytes
//!   bp_muladd_128 : 32 bytes  (2 * 128/8)
//!   bp_muladd_256 : 96 bytes  (2 * 256/8 + scratch_for_muladd_128 = 64 + 32)
//!   bp_muladd_512 : 224 bytes (2 * 512/8 + scratch_for_muladd_256 = 128 + 96)
//!   bp_mulmod_512 : 160 bytes (512/8 + scratch_for_muladd_256 = 64 + 96)
//!   basis_m2_mul  : 224 bytes (64 `[w1]` + 160 `[w2 for bp_mulmod_512]`)

// ---------------------------------------------------------------------------
// XOR helpers (port of bp_xor_N macros from hawk_sign.c:71-149)
// ---------------------------------------------------------------------------

/// Port of `bp_xor_N` for completeness of the C-macro surface. HAWK-512's
/// actual code path XORs inline, so this free function is currently unused —
/// kept for parity with hawk_sign.c's bp_xor_N family and for future callers.
#[allow(dead_code)]
#[inline]
pub(crate) fn bp_xor(d: &mut [u8], a: &[u8], b: &[u8]) {
    for i in 0..d.len() {
        d[i] = a[i] ^ b[i];
    }
}

// ---------------------------------------------------------------------------
// Carryless 32×32 → 64-bit multiply (port of bp_mul_32, hawk_sign.c:38-69)
//
// Classic "4-bit interleaved" technique: split x,y by nibble-lane (bits
// 0,4,8,... / 1,5,9,... / 2,6,10,... / 3,7,11,...).  Each lane's carry is
// absorbed by the gap to the next lane, so ordinary 64-bit integer multiplies
// give the GF(2) lane product without cross-carry.  Masking then keeps only
// the desired bits.
// ---------------------------------------------------------------------------

pub(crate) fn bp_mul_32(x: u32, y: u32) -> u64 {
    let x0 = x & 0x1111_1111u32;
    let x1 = x & 0x2222_2222u32;
    let x2 = x & 0x4444_4444u32;
    let x3 = x & 0x8888_8888u32;

    let y0 = y & 0x1111_1111u32;
    let y1 = y & 0x2222_2222u32;
    let y2 = y & 0x4444_4444u32;
    let y3 = y & 0x8888_8888u32;

    let m = |a: u32, b: u32| -> u64 { (a as u64).wrapping_mul(b as u64) };

    let z0 = m(x0, y0) ^ m(x1, y3) ^ m(x2, y2) ^ m(x3, y1);
    let z1 = m(x0, y1) ^ m(x1, y0) ^ m(x2, y3) ^ m(x3, y2);
    let z2 = m(x0, y2) ^ m(x1, y1) ^ m(x2, y0) ^ m(x3, y3);
    let z3 = m(x0, y3) ^ m(x1, y2) ^ m(x2, y1) ^ m(x3, y0);

    let z0 = z0 & 0x1111_1111_1111_1111u64;
    let z1 = z1 & 0x2222_2222_2222_2222u64;
    let z2 = z2 & 0x4444_4444_4444_4444u64;
    let z3 = z3 & 0x8888_8888_8888_8888u64;

    z0 | z1 | z2 | z3
}

// ---------------------------------------------------------------------------
// Little-endian byte-array <-> integer helpers
// (mirrors dec32le / dec64le / enc64le in hawk_inner.h)
// ---------------------------------------------------------------------------

#[inline]
fn dec32le(src: &[u8]) -> u32 {
    u32::from_le_bytes(src[..4].try_into().unwrap())
}

#[inline]
fn dec64le(src: &[u8]) -> u64 {
    u64::from_le_bytes(src[..8].try_into().unwrap())
}

#[inline]
fn enc64le(dst: &mut [u8], x: u64) {
    dst[..8].copy_from_slice(&x.to_le_bytes());
}

// ---------------------------------------------------------------------------
// bp_muladd_64(d, a, b, _tmp)
//
// d  : 16-byte accumulator (128-bit running product)
// a,b: 8-byte inputs       (64-bit polynomials)
// Port of hawk_sign.c:151-168.
// ---------------------------------------------------------------------------

pub(crate) fn bp_muladd_64(d: &mut [u8], a: &[u8], b: &[u8], _tmp: &mut [u8]) {
    debug_assert!(d.len() >= 16);
    debug_assert!(a.len() >= 8);
    debug_assert!(b.len() >= 8);

    let a0 = dec32le(a);
    let a1 = dec32le(&a[4..]);
    let b0 = dec32le(b);
    let b1 = dec32le(&b[4..]);

    let c0 = bp_mul_32(a0, b0);
    let c1 = bp_mul_32(a1, b1);
    let c2 = bp_mul_32(a0 ^ a1, b0 ^ b1) ^ c0 ^ c1;

    let lo = dec64le(d) ^ c0 ^ (c2 << 32);
    enc64le(&mut d[..8], lo);
    let hi = dec64le(&d[8..]) ^ c1 ^ (c2 >> 32);
    enc64le(&mut d[8..16], hi);
}

// ---------------------------------------------------------------------------
// Karatsuba muladd helpers — generic implementation of MKBP_MULADD.
//
// bp_muladd_N(d, a, b, tmp):
//   d  : 2*(N/8)-byte accumulator
//   a,b: N/8-byte inputs
//   tmp: scratch (see sizing comment at top of file)
//
// For each N we implement the exact macro expansion from hawk_sign.c.
// ---------------------------------------------------------------------------

// Helper: XOR n bytes of src into dst, both starting at offset 0.
#[inline]
fn xor_bytes(dst: &mut [u8], src: &[u8], n: usize) {
    for i in 0..n {
        dst[i] ^= src[i];
    }
}

/// bp_muladd_128: d[32] += a[16] × b[16], tmp ≥ 32 bytes.
pub(crate) fn bp_muladd_128(d: &mut [u8], a: &[u8], b: &[u8], tmp: &mut [u8]) {
    // n=128, hn=64 → N8=16, HN8=8
    const N8: usize = 16;
    const HN8: usize = 8;
    debug_assert!(d.len() >= 2 * N8);
    debug_assert!(a.len() >= N8);
    debug_assert!(b.len() >= N8);
    debug_assert!(tmp.len() >= 2 * N8); // t1[16] + t2[16]; t3 unused by muladd_64

    // t2[0..8]  = a0 + a1
    // t2[8..16] = b0 + b1
    for i in 0..HN8 {
        tmp[N8 + i] = a[i] ^ a[HN8 + i];
        tmp[N8 + HN8 + i] = b[i] ^ b[HN8 + i];
    }
    // t1[0..16] = d[0..16] ^ d[16..32]
    for i in 0..N8 {
        tmp[i] = d[i] ^ d[N8 + i];
    }
    // t1 += (a0+a1)*(b0+b1)
    // bp_muladd_64(t1, t2[0..8], t2[8..16], t3)
    // We need disjoint refs: t1=[0..16], t2=[16..32]; pass empty slice as t3.
    {
        // Copy t2 halves to avoid aliasing with t1 during the call
        let ta = [
            tmp[N8],
            tmp[N8 + 1],
            tmp[N8 + 2],
            tmp[N8 + 3],
            tmp[N8 + 4],
            tmp[N8 + 5],
            tmp[N8 + 6],
            tmp[N8 + 7],
        ];
        let tb = [
            tmp[N8 + HN8],
            tmp[N8 + HN8 + 1],
            tmp[N8 + HN8 + 2],
            tmp[N8 + HN8 + 3],
            tmp[N8 + HN8 + 4],
            tmp[N8 + HN8 + 5],
            tmp[N8 + HN8 + 6],
            tmp[N8 + HN8 + 7],
        ];
        let (t1, _rest) = tmp.split_at_mut(N8);
        bp_muladd_64(t1, &ta, &tb, &mut []);
    }
    // d[0..16] += a[0..8]*b[0..8]
    {
        let (d_lo, _) = d.split_at_mut(N8);
        let a0 = &a[..HN8];
        let b0 = &b[..HN8];
        bp_muladd_64(d_lo, a0, b0, &mut []);
    }
    // d[16..32] += a[8..16]*b[8..16]
    {
        let (_, d_rest) = d.split_at_mut(N8);
        let a1 = &a[HN8..N8];
        let b1 = &b[HN8..N8];
        bp_muladd_64(d_rest, a1, b1, &mut []);
    }
    // t1 ^= d[0..16] ^ d[16..32]
    for i in 0..N8 {
        tmp[i] ^= d[i];
    }
    for i in 0..N8 {
        tmp[i] ^= d[N8 + i];
    }
    // d[8..24] ^= t1
    for i in 0..N8 {
        d[HN8 + i] ^= tmp[i];
    }
}

/// bp_muladd_256: d[64] += a[32] × b[32], tmp ≥ 96 bytes.
pub(crate) fn bp_muladd_256(d: &mut [u8], a: &[u8], b: &[u8], tmp: &mut [u8]) {
    // n=256, hn=128 → N8=32, HN8=16
    const N8: usize = 32;
    const HN8: usize = 16;
    debug_assert!(d.len() >= 2 * N8);
    debug_assert!(a.len() >= N8);
    debug_assert!(b.len() >= N8);
    debug_assert!(tmp.len() >= 2 * N8 + 32); // t1[32]+t2[32]+t3≥32 for muladd_128

    // t2[0..16]  = a0+a1
    // t2[16..32] = b0+b1
    for i in 0..HN8 {
        tmp[N8 + i] = a[i] ^ a[HN8 + i];
        tmp[N8 + HN8 + i] = b[i] ^ b[HN8 + i];
    }
    // t1[0..32] = d[0..32] ^ d[32..64]
    for i in 0..N8 {
        tmp[i] = d[i] ^ d[N8 + i];
    }
    // t1 += (a0+a1)*(b0+b1)
    {
        let t2a: [u8; 16] = tmp[N8..N8 + HN8].try_into().unwrap();
        let t2b: [u8; 16] = tmp[N8 + HN8..N8 + N8].try_into().unwrap();
        let (t1, rest) = tmp.split_at_mut(N8);
        let (_, t3) = rest.split_at_mut(N8);
        bp_muladd_128(t1, &t2a, &t2b, t3);
    }
    // d[0..32] += a[0..16]*b[0..16]
    {
        let a0: [u8; 16] = a[..HN8].try_into().unwrap();
        let b0: [u8; 16] = b[..HN8].try_into().unwrap();
        let (d_lo, _) = d.split_at_mut(N8);
        let (_, t3) = tmp.split_at_mut(N8);
        bp_muladd_128(d_lo, &a0, &b0, t3);
    }
    // d[32..64] += a[16..32]*b[16..32]
    {
        let a1: [u8; 16] = a[HN8..N8].try_into().unwrap();
        let b1: [u8; 16] = b[HN8..N8].try_into().unwrap();
        let (_, d_rest) = d.split_at_mut(N8);
        let (_, t3) = tmp.split_at_mut(N8);
        bp_muladd_128(d_rest, &a1, &b1, t3);
    }
    // t1 ^= d[0..32] ^ d[32..64]
    for i in 0..N8 {
        tmp[i] ^= d[i];
    }
    for i in 0..N8 {
        tmp[i] ^= d[N8 + i];
    }
    // d[16..48] ^= t1
    for i in 0..N8 {
        d[HN8 + i] ^= tmp[i];
    }
}

/// bp_muladd_512: d[128] += a[64] × b[64], tmp ≥ 224 bytes.
///
/// Part of the MKBP_MULADD macro expansion in hawk_sign.c:197. HAWK-512's
/// bp_mulmod_512 uses bp_muladd_256 internally, not _512, so this is kept
/// for completeness with the C macro surface.
#[allow(dead_code)]
pub(crate) fn bp_muladd_512(d: &mut [u8], a: &[u8], b: &[u8], tmp: &mut [u8]) {
    // n=512, hn=256 → N8=64, HN8=32
    const N8: usize = 64;
    const HN8: usize = 32;
    debug_assert!(d.len() >= 2 * N8);
    debug_assert!(a.len() >= N8);
    debug_assert!(b.len() >= N8);
    debug_assert!(tmp.len() >= 2 * N8 + 96); // t1[64]+t2[64]+t3≥96 for muladd_256

    // t2[0..32]  = a0+a1
    // t2[32..64] = b0+b1
    for i in 0..HN8 {
        tmp[N8 + i] = a[i] ^ a[HN8 + i];
        tmp[N8 + HN8 + i] = b[i] ^ b[HN8 + i];
    }
    // t1[0..64] = d[0..64] ^ d[64..128]
    for i in 0..N8 {
        tmp[i] = d[i] ^ d[N8 + i];
    }
    // t1 += (a0+a1)*(b0+b1)
    {
        let mut t2a = [0u8; 32];
        let mut t2b = [0u8; 32];
        t2a.copy_from_slice(&tmp[N8..N8 + HN8]);
        t2b.copy_from_slice(&tmp[N8 + HN8..N8 + N8]);
        let (t1, rest) = tmp.split_at_mut(N8);
        let (_, t3) = rest.split_at_mut(N8);
        bp_muladd_256(t1, &t2a, &t2b, t3);
    }
    // d[0..64] += a[0..32]*b[0..32]
    {
        let mut a0 = [0u8; 32];
        let mut b0 = [0u8; 32];
        a0.copy_from_slice(&a[..HN8]);
        b0.copy_from_slice(&b[..HN8]);
        let (d_lo, _) = d.split_at_mut(N8);
        let (_, t3) = tmp.split_at_mut(N8);
        bp_muladd_256(d_lo, &a0, &b0, t3);
    }
    // d[64..128] += a[32..64]*b[32..64]
    {
        let mut a1 = [0u8; 32];
        let mut b1 = [0u8; 32];
        a1.copy_from_slice(&a[HN8..N8]);
        b1.copy_from_slice(&b[HN8..N8]);
        let (_, d_rest) = d.split_at_mut(N8);
        let (_, t3) = tmp.split_at_mut(N8);
        bp_muladd_256(d_rest, &a1, &b1, t3);
    }
    // t1 ^= d[0..64] ^ d[64..128]
    for i in 0..N8 {
        tmp[i] ^= d[i];
    }
    for i in 0..N8 {
        tmp[i] ^= d[N8 + i];
    }
    // d[32..96] ^= t1
    for i in 0..N8 {
        d[HN8 + i] ^= tmp[i];
    }
}

// ---------------------------------------------------------------------------
// bp_mulmod_512: d[64] = a[64]*b[64] mod X^512+1
//
// Port of MKBP_MULMOD(512, 256) (hawk_sign.c:199-225).
// tmp must be ≥ 160 bytes: t1[64] + t2[96] for bp_muladd_256.
//
// Algorithm (n=512, hn=256, n/8=64, hn/8=32):
//   a = a0 + a1*X^256,  b = b0 + b1*X^256
//   t1  = (a0+a1)*(b0+b1)          [512-bit product, stored in t1[64]]
//   d   = a0*b0 + a1*b1            [mod X^512+1: both halves fold back]
//   t1 ^= d                        → t1 = a0*b1 + a1*b0  (cross term)
//   Fold by X^256 rotation mod X^512+1 (X^512 ≡ 1):
//   d[0..32]  ^= t1[32..64]
//   d[32..64] ^= t1[0..32]
// ---------------------------------------------------------------------------

/// Multiplication in `GF(2)[X]/(X^512+1)`.
/// `d` = `a`×`b` mod `X^512+1`.
/// All three are 64-byte arrays (512 bits); `tmp` ≥ 160 bytes scratch.
pub fn bp_mulmod_512(d: &mut [u8], a: &[u8], b: &[u8], tmp: &mut [u8]) {
    // n=512, hn=256 → N8=64, HN8=32
    const N8: usize = 64;
    const HN8: usize = 32;

    debug_assert!(d.len() >= N8);
    debug_assert!(a.len() >= N8);
    debug_assert!(b.len() >= N8);
    debug_assert!(tmp.len() >= N8 + 96); // t1[64] + t2[96]

    // Step: d[0..32] = a0+a1, d[32..64] = b0+b1  (use d as scratch)
    for i in 0..HN8 {
        d[i] = a[i] ^ a[HN8 + i];
        d[HN8 + i] = b[i] ^ b[HN8 + i];
    }
    // t1[0..64] = 0
    for i in 0..N8 {
        tmp[i] = 0;
    }
    // t1 += (a0+a1)*(b0+b1)  via bp_muladd_256(t1, d[0..32], d[32..64], t2)
    {
        let da: [u8; 32] = d[..HN8].try_into().unwrap();
        let db: [u8; 32] = d[HN8..N8].try_into().unwrap();
        let (t1, t2) = tmp.split_at_mut(N8);
        bp_muladd_256(t1, &da, &db, t2);
    }
    // d = 0; then d += a0*b0 + a1*b1
    for i in 0..N8 {
        d[i] = 0;
    }
    {
        let a0: [u8; 32] = a[..HN8].try_into().unwrap();
        let b0: [u8; 32] = b[..HN8].try_into().unwrap();
        let a1: [u8; 32] = a[HN8..N8].try_into().unwrap();
        let b1: [u8; 32] = b[HN8..N8].try_into().unwrap();
        let (_, t2) = tmp.split_at_mut(N8);
        bp_muladd_256(d, &a0, &b0, t2);
        bp_muladd_256(d, &a1, &b1, t2);
    }
    // t1 ^= d  → t1 = a0*b1 + a1*b0
    xor_bytes(tmp, d, N8);
    // d ^= rotate_256(t1): d[0..32] ^= t1[32..64], d[32..64] ^= t1[0..32]
    // Note: the cross-term t1 has 512 bits but multiplying two 256-bit halves
    // gives at most a 511-bit result, so t1 fits in 64 bytes with possible
    // overflow into the full 64 bytes.  The rotation XORs:
    //   d[lo half] ^= t1[hi half]   (high bits of cross-term fold by X^512≡1)
    //   d[hi half] ^= t1[lo half]   (low bits of cross-term shift up)
    for i in 0..HN8 {
        d[i] ^= tmp[HN8 + i];
        d[HN8 + i] ^= tmp[i];
    }
}

// ---------------------------------------------------------------------------
// basis_m2_mul for HAWK-512 (logn=9)
//
// (t0, t1) = B * h  mod 2, where B = [[f F],[g G]].
// Port of basis_m2_mul case 9 (hawk_sign.c:279-286).
//
// All byte arrays are n/8 = 64 bytes.
// tmp must be ≥ 224 bytes: w1[64] + w2[160] (for bp_mulmod_512).
// ---------------------------------------------------------------------------

/// Compute `(t0, t1) = B × h mod 2` for HAWK-512.
///
/// - `t0 = f2·h0 ⊕ f_cap2·h1`
/// - `t1 = g2·h0 ⊕ g_cap2·h1`
///
/// All slice arguments must be exactly 64 bytes.  `tmp` must be ≥ 224 bytes.
pub fn basis_m2_mul(
    t0: &mut [u8],
    t1: &mut [u8],
    h0: &[u8],
    h1: &[u8],
    f2: &[u8],
    g2: &[u8],
    f_cap2: &[u8],
    g_cap2: &[u8],
    tmp: &mut [u8],
) {
    assert_eq!(t0.len(), 64);
    assert_eq!(t1.len(), 64);
    assert_eq!(h0.len(), 64);
    assert_eq!(h1.len(), 64);
    assert_eq!(f2.len(), 64);
    assert_eq!(g2.len(), 64);
    assert_eq!(f_cap2.len(), 64);
    assert_eq!(g_cap2.len(), 64);
    assert!(tmp.len() >= 224);

    // w1 = tmp[0..64], w2 = tmp[64..224]
    let (w1_slice, w2_slice) = tmp.split_at_mut(64);

    // t0 = h0 * f2
    bp_mulmod_512(t0, h0, f2, w2_slice);
    // w1 = h1 * F2
    bp_mulmod_512(w1_slice, h1, f_cap2, w2_slice);
    // t0 ^= w1
    for i in 0..64 {
        t0[i] ^= w1_slice[i];
    }

    // t1 = h0 * g2
    bp_mulmod_512(t1, h0, g2, w2_slice);
    // w1 = h1 * G2
    bp_mulmod_512(w1_slice, h1, g_cap2, w2_slice);
    // t1 ^= w1
    for i in 0..64 {
        t1[i] ^= w1_slice[i];
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bp_xor_basic() {
        let a = [0xAAu8, 0x55, 0xFF, 0x00];
        let b = [0xFFu8, 0xFF, 0x0F, 0xF0];
        let mut d = [0u8; 4];
        bp_xor(&mut d, &a, &b);
        assert_eq!(d, [0x55, 0xAA, 0xF0, 0xF0]);
    }

    #[test]
    fn bp_mul_32_identity() {
        // 1 * x = x in GF(2)[X]. bp_mul_32(1, x) = x (with high bits 0).
        let r = bp_mul_32(1, 0xDEAD_BEEF);
        assert_eq!(r, 0xDEAD_BEEF);
    }

    #[test]
    fn bp_mul_32_commutes() {
        let a: u32 = 0x1234_5678;
        let b: u32 = 0x9ABC_DEF0;
        assert_eq!(bp_mul_32(a, b), bp_mul_32(b, a));
    }

    #[test]
    fn bp_mul_32_known_value() {
        // (X+1)^2 = X^2+1 in GF(2)[X].  (X+1) = 0b11 = 3.
        // So bp_mul_32(3, 3) should equal 5 (X^2 + 1 = 0b101).
        assert_eq!(bp_mul_32(3, 3), 5);
    }

    #[test]
    fn bp_mul_32_degree_product() {
        // X^3 * X^4 = X^7. Bits: 0b1000 * 0b10000 = 0b10000000
        assert_eq!(bp_mul_32(1 << 3, 1 << 4), 1u64 << 7);
    }

    #[test]
    fn bp_muladd_64_basic() {
        // 1 * 1 = 1.  After accumulate from zero: d[0] = 1.
        let mut d = [0u8; 16];
        let a = [1u8, 0, 0, 0, 0, 0, 0, 0];
        let b = [1u8, 0, 0, 0, 0, 0, 0, 0];
        bp_muladd_64(&mut d, &a, &b, &mut []);
        assert_eq!(d[0], 1);
        assert_eq!(d[1..], [0u8; 15]);
    }

    #[test]
    fn bp_muladd_128_identity_poly() {
        // a * 1 = a: muladd_128(d=0, a, b=[1,...]) → d = a in first 16 bytes.
        let mut a = [0u8; 16];
        for i in 0..16 {
            a[i] = (i as u8).wrapping_mul(7).wrapping_add(3);
        }
        let mut one = [0u8; 16];
        one[0] = 1;
        let mut d = [0u8; 32];
        let mut tmp = [0u8; 32];
        bp_muladd_128(&mut d, &a, &one, &mut tmp);
        assert_eq!(&d[..16], &a[..]);
        assert_eq!(&d[16..], &[0u8; 16]);
    }

    #[test]
    fn bp_muladd_256_identity_poly() {
        // a * 1 = a: muladd_256(d=0, a, b=[1,...]) → d[0..32] = a.
        let mut a = [0u8; 32];
        for i in 0..32 {
            a[i] = (i as u8).wrapping_mul(11).wrapping_add(5);
        }
        let mut one = [0u8; 32];
        one[0] = 1;
        let mut d = [0u8; 64];
        let mut tmp = [0u8; 96];
        bp_muladd_256(&mut d, &a, &one, &mut tmp);
        assert_eq!(&d[..32], &a[..]);
        assert_eq!(&d[32..], &[0u8; 32]);
    }

    #[test]
    fn bp_mulmod_512_zero_input() {
        // a * 0 = 0 in GF(2)[X]/(X^512+1). Verified against C reference.
        let a = [0xFFu8; 64];
        let zero = [0u8; 64];
        let mut d = [0u8; 64];
        let mut tmp = [0u8; 160];
        bp_mulmod_512(&mut d, &a, &zero, &mut tmp);
        assert_eq!(d, [0u8; 64]);
    }

    #[test]
    fn bp_mulmod_512_times_one() {
        // b * 1 = b for a sparse monomial input (a[0]=1, rest zero).
        // This verifies the port handles sparse inputs correctly — the
        // "nonzero halves" claim from the original implementer was wrong;
        // the bug was an incorrect Karatsuba offset, now fixed.
        let mut b = [0u8; 64];
        for i in 0..64 {
            b[i] = (i as u8).wrapping_mul(17).wrapping_add(3);
        }
        let mut one = [0u8; 64];
        one[0] = 1; // polynomial 1
        let mut d = [0u8; 64];
        let mut tmp = [0u8; 160];
        bp_mulmod_512(&mut d, &b, &one, &mut tmp);
        assert_eq!(d, b);
    }

    #[test]
    fn bp_mulmod_512_commutes() {
        // a*b = b*a in GF(2)[X]/(X^512+1).
        let mut a = [0u8; 64];
        let mut b = [0u8; 64];
        for i in 0..64 {
            a[i] = (i as u8).wrapping_mul(17).wrapping_add(3);
            b[i] = (i as u8).wrapping_mul(31).wrapping_add(7);
        }
        let mut d1 = [0u8; 64];
        let mut d2 = [0u8; 64];
        let mut tmp = [0u8; 160];
        bp_mulmod_512(&mut d1, &a, &b, &mut tmp);
        bp_mulmod_512(&mut d2, &b, &a, &mut tmp);
        assert_eq!(d1, d2);
    }

    #[test]
    fn bp_mulmod_512_distributive() {
        // (a+b)*c = a*c + b*c for dense polynomials.
        let mut a = [0u8; 64];
        let mut b = [0u8; 64];
        let mut c = [0u8; 64];
        for i in 0..64 {
            a[i] = (i as u8).wrapping_mul(17).wrapping_add(3);
            b[i] = (i as u8).wrapping_mul(23).wrapping_add(7);
            c[i] = (i as u8).wrapping_mul(31).wrapping_add(11);
        }
        let mut ab = [0u8; 64];
        for i in 0..64 {
            ab[i] = a[i] ^ b[i];
        }
        let mut d1 = [0u8; 64]; // (a^b)*c
        let mut d2 = [0u8; 64]; // a*c
        let mut d3 = [0u8; 64]; // b*c
        let mut tmp = [0u8; 160];
        bp_mulmod_512(&mut d1, &ab, &c, &mut tmp);
        bp_mulmod_512(&mut d2, &a, &c, &mut tmp);
        bp_mulmod_512(&mut d3, &b, &c, &mut tmp);
        // (a^b)*c should equal a*c ^ b*c
        let mut expected = [0u8; 64];
        for i in 0..64 {
            expected[i] = d2[i] ^ d3[i];
        }
        assert_eq!(d1, expected);
    }

    #[test]
    fn bp_mulmod_512_x_times_x511_equals_1() {
        // X * X^511 = X^512 = 1 mod X^512+1 (since X^512 ≡ 1).
        // Verified against C reference.
        let mut a = [0u8; 64];
        a[0] = 2; // polynomial X (bit 1)
        let mut b = [0u8; 64];
        b[63] = 0x80; // polynomial X^511 (bit 511)
        let mut d = [0u8; 64];
        let mut tmp = [0u8; 160];
        bp_mulmod_512(&mut d, &a, &b, &mut tmp);
        // result should be polynomial 1 (bit 0)
        assert_eq!(d[0], 1);
        assert_eq!(&d[1..], &[0u8; 63]);
    }

    #[test]
    fn basis_m2_mul_with_zero_matrix_gives_zero() {
        // B = [[0 0],[0 0]] → t0 = t1 = 0 regardless of h.
        let zero = [0u8; 64];
        let h0 = [0xFFu8; 64];
        let h1 = [0xAAu8; 64];
        let mut t0 = [0u8; 64];
        let mut t1 = [0u8; 64];
        let mut tmp = [0u8; 224];
        basis_m2_mul(
            &mut t0, &mut t1, &h0, &h1, &zero, &zero, &zero, &zero, &mut tmp,
        );
        assert_eq!(t0, [0u8; 64]);
        assert_eq!(t1, [0u8; 64]);
    }

    #[test]
    fn basis_m2_mul_linearity() {
        // Verify basis_m2_mul is linear in h:
        //   B*(h+h') = B*h + B*h'
        // Use dense polynomials throughout.
        let mut f2 = [0u8; 64];
        let mut g2 = [0u8; 64];
        let mut f_cap2 = [0u8; 64];
        let mut g_cap2 = [0u8; 64];
        let mut h0a = [0u8; 64];
        let mut h1a = [0u8; 64];
        let mut h0b = [0u8; 64];
        let mut h1b = [0u8; 64];
        for i in 0..64 {
            f2[i] = (i as u8).wrapping_mul(17).wrapping_add(3);
            g2[i] = (i as u8).wrapping_mul(19).wrapping_add(5);
            f_cap2[i] = (i as u8).wrapping_mul(23).wrapping_add(7);
            g_cap2[i] = (i as u8).wrapping_mul(29).wrapping_add(11);
            h0a[i] = (i as u8).wrapping_mul(31).wrapping_add(13);
            h1a[i] = (i as u8).wrapping_mul(37).wrapping_add(17);
            h0b[i] = (i as u8).wrapping_mul(41).wrapping_add(19);
            h1b[i] = (i as u8).wrapping_mul(43).wrapping_add(23);
        }
        let mut h0ab = [0u8; 64];
        let mut h1ab = [0u8; 64];
        for i in 0..64 {
            h0ab[i] = h0a[i] ^ h0b[i];
            h1ab[i] = h1a[i] ^ h1b[i];
        }
        let mut ta0 = [0u8; 64];
        let mut ta1 = [0u8; 64];
        let mut tb0 = [0u8; 64];
        let mut tb1 = [0u8; 64];
        let mut tab0 = [0u8; 64];
        let mut tab1 = [0u8; 64];
        let mut tmp = [0u8; 224];
        basis_m2_mul(
            &mut ta0, &mut ta1, &h0a, &h1a, &f2, &g2, &f_cap2, &g_cap2, &mut tmp,
        );
        basis_m2_mul(
            &mut tb0, &mut tb1, &h0b, &h1b, &f2, &g2, &f_cap2, &g_cap2, &mut tmp,
        );
        basis_m2_mul(
            &mut tab0, &mut tab1, &h0ab, &h1ab, &f2, &g2, &f_cap2, &g_cap2, &mut tmp,
        );
        let mut expected0 = [0u8; 64];
        let mut expected1 = [0u8; 64];
        for i in 0..64 {
            expected0[i] = ta0[i] ^ tb0[i];
            expected1[i] = ta1[i] ^ tb1[i];
        }
        assert_eq!(tab0, expected0);
        assert_eq!(tab1, expected1);
    }
}
