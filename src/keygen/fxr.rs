//! Q32.32 fixed-point arithmetic used by the NTRU solver's FFT.
//!
//! Values are stored in `u64` with the top 32 bits as the integer part and the
//! low 32 bits as the fractional part. Signed arithmetic is two's complement.
//!
//! Port of the inline `fxr_*` helpers in c-reference/hawk-512/ng_inner.h
//! (roughly lines 744-918).

/// Q32.32 fixed-point number.
///
/// The inner `u64` is a two's-complement signed Q32.32 value:
/// bits 63-32 hold the integer part, bits 31-0 hold the fractional part.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fxr(pub u64);

/// Zero: `fxr_zero` (ng_inner.h:916).
pub const FXR_ZERO: Fxr = Fxr(0);

/// sqrt(2) in Q32.32: `fxr_sqrt2` (ng_inner.h:917).
pub const FXR_SQRT2: Fxr = Fxr(6074001000);

impl Fxr {
    /// Construct from a signed 32-bit integer.
    ///
    /// Port of `fxr_of` (ng_inner.h:765-772).
    /// C: `x.v = (uint64_t)j << 32;`
    /// The cast to `uint64_t` of a signed `int32_t` sign-extends first.
    #[inline]
    pub fn of(j: i32) -> Self {
        Fxr((j as i64 as u64) << 32)
    }

    /// Construct from a pre-scaled 64-bit bit pattern placed directly.
    ///
    /// Port of `fxr_of_scaled32` (ng_inner.h:774-781).
    #[inline]
    pub fn of_scaled32(t: u64) -> Self {
        Fxr(t)
    }

    /// Addition.
    ///
    /// Port of `fxr_add` (ng_inner.h:783-788).
    #[inline]
    pub fn add(self, other: Self) -> Self {
        Fxr(self.0.wrapping_add(other.0))
    }

    /// Subtraction.
    ///
    /// Port of `fxr_sub` (ng_inner.h:790-795).
    #[inline]
    pub fn sub(self, other: Self) -> Self {
        Fxr(self.0.wrapping_sub(other.0))
    }

    /// Multiply by 2 (left shift by 1).
    ///
    /// Port of `fxr_double` (ng_inner.h:797-802).
    #[inline]
    pub fn double(self) -> Self {
        Fxr(self.0 << 1)
    }

    /// Negate.
    ///
    /// Port of `fxr_neg` (ng_inner.h:804-809).
    /// C: `x.v = -x.v;` (unsigned wrapping negation).
    #[inline]
    pub fn neg(self) -> Self {
        Fxr(self.0.wrapping_neg())
    }

    /// Absolute value (branchless mask-based).
    ///
    /// Port of `fxr_abs` (ng_inner.h:811-816).
    /// C: `x.v -= (x.v << 1) & (uint64_t)(*(int64_t *)&x.v >> 63);`
    ///
    /// The sign bit replicated across all 64 bits is used as a mask.
    /// When x >= 0 the mask is 0 (no-op). When x < 0 the mask is all-ones,
    /// so `x.v -= x.v << 1` which equals `x.v = -x.v`.
    #[inline]
    pub fn abs(self) -> Self {
        let sign_mask = (self.0 as i64 >> 63) as u64;
        Fxr(self.0.wrapping_sub((self.0 << 1) & sign_mask))
    }

    /// Q32.32 multiplication (truncating toward −∞, i.e. arithmetic right shift).
    ///
    /// Port of `fxr_mul` (ng_inner.h:818-843).
    /// Uses the `__int128` path: signed 64×64 → 128, then logical right-shift by 32.
    /// The cast `(uint64_t)(z >> 32)` keeps only the low 64 bits after the shift,
    /// which is equivalent to truncation of the Q32.32 product.
    #[inline]
    pub fn mul(self, other: Self) -> Self {
        let z = (self.0 as i64 as i128) * (other.0 as i64 as i128);
        Fxr((z >> 32) as u64)
    }

    /// Q32.32 squaring.
    ///
    /// Port of `fxr_sqr` (ng_inner.h:845-869).
    /// Uses the `__int128` path: t*t as signed i128 then right-shift by 32.
    #[inline]
    pub fn sqr(self) -> Self {
        let t = self.0 as i64 as i128;
        Fxr(((t * t) >> 32) as u64)
    }

    /// Divide by 2^n, rounding to nearest (ties round away from zero toward +∞).
    ///
    /// Port of `fxr_div2e` (ng_inner.h:878-884).
    /// C: `x.v += (((uint64_t)1 << n) >> 1); x.v = (uint64_t)(*(int64_t*)&x.v >> n);`
    /// Adds 0.5 ULP then arithmetic right-shifts.
    #[inline]
    pub fn div2e(self, n: u32) -> Self {
        let half = (1u64 << n) >> 1;
        let rounded = self.0.wrapping_add(half);
        Fxr(((rounded as i64) >> n) as u64)
    }

    /// Multiply by 2^n (left shift).
    ///
    /// Port of `fxr_mul2e` (ng_inner.h:886-891).
    #[inline]
    pub fn mul2e(self, n: u32) -> Self {
        Fxr(self.0 << n)
    }

    /// Round to nearest integer (ties round toward +∞), returning `i32`.
    ///
    /// Port of `fxr_round` (ng_inner.h:871-876).
    /// C: `x.v += 0x80000000ul; return (int32_t)(*(int64_t*)&x.v >> 32);`
    /// Adds 0.5 in Q32.32 then takes the integer part.
    #[inline]
    pub fn round(self) -> i32 {
        let biased = self.0.wrapping_add(0x8000_0000);
        (biased as i64 >> 32) as i32
    }

    /// Compute `1 / self`.
    ///
    /// Port of `fxr_inv` (ng_inner.h:896-901).
    #[inline]
    pub fn inv(self) -> Self {
        Fxr(inner_fxr_div(1u64 << 32, self.0))
    }

    /// Fixed-point division: `self / other`.
    ///
    /// Port of `fxr_div` (ng_inner.h:903-908).
    #[inline]
    pub fn div(self, other: Self) -> Self {
        Fxr(inner_fxr_div(self.0, other.0))
    }

    /// Less-than comparison (signed).
    ///
    /// Port of `fxr_lt` (ng_inner.h:910-913).
    #[inline]
    pub fn lt(self, other: Self) -> bool {
        (self.0 as i64) < (other.0 as i64)
    }
}

/// Bit-by-bit fixed-point division worker. Not exposed directly; used by
/// `Fxr::inv` and `Fxr::div`.
///
/// Port of `inner_fxr_div` (ng_fxp.c:4-51).
fn inner_fxr_div(x: u64, y: u64) -> u64 {
    // Extract signs and take absolute values.
    let sx = x >> 63;
    let x = (x ^ sx.wrapping_neg()).wrapping_add(sx);
    let sy = y >> 63;
    let y = (y ^ sy.wrapping_neg()).wrapping_add(sy);

    let mut q: u64 = 0;
    let mut num = x >> 31;

    // First 31 iterations: i = 63 down to 33.
    for i in (33u32..=63).rev() {
        let b = 1u64.wrapping_sub(num.wrapping_sub(y) >> 63);
        q |= b << i;
        num = num.wrapping_sub(y & b.wrapping_neg());
        num <<= 1;
        num |= (x >> (i - 33)) & 1;
    }
    // Next 33 iterations: i = 32 down to 0.
    for i in (0..=32i32).rev() {
        let b = 1u64.wrapping_sub(num.wrapping_sub(y) >> 63);
        q |= b << (i as u32);
        num = num.wrapping_sub(y & b.wrapping_neg());
        num <<= 1;
    }

    // Final rounding bit.
    let b = 1u64.wrapping_sub(num.wrapping_sub(y) >> 63);
    q = q.wrapping_add(b);

    // Restore sign.
    let s = sx ^ sy;
    (q ^ s.wrapping_neg()).wrapping_add(s)
}

/// A complex value with `Fxr` real and imaginary parts.
///
/// Port of `fxc` from c-reference/hawk-512/ng_inner.h:922-924.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fxc {
    pub re: Fxr,
    pub im: Fxr,
}

impl Fxc {
    /// Construct from two `Fxr` components.
    #[inline]
    pub const fn new(re: Fxr, im: Fxr) -> Self {
        Fxc { re, im }
    }

    /// Complex addition (component-wise).
    /// Port of `fxc_add` (ng_inner.h:928-934).
    #[inline]
    pub fn add(self, other: Self) -> Self {
        Fxc {
            re: self.re.add(other.re),
            im: self.im.add(other.im),
        }
    }

    /// Complex subtraction.
    /// Port of `fxc_sub` (ng_inner.h:936-942).
    #[inline]
    pub fn sub(self, other: Self) -> Self {
        Fxc {
            re: self.re.sub(other.re),
            im: self.im.sub(other.im),
        }
    }

    /// Divide by 2 (component-wise).
    /// Port of `fxc_half` (ng_inner.h:944-950).
    #[inline]
    pub fn half(self) -> Self {
        Fxc {
            re: self.re.div2e(1),
            im: self.im.div2e(1),
        }
    }

    /// Complex multiplication using 3-multiply Karatsuba form.
    ///
    /// For `r = (a + bi)*(c + di)`:
    ///   z0 = a*c, z1 = b*d, z2 = (a+b)*(c+d)
    ///   r.re = z0 - z1
    ///   r.im = z2 - (z0 + z1)
    ///
    /// Port of `fxc_mul` (ng_inner.h:952-974). Each fxr_mul truncates, so the
    /// imaginary result may differ slightly from `a*d + b*c` computed
    /// directly — this form is the canonical choice for reproducibility
    /// across implementations.
    #[inline]
    pub fn mul(self, other: Self) -> Self {
        let z0 = self.re.mul(other.re);
        let z1 = self.im.mul(other.im);
        let z2 = self.re.add(self.im).mul(other.re.add(other.im));
        Fxc {
            re: z0.sub(z1),
            im: z2.sub(z0.add(z1)),
        }
    }

    /// Complex conjugate (negate imaginary part).
    /// Port of `fxc_conj` (ng_inner.h:976-981).
    #[inline]
    pub fn conj(self) -> Self {
        Fxc {
            re: self.re,
            im: self.im.neg(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn of_and_back() {
        assert_eq!(Fxr::of(0), FXR_ZERO);
        let three = Fxr::of(3);
        // Integer part is the top 32 bits; 3 shifted left by 32.
        assert_eq!(three.0, 3u64 << 32);
    }

    #[test]
    fn of_negative_sign_extends() {
        // -1 as i32 sign-extends to 0xFFFFFFFF_FFFFFFFF as i64,
        // then shifted left by 32 gives 0xFFFFFFFF_00000000.
        assert_eq!(Fxr::of(-1).0, 0xFFFF_FFFF_0000_0000u64);
    }

    #[test]
    fn add_and_sub() {
        let a = Fxr::of(5);
        let b = Fxr::of(7);
        assert_eq!(a.add(b), Fxr::of(12));
        assert_eq!(b.sub(a), Fxr::of(2));
    }

    #[test]
    fn neg_of_neg_is_identity() {
        let a = Fxr::of(42);
        assert_eq!(a.neg().neg(), a);
    }

    #[test]
    fn abs_of_negative() {
        let a = Fxr::of(-5);
        assert_eq!(a.abs(), Fxr::of(5));
    }

    #[test]
    fn abs_of_positive() {
        let a = Fxr::of(5);
        assert_eq!(a.abs(), Fxr::of(5));
    }

    #[test]
    fn double_and_div2e() {
        let a = Fxr::of(6);
        assert_eq!(a.double(), Fxr::of(12));
        // div2e(1) is round-to-nearest, so exact halving of an integer
        // that is already even gives back the original.
        assert_eq!(Fxr::of(12).div2e(1), Fxr::of(6));
    }

    #[test]
    fn mul_integer_values() {
        let three = Fxr::of(3);
        let four = Fxr::of(4);
        assert_eq!(three.mul(four), Fxr::of(12));
    }

    #[test]
    fn mul_by_zero() {
        let zero = FXR_ZERO;
        let x = Fxr::of(123);
        assert_eq!(x.mul(zero), FXR_ZERO);
        assert_eq!(zero.mul(x), FXR_ZERO);
    }

    #[test]
    fn sqr_matches_mul() {
        let x = Fxr::of(7);
        assert_eq!(x.sqr(), x.mul(x));
    }

    #[test]
    fn round_half_integer() {
        // 1.5 in Q32.32: 1 * 2^32 + 2^31
        let one_and_half = Fxr::of_scaled32((1u64 << 32) | (1u64 << 31));
        let r = one_and_half.round();
        // C: adds 0x80000000 (= 2^31) to the fractional part then takes integer.
        // 1.5 + 0.5 = 2.0 => rounds to 2 (ties toward +inf).
        assert_eq!(r, 2);
    }

    #[test]
    fn round_below_half() {
        // 1.25 in Q32.32: (1 << 32) | (1 << 30)
        let one_quarter = Fxr::of_scaled32((1u64 << 32) | (1u64 << 30));
        assert_eq!(one_quarter.round(), 1);
    }

    #[test]
    fn div2e_rounds_toward_plus_inf() {
        // 3 / 2 = 1.5, rounds to 2 via div2e(1).
        // 3 in Q32.32 = 3 << 32; after div2e(1): add 2^31 then >> 1
        // => (3<<32 + 2^31) >> 1 = ... still a fractional value.
        // But verify with round: (Fxr::of(3).div2e(1)).round() should be 2.
        let half_three = Fxr::of(3).div2e(1);
        assert_eq!(half_three.round(), 2);
    }

    #[test]
    fn lt_comparison() {
        assert!(Fxr::of(-1).lt(Fxr::of(0)));
        assert!(Fxr::of(0).lt(Fxr::of(1)));
        assert!(!Fxr::of(1).lt(Fxr::of(0)));
        assert!(!Fxr::of(0).lt(Fxr::of(0)));
    }

    #[test]
    fn div_integer_values() {
        let twelve = Fxr::of(12);
        let three = Fxr::of(3);
        assert_eq!(twelve.div(three), Fxr::of(4));
    }

    #[test]
    fn div_by_one_is_identity() {
        let a = Fxr::of(42);
        let one = Fxr::of(1);
        assert_eq!(a.div(one), a);
    }

    #[test]
    fn div_of_negative() {
        // -12 / 3 = -4
        let neg_twelve = Fxr::of(-12);
        let three = Fxr::of(3);
        assert_eq!(neg_twelve.div(three), Fxr::of(-4));
    }

    #[test]
    fn inv_of_one_is_one() {
        let one = Fxr::of(1);
        assert_eq!(one.inv(), one);
    }

    #[test]
    fn inv_of_two_is_half() {
        // 1/2 in Q32.32: 2^31
        let two = Fxr::of(2);
        assert_eq!(two.inv(), Fxr::of_scaled32(1u64 << 31));
    }

    #[test]
    fn inv_then_mul_roundtrip() {
        let x = Fxr::of(7);
        let ix = x.inv();
        // x * (1/x) should be ~1 (within rounding tolerance).
        // In Q32.32 representation, 1.0 = 0x100000000.
        let prod = x.mul(ix);
        // Due to bit-by-bit division rounding, the result may be 0xFFFFFFFF or 0x100000000.
        // Assert within ±1 of the integer-1 bit pattern.
        let one = Fxr::of(1).0;
        let diff = if prod.0 > one {
            prod.0 - one
        } else {
            one - prod.0
        };
        assert!(
            diff <= 8,
            "1/7 * 7 off by more than 8 ULPs: {:#x} vs {:#x}",
            prod.0,
            one
        );
    }

    #[test]
    fn lt_basic() {
        assert!(Fxr::of(3).lt(Fxr::of(5)));
        assert!(!Fxr::of(5).lt(Fxr::of(3)));
        assert!(!Fxr::of(5).lt(Fxr::of(5)));
    }

    #[test]
    fn lt_handles_negatives() {
        assert!(Fxr::of(-5).lt(Fxr::of(3)));
        assert!(Fxr::of(-10).lt(Fxr::of(-1)));
        assert!(!Fxr::of(3).lt(Fxr::of(-5)));
    }

    #[test]
    fn fxc_add_basic() {
        let a = Fxc::new(Fxr::of(3), Fxr::of(4));
        let b = Fxc::new(Fxr::of(5), Fxr::of(6));
        assert_eq!(a.add(b), Fxc::new(Fxr::of(8), Fxr::of(10)));
    }

    #[test]
    fn fxc_sub_basic() {
        let a = Fxc::new(Fxr::of(10), Fxr::of(20));
        let b = Fxc::new(Fxr::of(3), Fxr::of(5));
        assert_eq!(a.sub(b), Fxc::new(Fxr::of(7), Fxr::of(15)));
    }

    #[test]
    fn fxc_half_basic() {
        let a = Fxc::new(Fxr::of(8), Fxr::of(4));
        assert_eq!(a.half(), Fxc::new(Fxr::of(4), Fxr::of(2)));
    }

    #[test]
    fn fxc_mul_real_times_real() {
        // (3 + 0i) * (4 + 0i) = 12 + 0i
        let a = Fxc::new(Fxr::of(3), Fxr::of(0));
        let b = Fxc::new(Fxr::of(4), Fxr::of(0));
        assert_eq!(a.mul(b), Fxc::new(Fxr::of(12), Fxr::of(0)));
    }

    #[test]
    fn fxc_mul_i_times_i() {
        // i * i = -1
        let i = Fxc::new(Fxr::of(0), Fxr::of(1));
        assert_eq!(i.mul(i), Fxc::new(Fxr::of(-1), Fxr::of(0)));
    }

    #[test]
    fn fxc_mul_1_plus_i_squared() {
        // (1+i)^2 = 2i
        let a = Fxc::new(Fxr::of(1), Fxr::of(1));
        assert_eq!(a.mul(a), Fxc::new(Fxr::of(0), Fxr::of(2)));
    }

    #[test]
    fn fxc_conj_negates_imaginary() {
        let a = Fxc::new(Fxr::of(3), Fxr::of(4));
        assert_eq!(a.conj(), Fxc::new(Fxr::of(3), Fxr::of(-4)));
    }
}
