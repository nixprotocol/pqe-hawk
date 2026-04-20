//! Signature sym-break check: sign of the first non-zero coefficient.
//!
//! Port of `poly_symbreak` (hawk_sign.c:658-671).

/// Returns the sign of the first non-zero coefficient of `s`:
///
/// - `+1` (as `u32` bit pattern for `1`) when the first non-zero coeff is
///   positive,
/// - `0xFFFF_FFFF` (bit pattern for `-1i32` cast to u32, i.e. `tbmask(x)|1`)
///   when the first non-zero coeff is negative,
/// - `0` when `s` is all zeros.
///
/// The C returns an `int32_t` that is one of {-1, 0, 1} (via the expression
/// `tbmask(x) | 1`).  We do the same.  Cast to `i32` to get the signed value.
///
/// C reference: hawk_sign.c:658-671.
pub(crate) fn poly_symbreak(s: &[i16]) -> i32 {
    let mut r: u32 = 0;
    let mut c: u32 = 0xFFFF_FFFF;
    for &coef in s.iter() {
        // C: uint32_t x = (uint32_t)s[u];
        // int16_t → int32_t (sign-extend) → uint32_t (reinterpret)
        let x = coef as i32 as u32;
        // nz = c & tbmask(x | -x): all-ones when c is still live AND x != 0
        // For x=0: x | -x = 0, tbmask(0) = 0, nz = 0.
        // For x≠0: x | wrapping_neg(x) has its top bit set, tbmask = 0xFFFFFFFF.
        let nz = c & (((x | x.wrapping_neg()) as i32 >> 31) as u32);
        c &= !nz;
        // tbmask(x): all-ones if x ≥ 0x80000000 (i.e. sign bit set → negative i16)
        let sign_mask = (x as i32 >> 31) as u32;
        r |= nz & (sign_mask | 1);
    }
    r as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_zeros_returns_zero() {
        let s = vec![0i16; 16];
        assert_eq!(poly_symbreak(&s), 0);
    }

    #[test]
    fn first_positive_returns_one() {
        let mut s = vec![0i16; 16];
        s[5] = 3;
        assert_eq!(poly_symbreak(&s), 1);
    }

    #[test]
    fn first_negative_returns_neg_one() {
        let mut s = vec![0i16; 16];
        s[5] = -3;
        assert_eq!(poly_symbreak(&s), -1);
    }

    #[test]
    fn later_nonzero_ignored() {
        let mut s = vec![0i16; 16];
        s[3] = 1;
        s[10] = -100;
        assert_eq!(poly_symbreak(&s), 1);
    }

    #[test]
    fn neg_before_pos() {
        let mut s = vec![0i16; 16];
        s[2] = -5;
        s[7] = 3;
        assert_eq!(poly_symbreak(&s), -1);
    }

    #[test]
    fn single_element_negative() {
        assert_eq!(poly_symbreak(&[-1i16]), -1);
    }

    #[test]
    fn single_element_positive() {
        assert_eq!(poly_symbreak(&[1i16]), 1);
    }

    #[test]
    fn min_i16() {
        // i16::MIN = -32768 is negative
        assert_eq!(poly_symbreak(&[i16::MIN]), -1);
    }
}
