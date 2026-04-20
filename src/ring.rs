//! Polynomials over Z_q[X]/(X^n+1) where n=HAWK_N, q=HAWK_Q.
//!
//! This is HAWK's ring — NOT the Frog ring used by pqe-ring. Do not reuse.
//! Reference: `c-reference/hawk-512/ng_inner.h` — `mp_add` / `mp_sub`
//! (coefficient-wise, standard [0, q) representation).

use crate::params::{HAWK_N, HAWK_Q};

/// A polynomial in Z_q[X]/(X^n+1) with coefficients stored little-endian
/// (coefficient at index i is the coefficient of X^i).
///
/// Coefficients are stored reduced mod q in [0, q).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Poly {
    /// Coefficients are public to allow FFI cross-check tests (Tasks 9–10) to
    /// read them directly without going through an accessor.
    pub coeffs: [u16; HAWK_N],
}

impl Poly {
    /// The additive identity: all coefficients zero.
    pub const ZERO: Poly = Poly {
        coeffs: [0; HAWK_N],
    };

    pub fn new(coeffs: [u16; HAWK_N]) -> Self {
        Poly { coeffs }
    }

    /// Coefficient-wise addition mod q.
    ///
    /// Ported from `ng_inner.h` — `mp_add(a, b, p)`:
    ///   d = a + b - p; return d + (p & tbmask(d));
    /// where `tbmask(x)` propagates the sign bit of x treated as a signed i32
    /// (all-ones mask if x >= 2^31, else zero). Works for inputs in [0, q).
    pub fn add(&self, other: &Self) -> Self {
        let q = HAWK_Q;
        let mut out = [0u16; HAWK_N];
        for i in 0..HAWK_N {
            let a = self.coeffs[i] as u32;
            let b = other.coeffs[i] as u32;
            let d = a.wrapping_add(b).wrapping_sub(q);
            // tbmask: propagate sign bit (all-ones if d >= 2^31, else 0)
            let mask = ((d as i32) >> 31) as u32;
            out[i] = d.wrapping_add(q & mask) as u16;
        }
        Poly { coeffs: out }
    }

    /// Coefficient-wise subtraction mod q.
    ///
    /// Ported from `ng_inner.h` — `mp_sub(a, b, p)`:
    ///   d = a - b; return d + (p & tbmask(d));
    /// where `tbmask(x)` propagates the sign bit of x treated as a signed i32
    /// (all-ones mask if x >= 2^31, else zero). Works for inputs in [0, q).
    pub fn sub(&self, other: &Self) -> Self {
        let q = HAWK_Q;
        let mut out = [0u16; HAWK_N];
        for i in 0..HAWK_N {
            let a = self.coeffs[i] as u32;
            let b = other.coeffs[i] as u32;
            let d = a.wrapping_sub(b);
            // tbmask: propagate sign bit (all-ones if d >= 2^31, else 0)
            let mask = ((d as i32) >> 31) as u32;
            out[i] = d.wrapping_add(q & mask) as u16;
        }
        Poly { coeffs: out }
    }

    /// Schoolbook polynomial multiplication in Z_q[X]/(X^n+1).
    ///
    /// O(n²) — used as a reference implementation; NTT (in `crate::ntt`) is faster.
    /// The negacyclic reduction (X^n = -1) is applied for indices i+j >= n.
    /// Inputs and outputs are coefficients in [0, q).
    pub fn mul_schoolbook(&self, other: &Self) -> Self {
        let q = HAWK_Q;
        let mut out = [0u32; HAWK_N];
        for i in 0..HAWK_N {
            for j in 0..HAWK_N {
                let prod = (self.coeffs[i] as u32) * (other.coeffs[j] as u32);
                if i + j < HAWK_N {
                    out[i + j] = (out[i + j] + prod) % q;
                } else {
                    // X^{i+j} = -X^{i+j-n} (negacyclic reduction, since X^n = -1)
                    let idx = i + j - HAWK_N;
                    out[idx] = (out[idx] + q - (prod % q)) % q;
                }
            }
        }
        let mut result = [0u16; HAWK_N];
        for i in 0..HAWK_N {
            result[i] = out[i] as u16;
        }
        Poly { coeffs: result }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_is_additive_identity() {
        let a = Poly::new([5; HAWK_N]);
        let z = Poly::ZERO;
        assert_eq!(a.add(&z), a);
        assert_eq!(z.add(&a), a);
    }

    #[test]
    fn self_sub_is_zero() {
        let a = Poly::new([5; HAWK_N]);
        assert_eq!(a.sub(&a), Poly::ZERO);
    }

    #[test]
    fn add_is_commutative() {
        let a = Poly::new([3; HAWK_N]);
        let b = Poly::new([7; HAWK_N]);
        assert_eq!(a.add(&b), b.add(&a));
    }

    #[test]
    fn mul_by_one_returns_self() {
        let mut one_coeffs = [0u16; HAWK_N];
        one_coeffs[0] = 1;
        let one = Poly::new(one_coeffs);
        let a = Poly::new([3; HAWK_N]);
        assert_eq!(a.mul_schoolbook(&one), a);
    }

    #[test]
    fn mul_x_shifts_coefficients() {
        // In Z[X]/(X^n+1), X · (a_0 + a_1 X + ... + a_{n-1} X^{n-1})
        // = a_0 X + ... + a_{n-2} X^{n-1} + a_{n-1} X^n
        // = a_0 X + ... + a_{n-2} X^{n-1} - a_{n-1}  (since X^n = -1)
        let mut a_coeffs = [0u16; HAWK_N];
        for i in 0..HAWK_N {
            a_coeffs[i] = (i as u16) % HAWK_Q as u16;
        }
        let a = Poly::new(a_coeffs);

        let mut x_coeffs = [0u16; HAWK_N];
        x_coeffs[1] = 1;
        let x = Poly::new(x_coeffs);

        let shifted = a.mul_schoolbook(&x);
        // shifted[0] = -a[n-1] mod q = q - a[n-1]
        let expected_0 = ((HAWK_Q as u16) - a_coeffs[HAWK_N - 1]) % (HAWK_Q as u16);
        assert_eq!(shifted.coeffs[0], expected_0);
        // shifted[1] = a[0], etc.
        for i in 1..HAWK_N {
            assert_eq!(shifted.coeffs[i], a_coeffs[i - 1]);
        }
    }
}
