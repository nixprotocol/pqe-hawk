//! HAWK-512 parameter constants.
//!
//! Values transcribed from the upstream reference at
//! `c-reference/hawk-512/hawk_config.h`, `hawk_inner.h`, `hawk.h`, and `api.h`.
//!
//! Each constant has a doc comment naming the source file and the `#define`
//! it corresponds to. If upstream parameters change, update both this file
//! and regenerate KATs.

/// Ring dimension. HAWK-512 uses degree-512 polynomials.
/// Source: `hawk.h` — degree is 2^logn, with logn = 9 for HAWK-512
/// (comment at line 101: "Degree logarithm is provided as logn (8, 9 or 10)").
pub const HAWK_N: usize = 512;

/// log₂ of the ring dimension.
/// Source: `hawk.h` — logn = 9 is the HAWK-512 variant (n = 2^9 = 512).
pub const HAWK_LOGN: usize = 9;

/// HAWK's small modulus for the public Q polynomials (used in signing NTT).
/// Source: `c-reference/hawk-512/hawk_sign.c` line 6:
///   `#define Q   18433`
/// with the comment at line 4: "We use computations with polynomials modulo
/// X^n+1 and modulo 18433."
/// 18433 < 2^15 = 32768, which comfortably fits in u16.
pub const HAWK_Q: u32 = 18433;

/// Public key byte length (serialized `HawkPublicKey::to_bytes()`).
/// Source: `c-reference/hawk-512/api.h` — `CRYPTO_PUBLICKEYBYTES = 1024`.
/// Cross-check: `HAWK_PUBKEY_SIZE(9)` from `hawk.h` evaluates to
///   450 + 574*(2>>(10-9)) + 842*(1>>(10-9)) = 450 + 574 + 0 = 1024.
pub const HAWK_PUBLIC_KEY_BYTES: usize = 1024;

/// Secret key byte length (serialized `HawkSecretKey::to_bytes()`).
/// Source: `c-reference/hawk-512/api.h` — `CRYPTO_SECRETKEYBYTES = 184`.
/// Cross-check: `HAWK_PRIVKEY_SIZE(9)` from `hawk.h` evaluates to
///   8 + (11 << (9-5)) = 8 + 11*16 = 8 + 176 = 184.
pub const HAWK_SECRET_KEY_BYTES: usize = 184;

/// Signature byte length (serialized `HawkSignature::to_bytes()`).
/// Source: `c-reference/hawk-512/api.h` — `CRYPTO_BYTES = 555`.
/// Cross-check: `HAWK_SIG_SIZE(9)` from `hawk.h` evaluates to
///   249 + 306*(2>>(10-9)) + 360*(1>>(10-9)) = 249 + 306 + 0 = 555.
pub const HAWK_SIGNATURE_BYTES: usize = 555;

/// Salt length in bytes: per-signature fresh randomness absorbed by the hash.
/// Source: `c-reference/hawk-512/hawk_sign.c` lines 803–826 (comment and code):
///   "salt_len (bits): 192" for logn = 9, i.e. 192/8 = 24 bytes.
///   ```c
///   case 9:  salt_len = 24; break;
///   ```
/// Also confirmed in `hawk_vrfy.c` line 1548:
///   ```c
///   case 9:  salt_len = 24; break;
///   ```
pub const HAWK_SALT_BYTES: usize = 24;

/// Rejection-sample retry budget.
///
/// The upstream C reference (`hawk_sign.c` line 904) uses an unbounded loop:
///   ```c
///   for (uint32_t attempt = 0;; attempt += 2) { … }
///   ```
/// There is no explicit retry-limit constant in the reference implementation.
/// The reference relies on the fact that the Gaussian sampler very rarely
/// produces a signature that exceeds the size bound (less than 1-in-a-million
/// per the comment in `hawk.h`). This Rust port imposes a finite budget for
/// fail-safe behaviour: if the budget is exhausted, [`crate::HawkSecretKey::sign`]
/// returns [`crate::HawkError::SamplingFailure`] instead of looping forever.
///
/// # Probability analysis
///
/// Each retry has an independent probability `p < 2^-20` of occurring
/// (from the HAWK paper's Gaussian-sampler overflow analysis; tighter
/// empirical bounds are far lower — typical keygen produces a valid signature
/// on the *first* attempt with probability > 99.9%). For `p = 2^-20` and a
/// retry budget of 1024 = 2^10, the probability that a single sign call
/// exhausts the budget is approximately `p^1024 < 2^-20480` — far below any
/// physical or cryptographic attack threshold.
///
/// # Consensus-critical contexts
///
/// Callers in adversarial or consensus-critical contexts (e.g., block
/// producers, threshold signing) should be aware that:
///
/// 1. A malicious RNG that can force the sampler to exceed the norm bound
///    on every attempt would still take `~2^20 * 1024 = ~2^30` rejected
///    samples per failed sign call — detectable and not a useful DoS vector.
/// 2. `sign()` is not constant-time with respect to retry count; a side
///    channel observer can distinguish "sign succeeded on attempt k" for
///    different k values. For HAWK this is considered acceptable because
///    the attempt count does not leak key material.
/// 3. On `HawkError::SamplingFailure`, callers should treat it as a
///    transient error and either retry (with fresh randomness) or log and
///    fail the operation. Infinite retry loops in user code negate the
///    benefit of the finite budget.
pub const HAWK_SAMPLER_RETRY_BUDGET: usize = 1024;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn n_is_power_of_two() {
        assert!(HAWK_N.is_power_of_two());
        assert_eq!(HAWK_N, 1 << HAWK_LOGN);
    }

    #[test]
    fn n_is_512() {
        assert_eq!(HAWK_N, 512);
        assert_eq!(HAWK_LOGN, 9);
    }

    #[test]
    fn byte_lengths_are_positive() {
        assert!(HAWK_PUBLIC_KEY_BYTES > 0);
        assert!(HAWK_SECRET_KEY_BYTES > 0);
        assert!(HAWK_SIGNATURE_BYTES > 0);
    }

    #[test]
    fn q_fits_in_u16() {
        // HAWK uses a small modulus suitable for u16 coefficients.
        // If this fails, the ring module's coefficient width needs widening.
        assert!(HAWK_Q <= u16::MAX as u32 + 1);
    }
}
