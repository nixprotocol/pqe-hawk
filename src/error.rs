//! HAWK-specific error types.

use thiserror::Error;

/// Errors produced by `pqe-hawk` at the public API boundary.
///
/// Every operation that can fail returns `Result<_, HawkError>`; all variants
/// are non-fatal and safe to surface to end users. No variant carries secret
/// material — the `String` payloads describe *why* a byte buffer was rejected
/// (e.g. "padding nonzero") and never echo untrusted input.
#[derive(Debug, Error)]
pub enum HawkError {
    /// Signature verification failed. Returned by
    /// [`crate::HawkPublicKey::verify`] when the signature does not
    /// authenticate the given message under the given public key. No further
    /// detail is exposed (doing so would leak which check failed, potentially
    /// aiding an attacker crafting malleability attacks).
    #[error("signature verification failed")]
    InvalidSignature,

    /// A public-key byte buffer was malformed and could not be decoded.
    /// The payload describes the specific decoder check that failed
    /// (e.g. Golomb-Rice overflow, trailing-bits nonzero, padding nonzero).
    #[error("malformed public key: {0}")]
    MalformedPublicKey(String),

    /// A secret-key byte buffer was malformed and could not be decoded.
    /// Currently only surfaced from internal validation paths — the public
    /// [`crate::HawkSecretKey::from_bytes`] is infallible on any well-sized
    /// input, since the 184-byte layout has no rejectable fields.
    #[error("malformed secret key: {0}")]
    MalformedSecretKey(String),

    /// A signature byte buffer was malformed and could not be decoded.
    /// The payload describes the specific decoder check that failed
    /// (e.g. Golomb-Rice overflow, padding nonzero).
    #[error("malformed signature: {0}")]
    MalformedSignature(String),

    /// The bimodal Gaussian rejection sampler inside `sign` exceeded its
    /// retry budget without producing a valid signature. The default budget
    /// ([`crate::params::HAWK_SAMPLER_RETRY_BUDGET`]) is large enough that
    /// this is astronomically unlikely for a correctly generated key, so the
    /// variant typically signals either (a) a bit-flipped secret key, or
    /// (b) an exhausted / compromised RNG. Callers should treat it as a
    /// hard failure, not retry.
    #[error("rejection sampler exceeded retry budget ({retries} attempts)")]
    SamplingFailure {
        /// Number of retry attempts performed before giving up.
        retries: usize,
    },
}
