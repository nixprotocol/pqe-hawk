//! HAWK-specific error types.

use thiserror::Error;

#[derive(Debug, Error)]
pub enum HawkError {
    #[error("signature verification failed")]
    InvalidSignature,

    #[error("malformed public key: {0}")]
    MalformedPublicKey(String),

    #[error("malformed secret key: {0}")]
    MalformedSecretKey(String),

    #[error("malformed signature: {0}")]
    MalformedSignature(String),

    #[error("rejection sampler exceeded retry budget ({retries} attempts)")]
    SamplingFailure { retries: usize },
}
