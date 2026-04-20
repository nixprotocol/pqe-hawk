//! HAWK-512 signature verification.
//!
//! Port of c-reference/hawk-512/hawk_vrfy.c. See
//! docs/superpowers/plans/2026-04-17-pqe-hawk-verify-decomposition.md
//! for the sub-task breakdown.

// Implementation modules — `pub` only so FFI cross-check tests can reach
// them. Hidden from rustdoc since they are not part of the stable API.
#[doc(hidden)]
pub mod consts;
#[doc(hidden)]
pub mod fx32;
#[doc(hidden)]
pub mod helpers;
#[doc(hidden)]
pub mod mp;
#[doc(hidden)]
pub mod verify_inner;

use crate::error::HawkError;
use crate::keygen::HawkPublicKey;
use crate::sign::HawkSignature;

impl HawkPublicKey {
    /// Verify a HAWK-512 signature on `msg` against this public key.
    ///
    /// Returns `Ok(())` if the signature is cryptographically valid.
    /// Returns `Err(HawkError::InvalidSignature)` if any verification check fails.
    pub fn verify(&self, msg: &[u8], sig: &HawkSignature) -> Result<(), HawkError> {
        // Encode the public key to derive hpub inside verify_inner.
        let pub_bytes = self.to_bytes()?;
        verify_inner::verify_inner(msg, &pub_bytes, &self.q00, &self.q01, &sig.salt, &sig.s1)
    }
}
