//! HAWK-512 post-quantum signature scheme, pure Rust port.
//!
//! Implements the HAWK signature scheme (eprint 2022/1155) at the NIST Level I
//! parameter set (n=512). This crate is a pure Rust port of the reference
//! implementation at `github.com/hawk-sign/dev` (MIT-licensed), with an
//! optional feature flag to enable a dev-time FFI cross-check harness that
//! byte-diffs our Rust output against the reference C on every operation.
//!
//! # Security status
//!
//! **POC-grade only.** This implementation has not been audited. The upstream
//! HAWK authors explicitly note "no security review yet — use at own risk."
//! See `SECURITY.md` for details. Constant-time posture is best-effort; no
//! formal side-channel analysis has been performed.
//!
//! # API
//!
//! ```ignore
//! use pqe_hawk::HawkKeypair;
//! use rand::rngs::OsRng;
//!
//! let kp = HawkKeypair::generate(&mut OsRng);
//! let sig = kp.secret.sign(b"hello", &mut OsRng).unwrap();
//! assert!(kp.public.verify(b"hello", &sig).is_ok());
//! ```

pub mod error;
pub mod params;

// User-facing modules (types, sign/verify entry points).
pub mod keygen;
pub mod sign;
pub mod verify;

// Implementation modules — `pub` only so FFI cross-check tests in
// tests/cross_check.rs can reach them. Hidden from rustdoc since they are
// not part of the stable public API.
#[doc(hidden)]
pub mod hash;
#[doc(hidden)]
pub mod ntt;
#[doc(hidden)]
pub mod ring;
#[doc(hidden)]
pub mod sample;
#[doc(hidden)]
pub mod serialize;

pub use error::HawkError;
pub use keygen::{HawkKeypair, HawkPublicKey, HawkSecretKey};
pub use sign::HawkSignature;
