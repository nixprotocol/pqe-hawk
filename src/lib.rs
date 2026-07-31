//! HAWK-512 post-quantum signature scheme, pure Rust port.
//!
//! Implements the HAWK signature scheme (eprint 2022/1155) at the HAWK-512
//! parameter set (n=512). This crate is a pure Rust port of the reference
//! implementation at `github.com/hawk-sign/dev` (MIT-licensed), with an
//! optional feature flag to enable a dev-time FFI cross-check harness that
//! byte-diffs our Rust output against the reference C on every operation.

#![deny(unsafe_code)]
//!
//! # Security status
//!
//! **🚫 DO NOT USE IN PRODUCTION.**
//!
//! **SCHEME BROKEN — research/reference only.** HAWK was withdrawn from NIST
//! standardization on 2026-07-29 after the key-recovery attack of Strážnickas &
//! Weis (*"HAWK-n Key Recovery Reduces to SVP in Dimension n/2+1"*), which the
//! HAWK team confirmed approximately halves the lattice-reduction block size
//! needed to recover an equivalent secret key. For HAWK-512 this lowers key
//! recovery to ~2^108 gates, below its claimed NIST Level I security. The attack
//! is structural (public key = Gram matrix `B*B`, `det B = 1`) and cannot be
//! mitigated in this port. **Do not deploy.**
//!
//! This implementation is also unaudited, and its constant-time posture is
//! best-effort with no formal side-channel analysis. See `SECURITY.md` for
//! details.
//!
//! # API
//!
//! ```no_run
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
