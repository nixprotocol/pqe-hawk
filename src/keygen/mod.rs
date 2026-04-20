//! HAWK key generation.
//!
//! Port of c-reference/hawk-512/hawk_kgen.c + ng_hawk.c + ng_ntru.c +
//! ng_fxp.c + ng_zint31.c + ng_mp31.c. All sub-modules are byte-exact
//! against the reference C, verified via FFI cross-check.
//!
//! The main public entry point is [`hawk_keygen::hawk_keygen_512`], which
//! generates a full HAWK-512 keypair from an RNG-supplying closure.

// Implementation modules — `pub` only so FFI cross-check tests in
// tests/cross_check.rs can reach them. Hidden from rustdoc since they are
// not part of the stable public API.
#[doc(hidden)]
pub mod fxp;
#[doc(hidden)]
pub mod fxr;
#[doc(hidden)]
pub mod hawk_keygen;
#[doc(hidden)]
pub mod mp31;
#[doc(hidden)]
pub mod ntru;
#[doc(hidden)]
pub mod ntru_profile;
#[doc(hidden)]
pub mod poly;
#[doc(hidden)]
pub mod primes;
#[doc(hidden)]
pub mod q_derive;
#[doc(hidden)]
pub mod regen_fg;
#[doc(hidden)]
pub mod zint31;

pub mod types;

pub use hawk_keygen::hawk_keygen_512;
pub use types::{HawkKeypair, HawkPublicKey, HawkSecretKey};
