//! SHAKE-256 hashing primitives for HAWK.
//!
//! HAWK uses two distinct SHAKE-256 absorb/squeeze operations during signing
//! and verification (see c-reference/hawk-512/hawk_sign.c lines 884-966):
//!
//! 1. `hm = SHAKE256(message || hpub)` — 64-byte message digest binding
//!    the message to the public key.
//! 2. `h = SHAKE256(hm || salt)` — 2n-bit target derived per signature
//!    attempt, conceptually two n-bit vectors `(h0, h1)`.
//!
//! These are exposed at the C's natural granularity. The bit-vector unpacking
//! happens in [`crate::sample`].

use crate::params::HAWK_N;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;

/// Compute `hm = SHAKE256(msg || hpub)` (64 bytes).
///
/// Port of c-reference/hawk-512/hawk_sign.c:893-897. The C version reuses
/// a partially-absorbed `shake_context` (already containing `msg`) and only
/// injects `hpub` before flipping; this Rust API takes the message bytes
/// directly for simplicity. The hashing semantics are identical.
pub fn compute_hm(msg: &[u8], hpub: &[u8]) -> [u8; 64] {
    let mut shake = Shake256::default();
    shake.update(msg);
    shake.update(hpub);
    let mut reader = shake.finalize_xof();
    let mut hm = [0u8; 64];
    reader.read(&mut hm);
    hm
}

/// Length of one packed bit-vector h0 / h1 in bytes (n bits).
pub const H_HALF_LEN: usize = HAWK_N / 8;

/// Compute `h = SHAKE256(hm || salt)` and return as two packed `n`-bit
/// vectors `(h0, h1)`, each `H_HALF_LEN` bytes.
///
/// Port of c-reference/hawk-512/hawk_sign.c:962-966. The C extracts a single
/// `2n / 8` byte buffer; this Rust API splits that into `(h0, h1)` since
/// the two halves are consumed independently downstream.
pub fn compute_h(hm: &[u8; 64], salt: &[u8]) -> ([u8; H_HALF_LEN], [u8; H_HALF_LEN]) {
    let mut shake = Shake256::default();
    shake.update(hm.as_ref());
    shake.update(salt);
    let mut reader = shake.finalize_xof();
    let mut h0 = [0u8; H_HALF_LEN];
    let mut h1 = [0u8; H_HALF_LEN];
    reader.read(&mut h0);
    reader.read(&mut h1);
    (h0, h1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hm_lengths() {
        let hm = compute_hm(&[], &[]);
        assert_eq!(hm.len(), 64);
    }

    #[test]
    fn h_lengths() {
        let (h0, h1) = compute_h(&[0u8; 64], &[]);
        assert_eq!(h0.len(), HAWK_N / 8);
        assert_eq!(h1.len(), HAWK_N / 8);
    }

    #[test]
    fn hm_changes_with_input() {
        let a = compute_hm(b"hello", b"world");
        let b = compute_hm(b"hello", b"World");
        assert_ne!(a, b);
        let c = compute_hm(b"Hello", b"world");
        assert_ne!(a, c);
    }

    #[test]
    fn h_changes_with_salt() {
        let hm = [0u8; 64];
        let (h0_a, h1_a) = compute_h(&hm, &[1u8; 24]);
        let (h0_b, h1_b) = compute_h(&hm, &[2u8; 24]);
        assert_ne!(h0_a, h0_b);
        assert_ne!(h1_a, h1_b);
    }

    #[test]
    fn h_halves_differ() {
        // h0 and h1 are consecutive squeezes of the same SHAKE state, so
        // overwhelmingly should not coincide for any reasonable input.
        let (h0, h1) = compute_h(&[0u8; 64], &[]);
        assert_ne!(h0, h1);
    }

    #[test]
    fn empty_inputs_produce_known_constants() {
        // SHAKE256(empty) — deterministic, so this test pins the output.
        // Computed by running compute_hm(b"", b"") locally; if this changes
        // we have a wrapper bug. Verified against `sha3` crate reference impl.
        let hm = compute_hm(&[], &[]);
        // SHAKE256("") first 8 bytes (well-known KAT):
        // 46 b9 dd 2b 0b a8 8d 13
        assert_eq!(&hm[..8], &[0x46, 0xb9, 0xdd, 0x2b, 0x0b, 0xa8, 0x8d, 0x13]);
    }
}
