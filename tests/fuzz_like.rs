//! Quasi-fuzz tests runnable under stable Rust via proptest.
//!
//! The full coverage-guided fuzz harness is under `fuzz/` and requires
//! nightly Rust + cargo-fuzz. These tests provide a lightweight version
//! that catches the most common malformed-input panic classes in CI
//! without extra tooling.

use pqe_hawk::params::{HAWK_PUBLIC_KEY_BYTES, HAWK_SIGNATURE_BYTES};
use pqe_hawk::{HawkKeypair, HawkPublicKey, HawkSignature};
use proptest::prelude::*;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 256,
        ..ProptestConfig::default()
    })]

    /// `HawkPublicKey::from_bytes` must never panic on any 1024-byte input.
    #[test]
    fn decode_pubkey_never_panics(bytes in prop::array::uniform32(any::<u8>())) {
        // Assemble a full 1024-byte buffer from a 32-byte seed via replication;
        // this gives varied patterns while keeping the generator simple.
        let mut buf = [0u8; HAWK_PUBLIC_KEY_BYTES];
        for (i, b) in buf.iter_mut().enumerate() {
            *b = bytes[i % 32];
        }
        // Must return Ok or Err, never panic/abort.
        let _ = HawkPublicKey::from_bytes(&buf);
    }

    /// `HawkSignature::from_bytes` must never panic on any 555-byte input.
    #[test]
    fn decode_signature_never_panics(bytes in prop::array::uniform32(any::<u8>())) {
        let mut buf = [0u8; HAWK_SIGNATURE_BYTES];
        for (i, b) in buf.iter_mut().enumerate() {
            *b = bytes[i % 32];
        }
        let _ = HawkSignature::from_bytes(&buf);
    }

    /// Full verify pipeline must never panic on any (sig, msg) pair against a
    /// fixed pinned pubkey. Any panic indicates a potential DoS.
    #[test]
    fn verify_never_panics(
        sig_seed in prop::array::uniform32(any::<u8>()),
        msg in prop::collection::vec(any::<u8>(), 0..=256),
    ) {
        // Pinned pubkey for reproducibility.
        let mut kp_rng = ChaCha20Rng::from_seed([1u8; 32]);
        let pk = HawkKeypair::generate(&mut kp_rng).public;

        // Build a signature buffer from the seed.
        let mut sig_buf = [0u8; HAWK_SIGNATURE_BYTES];
        for (i, b) in sig_buf.iter_mut().enumerate() {
            *b = sig_seed[i % 32];
        }

        // Decode may fail (Err); if it succeeds, verify must also succeed or
        // fail gracefully, but never panic.
        if let Ok(sig) = HawkSignature::from_bytes(&sig_buf) {
            let _ = pk.verify(&msg, &sig);
        }
    }

    /// Mutating a single byte of a valid signature must not cause verify to
    /// panic (must return Err or — in the extraordinarily unlikely case of
    /// the mutation landing on a bit-pattern that still satisfies the
    /// verification equation — Ok).
    ///
    /// Slow (keygen + sign per case); run with fewer cases than the decoder
    /// tests above.
    #[test]
    #[cfg_attr(not(feature = "slow-tests"), ignore)]
    fn mutated_signature_verify_never_panics(
        key_seed in prop::array::uniform32(any::<u8>()),
        sign_seed in prop::array::uniform32(any::<u8>()),
        msg in prop::collection::vec(any::<u8>(), 0..=128),
        byte_idx in 0usize..HAWK_SIGNATURE_BYTES,
        bit_idx in 0u8..8,
    ) {
        let mut kp_rng = ChaCha20Rng::from_seed(key_seed);
        let kp = HawkKeypair::generate(&mut kp_rng);
        let mut sign_rng = ChaCha20Rng::from_seed(sign_seed);
        let Ok(sig) = kp.secret.sign(&msg, &mut sign_rng) else { return Ok(()); };
        let Ok(mut bytes) = sig.to_bytes() else { return Ok(()); };
        bytes[byte_idx] ^= 1 << bit_idx;
        if let Ok(mutated) = HawkSignature::from_bytes(&bytes) {
            let _ = kp.public.verify(&msg, &mutated);
        }
    }

    /// Corrupting the trailing zero-padding region of a valid pubkey must
    /// cause `from_bytes` to reject with `Err(MalformedPublicKey)` and never
    /// panic. Guards the `"padding nonzero"` check against accidental removal
    /// or off-by-one in the padding scan.
    ///
    /// The HAWK-512 pubkey Golomb-Rice encoding leaves a run of zero bytes at
    /// the tail of the 1024-byte buffer; flipping any of those bytes to a
    /// nonzero value MUST be rejected as malformed (otherwise an attacker
    /// could smuggle data in the padding region).
    #[test]
    #[cfg_attr(not(feature = "slow-tests"), ignore)]
    fn corrupted_pubkey_padding_rejected(
        key_seed in prop::array::uniform32(any::<u8>()),
        byte_offset_from_end in 1usize..=16,
        corruption in 1u8..=255,
    ) {
        let mut kp_rng = ChaCha20Rng::from_seed(key_seed);
        let pk = HawkKeypair::generate(&mut kp_rng).public;
        let mut bytes = pk.to_bytes().expect("pk encode");
        // Overwrite a byte near the end. The last few bytes of a freshly
        // encoded pubkey are (with overwhelming probability) part of the
        // zero-padding tail, not the encoded variable-length region.
        let idx = HAWK_PUBLIC_KEY_BYTES - byte_offset_from_end;
        // If the byte is already zero, corrupting to nonzero must be rejected.
        // If the byte is part of the encoded region (unlikely this close to
        // the end) the decoder still must not panic.
        let original = bytes[idx];
        bytes[idx] = corruption;
        if original == 0 {
            prop_assert!(
                HawkPublicKey::from_bytes(&bytes).is_err(),
                "corrupted padding byte at offset {idx} should be rejected"
            );
        } else {
            // Not padding — just assert no panic.
            let _ = HawkPublicKey::from_bytes(&bytes);
        }
    }

    /// Same guarantee for signature padding: overwriting the zero-tail of a
    /// valid signature must be rejected. Guards the signature `"padding
    /// nonzero"` check.
    #[test]
    #[cfg_attr(not(feature = "slow-tests"), ignore)]
    fn corrupted_signature_padding_rejected(
        key_seed in prop::array::uniform32(any::<u8>()),
        sign_seed in prop::array::uniform32(any::<u8>()),
        msg in prop::collection::vec(any::<u8>(), 0..=64),
        byte_offset_from_end in 1usize..=16,
        corruption in 1u8..=255,
    ) {
        let mut kp_rng = ChaCha20Rng::from_seed(key_seed);
        let kp = HawkKeypair::generate(&mut kp_rng);
        let mut sign_rng = ChaCha20Rng::from_seed(sign_seed);
        let Ok(sig) = kp.secret.sign(&msg, &mut sign_rng) else { return Ok(()); };
        let Ok(mut bytes) = sig.to_bytes() else { return Ok(()); };
        let idx = HAWK_SIGNATURE_BYTES - byte_offset_from_end;
        let original = bytes[idx];
        bytes[idx] = corruption;
        if original == 0 {
            prop_assert!(
                HawkSignature::from_bytes(&bytes).is_err(),
                "corrupted signature padding byte at offset {idx} should be rejected"
            );
        } else {
            let _ = HawkSignature::from_bytes(&bytes);
        }
    }

    /// Roundtrip: encoding a freshly decoded valid pubkey must yield the same
    /// bytes. Ensures `from_bytes ∘ to_bytes` is idempotent on valid inputs,
    /// ruling out lossy or nondeterministic decode paths.
    #[test]
    fn pubkey_roundtrip_idempotent(key_seed in prop::array::uniform32(any::<u8>())) {
        let mut kp_rng = ChaCha20Rng::from_seed(key_seed);
        let pk1 = HawkKeypair::generate(&mut kp_rng).public;
        let bytes1 = pk1.to_bytes().expect("encode 1");
        let pk2 = HawkPublicKey::from_bytes(&bytes1).expect("decode valid pk");
        let bytes2 = pk2.to_bytes().expect("encode 2");
        prop_assert_eq!(bytes1, bytes2);
    }

    /// Pathological all-zeros pubkey: decoder must reject without panicking.
    /// Under the Golomb-Rice wire format an all-zero buffer decodes to an
    /// all-zero polynomial for q00/q01, which may or may not pass the
    /// decoder's internal checks — but must never panic.
    #[test]
    fn all_zeros_pubkey_never_panics(_ in Just(())) {
        let bytes = [0u8; HAWK_PUBLIC_KEY_BYTES];
        let _ = HawkPublicKey::from_bytes(&bytes);
    }

    /// Pathological all-ones pubkey (every byte 0xFF): decoder must handle
    /// without panicking. 0xFF bytes encode maximal unary runs in the
    /// variable-length Golomb-Rice tail, which historically has been a rich
    /// source of integer-overflow bugs in C implementations.
    #[test]
    fn all_ones_pubkey_never_panics(_ in Just(())) {
        let bytes = [0xFFu8; HAWK_PUBLIC_KEY_BYTES];
        let _ = HawkPublicKey::from_bytes(&bytes);
    }

    /// Pathological all-ones signature: same rationale as the pubkey case.
    #[test]
    fn all_ones_signature_never_panics(_ in Just(())) {
        let bytes = [0xFFu8; HAWK_SIGNATURE_BYTES];
        let _ = HawkSignature::from_bytes(&bytes);
    }
}
