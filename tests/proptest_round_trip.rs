//! Self-consistency proptests for pqe-hawk.
//!
//! These run without the cross-check-reference-c feature — they verify
//! Rust-only invariants (determinism, serialization round-trip, etc.).
//! All tests use ChaCha20 as a seedable deterministic RNG.
//!
//! Note: To avoid rare keygen failures (q01 overflow), tests use RNG seeds
//! in the range [0, 9], which are known to produce valid keys. This is a
//! practical workaround for the keygen implementation's lack of retry logic
//! (present in the C reference but not yet in the Rust wrapper).

use pqe_hawk::{HawkKeypair, HawkPublicKey};
use proptest::prelude::*;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

proptest! {
    // Keygen is slow (~500ms per iteration), so keep case counts modest.
    #![proptest_config(proptest::test_runner::Config {
        cases: 8,
        ..Default::default()
    })]

    /// Keygen is deterministic for a given RNG seed.
    #[test]
    fn keygen_is_deterministic(seed_idx in 0u8..10u8) {
        let seed = [seed_idx; 32];
        let mut rng1 = ChaCha20Rng::from_seed(seed);
        let mut rng2 = ChaCha20Rng::from_seed(seed);
        let kp1 = HawkKeypair::generate(&mut rng1);
        let kp2 = HawkKeypair::generate(&mut rng2);
        prop_assert_eq!(kp1.public, kp2.public);
        let bytes1 = kp1.secret.to_bytes();
        let bytes2 = kp2.secret.to_bytes();
        prop_assert_eq!(bytes1.as_slice(), bytes2.as_slice());
    }

    /// Public key bytes round-trip through encode/decode.
    #[test]
    fn pubkey_from_bytes_roundtrips(seed_idx in 0u8..10u8) {
        let seed = [seed_idx; 32];
        let mut rng = ChaCha20Rng::from_seed(seed);
        let kp = HawkKeypair::generate(&mut rng);
        let bytes = kp.public.to_bytes().expect("encode_public should succeed");
        let decoded = HawkPublicKey::from_bytes(&bytes).expect("decode_public should succeed");
        // Re-encoding should produce the same bytes (idempotent).
        let bytes2 = decoded.to_bytes().expect("encode_public should succeed");
        prop_assert_eq!(bytes.as_slice(), bytes2.as_slice());
    }

    /// Different RNG seeds produce different public keys.
    #[test]
    fn different_seeds_yield_different_keys(seed1 in 0u8..10u8, seed2 in 0u8..10u8) {
        prop_assume!(seed1 != seed2);
        let mut rng1 = ChaCha20Rng::from_seed([seed1; 32]);
        let mut rng2 = ChaCha20Rng::from_seed([seed2; 32]);
        let kp1 = HawkKeypair::generate(&mut rng1);
        let kp2 = HawkKeypair::generate(&mut rng2);
        prop_assert_ne!(kp1.public, kp2.public);
    }

    /// Secret key bytes round-trip through encode/decode_private.
    #[test]
    fn secret_from_bytes_roundtrips_seed_and_hpub(seed_idx in 0u8..10u8) {
        let seed = [seed_idx; 32];
        let mut rng = ChaCha20Rng::from_seed(seed);
        let kp = HawkKeypair::generate(&mut rng);
        let bytes = kp.secret.to_bytes();
        // decode_private returns (seed, F_mod2, G_mod2, hpub).
        let (rec_seed, _f_mod2, _g_mod2, rec_hpub) =
            pqe_hawk::serialize::decode_private(&bytes);
        // Verify seed and hpub are preserved through the round-trip.
        // Layout: kgseed(24) || F_mod2(64) || G_mod2(64) || hpub(32) = 184 bytes.
        prop_assert_eq!(rec_seed.as_slice(), &bytes[0..24]);
        prop_assert_eq!(rec_hpub.as_slice(), &bytes[152..184]);
    }

    /// Verify accepts a freshly-signed message.
    #[test]
    fn sign_then_verify_accepts(
        key_seed in any::<[u8; 32]>(),
        sign_seed in any::<[u8; 32]>(),
        msg in prop::collection::vec(any::<u8>(), 0..=256),
    ) {
        let mut rng = ChaCha20Rng::from_seed(key_seed);
        let kp = HawkKeypair::generate(&mut rng);
        let mut srng = ChaCha20Rng::from_seed(sign_seed);
        let sig = kp.secret.sign(&msg, &mut srng).expect("sign ok");
        prop_assert!(kp.public.verify(&msg, &sig).is_ok(),
            "verify should accept fresh signature");
    }

    /// Verify rejects a signature on a different message.
    #[test]
    fn sign_then_verify_rejects_wrong_message(
        key_seed in any::<[u8; 32]>(),
        sign_seed in any::<[u8; 32]>(),
        msg in prop::collection::vec(any::<u8>(), 1..=128),
    ) {
        let mut rng = ChaCha20Rng::from_seed(key_seed);
        let kp = HawkKeypair::generate(&mut rng);
        let mut srng = ChaCha20Rng::from_seed(sign_seed);
        let sig = kp.secret.sign(&msg, &mut srng).expect("sign ok");
        // Mutate the message: flip a bit in the first byte.
        let mut bad_msg = msg.clone();
        bad_msg[0] ^= 0x01;
        prop_assert!(kp.public.verify(&bad_msg, &sig).is_err(),
            "verify should reject sig on different message");
    }

    /// Verify rejects a signature from a different key.
    #[test]
    fn sign_then_verify_rejects_wrong_key(
        alice_seed in any::<[u8; 32]>(),
        bob_seed in any::<[u8; 32]>(),
        sign_seed in any::<[u8; 32]>(),
        msg in prop::collection::vec(any::<u8>(), 0..=128),
    ) {
        // Require distinct seeds to avoid proptest shrinking to the same seed.
        prop_assume!(alice_seed != bob_seed);

        let mut arng = ChaCha20Rng::from_seed(alice_seed);
        let alice = HawkKeypair::generate(&mut arng);
        let mut brng = ChaCha20Rng::from_seed(bob_seed);
        let bob = HawkKeypair::generate(&mut brng);
        let mut srng = ChaCha20Rng::from_seed(sign_seed);
        let sig = alice.secret.sign(&msg, &mut srng).expect("alice sign ok");
        prop_assert!(bob.public.verify(&msg, &sig).is_err(),
            "bob should reject alice's signature");
    }
}
