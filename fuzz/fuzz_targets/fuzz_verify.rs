//! Fuzz target: full verify pipeline must never panic on any input.
//!
//! Uses a fixed known-good keypair and message, fuzzes the signature bytes.
//! Any panic indicates a potential DoS / memory-safety vulnerability in the
//! verify path — verify must reject gracefully via `Err(...)` for any input.

#![no_main]

use libfuzzer_sys::fuzz_target;
use pqe_hawk::{HawkKeypair, HawkPublicKey, HawkSignature};
use pqe_hawk::params::{HAWK_PUBLIC_KEY_BYTES, HAWK_SIGNATURE_BYTES};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::sync::OnceLock;

static PINNED_PUBKEY: OnceLock<HawkPublicKey> = OnceLock::new();

fn pinned_pubkey() -> &'static HawkPublicKey {
    PINNED_PUBKEY.get_or_init(|| {
        let mut rng = ChaCha20Rng::from_seed([1u8; 32]);
        let kp = HawkKeypair::generate(&mut rng);
        kp.public
    })
}

fuzz_target!(|data: &[u8]| {
    // Split the fuzzer-provided bytes into (signature, message).
    // Use the first HAWK_SIGNATURE_BYTES (or fewer, zero-padded) as the sig;
    // the remainder is the message.
    if data.is_empty() {
        return;
    }

    let sig_len = HAWK_SIGNATURE_BYTES.min(data.len());
    let mut sig_buf = [0u8; HAWK_SIGNATURE_BYTES];
    sig_buf[..sig_len].copy_from_slice(&data[..sig_len]);

    let msg: &[u8] = if data.len() > sig_len { &data[sig_len..] } else { &[] };

    // Decode the signature — may fail (Err).
    let sig = match HawkSignature::from_bytes(&sig_buf) {
        Ok(s) => s,
        Err(_) => return,
    };

    // Verify against a fixed pinned pubkey. Must return Ok or Err — never panic.
    let _ = pinned_pubkey().verify(msg, &sig);
});

#[allow(dead_code)]
fn size_hint_constants() {
    // Compile-time reference to keep the params public surface exercised.
    let _ = HAWK_PUBLIC_KEY_BYTES;
}
