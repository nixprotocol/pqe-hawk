//! Fuzz target: HawkPublicKey::from_bytes must never panic on any input.

#![no_main]

use libfuzzer_sys::fuzz_target;
use pqe_hawk::HawkPublicKey;
use pqe_hawk::params::HAWK_PUBLIC_KEY_BYTES;

fuzz_target!(|data: &[u8]| {
    // Accept any length — if the slice is exactly the expected size, try
    // decoding; otherwise ignore. Also test truncated / extended buffers as
    // arrays of exactly the right size but filled from arbitrary prefixes.
    if data.len() >= HAWK_PUBLIC_KEY_BYTES {
        let mut buf = [0u8; HAWK_PUBLIC_KEY_BYTES];
        buf.copy_from_slice(&data[..HAWK_PUBLIC_KEY_BYTES]);
        let _ = HawkPublicKey::from_bytes(&buf);
    }
    // Also try zero-padded short input — decode_public must still return Err
    // cleanly for any malformed fixed-size buffer.
    else if !data.is_empty() {
        let mut buf = [0u8; HAWK_PUBLIC_KEY_BYTES];
        buf[..data.len()].copy_from_slice(data);
        let _ = HawkPublicKey::from_bytes(&buf);
    }
});
