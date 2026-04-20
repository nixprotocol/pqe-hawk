//! Fuzz target: HawkSignature::from_bytes must never panic on any input.

#![no_main]

use libfuzzer_sys::fuzz_target;
use pqe_hawk::HawkSignature;
use pqe_hawk::params::HAWK_SIGNATURE_BYTES;

fuzz_target!(|data: &[u8]| {
    if data.len() >= HAWK_SIGNATURE_BYTES {
        let mut buf = [0u8; HAWK_SIGNATURE_BYTES];
        buf.copy_from_slice(&data[..HAWK_SIGNATURE_BYTES]);
        let _ = HawkSignature::from_bytes(&buf);
    } else if !data.is_empty() {
        let mut buf = [0u8; HAWK_SIGNATURE_BYTES];
        buf[..data.len()].copy_from_slice(data);
        let _ = HawkSignature::from_bytes(&buf);
    }
});
