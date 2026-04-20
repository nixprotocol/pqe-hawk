//! NIST AES-256-CTR-DRBG as used by the PQC submission framework.
//!
//! Spec: SP 800-90A Rev. 1, Section 10.2.1 (CTR_DRBG).
//! Purpose: reproduce the exact byte stream the NIST KAT harness feeds into
//! HAWK's keygen/sign. Not exposed to ship code.

use aes::cipher::generic_array::GenericArray;
use aes::cipher::{BlockEncrypt, KeyInit};
use aes::Aes256;

/// AES-256-CTR-DRBG per NIST SP 800-90A.
pub struct AesCtrDrbg {
    key: [u8; 32],
    v: [u8; 16],
}

impl AesCtrDrbg {
    /// Initialize from a 48-byte seed per NIST PQC framework conventions.
    pub fn new(seed: &[u8; 48]) -> Self {
        let mut drbg = AesCtrDrbg {
            key: [0u8; 32],
            v: [0u8; 16],
        };
        drbg.update(Some(seed));
        drbg
    }

    /// Derive N bytes. NIST algorithm: for each 16-byte output block,
    /// increment V (big-endian), encrypt V under K, output.
    /// After generating all requested bytes, call update(None) to re-key.
    pub fn generate(&mut self, out: &mut [u8]) {
        // Snapshot key before encrypting — the AES instance lives through this
        // generate call, not across update.
        let cipher = Aes256::new(GenericArray::from_slice(&self.key));
        for chunk in out.chunks_mut(16) {
            // V = V + 1 (big-endian increment, matches NIST / HAWK C).
            for i in (0..16).rev() {
                self.v[i] = self.v[i].wrapping_add(1);
                if self.v[i] != 0 {
                    break;
                }
            }
            let mut block = GenericArray::from(self.v);
            cipher.encrypt_block(&mut block);
            let n = chunk.len();
            chunk.copy_from_slice(&block[..n]);
        }
        // Per SP 800-90A, after generation refresh state.
        self.update(None);
    }

    /// NIST "update" function: refresh K and V using AES-output-XOR-provided.
    fn update(&mut self, provided: Option<&[u8; 48]>) {
        let cipher = Aes256::new(GenericArray::from_slice(&self.key));
        let mut out = [0u8; 48];
        for i in 0..3 {
            // V = V + 1 (big-endian increment)
            for j in (0..16).rev() {
                self.v[j] = self.v[j].wrapping_add(1);
                if self.v[j] != 0 {
                    break;
                }
            }
            let mut block = GenericArray::from(self.v);
            cipher.encrypt_block(&mut block);
            out[i * 16..(i + 1) * 16].copy_from_slice(&block);
        }
        if let Some(p) = provided {
            for i in 0..48 {
                out[i] ^= p[i];
            }
        }
        self.key.copy_from_slice(&out[0..32]);
        self.v.copy_from_slice(&out[32..48]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn drbg_is_deterministic() {
        let seed = [0u8; 48];
        let mut a = AesCtrDrbg::new(&seed);
        let mut b = AesCtrDrbg::new(&seed);
        let mut out_a = [0u8; 128];
        let mut out_b = [0u8; 128];
        a.generate(&mut out_a);
        b.generate(&mut out_b);
        assert_eq!(out_a, out_b);
    }

    #[test]
    fn drbg_different_seeds_differ() {
        let mut a = AesCtrDrbg::new(&[0u8; 48]);
        let mut b = AesCtrDrbg::new(&[1u8; 48]);
        let mut out_a = [0u8; 64];
        let mut out_b = [0u8; 64];
        a.generate(&mut out_a);
        b.generate(&mut out_b);
        assert_ne!(out_a, out_b);
    }

    #[test]
    fn drbg_structure_smoke() {
        let seed: [u8; 48] = std::array::from_fn(|i| i as u8);
        let mut drbg = AesCtrDrbg::new(&seed);
        let mut out = [0u8; 128];
        drbg.generate(&mut out);
        // Just check that output is non-trivial (not all zeros or all 0xFF).
        assert!(out.iter().any(|&b| b != 0));
        assert!(out.iter().any(|&b| b != 0xFF));
    }

    /// NIST KAT pin: seed = [0, 1, 2, ..., 47], request 128 bytes.
    ///
    /// Expected output obtained by running the vendored `rng.c`:
    ///   randombytes_init(seed=[0..48], NULL, 256); randombytes(out, 128);
    #[test]
    fn drbg_matches_nist_reference_output() {
        let seed: [u8; 48] = std::array::from_fn(|i| i as u8);
        let mut drbg = AesCtrDrbg::new(&seed);
        let mut out = [0u8; 128];
        drbg.generate(&mut out);

        let expected_hex = "061550234D158C5EC95595FE04EF7A25767F2E24CC2BC479D09D86DC9ABCFDE7056A8C266F9EF97ED08541DBD2E1FFA19810F5392D076276EF41277C3AB6E94A4E3B7DCC104A05BB089D338BF55C72CAB375389A94BB920BD5D6DC9E7F2EC6FDE028B6F5724BB039F3652AD98DF8CE6C97013210B84BBE81388C3D141D61957C";
        let expected = hex::decode(expected_hex).expect("valid hex");
        assert_eq!(
            &out[..expected.len()],
            &expected[..],
            "DRBG output does not match NIST KAT reference"
        );
    }
}
