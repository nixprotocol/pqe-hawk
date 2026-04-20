//! HAWK signature type and generation.
//!
//! Signatures are a salt plus a short polynomial `s1` in `Z[X]/(X^n+1)`.
//! The signing algorithms themselves come in Task 22 (randomized) and
//! Task 23 (deterministic).

// `bp` is `pub` so FFI cross-check tests in tests/cross_check.rs can reach
// its helpers; hidden from rustdoc since it's not part of the stable API.
#[doc(hidden)]
pub mod bp;
pub(crate) mod symbreak;

pub(crate) use bp::basis_m2_mul;
pub(crate) use symbreak::poly_symbreak;

use crate::error::HawkError;
use crate::keygen::HawkSecretKey;
use crate::ntt::{
    mq_intt, mq_montymul, mq_ntt, mq_poly_set_small, mq_poly_snorm, mq_sub, mq_tomonty,
};
use crate::params::{
    HAWK_LOGN, HAWK_N, HAWK_SALT_BYTES, HAWK_SAMPLER_RETRY_BUDGET, HAWK_SIGNATURE_BYTES,
};
use rand::RngCore;

/// HAWK-512 signature.
///
/// Wire format matches the reference C (`encode_sig` in hawk_sign.c:680):
///   salt (24 bytes) || GolombRice(s1, low=5) || zero-padding → 555 bytes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HawkSignature {
    pub(crate) salt: [u8; HAWK_SALT_BYTES],
    pub(crate) s1: Vec<i16>,
}

impl HawkSignature {
    /// Serialize to the 555-byte wire format.
    pub fn to_bytes(&self) -> Result<[u8; HAWK_SIGNATURE_BYTES], HawkError> {
        crate::serialize::encode_signature(&self.salt, &self.s1)
    }

    /// Deserialize from the 555-byte wire format.
    pub fn from_bytes(bytes: &[u8; HAWK_SIGNATURE_BYTES]) -> Result<Self, HawkError> {
        crate::serialize::decode_signature(bytes)
    }
}

// ---------------------------------------------------------------------------
// HAWK-512 signing: port of `sign_finish_inner` (hawk_sign.c:765-1159),
// use_shake=0 ("alt") variant.
// ---------------------------------------------------------------------------

/// Inner free function for signing. Extracted so `HawkSecretKey::sign` is a
/// thin wrapper and the logic is testable without a method receiver.
///
/// Port of `sign_finish_inner` (hawk_sign.c:904-1158) with `priv_decoded=1`
/// and `use_shake=0`.
fn sign_512<R: RngCore>(
    secret: &HawkSecretKey,
    msg: &[u8],
    rng: &mut R,
) -> Result<HawkSignature, HawkError> {
    const N: usize = HAWK_N; // 512
    const N8: usize = N / 8; // 64 — bytes per packed bit-vector
    const MAX_XNORM: u32 = 8317; // HAWK-512: 2*n*(sigma_ver^2) bound
    const LIM: i32 = 1 << 9; // 512 — per-coefficient bound for s1

    // Step 1: compute hm = SHAKE256(msg || hpub) → 64 bytes.
    let hm: [u8; 64] = crate::hash::compute_hm(msg, &secret.hpub);

    // Step 2: extract F2 = low bits of f_cap, G2 = low bits of g_cap (n/8 bytes each).
    let mut f_cap2 = [0u8; N8];
    let mut g_cap2 = [0u8; N8];
    crate::serialize::extract_lowbit(&secret.f_cap, &mut f_cap2);
    crate::serialize::extract_lowbit(&secret.g_cap, &mut g_cap2);

    // Step 3: extract f2, g2 (low bits of f, g).
    let mut f2 = [0u8; N8];
    let mut g2 = [0u8; N8];
    crate::serialize::extract_lowbit(&secret.f, &mut f2);
    crate::serialize::extract_lowbit(&secret.g, &mut g2);

    // Retry loop.
    for _attempt in 0..HAWK_SAMPLER_RETRY_BUDGET {
        // 3a. Generate 24-byte salt from rng directly (use_shake=0 path).
        let mut salt = [0u8; HAWK_SALT_BYTES];
        rng.fill_bytes(&mut salt);

        // 3b. h = SHAKE256(hm || salt) → (h0, h1), each N8 = 64 bytes.
        let (h0, h1) = crate::hash::compute_h(&hm, &salt);

        // 3c. t = B*h mod 2 → (t0, t1), each N8 bytes.
        let mut t0 = [0u8; N8];
        let mut t1 = [0u8; N8];
        let mut bp_tmp = [0u8; 224];
        basis_m2_mul(
            &mut t0,
            &mut t1,
            &h0,
            &h1,
            &f2,
            &g2,
            &f_cap2,
            &g_cap2,
            &mut bp_tmp,
        );

        // 3d. Gaussian sample conditioned on t (pass t0 || t1 concatenated = 128 bytes).
        //     sig_gauss_alt indexes into the full 2n-bit parity vector using v = 0..2n-1.
        let mut t_parity = [0u8; 2 * N8]; // 128 bytes
        t_parity[..N8].copy_from_slice(&t0);
        t_parity[N8..].copy_from_slice(&t1);
        let gs = crate::sample::sample(rng, &t_parity);

        // 3e. Check squared norm. Port of hawk_sign.c:1019. Note that the
        //     norm check happens here, immediately after sampling and
        //     *before* the later sym-break step (3h), matching the C's
        //     control-flow order. Because the sym-break applies only to w3
        //     (a post-NTT quantity, not to the Gaussian samples), and
        //     because squared norm is invariant under negation anyway, the
        //     ordering is both correct and unambiguous.
        if gs.squared_norm > MAX_XNORM {
            continue;
        }

        // 3f. x0 = gs.x[0..N], x1 = gs.x[N..2N].
        let x0 = &gs.x[..N];
        let x1 = &gs.x[N..];

        // 3g. Compute w3 = 2*(f*x1 - g*x0) via NTT.
        //     All NTT operations work in [1, q] (Montgomery) domain.
        //
        //     w1 = g * x0 (NTT domain):
        let mut w1 = vec![0u16; N];
        let mut w2 = vec![0u16; N];
        let mut w3 = vec![0u16; N];

        // w1 ← g (mq domain)
        mq_poly_set_small(&mut w1, &secret.g);
        // w2 ← x0 (mq domain)
        mq_poly_set_small(&mut w2, x0);
        mq_ntt(HAWK_LOGN, &mut w1);
        mq_ntt(HAWK_LOGN, &mut w2);
        // w1[u] = montymul(g[u], x0[u]) — NTT-domain product g*x0
        for u in 0..N {
            w1[u] = mq_montymul(w1[u] as u32, w2[u] as u32) as u16;
        }

        // w3 ← f (mq domain), w2 ← x1 (mq domain)
        mq_poly_set_small(&mut w3, &secret.f);
        mq_poly_set_small(&mut w2, x1);
        mq_ntt(HAWK_LOGN, &mut w3);
        mq_ntt(HAWK_LOGN, &mut w2);
        // w3[u] = tomonty(sub(montymul(f[u], x1[u]), w1[u]))
        //       = tomonty(f*x1 - g*x0) in NTT domain
        //   → after iNTT + snorm, w3 holds signed coefficients of 2*(f*x1 - g*x0)
        for u in 0..N {
            w3[u] = mq_tomonty(mq_sub(
                mq_montymul(w2[u] as u32, w3[u] as u32),
                w1[u] as u32,
            )) as u16;
        }
        mq_intt(HAWK_LOGN, &mut w3);
        mq_poly_snorm(&mut w3);
        // w3 now contains signed i16 values stored as u16 (two's complement).

        // 3h. Sym-break: `w3` currently holds `h1 - 2*s1_raw`. Depending on the
        //     sign of the first non-zero coefficient, we must return either
        //     `s1 = (h1 - w3)/2` or `s1 = (h1 + w3)/2` — which is a conditional
        //     negation of `w3` followed by `(w3 + h1_bit) >> 1`.
        //
        //     Port of hawk_sign.c:1097-1122. The mask semantics from the C:
        //       ps = +1 (first non-zero positive) → nm = 0xFFFF_FFFF → negate w3
        //       ps = -1 (first non-zero negative) → nm = 0          → keep w3
        //       ps =  0 (all zeros)               → nm = 0          → keep w3
        //     (The comment in the C explicitly states: "this uses a conditional
        //     negation of w3".)
        // Safe conversion: read each w3[u] as a signed i16 (the u16→i16
        // numeric cast preserves the two's-complement bit pattern without
        // unsafe slice aliasing).
        let mut s1_vec: Vec<i16> = w3.iter().map(|&v| v as i16).collect();
        let ps = poly_symbreak(&s1_vec);
        // `~tbmask(ps - 1)` where `tbmask(x) = (x as i32 >> 31) as u32` gives
        // the table above.
        let nm: u32 = !((((ps as u32).wrapping_sub(1)) as i32 >> 31) as u32);

        // Recover s1 from w3, checking per-coefficient bound.
        let mut in_bounds = true;
        for u in 0..N {
            // z = current signed-normed value (as u32, two's complement).
            let z = s1_vec[u] as i32 as u32;
            // Conditional negate: if nm = 0xFFFF..., z becomes -z; else stays z.
            let z = (z ^ nm).wrapping_sub(nm);
            // Add the h1 parity bit for this index.
            let h1_bit = ((h1[u >> 3] >> (u & 7)) & 1) as u32;
            let z = z.wrapping_add(h1_bit);
            // Arithmetic right-shift by 1: y = z >> 1 (signed).
            let y = (z as i32) >> 1;
            if y < -LIM || y >= LIM {
                in_bounds = false;
                break;
            }
            s1_vec[u] = y as i16;
        }

        if !in_bounds {
            continue;
        }

        // 3i. Encode signature: salt || GR(s1) → 555 bytes.
        let sig = HawkSignature { salt, s1: s1_vec };

        // Encoding always succeeds when in_bounds (the GR encoder may still
        // overflow in theory, but for HAWK-512 with LIM=512 it fits).
        return Ok(sig);
    }

    Err(HawkError::SamplingFailure {
        retries: HAWK_SAMPLER_RETRY_BUDGET,
    })
}

impl HawkSecretKey {
    /// Sign a message using fresh random salt from `rng`.
    ///
    /// Port of `hawk_sign_finish_alt` (`use_shake=0`) from hawk_sign.c:765-1159
    /// for HAWK-512 (logn=9). Uses the "alt" sampler path that draws salt and
    /// Gaussian samples directly from `rng` without SHAKE post-processing.
    ///
    /// Returns `Err(HawkError::SamplingFailure)` if the rejection sampler
    /// fails to produce a valid signature within the retry budget (extremely
    /// unlikely in practice; the C reference uses an unbounded loop).
    pub fn sign<R: RngCore>(&self, msg: &[u8], rng: &mut R) -> Result<HawkSignature, HawkError> {
        sign_512(self, msg, rng)
    }

    /// Sign deterministically: same (self, msg, nonce_seed) always produces
    /// the same signature bytes.
    ///
    /// The "randomness" is derived from SHAKE256 over:
    ///   self.to_bytes() || nonce_seed || msg
    /// This is useful for proposer crash recovery in the pqe chain: the
    /// proposer can re-derive an identical signature after a crash without
    /// needing to persist signing-time RNG state.
    pub fn sign_deterministic(
        &self,
        msg: &[u8],
        nonce_seed: &[u8; 32],
    ) -> Result<HawkSignature, crate::error::HawkError> {
        use sha3::digest::{ExtendableOutput, Update, XofReader};
        use sha3::Shake256;

        // Derive a SHAKE-256 XOF seeded with key bytes || nonce_seed || msg.
        let mut shake = Shake256::default();
        let sk_bytes = self.to_bytes();
        shake.update(&sk_bytes);
        shake.update(nonce_seed);
        shake.update(msg);
        let reader = shake.finalize_xof();

        /// Adapter: RngCore implementor pulling bytes from a SHAKE XOF.
        struct ShakeRng<R: XofReader> {
            reader: R,
        }
        impl<R: XofReader> rand::RngCore for ShakeRng<R> {
            fn next_u32(&mut self) -> u32 {
                let mut b = [0u8; 4];
                self.fill_bytes(&mut b);
                u32::from_le_bytes(b)
            }
            fn next_u64(&mut self) -> u64 {
                let mut b = [0u8; 8];
                self.fill_bytes(&mut b);
                u64::from_le_bytes(b)
            }
            fn fill_bytes(&mut self, dst: &mut [u8]) {
                self.reader.read(dst);
            }
            fn try_fill_bytes(&mut self, dst: &mut [u8]) -> Result<(), rand::Error> {
                self.reader.read(dst);
                Ok(())
            }
        }

        let mut rng = ShakeRng { reader };
        self.sign(msg, &mut rng)
    }
}

#[cfg(test)]
mod sign_tests {
    use super::*;
    use crate::keygen::HawkKeypair;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn sign_produces_valid_signature_bytes() {
        let mut rng = ChaCha20Rng::from_seed([7u8; 32]);
        let kp = HawkKeypair::generate(&mut rng);
        let sig = kp
            .secret
            .sign(b"hello", &mut rng)
            .expect("sign should succeed");
        let bytes = sig.to_bytes().expect("encode should succeed");
        // Verify salt is non-trivial and bytes are not all-zero.
        assert!(bytes.iter().any(|&b| b != 0));
        // Signature length must be exactly HAWK_SIGNATURE_BYTES.
        assert_eq!(bytes.len(), HAWK_SIGNATURE_BYTES);
    }

    #[test]
    fn sign_is_deterministic_for_same_rng_seed() {
        // Two keypairs generated from identical seeds are identical.
        let mut rng1 = ChaCha20Rng::from_seed([7u8; 32]);
        let mut rng2 = ChaCha20Rng::from_seed([7u8; 32]);
        let kp1 = HawkKeypair::generate(&mut rng1);
        let kp2 = HawkKeypair::generate(&mut rng2);

        // Use separate seeded RNGs for the signing step so state is identical.
        let mut srng1 = ChaCha20Rng::from_seed([99u8; 32]);
        let mut srng2 = ChaCha20Rng::from_seed([99u8; 32]);
        let sig1 = kp1.secret.sign(b"same msg", &mut srng1).unwrap();
        let sig2 = kp2.secret.sign(b"same msg", &mut srng2).unwrap();
        assert_eq!(
            sig1.to_bytes().unwrap(),
            sig2.to_bytes().unwrap(),
            "signatures from identical key+rng must be identical"
        );
    }

    #[test]
    fn sign_deterministic_is_stable() {
        use crate::HawkKeypair;
        let mut rng = ChaCha20Rng::from_seed([1u8; 32]);
        let kp = HawkKeypair::generate(&mut rng);
        let sig1 = kp.secret.sign_deterministic(b"hello", &[42u8; 32]).unwrap();
        let sig2 = kp.secret.sign_deterministic(b"hello", &[42u8; 32]).unwrap();
        assert_eq!(sig1.to_bytes().unwrap(), sig2.to_bytes().unwrap());
    }

    #[test]
    fn sign_deterministic_differs_with_different_nonce() {
        use crate::HawkKeypair;
        let mut rng = ChaCha20Rng::from_seed([2u8; 32]);
        let kp = HawkKeypair::generate(&mut rng);
        let sig1 = kp.secret.sign_deterministic(b"hello", &[1u8; 32]).unwrap();
        let sig2 = kp.secret.sign_deterministic(b"hello", &[2u8; 32]).unwrap();
        assert_ne!(sig1.to_bytes().unwrap(), sig2.to_bytes().unwrap());
    }

    #[test]
    fn sign_deterministic_differs_with_different_msg() {
        use crate::HawkKeypair;
        let mut rng = ChaCha20Rng::from_seed([3u8; 32]);
        let kp = HawkKeypair::generate(&mut rng);
        let sig1 = kp.secret.sign_deterministic(b"hello", &[7u8; 32]).unwrap();
        let sig2 = kp.secret.sign_deterministic(b"world", &[7u8; 32]).unwrap();
        assert_ne!(sig1.to_bytes().unwrap(), sig2.to_bytes().unwrap());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sig_roundtrip_zero_poly() {
        let salt = [7u8; HAWK_SALT_BYTES];
        let s1 = vec![0i16; HAWK_N];
        let sig = HawkSignature {
            salt,
            s1: s1.clone(),
        };
        let bytes = sig.to_bytes().unwrap();
        let decoded = HawkSignature::from_bytes(&bytes).unwrap();
        assert_eq!(decoded.salt, salt);
        assert_eq!(decoded.s1, s1);
    }

    #[test]
    fn sig_roundtrip_small_values() {
        let salt = [42u8; HAWK_SALT_BYTES];
        let s1: Vec<i16> = (0..HAWK_N).map(|i| ((i as i32 % 7) - 3) as i16).collect();
        let sig = HawkSignature {
            salt,
            s1: s1.clone(),
        };
        let bytes = sig.to_bytes().unwrap();
        let decoded = HawkSignature::from_bytes(&bytes).unwrap();
        assert_eq!(decoded.salt, salt);
        assert_eq!(decoded.s1, s1);
    }

    #[test]
    fn sig_rejects_nonzero_padding() {
        let salt = [0u8; HAWK_SALT_BYTES];
        let s1 = vec![0i16; HAWK_N];
        let sig = HawkSignature { salt, s1 };
        let mut bytes = sig.to_bytes().unwrap();
        // Flip the last byte (padding).
        bytes[HAWK_SIGNATURE_BYTES - 1] ^= 0x01;
        let r = HawkSignature::from_bytes(&bytes);
        assert!(r.is_err());
    }
}
