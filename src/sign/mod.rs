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
    // Production path: no-op hook (see `sign_512_hooked`).
    sign_512_hooked(secret, msg, rng, &mut |_attempt, _s1| {})
}

/// `sign_512` with a test seam.
///
/// `s1_hook` is invoked once per retry-loop iteration, immediately before the
/// encode step, with the attempt index and a mutable `s1`. Production passes a
/// no-op, so this is byte-for-byte identical to the un-hooked path (pinned by
/// `hooked_signing_with_noop_equals_production`). Tests use it to force the
/// rare encode-overflow path deterministically and verify the loop recovers by
/// resigning.
fn sign_512_hooked<R: RngCore>(
    secret: &HawkSecretKey,
    msg: &[u8],
    rng: &mut R,
    s1_hook: &mut dyn FnMut(usize, &mut Vec<i16>),
) -> Result<HawkSignature, HawkError> {
    const N: usize = HAWK_N; // 512
    const N8: usize = N / 8; // 64 — bytes per packed bit-vector
    const MAX_XNORM: u32 = 8317; // HAWK-512: 2*n*(sigma_ver^2) bound
    const LIM: i32 = 1 << 9; // 512 — per-coefficient bound for s1

    // Step 1: compute hm = SHAKE256(msg || hpub) → 64 bytes.
    //
    // Note: hm is NOT wrapped in `Zeroizing`. It is derived entirely from
    // public inputs — `msg` is chosen by the caller, and `hpub` is
    // `SHAKE256(encoded_public_key)[..32]` which is a function of the
    // public key alone. Zeroizing it would be defense-theatre, not
    // defense-in-depth. The same reasoning applies to `salt` (public by
    // construction; part of the output signature) and `h0, h1`
    // (SHAKE256(hm || salt) — both inputs public).
    let hm: [u8; 64] = crate::hash::compute_hm(msg, &secret.hpub);

    // Step 2: extract F2 = low bits of f_cap, G2 = low bits of g_cap (n/8
    // bytes each). These 64-byte masks are the `mod 2` fingerprint of the
    // secret key's capital-F, capital-G polynomials; treat them as secret and
    // wrap in Zeroizing so the stack buffers are cleared on scope exit.
    let mut f_cap2 = zeroize::Zeroizing::new([0u8; N8]);
    let mut g_cap2 = zeroize::Zeroizing::new([0u8; N8]);
    crate::serialize::extract_lowbit(&secret.f_cap, &mut *f_cap2);
    crate::serialize::extract_lowbit(&secret.g_cap, &mut *g_cap2);

    // Step 3: extract f2, g2 (low bits of f, g). Same secret-fingerprint
    // treatment — zeroize on drop.
    let mut f2 = zeroize::Zeroizing::new([0u8; N8]);
    let mut g2 = zeroize::Zeroizing::new([0u8; N8]);
    crate::serialize::extract_lowbit(&secret.f, &mut *f2);
    crate::serialize::extract_lowbit(&secret.g, &mut *g2);

    // Retry loop.
    for attempt in 0..HAWK_SAMPLER_RETRY_BUDGET {
        // 3a. Generate 24-byte salt from rng directly (use_shake=0 path).
        let mut salt = [0u8; HAWK_SALT_BYTES];
        rng.fill_bytes(&mut salt);

        // 3b. h = SHAKE256(hm || salt) → (h0, h1), each N8 = 64 bytes.
        let (h0, h1) = crate::hash::compute_h(&hm, &salt);

        // 3c. t = B*h mod 2 → (t0, t1), each N8 bytes. (t0, t1) and the
        // scratch buffer `bp_tmp` hold secret-derived values (linear
        // combinations of the secret mod-2 basis with the public h vector),
        // so zeroize on retry-loop-iteration scope exit.
        let mut t0 = zeroize::Zeroizing::new([0u8; N8]);
        let mut t1 = zeroize::Zeroizing::new([0u8; N8]);
        let mut bp_tmp = zeroize::Zeroizing::new([0u8; 224]);
        basis_m2_mul(
            &mut *t0,
            &mut *t1,
            &h0,
            &h1,
            &*f2,
            &*g2,
            &*f_cap2,
            &*g_cap2,
            &mut *bp_tmp,
        );

        // 3d. Gaussian sample conditioned on t (pass t0 || t1 concatenated = 128 bytes).
        //     sig_gauss_alt indexes into the full 2n-bit parity vector using v = 0..2n-1.
        //     t_parity is a re-packing of the secret-derived (t0, t1);
        //     zeroize on scope exit.
        let mut t_parity = zeroize::Zeroizing::new([0u8; 2 * N8]); // 128 bytes
        t_parity[..N8].copy_from_slice(&*t0);
        t_parity[N8..].copy_from_slice(&*t1);
        let gs = crate::sample::sample(rng, &*t_parity);

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
        //     w1, w2, w3 hold NTT-domain intermediates derived from secret
        //     (f, g). `Zeroizing` wraps them so the backing `Vec<u16>` is
        //     zeroed when it goes out of scope (including via early `continue`
        //     in the retry loop below), matching the treatment of the secret
        //     key itself (`ZeroizeOnDrop`).
        //
        //     w1 = g * x0 (NTT domain):
        let mut w1 = zeroize::Zeroizing::new(vec![0u16; N]);
        let mut w2 = zeroize::Zeroizing::new(vec![0u16; N]);
        let mut w3 = zeroize::Zeroizing::new(vec![0u16; N]);

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
        //
        // Wrap in `Zeroizing` so a failed-bounds `continue` path below still
        // zeroes the intermediate (it holds secret-derived values until the
        // final normalization loop completes). The success path moves out of
        // the wrapper into the public signature struct.
        let mut s1_vec =
            zeroize::Zeroizing::new(w3.iter().map(|&v| v as i16).collect::<Vec<i16>>());
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
            if !(-LIM..LIM).contains(&y) {
                in_bounds = false;
                break;
            }
            s1_vec[u] = y as i16;
        }

        if !in_bounds {
            continue;
        }

        // 3i. Encode signature: salt || GR(s1) → 555 bytes. Take the inner
        // Vec out of the Zeroizing wrapper (the signature itself is public
        // output — no need to zeroize it). Replacing with an empty Vec
        // leaves a zero-length buffer for Zeroizing to handle on drop.
        //
        // The per-coefficient bound (3h) does NOT guarantee the Golomb-Rice
        // encoding fits the fixed 555-byte buffer: its variable section grows
        // with sum(|s1[i]| >> 5), which can overflow even when every
        // coefficient is in-bounds. If it overflows, resign with a fresh salt
        // rather than return a signature that cannot be serialized. This
        // matches the reference C, which restarts on encode failure
        // (hawk_sign.c:1142-1158).
        let mut s1 = std::mem::take(&mut *s1_vec);
        // Test seam: production passes a no-op (see `sign_512`). Tests use it to
        // force `s1` to an overflowing value on a chosen attempt, exercising the
        // resign path below.
        s1_hook(attempt, &mut s1);
        match try_encode_attempt(salt, s1) {
            Some(sig) => return Ok(sig),
            None => continue,
        }
    }

    Err(HawkError::SamplingFailure {
        retries: HAWK_SAMPLER_RETRY_BUDGET,
    })
}

/// Build a `HawkSignature` from `(salt, s1)` only if it serializes into the
/// fixed 555-byte wire format. Returns `None` on Golomb-Rice buffer overflow,
/// signalling the signing loop to resign with a fresh salt (port of the
/// `encode_sig` failure path in hawk_sign.c:1142-1158).
///
/// Separated from `sign_512` so the encode-or-resign decision is unit-testable
/// with crafted `s1` values without driving the full sampler.
fn try_encode_attempt(salt: [u8; HAWK_SALT_BYTES], s1: Vec<i16>) -> Option<HawkSignature> {
    let sig = HawkSignature { salt, s1 };
    // `to_bytes` overflow is the only expected error here; treat any encode
    // failure as "does not fit, resign".
    sig.to_bytes().ok().map(|_| sig)
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

    // --- resign-on-encode-overflow (regression guard) ---
    //
    // The Golomb-Rice signature encoder writes a fixed 555-byte buffer. Its
    // variable section grows with sum(|s1[i]| >> 5), so a signature can be
    // in-bounds per-coefficient (|s1[i]| < LIM = 512) yet still overflow in
    // aggregate. `sign_512` must NOT return such a signature; it must resign.
    // The encode-or-resign decision is `try_encode_attempt`.

    /// An s1 that is in-bounds per-coefficient but overflows the buffer in
    /// aggregate. From the empirical boundary, sum(|c|>>5) > ~664 overflows;
    /// 50 coefficients at 511 give sum = 50*15 = 750, comfortably over. Every
    /// coefficient is < LIM (512), so the signing bounds check would pass it.
    fn overflowing_s1() -> Vec<i16> {
        let mut s1 = vec![0i16; HAWK_N];
        for s in s1.iter_mut().take(50) {
            *s = 511;
        }
        s1
    }

    #[test]
    fn try_encode_attempt_returns_none_on_aggregate_overflow() {
        let salt = [0u8; HAWK_SALT_BYTES];
        let s1 = overflowing_s1();
        // Every coefficient is within the per-coefficient signing bound...
        assert!(s1.iter().all(|&c| (-512..512).contains(&(c as i32))));
        // ...yet the encoding overflows, so the attempt must be rejected.
        assert!(
            try_encode_attempt(salt, s1).is_none(),
            "an s1 that overflows the 555-byte buffer must not yield a signature"
        );
    }

    #[test]
    fn try_encode_attempt_returns_some_for_encodable_s1() {
        let salt = [9u8; HAWK_SALT_BYTES];
        let s1: Vec<i16> = (0..HAWK_N).map(|i| ((i as i32 % 7) - 3) as i16).collect();
        let sig = try_encode_attempt(salt, s1.clone())
            .expect("an in-budget s1 must encode to a signature");
        assert_eq!(sig.salt, salt);
        assert_eq!(sig.s1, s1);
        // And the returned signature is guaranteed serializable.
        assert!(sig.to_bytes().is_ok());
    }

    #[test]
    fn signed_signatures_always_serialize() {
        use crate::keygen::HawkKeypair;
        use rand::SeedableRng;
        use rand_chacha::ChaCha20Rng;

        // The invariant the bug violated: every signature `sign` returns must
        // serialize. Drive many signings; assert each output round-trips.
        // NOTE: overflow is ~1-in-thousands, so these signings almost certainly
        // never hit the resign path. This guards the invariant but does NOT
        // exercise the retry; `sign_recovers_from_forced_encode_overflow` does.
        let mut rng = ChaCha20Rng::from_seed([5u8; 32]);
        let kp = HawkKeypair::generate(&mut rng);
        for i in 0..200u32 {
            let msg = format!("msg-{i}");
            let sig = kp
                .secret
                .sign(msg.as_bytes(), &mut rng)
                .expect("signing within budget");
            assert!(
                sig.to_bytes().is_ok(),
                "sign returned a signature that fails to serialize (msg {i})"
            );
        }
    }

    // --- integration: the loop actually RECOVERS from an overflow ---
    //
    // The tests above prove `try_encode_attempt` rejects a bad s1, but not that
    // `sign_512` responds to that rejection by looping and ultimately returning
    // a VALID signature. Overflow is too rare to hit by chance, so we use the
    // `s1_hook` seam to force the first attempt's s1 to overflow, then assert
    // the loop resigns and yields a serializable signature.

    /// Drive `sign_512_hooked` while forcing the first attempt's s1 to a value
    /// that overflows the wire buffer, so the loop is made to traverse the
    /// overflow -> resign -> success path deterministically.
    #[test]
    fn sign_recovers_from_forced_encode_overflow() {
        use crate::keygen::HawkKeypair;
        use rand::SeedableRng;
        use rand_chacha::ChaCha20Rng;
        use std::cell::Cell;

        let mut rng = ChaCha20Rng::from_seed([11u8; 32]);
        let kp = HawkKeypair::generate(&mut rng);

        let attempts = Cell::new(0usize);
        let sig = sign_512_hooked(&kp.secret, b"recover", &mut rng, &mut |attempt, s1| {
            attempts.set(attempt + 1);
            // On the FIRST attempt only, clobber s1 so it overflows on encode.
            // Later attempts are left untouched so the real sampler output is
            // used and the loop can succeed.
            if attempt == 0 {
                for s in s1.iter_mut().take(50) {
                    *s = 511;
                }
            }
        })
        .expect("signing must recover by resigning, not error out");

        // The loop must have resigned at least once (the forced first attempt
        // overflowed, so a real second attempt was required).
        assert!(
            attempts.get() >= 2,
            "expected a resign (>=2 attempts), got {}",
            attempts.get()
        );
        // The recovered signature must be valid and serializable.
        assert!(
            sig.to_bytes().is_ok(),
            "recovered signature must serialize to 555 bytes"
        );
        // And it must verify against the public key (a real signature, not junk).
        assert!(
            kp.public.verify(b"recover", &sig).is_ok(),
            "recovered signature must verify"
        );
    }

    /// The unmodified loop (no-op hook) must equal the production `sign` path
    /// byte-for-byte, proving the seam changes nothing in production.
    #[test]
    fn hooked_signing_with_noop_equals_production() {
        use crate::keygen::HawkKeypair;
        use rand::SeedableRng;
        use rand_chacha::ChaCha20Rng;

        let mut kgen = ChaCha20Rng::from_seed([21u8; 32]);
        let kp = HawkKeypair::generate(&mut kgen);

        let mut rng_a = ChaCha20Rng::from_seed([22u8; 32]);
        let mut rng_b = ChaCha20Rng::from_seed([22u8; 32]);
        let via_prod = kp.secret.sign(b"same", &mut rng_a).unwrap();
        let via_hook = sign_512_hooked(&kp.secret, b"same", &mut rng_b, &mut |_, _| {}).unwrap();
        assert_eq!(via_prod.to_bytes().unwrap(), via_hook.to_bytes().unwrap());
    }

    /// Encoder overflow boundary: pin the threshold so a future encoder change
    /// that shifts it is caught. sum(|c|>>5) must be < the budget to fit.
    #[test]
    fn encode_overflow_boundary() {
        let salt = [0u8; HAWK_SALT_BYTES];
        // k coefficients at 511 contribute k*(511>>5) = k*15 to the unary sum.
        // From the empirical boundary the crossover is between 42 and 45.
        let mk = |k: usize| {
            let mut s1 = vec![0i16; HAWK_N];
            for s in s1.iter_mut().take(k) {
                *s = 511;
            }
            s1
        };
        // 42*15 = 630 fits; 45*15 = 675 overflows.
        assert!(
            crate::serialize::encode_signature(&salt, &mk(42)).is_ok(),
            "sum=630 must fit the 555-byte buffer"
        );
        assert!(
            crate::serialize::encode_signature(&salt, &mk(45)).is_err(),
            "sum=675 must overflow the 555-byte buffer"
        );
    }
}
