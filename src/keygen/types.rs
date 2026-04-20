//! HAWK keypair types: [`HawkPublicKey`], [`HawkSecretKey`], [`HawkKeypair`].

use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use subtle::ConstantTimeEq;
use zeroize::{Zeroize, ZeroizeOnDrop};

/// HAWK-512 public key: the Q-polynomials used for signature verification.
///
/// `q00` and `q01` are the two polynomials of the HAWK public key matrix.
/// `q11` is not stored here; it is derived from `q00` and `q01` during
/// verification.
///
/// Equality uses a data-dependent check (public-key comparison is not
/// security-sensitive). For constant-time comparison of encoded bytes,
/// call [`HawkPublicKey::ct_eq`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HawkPublicKey {
    pub(crate) q00: Vec<i16>,
    pub(crate) q01: Vec<i16>,
}

impl HawkPublicKey {
    /// Constant-time comparison of two public keys.
    ///
    /// Useful when comparing against an expected pinned value — avoids
    /// leaking the first-differing-coefficient position via timing.
    pub fn ct_eq(&self, other: &Self) -> subtle::Choice {
        // Compare coefficient-wise. Both q00 and q01 are length HAWK_N,
        // i.e. 512 for HAWK-512; identical sizes are guaranteed by the type.
        // Serialize each i16 to little-endian bytes so the subtle crate can
        // do byte-wise constant-time comparison.
        let mut a = Vec::with_capacity(4 * crate::params::HAWK_N);
        let mut b = Vec::with_capacity(4 * crate::params::HAWK_N);
        for (x, y) in self.q00.iter().zip(other.q00.iter()) {
            a.extend_from_slice(&x.to_le_bytes());
            b.extend_from_slice(&y.to_le_bytes());
        }
        for (x, y) in self.q01.iter().zip(other.q01.iter()) {
            a.extend_from_slice(&x.to_le_bytes());
            b.extend_from_slice(&y.to_le_bytes());
        }
        a.ct_eq(&b)
    }
}

/// HAWK-512 secret key: the NTRU lattice polynomials plus salt-derivation seed.
///
/// Fields match what the C reference `encode_private` serialises:
/// `kgseed` (24 bytes) || `F mod 2` (64 bytes) || `G mod 2` (64 bytes) ||
/// `hpub` (32 bytes).  The full-precision `f_cap` / `g_cap` vectors are also
/// retained in-memory to avoid re-expanding during signing.
///
/// `hpub` is computed as `SHAKE256(encode_public(q00, q01))[..32]`,
/// matching the C reference `encode_private` (hawk_kgen.c:56-58).
///
/// # Memory safety
///
/// All secret material is zeroed on drop via [`zeroize::ZeroizeOnDrop`].
/// This helps defeat casual memory-scraping attacks where a dropped secret
/// key's backing allocations remain readable (e.g. in core dumps, swap,
/// `/proc/*/mem`). It does **not** defend against attackers with real-time
/// memory access, side-channel capabilities, or compromised allocators.
#[derive(Clone, Debug, Zeroize, ZeroizeOnDrop)]
pub struct HawkSecretKey {
    pub(crate) f: Vec<i8>,
    pub(crate) g: Vec<i8>,
    /// F (capital-F NTRU polynomial).
    pub(crate) f_cap: Vec<i8>,
    /// G (capital-G NTRU polynomial).
    pub(crate) g_cap: Vec<i8>,
    /// kgseed: 24-byte seed used to derive (f, g) via HAWK_REGEN_FG.
    pub(crate) seed: Vec<u8>,
    /// SHAKE-256 hash of the encoded public key, truncated to 32 bytes.
    /// Mixed into the signing hash (hm) during signature generation.
    pub(crate) hpub: [u8; 32],
}

/// A HAWK-512 keypair.
///
/// Generate with [`HawkKeypair::generate`]:
///
/// ```no_run
/// use pqe_hawk::HawkKeypair;
/// use rand::rngs::OsRng;
///
/// let kp = HawkKeypair::generate(&mut OsRng);
/// ```
///
/// The inner [`HawkSecretKey`] zeroes its secret material on drop.
#[derive(Clone, Debug)]
pub struct HawkKeypair {
    pub public: HawkPublicKey,
    pub secret: HawkSecretKey,
}

impl HawkKeypair {
    /// Generate a fresh HAWK-512 keypair using `rng`.
    ///
    /// This is a thin wrapper over [`crate::keygen::hawk_keygen_512`] that
    /// accepts any [`rand::RngCore`] implementor and bundles the output into
    /// the ergonomic [`HawkKeypair`] type.
    pub fn generate<R: rand::RngCore>(rng: &mut R) -> Self {
        // Outer retry loop mirrors the C public keygen (hawk_kgen.c:276-308):
        // after the internal Hawk_keygen candidate passes all algebraic checks,
        // the C also requires the pubkey to Golomb-Rice-encode within the fixed
        // 1024-byte budget. A rare (f,g) distribution produces q00/q01
        // coefficients whose unary tails overflow — those candidates must be
        // discarded. The internal `hawk_keygen_512` stays byte-exact with the
        // C's internal `Hawk_keygen`; the encode-budget check lives here.
        loop {
            let out = crate::keygen::hawk_keygen::hawk_keygen_512(|buf| {
                rng.fill_bytes(buf);
            });

            let pub_bytes = match crate::serialize::encode_public(&out.q00, &out.q01) {
                Ok(b) => b,
                Err(_) => continue,
            };
            // hpub = SHAKE256(encoded_pub_bytes)[..32], matching
            // `encode_private` (hawk_kgen.c:56-58):
            //   shake_init(&sc, 256); shake_inject(&sc, pub, pub_len); shake_flip(&sc);
            //   shake_extract(&sc, buf, hpub_len);
            let mut shake = Shake256::default();
            shake.update(&pub_bytes);
            let mut reader = shake.finalize_xof();
            let mut hpub = [0u8; 32];
            reader.read(&mut hpub);

            return HawkKeypair {
                public: HawkPublicKey {
                    q00: out.q00,
                    q01: out.q01,
                },
                secret: HawkSecretKey {
                    f: out.f,
                    g: out.g,
                    f_cap: out.f_cap,
                    g_cap: out.g_cap,
                    seed: out.seed,
                    hpub,
                },
            };
        }
    }
}

impl HawkPublicKey {
    /// Serialize this public key to its 1024-byte wire format.
    ///
    /// The encoding is Golomb-Rice for q00 (half-size) followed by
    /// Golomb-Rice for q01 (full-size), zero-padded to 1024 bytes.
    pub fn to_bytes(
        &self,
    ) -> Result<[u8; crate::params::HAWK_PUBLIC_KEY_BYTES], crate::error::HawkError> {
        crate::serialize::encode_public(&self.q00, &self.q01)
    }

    /// Deserialize a public key from its 1024-byte wire format.
    pub fn from_bytes(
        bytes: &[u8; crate::params::HAWK_PUBLIC_KEY_BYTES],
    ) -> Result<Self, crate::error::HawkError> {
        let (q00, q01) = crate::serialize::decode_public(bytes)?;
        Ok(HawkPublicKey { q00, q01 })
    }
}

impl HawkSecretKey {
    /// Serialize this secret key to its 184-byte wire format.
    ///
    /// Layout: kgseed(24) || F_mod2(64) || G_mod2(64) || hpub(32).
    pub fn to_bytes(&self) -> [u8; crate::params::HAWK_SECRET_KEY_BYTES] {
        crate::serialize::encode_private(&self.seed, &self.f_cap, &self.g_cap, &self.hpub)
    }

    /// Deserialize a secret key from its 184-byte wire format.
    ///
    /// Reconstructs the in-memory representation:
    /// - `(f, g)` are regenerated from the 24-byte `kgseed` via
    ///   `Hawk_regen_fg` (identical to the keygen derivation).
    /// - `(F, G)` are lifted from the stored `F_mod2 / G_mod2` bitmasks into
    ///   `Vec<i8>` so that their low bits match the mod-2 bytes. The signing
    ///   path only consumes `F_cap & 1` and `G_cap & 1` (via
    ///   `serialize::extract_lowbit`), so a 0/1 lift is sufficient for
    ///   signing; full-precision `F, G` are not needed post-keygen.
    /// - `hpub` is copied verbatim.
    ///
    /// Decoding is infallible on a well-sized input: the 184-byte layout is
    /// fixed, every field has a known position, and no validation is
    /// performed on the kgseed or the mod-2 bits (both are free-form bitmasks
    /// in the HAWK encoding).
    pub fn from_bytes(bytes: &[u8; crate::params::HAWK_SECRET_KEY_BYTES]) -> Self {
        let (seed, f_mod2, g_mod2, hpub_vec) = crate::serialize::decode_private(bytes);

        let mut f = vec![0i8; crate::params::HAWK_N];
        let mut g = vec![0i8; crate::params::HAWK_N];
        crate::keygen::regen_fg::hawk_regen_fg(&mut f, &mut g, &seed);

        let f_cap = lift_mod2_bits(&f_mod2);
        let g_cap = lift_mod2_bits(&g_mod2);

        let mut hpub = [0u8; 32];
        hpub.copy_from_slice(&hpub_vec);

        HawkSecretKey {
            f,
            g,
            f_cap,
            g_cap,
            seed,
            hpub,
        }
    }
}

/// Lift a 64-byte mod-2 bitmask (one bit per coefficient, packed 8-per-byte,
/// low-bit-first) back into a 512-element `Vec<i8>` whose `x & 1` matches
/// the original bit. Inverse of `serialize::extract_lowbit`.
fn lift_mod2_bits(packed: &[u8]) -> Vec<i8> {
    debug_assert_eq!(packed.len(), crate::params::HAWK_N / 8);
    let mut out = vec![0i8; crate::params::HAWK_N];
    for (u, &byte) in packed.iter().enumerate() {
        for i in 0..8 {
            out[u * 8 + i] = ((byte >> i) & 1) as i8;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn hawk_keypair_generate_produces_valid_sizes() {
        let mut rng = ChaCha20Rng::from_seed([7u8; 32]);
        let kp = HawkKeypair::generate(&mut rng);
        assert_eq!(kp.public.q00.len(), crate::params::HAWK_N);
        assert_eq!(kp.public.q01.len(), crate::params::HAWK_N);
        assert_eq!(kp.secret.f.len(), crate::params::HAWK_N);
        assert_eq!(kp.secret.g.len(), crate::params::HAWK_N);
        assert_eq!(kp.secret.f_cap.len(), crate::params::HAWK_N);
        assert_eq!(kp.secret.g_cap.len(), crate::params::HAWK_N);
        assert_eq!(kp.secret.seed.len(), 24);
        assert_eq!(kp.secret.hpub.len(), 32);
    }

    #[test]
    fn hawk_keypair_generate_is_deterministic_for_same_seed() {
        let mut rng1 = ChaCha20Rng::from_seed([42u8; 32]);
        let mut rng2 = ChaCha20Rng::from_seed([42u8; 32]);
        let kp1 = HawkKeypair::generate(&mut rng1);
        let kp2 = HawkKeypair::generate(&mut rng2);
        assert_eq!(kp1.public, kp2.public);
        assert_eq!(kp1.secret.f, kp2.secret.f);
        assert_eq!(kp1.secret.g, kp2.secret.g);
        assert_eq!(kp1.secret.f_cap, kp2.secret.f_cap);
        assert_eq!(kp1.secret.g_cap, kp2.secret.g_cap);
        assert_eq!(kp1.secret.seed, kp2.secret.seed);
        assert_eq!(kp1.secret.hpub, kp2.secret.hpub);
    }

    #[test]
    fn hawk_keypair_different_seeds_differ() {
        let mut rng1 = ChaCha20Rng::from_seed([1u8; 32]);
        let mut rng2 = ChaCha20Rng::from_seed([2u8; 32]);
        let kp1 = HawkKeypair::generate(&mut rng1);
        let kp2 = HawkKeypair::generate(&mut rng2);
        assert_ne!(kp1.public, kp2.public);
    }

    #[test]
    fn pubkey_roundtrip_from_keygen() {
        let mut rng = ChaCha20Rng::from_seed([7u8; 32]);
        let kp = HawkKeypair::generate(&mut rng);
        let bytes = kp.public.to_bytes().unwrap();
        let decoded = HawkPublicKey::from_bytes(&bytes).unwrap();
        // q00[0..n/2] and q01 must round-trip exactly.
        for i in 0..(crate::params::HAWK_N / 2) {
            assert_eq!(decoded.q00[i], kp.public.q00[i], "q00[{}] mismatch", i);
        }
        assert_eq!(decoded.q01, kp.public.q01);
        // q00[n/2..n] is not encoded (auto-adjoint property);
        // decoded slots are zero.
        for i in (crate::params::HAWK_N / 2)..crate::params::HAWK_N {
            assert_eq!(decoded.q00[i], 0, "q00[{}] should be zero after decode", i);
        }
    }

    #[test]
    fn secret_roundtrip_from_keygen() {
        let mut rng = ChaCha20Rng::from_seed([7u8; 32]);
        let kp = HawkKeypair::generate(&mut rng);
        let bytes = kp.secret.to_bytes();
        let (seed, f_mod2, g_mod2, hpub) = crate::serialize::decode_private(&bytes);
        assert_eq!(seed, kp.secret.seed);
        assert_eq!(hpub.as_slice(), kp.secret.hpub.as_slice());
        assert_eq!(f_mod2.len(), 64);
        assert_eq!(g_mod2.len(), 64);
        // Verify F_mod2 bits match the low bits of f_cap.
        for i in 0..crate::params::HAWK_N {
            let bit = (f_mod2[i / 8] >> (i % 8)) & 1;
            assert_eq!(
                bit as i8,
                kp.secret.f_cap[i] & 1,
                "F_mod2 bit {} mismatch",
                i
            );
        }
    }

    #[test]
    fn zeroize_on_drop_clears_secret_material() {
        use zeroize::Zeroize;
        let mut rng = ChaCha20Rng::from_seed([99u8; 32]);
        let kp = HawkKeypair::generate(&mut rng);
        let mut sk = kp.secret;
        // Populated with non-trivial material.
        assert!(sk.f.iter().any(|&v| v != 0));
        assert!(sk.seed.iter().any(|&v| v != 0));
        assert!(sk.hpub.iter().any(|&v| v != 0));

        // Explicit zeroize (the ZeroizeOnDrop derive would do this at drop).
        sk.zeroize();
        assert!(sk.f.iter().all(|&v| v == 0), "f not zeroed");
        assert!(sk.g.iter().all(|&v| v == 0), "g not zeroed");
        assert!(sk.f_cap.iter().all(|&v| v == 0), "f_cap not zeroed");
        assert!(sk.g_cap.iter().all(|&v| v == 0), "g_cap not zeroed");
        assert!(sk.seed.iter().all(|&v| v == 0), "seed not zeroed");
        assert!(sk.hpub.iter().all(|&v| v == 0), "hpub not zeroed");
    }

    #[test]
    fn pubkey_ct_eq_agrees_with_eq() {
        let mut rng = ChaCha20Rng::from_seed([77u8; 32]);
        let kp1 = HawkKeypair::generate(&mut rng);
        let mut rng2 = ChaCha20Rng::from_seed([77u8; 32]);
        let kp2 = HawkKeypair::generate(&mut rng2);
        let mut rng3 = ChaCha20Rng::from_seed([88u8; 32]);
        let kp3 = HawkKeypair::generate(&mut rng3);
        assert!(bool::from(kp1.public.ct_eq(&kp2.public)));
        assert!(!bool::from(kp1.public.ct_eq(&kp3.public)));
    }

    #[test]
    fn secret_from_bytes_roundtrip_fields() {
        let mut rng = ChaCha20Rng::from_seed([13u8; 32]);
        let kp = HawkKeypair::generate(&mut rng);
        let bytes = kp.secret.to_bytes();
        let sk2 = HawkSecretKey::from_bytes(&bytes);
        assert_eq!(sk2.seed, kp.secret.seed);
        assert_eq!(sk2.hpub, kp.secret.hpub);
        assert_eq!(sk2.f, kp.secret.f);
        assert_eq!(sk2.g, kp.secret.g);
        for i in 0..crate::params::HAWK_N {
            assert_eq!(
                sk2.f_cap[i] & 1,
                kp.secret.f_cap[i] & 1,
                "F low-bit mismatch at {i}"
            );
            assert_eq!(
                sk2.g_cap[i] & 1,
                kp.secret.g_cap[i] & 1,
                "G low-bit mismatch at {i}"
            );
        }
    }

    #[test]
    fn secret_from_bytes_produces_usable_signing_key() {
        let mut rng = ChaCha20Rng::from_seed([21u8; 32]);
        let kp = HawkKeypair::generate(&mut rng);
        let bytes = kp.secret.to_bytes();
        let sk2 = HawkSecretKey::from_bytes(&bytes);

        let msg = b"restored-key signing test";
        let mut sig_rng = ChaCha20Rng::from_seed([22u8; 32]);
        let sig = sk2.sign(msg, &mut sig_rng).expect("restored sk must sign");
        // The associated public key (unchanged across serialize/deserialize)
        // must accept the signature.
        assert!(kp.public.verify(msg, &sig).is_ok());
    }

    #[test]
    fn secret_to_bytes_roundtrip_is_idempotent() {
        let mut rng = ChaCha20Rng::from_seed([33u8; 32]);
        let kp = HawkKeypair::generate(&mut rng);
        let bytes1 = kp.secret.to_bytes();
        let sk2 = HawkSecretKey::from_bytes(&bytes1);
        let bytes2 = sk2.to_bytes();
        assert_eq!(bytes1, bytes2);
    }
}
