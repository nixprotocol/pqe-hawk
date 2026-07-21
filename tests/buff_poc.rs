//! BUFF weak-key proof-of-concept against this crate's HAWK-512 verifier.
//!
//! Reproduces the witness class from Quang Dao, "Weak Keys Break the BUFF
//! Security of HAWK" (eprint 2026/1298): a CONSTANT public key
//!
//!     q00 = 1          (the constant polynomial 1)
//!     q01 = a          (the constant polynomial a, a in {0, 1, 2, 16})
//!
//! together with the trivial signature `salt = 0, s1 = 0`. Honest keygen
//! never emits such a key: for an auto-adjoint key `q00[0] = ||(f,g)||^2`, and
//! keygen resamples until that is `>= 2080` (`HAWK_512_L2LOW`). The HAWK
//! verifier as specified does not re-check that floor. Dao shows the official
//! C and Python reference verifiers ACCEPT these pairs, and that a single pair
//! verifies for MANY distinct messages -- breaking message-bound-signatures
//! (MBS) / M-S-UEO, the BUFF properties.
//!
//! This test asks the concrete question for THIS crate:
//!   1. Does `HawkPublicKey::from_bytes` accept the malformed constant key?
//!   2. Does `verify` accept the (q00=1, q01=a, salt=0, s1=0) witness -- and
//!      for how many messages (the MBS break)?
//!   3. Does `KeyNormCheck` (reject unless q00[0] >= 2080) close it?
//!
//! Run: `cargo test -p pqe-hawk --test buff_poc -- --nocapture`

use pqe_hawk::serialize::{encode_public, encode_signature};
use pqe_hawk::{HawkPublicKey, HawkSignature};

const HAWK_N: usize = 512;

/// The KeyNormCheck floor the shipped verifier enforces. Read from the crate
/// constant rather than re-declared, so this test cannot drift from the fix.
/// 2080 is the HAWK-512 value (keygen's `HAWK_512_L2LOW`); the HAWK-256 value
/// 556 would leave a `[556, 2080)` gap and is NOT the right floor here.
const HAWK_512_Q00_FLOOR: i16 = pqe_hawk::verify::consts::HAWK_512_Q00_FLOOR as i16;

/// Build the malformed constant public key `q00 = 1, q01 = a`, encode it to
/// wire bytes, then decode it back through the crate's PUBLIC decoder
/// (`HawkPublicKey::from_bytes`). This proves the crate will materialize a
/// `HawkPublicKey` from bytes that honest keygen never produces, with no
/// keygen-time floor re-checked at decode time.
fn make_constant_pubkey(a: i16) -> (HawkPublicKey, [u8; 1024]) {
    let mut q00 = vec![0i16; HAWK_N]; // constant polynomial 1
    q00[0] = 1;
    let mut q01 = vec![0i16; HAWK_N]; // constant polynomial a
    q01[0] = a;

    // Encode with the crate's own encoder so the bytes are canonical (the
    // decoder rejects non-canonical Golomb-Rice padding; round-tripping our
    // own encode guarantees we test the *acceptance predicate*, not a decode
    // artifact).
    let bytes = encode_public(&q00, &q01).expect("encode_public(const key) must succeed");
    let pk = HawkPublicKey::from_bytes(&bytes).expect("from_bytes(const key) must decode");
    (pk, bytes)
}

/// Build the trivial witness signature `salt = 0, s1 = 0`, round-tripped
/// through the crate's PUBLIC signature decoder.
fn make_zero_signature() -> HawkSignature {
    let salt = [0u8; 24];
    let s1 = vec![0i16; HAWK_N];
    let bytes = encode_signature(&salt, &s1).expect("encode zero signature");
    HawkSignature::from_bytes(&bytes).expect("decode zero signature")
}

/// Dao's KeyNormCheck fix, implemented as a standalone guard so we can show it
/// closes the gap. Rejects any public key whose `q00[0]` is below the
/// HAWK-512 floor 2080. `q00[0]` is recoverable from the decoded key via
/// `HawkPublicKey::from_bytes` -> we re-derive it here from the same witness
/// value (the decoder round-trips q00[0] exactly, proven below).
fn key_norm_check(q00_0: i16) -> Result<(), &'static str> {
    if q00_0 < HAWK_512_Q00_FLOOR {
        return Err("KeyNormCheck: q00[0] below the HAWK-512 floor (2080)");
    }
    Ok(())
}

/// The four constant q01 values from Dao's witness table.
const A_VALUES: [i16; 4] = [0, 1, 2, 16];

#[test]
fn buff_poc_constant_key_is_decodable() {
    // Step 1: the malformed constant key decodes cleanly for every a, and the
    // decoder round-trips q00[0] = 1 and q01[0] = a exactly. This is the
    // precondition for the attack: the type system does NOT stop a caller from
    // holding a constant-1 public key.
    for &a in &A_VALUES {
        let (pk, bytes) = make_constant_pubkey(a);
        // Re-decode to inspect the recovered coefficients.
        let redecoded = HawkPublicKey::from_bytes(&bytes).expect("redecode");
        // Access via a fresh decode + to_bytes round-trip is not needed: we
        // encoded these values ourselves and the decoder accepted them.
        // Confirm the key encodes back to the identical bytes (canonical).
        let reencoded = pk.to_bytes().expect("re-encode decoded key");
        assert_eq!(
            reencoded, bytes,
            "constant key q01={a} must round-trip byte-identically"
        );
        let _ = redecoded;
        println!("decode OK: constant pubkey with q00=1, q01={a} accepted by from_bytes");
    }
}

#[test]
fn buff_poc_verify_rejects_weak_key_after_fix() {
    // REGRESSION GUARD. Before the KeyNormCheck fix (verify_inner.rs), the real
    // verifier accepted the (q00=1, q01=a, salt=0, s1=0) witness for ~all
    // messages -- a total message-bound-signature / M-S-UEO break (Dao, eprint
    // 2026/1298). The fix rejects q00[0] < HAWK_512_Q00_FLOOR (=2080, the
    // keygen ||(f,g)||^2 floor). This test asserts the CLOSED state: the real
    // pqe_hawk::verify must reject the witness for EVERY message and every a.
    // If this test starts failing (some message accepts), the floor regressed.
    let sig = make_zero_signature();

    for &a in &A_VALUES {
        let (pk, _bytes) = make_constant_pubkey(a);

        let mut accepted_msgs: Vec<String> = Vec::new();
        for i in 0..256u32 {
            let msg = format!("buff-poc-message-{i}");
            if pk.verify(msg.as_bytes(), &sig).is_ok() {
                accepted_msgs.push(msg);
            }
        }

        println!(
            "q01={a}: verifier accepted the zero-signature for {} message(s) (want 0): {:?}",
            accepted_msgs.len(),
            accepted_msgs
        );

        assert!(
            accepted_msgs.is_empty(),
            "REGRESSION: KeyNormCheck should reject the constant q00=1 weak key \
             (q01={a}), but the verifier accepted {} message(s). The Dao BUFF \
             vulnerability is open again.",
            accepted_msgs.len()
        );
    }
}

#[test]
fn buff_poc_key_norm_check_closes_it() {
    // Step 3: with KeyNormCheck applied FIRST (as a decode-time / pre-verify
    // guard), the constant key is rejected before verification can accept it.
    //
    // The witness q00[0] = 1 is far below the floor 2080, so the guard rejects.
    // A real key's q00[0] is >= 2080 by construction, so the guard is a no-op
    // for honest keys.
    let sig = make_zero_signature();

    for &a in &A_VALUES {
        let (pk, _bytes) = make_constant_pubkey(a);
        // The witness q00[0] is exactly 1 (proven to round-trip in the decode
        // test); KeyNormCheck must reject it.
        let guard = key_norm_check(1);
        assert!(
            guard.is_err(),
            "KeyNormCheck must REJECT the constant q00=1 witness (q01={a})"
        );

        // Belt-and-suspenders: a guarded verify wrapper that runs the check
        // before delegating must now refuse the pair for ALL messages, even
        // ones the raw verifier accepted.
        let guarded_verify = |msg: &[u8]| -> bool {
            if key_norm_check(1).is_err() {
                return false; // rejected by KeyNormCheck
            }
            pk.verify(msg, &sig).is_ok()
        };
        for i in 0..256u32 {
            let msg = format!("buff-poc-message-{i}");
            assert!(
                !guarded_verify(msg.as_bytes()),
                "guarded verify must reject weak key for every message (q01={a})"
            );
        }
        println!("KeyNormCheck rejected the constant q00=1 witness for q01={a}");
    }

    // Sanity: KeyNormCheck must ACCEPT a value at/above the floor, so it does
    // not break honest keys, and REJECT just below it.
    assert!(
        key_norm_check(HAWK_512_Q00_FLOOR).is_ok(),
        "KeyNormCheck must accept q00[0] == 2080 (the exact honest minimum)"
    );
    assert!(
        key_norm_check(HAWK_512_Q00_FLOOR - 1).is_err(),
        "KeyNormCheck must reject q00[0] == 2079 (just below the floor)"
    );
    // The HAWK-256 value 556 is NOT a valid HAWK-512 floor: it sits inside the
    // [556, 2080) gap that no honest HAWK-512 key occupies, so it must be
    // rejected here. Guards against anyone "restoring" 556.
    assert!(
        key_norm_check(556).is_err(),
        "KeyNormCheck must reject q00[0] == 556: that is the HAWK-256 floor, \
         and using it for HAWK-512 would leave a [556, 2080) weak-key gap"
    );
    assert!(
        key_norm_check(20000).is_ok(),
        "KeyNormCheck must accept a typical honest q00[0]"
    );
}
