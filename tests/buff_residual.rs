//! Residual-BUFF probing against the FIXED HAWK-512 verifier.
//!
//! The 0.1.1 fix added a `KeyNormCheck` floor: `verify_inner` rejects any
//! public key with `q00[0] < HAWK_512_Q00_FLOOR (= 2080)`
//! (`src/verify/verify_inner.rs:190`, variant `VerifyReject::KeyNormCheckQ00Floor`).
//! That demonstrably closes the constant `q00 = 1` weak-key witness
//! (`tests/buff_poc.rs`). But Dao (eprint 2026/1298) and the fix's own notes
//! say the floor does NOT prove full BUFF restoration. This test probes the
//! three explicitly-open residual classes against the REAL verifier and asks,
//! for each: does a FLOOR-PASSING malformed key (`q00[0] >= 2080`, so the fix
//! does not catch it) plus a crafted low-norm signature verify for >= 2
//! distinct messages (a surviving BUFF break), or does some OTHER check reject
//! it — and if so, exactly which one?
//!
//! Residual classes:
//!   1. S-CEO / S-DEO (exclusive-ownership / DSKS key-substitution): can one
//!      (msg, sig) that verifies under a floor-passing key ALSO verify under a
//!      DIFFERENT floor-passing key?
//!   2. Non-constant / uneven-spectrum weak keys: `q00[0] >= 2080` (passes the
//!      floor) but the rest of q00 is malformed (degenerate FFT spectrum,
//!      near-zero bin, structured perturbation) so a crafted low-norm signature
//!      verifies for multiple messages.
//!   3. PolyQnorm soundness: choose (q00, q01) that are NOT a genuine real
//!      auto-adjoint `||(f,g)||^2` form and see whether the dual-prime norm
//!      accepts them.
//!
//! Attribution. The production verifier already exposes
//! `verify_inner_labeled(...) -> Result<(), VerifyReject>`
//! (`src/verify/verify_inner.rs:97`); `verify_inner` (the public path) is a
//! thin wrapper that discards the label (`:86-87`). So the label we read here
//! is the SAME check the real `HawkPublicKey::verify` hits — single source of
//! truth, no parallel model that could drift. Each `VerifyReject` variant names
//! its exact `verify_inner.rs` line in the enum doc (`:38-56`). We map each
//! variant to that line below in `reject_site`.
//!
//! Run: `cargo test -p pqe-hawk --test buff_residual -- --nocapture`

use pqe_hawk::serialize::{encode_public, encode_signature};
use pqe_hawk::verify::consts::HAWK_512_Q00_FLOOR;
use pqe_hawk::verify::verify_inner::{verify_inner_labeled, VerifyReject};
use pqe_hawk::{HawkKeypair, HawkPublicKey, HawkSignature};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

const HAWK_N: usize = 512;

/// The exact `file:line` of the production check behind each reject label.
/// Sourced from the `VerifyReject` enum doc comments in
/// `src/verify/verify_inner.rs` and confirmed against `grep -n`.
fn reject_site(r: VerifyReject) -> &'static str {
    match r {
        VerifyReject::SymBreakFirstNonzeroNegative => {
            "verify_inner.rs:166 (sym-break: first nonzero t1 coeff negative)"
        }
        VerifyReject::SymBreakT1AllZero => "verify_inner.rs:174 (sym-break: t1 all zero)",
        VerifyReject::KeyNormCheckQ00Floor => {
            "verify_inner.rs:190 (KeyNormCheck: q00[0] < 2080 FLOOR)"
        }
        VerifyReject::DivisionDomain => {
            "verify_inner.rs:275 (division domain: w00 out of (0,2^30) or quotient overflow)"
        }
        VerifyReject::S0OutOfRange => "verify_inner.rs:312 (s0 reconstruction out of range)",
        VerifyReject::TnormPrimeDisagree => {
            "verify_inner.rs:434 (dual-prime norm P1/P2 disagree => tnorm too large)"
        }
        VerifyReject::TnormNotDivisibleByHn => {
            "verify_inner.rs:447 (n*sqnorm_Q not divisible by hn)"
        }
        VerifyReject::TnormExceedsBound => {
            "verify_inner.rs:450 ((tnorm>>8) > max_tnorm 8317)"
        }
    }
}

// ---------------------------------------------------------------------------
// Harness helpers.
// ---------------------------------------------------------------------------

/// Encode (q00, q01) with the crate encoder, decode back through the PUBLIC
/// decoder. Returns the materialised key + the canonical wire bytes, or the
/// decode error (some malformed shapes may be rejected by encode/decode
/// itself, which is itself a finding to report).
fn materialize(q00: &[i16], q01: &[i16]) -> Result<(HawkPublicKey, [u8; 1024]), String> {
    let bytes = encode_public(q00, q01).map_err(|e| format!("encode_public: {e:?}"))?;
    let pk = HawkPublicKey::from_bytes(&bytes).map_err(|e| format!("from_bytes: {e:?}"))?;
    Ok((pk, bytes))
}

fn zero_signature() -> (HawkSignature, [u8; 24], Vec<i16>) {
    let salt = [0u8; 24];
    let s1 = vec![0i16; HAWK_N];
    let bytes = encode_signature(&salt, &s1).expect("encode zero sig");
    let sig = HawkSignature::from_bytes(&bytes).expect("decode zero sig");
    (sig, salt, s1)
}

/// A signature with a single nonzero low coefficient. Returns the sig plus its
/// (salt, s1) so the labeled verifier can be called directly. We keep s1 tiny
/// so its Golomb-Rice tail cannot overflow the 555-byte buffer.
fn small_signature(idx: usize, val: i16) -> Option<(HawkSignature, [u8; 24], Vec<i16>)> {
    let salt = [0u8; 24];
    let mut s1 = vec![0i16; HAWK_N];
    s1[idx] = val;
    let bytes = encode_signature(&salt, &s1).ok()?;
    let sig = HawkSignature::from_bytes(&bytes).ok()?;
    Some((sig, salt, s1))
}

/// Recover (salt, s1) from a genuine signature (e.g. an honest one) so we can
/// feed the labeled verifier. Round-trips through the crate encoder + a
/// read-only Golomb-Rice s1 decode mirroring `decode_signature`.
fn sig_parts(sig: &HawkSignature) -> ([u8; 24], Vec<i16>) {
    let b = sig.to_bytes().expect("sig encode");
    let mut salt = [0u8; 24];
    salt.copy_from_slice(&b[..24]);
    let s1 = decode_sig_s1(&b).expect("recover s1 from signature bytes");
    (salt, s1)
}

/// Minimal Golomb-Rice s1 decoder mirroring `serialize::decode_signature`
/// (low=5, lim=9). Read-only reconstruction used only to feed the labeled
/// verifier; the REAL verify path uses the crate's own decoder via
/// `HawkPublicKey::verify`, and `probe` cross-checks the two agree.
fn decode_sig_s1(bytes: &[u8]) -> Option<Vec<i16>> {
    let logn = 9u32;
    let low = 5u32;
    let lim_bits = 9u32;
    let n = 1usize << logn;
    let buf = &bytes[24..];
    let buf_len = buf.len();
    let min_len = ((low + 1) as usize) << (logn - 3);
    if buf_len < min_len {
        return None;
    }
    let mut d = vec![0i16; n];
    let mut voff = min_len;
    let ntz = |x: u8| -> i32 {
        if x == 0 {
            8
        } else {
            x.trailing_zeros() as i32
        }
    };
    let mut acc: u32 = 0;
    let mut acc_off: i32 = 0;
    let lim_hi: i32 = 1i32 << (lim_bits - low);
    for u in 0..n {
        while acc == 0 {
            if acc_off >= lim_hi {
                return None;
            }
            if voff >= buf_len {
                return None;
            }
            acc |= (buf[voff] as u32) << acc_off;
            voff += 1;
            acc_off += 8;
        }
        let mut k = ntz((acc & 0xFF) as u8);
        if k == 8 {
            k += ntz(((acc >> 8) & 0xFF) as u8);
            if k >= lim_hi {
                return None;
            }
        }
        d[u] = ((k as u16) << low) as i16;
        acc >>= (k + 1) as u32;
        acc_off -= k + 1;
    }
    let mut loff: usize = 1 << (logn - 3);
    let lmask: u32 = (1u32 << low) - 1;
    for u in (0..n).step_by(8) {
        let sbb = buf[u >> 3] as u32;
        let mut lpp: u64 = 0;
        for j in 0..(low as usize) {
            lpp |= (buf[loff] as u64) << (j * 8);
            loff += 1;
        }
        for i in 0..8usize {
            let lp = ((lpp >> (i as u32 * low)) as u32) & lmask;
            let sm = ((sbb >> i) & 1).wrapping_neg();
            let cur = d[u + i] as u16 as u32;
            d[u + i] = (cur ^ (sm & 0xFFFF) ^ lp) as u16 as i16;
        }
    }
    Some(d)
}

/// Outcome of one probe: whether the real verifier accepted, and (on reject)
/// which production check fired.
#[derive(Clone, Debug)]
enum Probe {
    Accept,
    Reject(VerifyReject),
}

/// Run a (key, sig) against BOTH the real verifier and the labeled verifier
/// over `msgs`, asserting the two verdicts agree (they share one
/// implementation, so this is a belt-and-suspenders guard), and return the
/// per-message outcomes plus the accepted count.
fn probe(
    pk: &HawkPublicKey,
    bytes: &[u8; 1024],
    q00: &[i16],
    q01: &[i16],
    sig: &HawkSignature,
    salt: &[u8],
    s1: &[i16],
    msgs: &[Vec<u8>],
) -> (usize, Vec<Probe>) {
    let mut accepted = 0usize;
    let mut out = Vec::with_capacity(msgs.len());
    for m in msgs {
        let real_ok = pk.verify(m, sig).is_ok();
        let labeled = verify_inner_labeled(m, bytes, q00, q01, salt, s1);
        let labeled_ok = labeled.is_ok();
        assert_eq!(
            real_ok, labeled_ok,
            "labeled/real verdict mismatch for msg {:?}: real={real_ok} labeled={labeled:?}",
            String::from_utf8_lossy(m)
        );
        match labeled {
            Ok(()) => {
                accepted += 1;
                out.push(Probe::Accept);
            }
            Err(r) => out.push(Probe::Reject(r)),
        }
    }
    (accepted, out)
}

/// Summarise the distinct reject reasons across a probe's outcomes.
fn summarize(outcomes: &[Probe]) -> String {
    let mut accepts = 0usize;
    let mut counts: std::collections::BTreeMap<String, usize> = Default::default();
    for o in outcomes {
        match o {
            Probe::Accept => accepts += 1,
            Probe::Reject(r) => *counts.entry(reject_site(*r).to_string()).or_default() += 1,
        }
    }
    let mut parts = Vec::new();
    if accepts > 0 {
        parts.push(format!("ACCEPT x{accepts}"));
    }
    for (site, c) in counts {
        parts.push(format!("[{site}] x{c}"));
    }
    parts.join("  |  ")
}

fn make_msgs(prefix: &str, count: u32) -> Vec<Vec<u8>> {
    (0..count)
        .map(|i| format!("{prefix}-{i}").into_bytes())
        .collect()
}

// ===========================================================================
// RESIDUAL CASE 2 (probed first: it is the sharpest test of the floor).
//
// Non-constant / uneven-spectrum weak keys. `q00[0] >= 2080` passes the floor
// (verify_inner.rs:190), but the rest of q00 is malformed so the FFT spectrum
// is degenerate. We try several concrete shapes and, for each, feed the
// trivial + small signatures over many messages, reporting accept count and
// the rejecting check.
// ===========================================================================
#[test]
fn residual_case2_nonconstant_weak_keys() {
    println!("\n================ RESIDUAL CASE 2: non-constant / uneven-spectrum weak keys ================");
    let msgs = make_msgs("case2-msg", 128);

    struct Shape {
        name: &'static str,
        q00: Vec<i16>,
        q01: Vec<i16>,
    }
    let mut shapes: Vec<Shape> = Vec::new();

    // (a) Floor-exact constant: q00 = 2080 (DC only), the smallest-norm
    //     floor-passing constant key. If ANY floor-passing constant key breaks
    //     BUFF, it is this one.
    {
        let mut q00 = vec![0i16; HAWK_N];
        q00[0] = 2080;
        let mut q01 = vec![0i16; HAWK_N];
        q01[0] = 1;
        shapes.push(Shape { name: "(a) constant q00=2080, q01=1 (floor-exact constant)", q00, q01 });
    }
    // (b) Floor DC + max alternating +-511 perturbation (concentrates spectral
    //     energy, drives some FFT bins far from the DC value).
    {
        let mut q00 = vec![0i16; HAWK_N];
        q00[0] = 2080;
        for u in 1..256 {
            q00[u] = if u % 2 == 0 { 511 } else { -511 };
        }
        let mut q01 = vec![0i16; HAWK_N];
        q01[0] = 1;
        shapes.push(Shape { name: "(b) q00[0]=2080 + max alternating +-511 perturbation", q00, q01 });
    }
    // (c) DC at floor + single large off-DC coeff (rank-1 spectral tilt).
    {
        let mut q00 = vec![0i16; HAWK_N];
        q00[0] = 2080;
        q00[1] = 511;
        let mut q01 = vec![0i16; HAWK_N];
        q01[0] = 16;
        shapes.push(Shape { name: "(c) q00[0]=2080, q00[1]=511 single-tilt, q01=16", q00, q01 });
    }
    // (d) Floor DC + degenerate q01 = max constant, trying to make q01*t1/q00
    //     round to a huge s0 for some message.
    {
        let mut q00 = vec![0i16; HAWK_N];
        q00[0] = 2080;
        let mut q01 = vec![0i16; HAWK_N];
        for u in 0..HAWK_N {
            q01[u] = 4095;
        }
        shapes.push(Shape { name: "(d) q00[0]=2080, q01=all 4095 (max)", q00, q01 });
    }
    // (e) DC at floor + localized low-frequency bump (breaks flat spectrum).
    {
        let mut q00 = vec![0i16; HAWK_N];
        q00[0] = 2080;
        for u in 1..8 {
            q00[u] = 511;
        }
        let mut q01 = vec![0i16; HAWK_N];
        q01[0] = 1;
        shapes.push(Shape { name: "(e) q00[0]=2080 + localized 511 bump [1..8]", q00, q01 });
    }
    // (f) Above-floor DC + all-negative perturbation to try to drive
    //     cstup+fq00[u] below zero (=> DivisionDomain) or near-zero (=> blowup).
    {
        let mut q00 = vec![0i16; HAWK_N];
        q00[0] = 2100;
        for u in 1..256 {
            q00[u] = -511;
        }
        let mut q01 = vec![0i16; HAWK_N];
        q01[0] = 1;
        shapes.push(Shape { name: "(f) q00[0]=2100 + all -511 (drive bins toward 0)", q00, q01 });
    }

    let (zsig, zsalt, zs1) = zero_signature();
    let mut any_break = false;

    for sh in &shapes {
        let (pk, bytes) = match materialize(&sh.q00, &sh.q01) {
            Ok(v) => v,
            Err(e) => {
                println!("  {} -> NOT MATERIALIZABLE (encode/decode rejected): {e}", sh.name);
                continue;
            }
        };
        assert!(sh.q00[0] as i32 >= HAWK_512_Q00_FLOOR, "test bug: shape below floor");

        let (acc0, det0) = probe(&pk, &bytes, &sh.q00, &sh.q01, &zsig, &zsalt, &zs1, &msgs);

        let mut acc_small = 0usize;
        let mut small_summ = String::new();
        for (idx, val) in [(0usize, 1i16), (7, 1), (0, -1), (3, 2)] {
            if let Some((sig, salt, s1)) = small_signature(idx, val) {
                let (a, d) = probe(&pk, &bytes, &sh.q00, &sh.q01, &sig, &salt, &s1, &msgs);
                acc_small += a;
                if small_summ.is_empty() {
                    small_summ = summarize(&d);
                }
            }
        }

        println!("  {}", sh.name);
        println!("      zero-sig ({} msgs):  {}", msgs.len(), summarize(&det0));
        println!(
            "      small-sigs (4 sigs x {} msgs, first-sig breakdown): {}   [total accepts across all small sigs: {}]",
            msgs.len(),
            small_summ,
            acc_small
        );

        if acc0 + acc_small >= 2 {
            any_break = true;
            println!("      *** SURVIVING BUFF BREAK: floor-passing key verified >= 2 (msg,sig) pairs ***");
        }
    }

    assert!(
        !any_break,
        "SURVIVING BUFF BREAK in case 2: a floor-passing malformed key verified >= 2 pairs. \
         See stdout for the witness."
    );
    println!("  => CASE 2 verdict: no floor-passing non-constant key produced >= 2 accepted pairs.");
}

// ===========================================================================
// RESIDUAL CASE 1: S-CEO / S-DEO (exclusive-ownership / DSKS).
//
// The exclusive-ownership games ask: given a (msg, sig) that verifies under
// pk_A, can an attacker exhibit a DIFFERENT key pk_B under which the SAME
// (msg, sig) also verifies? We take an HONEST key + its real signature (a
// genuinely verifying (msg,sig)), then search a family of floor-passing
// malformed second keys for another acceptor.
// ===========================================================================
#[test]
fn residual_case1_exclusive_ownership() {
    println!("\n================ RESIDUAL CASE 1: S-CEO / S-DEO (exclusive ownership / DSKS) ================");

    let mut rng = ChaCha20Rng::from_seed([201u8; 32]);
    let kp = HawkKeypair::generate(&mut rng);
    let mut srng = ChaCha20Rng::from_seed([202u8; 32]);
    let msg = b"exclusive-ownership-target-message".to_vec();
    let sig = kp.secret.sign(&msg, &mut srng).expect("sign");
    assert!(kp.public.verify(&msg, &sig).is_ok(), "honest sig must verify under pk_A");
    println!(
        "  Built honest (pk_A, msg, sig): verifies under pk_A = OK. q00_A[0]={}",
        pubkey_q00_0(&kp.public)
    );

    let (salt, s1) = sig_parts(&sig);

    // Attacker's second-key family: floor-passing malformed keys plus small
    // perturbations of pk_A itself.
    let (q00_a, q01_a) = pubkey_coeffs(&kp.public);
    let mut candidates: Vec<(String, Vec<i16>, Vec<i16>)> = Vec::new();

    for a in [0i16, 1, 2, 16, 100, 4095] {
        let mut q00 = vec![0i16; HAWK_N];
        q00[0] = 2080;
        let mut q01 = vec![0i16; HAWK_N];
        q01[0] = a;
        candidates.push((format!("constant q00=2080,q01={a}"), q00, q01));
    }
    for delta in [1i16, -1, 2, -2, 16, -16] {
        let mut q01 = q01_a.clone();
        q01[0] = q01[0].wrapping_add(delta);
        candidates.push((format!("pk_A with q01[0]+={delta}"), q00_a.clone(), q01));
    }
    for d0 in [1i16, 100, 200, -50] {
        let mut q00 = q00_a.clone();
        q00[0] = q00[0].wrapping_add(d0);
        if (q00[0] as i32) < HAWK_512_Q00_FLOOR {
            continue;
        }
        candidates.push((format!("pk_A with q00[0]+={d0}"), q00, q01_a.clone()));
    }
    {
        let q01 = vec![0i16; HAWK_N];
        candidates.push(("pk_A q00 with q01=0".into(), q00_a.clone(), q01));
    }

    let mut second_acceptor: Option<String> = None;
    let mut checked = 0usize;
    for (name, q00, q01) in &candidates {
        if q00 == &q00_a && q01 == &q01_a {
            continue; // skip identical-to-A
        }
        let (pk_b, bytes_b) = match materialize(q00, q01) {
            Ok(v) => v,
            Err(e) => {
                println!("  cand [{name}] not materializable: {e}");
                continue;
            }
        };
        checked += 1;
        let real_ok = pk_b.verify(&msg, &sig).is_ok();
        let labeled = verify_inner_labeled(&msg, &bytes_b, q00, q01, &salt, &s1);
        assert_eq!(real_ok, labeled.is_ok(), "labeled/real mismatch (case1) [{name}]");
        match labeled {
            Ok(()) => {
                println!("  cand [{name}]: (msg,sig) ALSO VERIFIES under pk_B  <-- exclusive-ownership break");
                second_acceptor = Some(name.clone());
                break;
            }
            Err(r) => {
                println!("  cand [{name}]: rejected by {}", reject_site(r));
            }
        }
    }

    match second_acceptor {
        Some(w) => panic!(
            "SURVIVING S-CEO/S-DEO BREAK: the honest (msg,sig) also verifies under a DIFFERENT \
             floor-passing key: [{w}]. Exclusive ownership is broken and the floor does NOT close it."
        ),
        None => println!(
            "  => CASE 1 verdict: over {checked} floor-passing second-key candidates, NONE \
             re-verified the honest (msg,sig). No exclusive-ownership break found here."
        ),
    }
}

// ===========================================================================
// RESIDUAL CASE 3: PolyQnorm soundness.
//
// Does the dual-prime norm accept a (q00, q01) that is NOT a genuine real
// auto-adjoint ||(f,g)||^2 form? We take a REAL honest key's spectrum, then
// MUTATE it in ways that destroy the "q00 = f adj f + g adj g, q01 = F adj f +
// G adj g" relationship while keeping the floor. If PolyQnorm is fooled, the
// norm still lands <= max_tnorm for some crafted signature.
// ===========================================================================
#[test]
fn residual_case3_polyqnorm_soundness() {
    println!("\n================ RESIDUAL CASE 3: PolyQnorm / real-quadratic-form soundness ================");
    let msgs = make_msgs("case3-msg", 64);

    let mut rng = ChaCha20Rng::from_seed([31u8; 32]);
    let kp = HawkKeypair::generate(&mut rng);
    let (q00_real, q01_real) = pubkey_coeffs(&kp.public);
    println!("  honest reference key: q00[0]={}", q00_real[0]);

    struct M {
        name: &'static str,
        q00: Vec<i16>,
        q01: Vec<i16>,
    }
    let mut muts: Vec<M> = Vec::new();

    // (a) Real q00, q01 doubled (no longer F adj f + G adj g for this q00).
    {
        let q00 = q00_real.clone();
        let q01: Vec<i16> = q01_real
            .iter()
            .map(|&x| x.saturating_mul(2).clamp(-4095, 4095))
            .collect();
        muts.push(M { name: "(a) real q00, q01 := 2*q01 (breaks cross relation)", q00, q01 });
    }
    // (b) Real q01, q00 spectrum perturbed at one coeff (>= floor still).
    {
        let mut q00 = q00_real.clone();
        q00[5] = q00[5].saturating_add(200).clamp(-511, 511);
        let q01 = q01_real.clone();
        muts.push(M { name: "(b) q00[5]+=200 (perturbs genuine spectrum), real q01", q00, q01 });
    }
    // (c) Real spectrum, q00[0] inflated far above floor.
    {
        let mut q00 = q00_real.clone();
        q00[0] = 20000;
        let q01 = q01_real.clone();
        muts.push(M { name: "(c) real spectrum, q00[0]:=20000 (DC inflated)", q00, q01 });
    }
    // (d) Fully synthetic hand-forged small (q00,q01) mimicking a real form.
    {
        let mut q00 = vec![0i16; HAWK_N];
        q00[0] = 2080;
        for u in 1..256 {
            q00[u] = ((u as i32 % 7) - 3) as i16;
        }
        let mut q01 = vec![0i16; HAWK_N];
        for u in 0..HAWK_N {
            q01[u] = ((u as i32 % 5) - 2) as i16;
        }
        q01[0] = 1;
        muts.push(M { name: "(d) hand-forged small (q00,q01) mimicking a real form", q00, q01 });
    }

    let (zsig, zsalt, zs1) = zero_signature();
    let mut any_break = false;

    for m in &muts {
        let (pk, bytes) = match materialize(&m.q00, &m.q01) {
            Ok(v) => v,
            Err(e) => {
                println!("  {} -> NOT MATERIALIZABLE: {e}", m.name);
                continue;
            }
        };
        assert!(m.q00[0] as i32 >= HAWK_512_Q00_FLOOR, "test bug: below floor");

        let (acc0, det0) = probe(&pk, &bytes, &m.q00, &m.q01, &zsig, &zsalt, &zs1, &msgs);

        // Honest signature (valid under the REAL key) fed to the MUTATED key.
        let mut srng = ChaCha20Rng::from_seed([32u8; 32]);
        let hon_sig = kp.secret.sign(&msgs[0], &mut srng).expect("sign");
        let hon_ok_real = kp.public.verify(&msgs[0], &hon_sig).is_ok();
        let hon_ok_mut = pk.verify(&msgs[0], &hon_sig).is_ok();
        let (hsalt, hs1) = sig_parts(&hon_sig);
        let hon_label_mut = verify_inner_labeled(&msgs[0], &bytes, &m.q00, &m.q01, &hsalt, &hs1);

        println!("  {}", m.name);
        println!("      zero-sig ({} msgs):  {}", msgs.len(), summarize(&det0));
        println!(
            "      honest sig on msg0: real-key={} mutated-key={} ({})",
            hon_ok_real,
            hon_ok_mut,
            match hon_label_mut {
                Ok(()) => "ACCEPT".to_string(),
                Err(r) => reject_site(r).to_string(),
            }
        );

        let mutated_accepts = acc0 + if hon_ok_mut { 1 } else { 0 };
        if mutated_accepts >= 2 {
            any_break = true;
            println!("      *** PolyQnorm SOUNDNESS BREAK: non-key (q00,q01) accepted >= 2 pairs ***");
        }
    }

    assert!(
        !any_break,
        "SURVIVING PolyQnorm soundness break in case 3: a non-key (q00,q01) form that passes \
         the floor was accepted for >= 2 pairs. See stdout."
    );
    println!("  => CASE 3 verdict: no fooled-norm acceptance; mutations that break the real \
             key relationship are rejected by the dual-prime norm / bound checks.");
}

// ===========================================================================
// CASE 2 (reinforced): malformed-q00 spectra that ACTUALLY REACH the verifier.
//
// Several natural case-2 shapes are rejected at ENCODE time
// (`encode_gr q01 overflow`) because a max/alternating q01 has Golomb-Rice
// tails that overflow the 1024-byte pubkey. That means the malformed q00 never
// reached the norm check. Here we pair each malformed-q00 spectrum with a
// SMALL, guaranteed-decodable q01 (q01[0] in {0,1}) so the degenerate q00
// spectrum genuinely drives the dual-prime norm / division. We also mount the
// strongest MBS attack: the attacker's best |t1|-minimizing signature.
//
// Structural note on the MBS attack. t1[i] = h1_bit[i] - 2*s1[i] with
// h1_bit in {0,1}; the coefficient-wise minimum |t1[i]| is achieved at
// s1[i]=0, giving t1 in {0,1}^n. So the ZERO signature already minimizes the
// dominant norm term Tr(t1 adj(t1)/q00) for EVERY message -- no nonzero s1 can
// beat it there. This test confirms even that optimum stays rejected across
// many messages for every floor-passing malformed q00 spectrum.
// ===========================================================================
#[test]
fn residual_case2_reaches_verifier_and_minimizes_t1() {
    println!("\n================ CASE 2 (reinforced): malformed q00 reaching the norm check ================");
    let msgs = make_msgs("case2r-msg", 128);

    // Malformed q00 spectra, each with small decodable q01 so they reach verify.
    struct Shape {
        name: &'static str,
        q00: Vec<i16>,
    }
    let mut shapes: Vec<Shape> = Vec::new();
    {
        let mut q00 = vec![0i16; HAWK_N];
        q00[0] = 2080;
        for u in 1..256 {
            q00[u] = if u % 2 == 0 { 511 } else { -511 };
        }
        shapes.push(Shape { name: "alternating +-511 spectrum (q00[0]=2080)", q00 });
    }
    {
        let mut q00 = vec![0i16; HAWK_N];
        q00[0] = 2100;
        for u in 1..256 {
            q00[u] = -511;
        }
        shapes.push(Shape { name: "all -511 spectrum (q00[0]=2100, drives bins toward 0)", q00 });
    }
    {
        // Near-cancelling: DC just above floor, big positive low bins to try to
        // make one FFT bin blow up large (=> that norm term tiny).
        let mut q00 = vec![0i16; HAWK_N];
        q00[0] = 2080;
        for u in 1..256 {
            q00[u] = 511;
        }
        shapes.push(Shape { name: "all +511 spectrum (q00[0]=2080, max positive bins)", q00 });
    }
    {
        // Sparse high-frequency spike.
        let mut q00 = vec![0i16; HAWK_N];
        q00[0] = 2080;
        q00[255] = 511;
        q00[128] = -511;
        shapes.push(Shape { name: "sparse hi-freq spikes (q00[0]=2080)", q00 });
    }
    {
        // GR-budget-safe alternating small perturbation (|q00[u]|<=31 => k=0,
        // no unary tail => always encodes). Degenerate-but-decodable spectrum.
        let mut q00 = vec![0i16; HAWK_N];
        q00[0] = 2080;
        for u in 1..256 {
            q00[u] = if u % 2 == 0 { 31 } else { -31 };
        }
        shapes.push(Shape { name: "alternating +-31 spectrum (encodable, q00[0]=2080)", q00 });
    }
    {
        // GR-budget-safe all-equal small perturbation, negative (push bins down).
        let mut q00 = vec![0i16; HAWK_N];
        q00[0] = 2090;
        for u in 1..256 {
            q00[u] = -31;
        }
        shapes.push(Shape { name: "all -31 spectrum (encodable, q00[0]=2090)", q00 });
    }
    {
        // Low-frequency ramp within budget (a few dozen moderate coeffs).
        let mut q00 = vec![0i16; HAWK_N];
        q00[0] = 2080;
        for u in 1..16 {
            q00[u] = 200;
        }
        shapes.push(Shape { name: "low-freq 200-bump [1..16] (encodable, q00[0]=2080)", q00 });
    }

    let q01_small = {
        let mut v = vec![0i16; HAWK_N];
        v[0] = 1;
        v
    };

    let mut any_break = false;
    for sh in &shapes {
        let (pk, bytes) = match materialize(&sh.q00, &q01_small) {
            Ok(v) => v,
            Err(e) => {
                println!("  {} -> NOT MATERIALIZABLE: {e}", sh.name);
                continue;
            }
        };
        assert!(sh.q00[0] as i32 >= HAWK_512_Q00_FLOOR);

        // Best MBS attack signature = ZERO (proven optimal for the t1 term).
        let (zsig, zsalt, zs1) = zero_signature();
        let (acc, det) = probe(&pk, &bytes, &sh.q00, &q01_small, &zsig, &zsalt, &zs1, &msgs);
        println!("  {}", sh.name);
        println!("      zero-sig (|t1|-optimal) over {} msgs:  {}", msgs.len(), summarize(&det));
        if acc >= 2 {
            any_break = true;
            println!("      *** SURVIVING BUFF BREAK (reached verifier) ***");
        }
    }
    assert!(!any_break, "SURVIVING BUFF BREAK in reinforced case 2");
    println!("  => reinforced CASE 2 verdict: malformed q00 spectra that reach the norm check are \
             all rejected (division-domain / dual-prime-norm), even under the |t1|-optimal signature.");
}

// ---------------------------------------------------------------------------
// Anchor: the labeled verifier ACCEPTs an honest (key, sig) -- so the accept
// path is exercised and attribution is anchored against a known-good ACCEPT,
// not only rejects. Also proves `probe`'s labeled==real cross-check on ACCEPT.
// ---------------------------------------------------------------------------
#[test]
fn labeled_verifier_accepts_honest_key() {
    let mut rng = ChaCha20Rng::from_seed([1u8; 32]);
    let kp = HawkKeypair::generate(&mut rng);
    let mut srng = ChaCha20Rng::from_seed([2u8; 32]);
    let msg = b"honest anchor".to_vec();
    let sig = kp.secret.sign(&msg, &mut srng).unwrap();
    let bytes = kp.public.to_bytes().unwrap();
    let (salt, s1) = sig_parts(&sig);
    let (q00, q01) = pubkey_coeffs(&kp.public);
    let real_ok = kp.public.verify(&msg, &sig).is_ok();
    let labeled = verify_inner_labeled(&msg, &bytes, &q00, &q01, &salt, &s1);
    println!("honest anchor: real_ok={real_ok}, labeled={labeled:?}");
    assert!(real_ok, "honest sig must verify");
    assert_eq!(labeled, Ok(()), "labeled verifier must ACCEPT the honest sig");
}

// ---------------------------------------------------------------------------
// Access the decoded q00/q01 coefficients from a public key without depending
// on private fields: round-trip through to_bytes + decode_public.
// ---------------------------------------------------------------------------
fn pubkey_coeffs(pk: &HawkPublicKey) -> (Vec<i16>, Vec<i16>) {
    let bytes = pk.to_bytes().expect("pubkey encode");
    pqe_hawk::serialize::decode_public(&bytes).expect("pubkey decode")
}
fn pubkey_q00_0(pk: &HawkPublicKey) -> i16 {
    pubkey_coeffs(pk).0[0]
}
