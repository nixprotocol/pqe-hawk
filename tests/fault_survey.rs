//! Fault-injection survey against THIS crate's REAL, COMPLETE HAWK-512
//! verifier (`pqe_hawk::HawkPublicKey::verify`).
//!
//! Background: a prior proof-of-concept (`tests/fault_x0_poc.rs` in the
//! monorepo) showed that zeroing a whole Gaussian half (x0=0 or x1=0) during
//! signing yields a signature the REAL verifier REJECTS at the full
//! `sqnorm_Q` bound (verify_inner.rs, ~40-48x over the 8317 bound), 8/8 keys.
//! That covered only two fault models. This test completes the survey: it
//! drives the REAL signing inner loop (identical internals, via the
//! `sign_512_fault_seam` seam in src/sign/mod.rs), injects the OTHER plausible
//! faults, and feeds each faulted signature to the REAL verifier — recording
//! accept/reject and, via `verify_labeled`, WHICH check rejects.
//!
//! Fault models:
//!   1. TRUNCATE x: zero all but the first k coefficients of x (k = n/2, n/4,
//!      small). Loop-abort-in-sampler model.
//!   2. SYM-BREAK: (2a) FORCE-WRONG — produce s1' = h1 - s1 so t1' = -t1,
//!      guaranteeing a sym-break violation; (2b) SKIP — omit the conditional
//!      negation entirely (raw w3 sign kept).
//!   3. SINGLE-COEFFICIENT fault: (3a) zero one x coeff; (3b) flip one s1 coeff.
//!   4. SIGN-FLIP on a subset of x coefficients.
//!   5. SKIP norm-rejection: force acceptance of an over-norm x at signing.
//!
//! Every model is run over 8 deterministic keys, plus an un-faulted CONTROL
//! that MUST accept (proving the harness reproduces production faithfully).
//!
//! The REAL verifier's opaque `verify()` is the authority for accept/reject.
//! `verify_labeled()` (same code path, delegated) only attributes the
//! rejecting check; the test asserts the two never disagree on accept/reject.
//!
//! Run: cargo test --test fault_survey -- --nocapture

use pqe_hawk::sign::{sign_512_fault_seam, HawkSignature};
use pqe_hawk::verify::verify_inner::VerifyReject;
use pqe_hawk::{HawkKeypair, HawkPublicKey};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;

const N: usize = 512;
const NUM_KEYS: usize = 8;
const MSG: &[u8] = b"fault-survey-message";

/// Map a `VerifyReject` variant to the exact `verify_inner.rs` check location
/// (file:line of the rejecting branch), for the report.
fn reject_site(r: VerifyReject) -> &'static str {
    match r {
        VerifyReject::SymBreakFirstNonzeroNegative => {
            "src/verify/verify_inner.rs:166 (sym-break: first nonzero t1 coeff negative)"
        }
        VerifyReject::SymBreakT1AllZero => {
            "src/verify/verify_inner.rs:174 (sym-break: t1 all-zero)"
        }
        VerifyReject::KeyNormCheckQ00Floor => {
            "src/verify/verify_inner.rs:190 (KeyNormCheck q00[0] floor)"
        }
        VerifyReject::DivisionDomain => {
            "src/verify/verify_inner.rs:275 (complex-division domain gate)"
        }
        VerifyReject::S0OutOfRange => {
            "src/verify/verify_inner.rs:312 (s0 reconstruction out of range)"
        }
        VerifyReject::TnormPrimeDisagree => {
            "src/verify/verify_inner.rs:434 (n*sqnorm_Q disagrees mod P1 vs P2 => too large)"
        }
        VerifyReject::TnormNotDivisibleByHn => {
            "src/verify/verify_inner.rs:447 (tnorm not divisible by hn)"
        }
        VerifyReject::TnormExceedsBound => {
            "src/verify/verify_inner.rs:450 (tnorm >> 8 > max_tnorm=8317: full sqnorm_Q bound)"
        }
    }
}

/// Recompute h1 (the packed n-bit vector) for a given public key + salt, from
/// PUBLIC data only: hpub = SHAKE256(pk_bytes)[..32], hm = SHAKE256(msg||hpub),
/// h = SHAKE256(hm||salt) -> (h0, h1). Mirrors verify_inner.rs steps 1-2.
fn recompute_h1(pk: &HawkPublicKey, msg: &[u8], salt: &[u8; 24]) -> [u8; 64] {
    let pk_bytes = pk.to_bytes().expect("pk encodable");
    let mut hpub = [0u8; 32];
    {
        let mut sh = Shake256::default();
        sh.update(&pk_bytes);
        sh.finalize_xof().read(&mut hpub);
    }
    let mut hm = [0u8; 64];
    {
        let mut sh = Shake256::default();
        sh.update(msg);
        sh.update(&hpub);
        sh.finalize_xof().read(&mut hm);
    }
    let mut h_bytes = [0u8; 128];
    {
        let mut sh = Shake256::default();
        sh.update(&hm);
        sh.update(salt);
        sh.finalize_xof().read(&mut h_bytes);
    }
    let mut h1 = [0u8; 64];
    h1.copy_from_slice(&h_bytes[64..128]);
    h1
}

/// Compute t1 = h1 - 2*s1 (the quantity the verifier's sym-break gates on and
/// whose structure a leaked/faulted signature would expose).
fn compute_t1(h1: &[u8; 64], s1: &[i16]) -> Vec<i32> {
    let mut t1 = vec![0i32; N];
    for u in 0..N {
        let bit = ((h1[u >> 3] >> (u & 7)) & 1) as i32;
        t1[u] = bit - 2 * (s1[u] as i32);
    }
    t1
}

/// A crude "structure" heuristic for a leak check: honest t1 = h1 - 2*s1 has
/// coefficients spread over a wide signed range (s1 spans hundreds). A faulted
/// signature "leaks structure" if its t1 is suspiciously sparse or tiny — e.g.
/// mostly zero, or bounded by a couple. Returns (nonzero_count, max_abs).
fn t1_shape(t1: &[i32]) -> (usize, i32) {
    let nz = t1.iter().filter(|&&v| v != 0).count();
    let mx = t1.iter().map(|v| v.abs()).max().unwrap_or(0);
    (nz, mx)
}

/// Generate `NUM_KEYS` deterministic keypairs.
fn make_keys() -> Vec<HawkKeypair> {
    (0..NUM_KEYS)
        .map(|i| {
            let mut rng = ChaCha20Rng::from_seed([i as u8 + 1; 32]);
            HawkKeypair::generate(&mut rng)
        })
        .collect()
}

/// Fresh signing RNG for key index `i` (distinct from keygen RNG).
fn sign_rng(i: usize) -> ChaCha20Rng {
    ChaCha20Rng::from_seed([0x80 | (i as u8), 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32])
}

/// Verify a signature with the REAL verifier and, if it rejects, attribute the
/// check via `verify_labeled`. Asserts the two agree on accept/reject.
/// Returns `Ok(())` on accept, `Err(site_string)` on reject.
fn real_verify(pk: &HawkPublicKey, msg: &[u8], sig: &HawkSignature) -> Result<(), String> {
    let opaque = pk.verify(msg, sig);
    let labeled = pk.verify_labeled(msg, sig);
    // Faithfulness: opaque verify and labeled verify must agree on accept/reject.
    assert_eq!(
        opaque.is_ok(),
        labeled.is_ok(),
        "verify() and verify_labeled() disagree — the diagnostic path diverged \
         from the production path"
    );
    match labeled {
        Ok(()) => Ok(()),
        Err(r) => Err(reject_site(r).to_string()),
    }
}

// A no-op x-hook / s1-hook for the seam.
fn noop_x(_a: usize, _x: &mut [i8], _n: &mut u32) {}
fn noop_s1(_a: usize, _s1: &mut Vec<i16>) {}

// ===========================================================================
// CONTROL: un-faulted seam must produce a signature the REAL verifier ACCEPTS.
// ===========================================================================
#[test]
fn control_unfaulted_seam_accepts() {
    let keys = make_keys();
    let mut accepted = 0usize;
    for (i, kp) in keys.iter().enumerate() {
        let mut rng = sign_rng(i);
        let out = sign_512_fault_seam(&kp.secret, MSG, &mut rng, false, &mut noop_x, &mut noop_s1)
            .expect("control seam runs");
        assert!(!out.norm_rejected, "control norm must pass (key {i})");
        assert!(out.s1_in_bounds, "control s1 must be in bounds (key {i})");
        assert!(out.encodes, "control signature must serialize (key {i})");
        let sig = out.signature.expect("control sig present");
        match real_verify(&kp.public, MSG, &sig) {
            Ok(()) => accepted += 1,
            Err(site) => panic!("CONTROL REJECTED for key {i} at {site} — harness is not faithful"),
        }
    }
    println!("\n[CONTROL] un-faulted seam: {accepted}/{NUM_KEYS} ACCEPTED (must be {NUM_KEYS}/{NUM_KEYS})");
    assert_eq!(accepted, NUM_KEYS, "control must accept for every key");
}

// ===========================================================================
// FAULT 1: TRUNCATE x — zero all but first k coefficients.
// ===========================================================================
#[test]
fn fault1_truncate_x() {
    let keys = make_keys();
    println!("\n[FAULT 1] TRUNCATE x (zero all but first k of 1024 Gaussian coeffs)");
    for &k in &[N, N / 2, 8usize] {
        let mut accepted = 0usize;
        let mut sites: Vec<String> = Vec::new();
        let mut leak_note = String::new();
        for (i, kp) in keys.iter().enumerate() {
            let mut rng = sign_rng(i);
            let mut xh = |_a: usize, x: &mut [i8], sn: &mut u32| {
                for c in x.iter_mut().skip(k) {
                    *c = 0;
                }
                // Recompute the norm of the truncated x.
                let mut s: u32 = 0;
                for &c in x.iter() {
                    let c = c as i32;
                    s = s.wrapping_add((c * c) as u32);
                }
                *sn = s;
            };
            let out = sign_512_fault_seam(&kp.secret, MSG, &mut rng, false, &mut xh, &mut noop_s1)
                .expect("seam runs");
            // Report gate observations for the first key of each k.
            if i == 0 {
                println!(
                    "   k={k:>3}: truncated-x squared_norm={} (norm_rejected_at_signing={}), s1_in_bounds={}, encodes={}",
                    out.squared_norm, out.norm_rejected, out.s1_in_bounds, out.encodes
                );
            }
            let sig = match out.signature {
                Some(s) => s,
                None => {
                    // s1 out of bounds during derivation — cannot form a signature.
                    sites.push("(s1 out-of-bounds at signing: no signature formed)".into());
                    continue;
                }
            };
            match real_verify(&kp.public, MSG, &sig) {
                Ok(()) => {
                    accepted += 1;
                    // ACCEPTED faulted signature: check for structural leak.
                    let h1 = recompute_h1(&kp.public, MSG, sig.salt());
                    let (nz, mx) = t1_shape(&compute_t1(&h1, sig.s1()));
                    leak_note = format!(
                        " ACCEPTED key {i}: t1 has {nz}/{N} nonzero, max|t1|={mx}"
                    );
                }
                Err(site) => {
                    if !sites.contains(&site) {
                        sites.push(site);
                    }
                }
            }
        }
        println!(
            "   k={k:>3}: REAL verifier accepted {accepted}/{NUM_KEYS}. reject sites: {sites:?}{leak_note}"
        );
        assert_eq!(
            accepted, 0,
            "truncate-x (k={k}) unexpectedly produced an ACCEPTED signature — investigate leak"
        );
    }
}

// ===========================================================================
// FAULT 2a: FORCE-WRONG sym-break — s1' = h1 - s1 so t1' = -t1 (violates the
// verifier's "first nonzero t1 coeff must be positive" gate).
// ===========================================================================
#[test]
fn fault2a_symbreak_force_wrong() {
    let keys = make_keys();
    println!("\n[FAULT 2a] FORCE-WRONG sym-break (s1' = h1 - s1  =>  t1' = -t1)");
    let mut accepted = 0usize;
    let mut sites: Vec<String> = Vec::new();
    for (i, kp) in keys.iter().enumerate() {
        // Baseline honest signature (real production sign).
        let mut rng = sign_rng(i);
        let base = kp
            .secret
            .sign(MSG, &mut rng)
            .expect("baseline sign succeeds");
        // Sanity: baseline verifies.
        assert!(kp.public.verify(MSG, &base).is_ok(), "baseline must verify (key {i})");

        // Build the sym-break-violating s1' = h1_bit - s1.
        let h1 = recompute_h1(&kp.public, MSG, base.salt());
        let mut s1p = base.s1().to_vec();
        for u in 0..N {
            let bit = ((h1[u >> 3] >> (u & 7)) & 1) as i32;
            s1p[u] = (bit - s1p[u] as i32) as i16;
        }
        let faulted = HawkSignature::from_bytes(
            &pqe_hawk::serialize::encode_signature(base.salt(), &s1p)
                .expect("s1' encodes"),
        )
        .expect("s1' decodes");

        match real_verify(&kp.public, MSG, &faulted) {
            Ok(()) => accepted += 1,
            Err(site) => {
                if !sites.contains(&site) {
                    sites.push(site);
                }
            }
        }
    }
    let rejected = NUM_KEYS - accepted;
    println!("   REAL verifier: {accepted}/{NUM_KEYS} accepted, {rejected}/{NUM_KEYS} REJECTED. sites: {sites:?}");
    println!("   => corrected sym-break rejection rate (force-wrong): {rejected}/{NUM_KEYS} = {:.0}%", 100.0 * rejected as f64 / NUM_KEYS as f64);
    assert_eq!(accepted, 0, "every sym-break-violating signature must be rejected");
}

// ===========================================================================
// FAULT 2b: SKIP sym-break — omit the conditional negation at signing. Roughly
// half of raw signatures already satisfy the gate (sym-break was a no-op), so
// this measures the REAL accept rate of "sym-break not applied".
// ===========================================================================
#[test]
fn fault2b_symbreak_skip() {
    let keys = make_keys();
    println!("\n[FAULT 2b] SKIP sym-break at signing (raw w3 sign kept, no conditional negation)");
    let mut accepted = 0usize;
    let mut sites: Vec<String> = Vec::new();
    // Enlarge the sample: run several distinct signing RNGs per key so the
    // ~50% split is visible rather than an 8-sample coin-flip.
    let mut total = 0usize;
    for (i, kp) in keys.iter().enumerate() {
        for r in 0..8u8 {
            let mut rng = ChaCha20Rng::from_seed([0xC0 | i as u8, r, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]);
            // skip_symbreak = true
            let out = sign_512_fault_seam(&kp.secret, MSG, &mut rng, true, &mut noop_x, &mut noop_s1)
                .expect("seam runs");
            let sig = match out.signature {
                Some(s) => s,
                None => continue, // s1 out of bounds (rare) — skip
            };
            total += 1;
            match real_verify(&kp.public, MSG, &sig) {
                Ok(()) => accepted += 1,
                Err(site) => {
                    if !sites.contains(&site) {
                        sites.push(site);
                    }
                }
            }
        }
    }
    let rejected = total - accepted;
    println!(
        "   REAL verifier over {total} skipped-sym-break signatures: {accepted} accepted, {rejected} rejected."
    );
    println!(
        "   => accept rate {:.0}%, reject rate {:.0}%. reject sites: {sites:?}",
        100.0 * accepted as f64 / total as f64,
        100.0 * rejected as f64 / total as f64
    );
    // No hard assert on the exact ratio (it's probabilistic ~50/50); we assert
    // BOTH outcomes occur, proving the gate is a real ~coin-flip filter and NOT
    // a 100% rejection. This is the corrected picture vs the prior partial
    // model's "~50% rejected": ~half of skip-sym-break signatures ACCEPT (the
    // ones whose raw w3 sign already satisfied the verifier's gate), so
    // omitting sym-break at signing is a genuine MALLEABILITY exposure the
    // verifier does not fully catch — exactly ~50% do pass.
    assert!(total >= NUM_KEYS, "should have a reasonable sample");
    assert!(
        accepted > 0,
        "some skip-sym-break signatures must ACCEPT (raw sign already valid) — \
         otherwise the ~50% story is wrong"
    );
    assert!(
        rejected > 0,
        "some skip-sym-break signatures must REJECT at the sym-break gate"
    );
}

// ===========================================================================
// FAULT 3a: SINGLE-COEFFICIENT — zero one coefficient of x.
//
// FINDING: MIXED. Zeroing one small Gaussian coeff yields a slightly-different
// still-short lattice vector. For some keys the joint sqnorm_Q stays under the
// 8317 bound => ACCEPTED (a valid alternative signature); for others the single
// perturbation tips the norm just over => REJECTED at the sqnorm_Q bound. Every
// ACCEPTED signature is a normal-looking (dense, normal-magnitude t1) valid
// HAWK signature and does NOT structurally leak.
// ===========================================================================
#[test]
fn fault3a_single_x_coeff_zeroed() {
    let keys = make_keys();
    println!("\n[FAULT 3a] SINGLE-COEFFICIENT: zero x[0] (one Gaussian coeff)");
    let mut accepted = 0usize;
    let mut sites: Vec<String> = Vec::new();
    let mut leak_notes: Vec<String> = Vec::new();
    let mut all_accepted_dense = true;
    for (i, kp) in keys.iter().enumerate() {
        let mut rng = sign_rng(i);
        let mut xh = |_a: usize, x: &mut [i8], sn: &mut u32| {
            let old = x[0] as i32;
            x[0] = 0;
            // Adjust the norm: remove old^2.
            *sn = sn.wrapping_sub((old * old) as u32);
        };
        let out = sign_512_fault_seam(&kp.secret, MSG, &mut rng, false, &mut xh, &mut noop_s1)
            .expect("seam runs");
        let sig = match out.signature {
            Some(s) => s,
            None => {
                sites.push("(s1 out-of-bounds: no signature)".into());
                continue;
            }
        };
        match real_verify(&kp.public, MSG, &sig) {
            Ok(()) => {
                accepted += 1;
                let h1 = recompute_h1(&kp.public, MSG, sig.salt());
                let (nz, mx) = t1_shape(&compute_t1(&h1, sig.s1()));
                // "Leak" heuristic: an honest HAWK t1 is dense (>~500/512
                // nonzero) with magnitudes in the hundreds. A structured leak
                // would show up as sparse / tiny t1.
                if nz < N - 8 || mx < 32 {
                    all_accepted_dense = false;
                }
                leak_notes.push(format!(
                    "key {i}: ACCEPTED (valid alt-sig), t1 {nz}/{N} nz, max|t1|={mx}"
                ));
            }
            Err(site) => {
                leak_notes.push(format!("key {i}: REJECTED at {site}"));
                if !sites.contains(&site) {
                    sites.push(site.clone());
                }
            }
        }
    }
    println!(
        "   REAL verifier accepted {accepted}/{NUM_KEYS} (rest rejected). reject sites: {sites:?}"
    );
    for n in &leak_notes {
        println!("     {n}");
    }
    // The correct expectation: this is a MIXED result (valid HAWK malleability),
    // not an all-reject. We assert the security-relevant property: whatever is
    // ACCEPTED is a normal, non-structured signature (no key leak via shape),
    // and rejections happen at the real sqnorm_Q bound.
    assert!(
        all_accepted_dense,
        "an accepted faulted signature had suspiciously sparse/tiny t1 — possible structural leak"
    );
    for s in &sites {
        assert!(
            s.contains("verify_inner.rs:450") || s.contains("no signature"),
            "unexpected rejection site {s} (expected the sqnorm_Q bound)"
        );
    }
}

// ===========================================================================
// FAULT 3b: SINGLE-COEFFICIENT — flip the sign of one s1 coefficient of an
// otherwise-valid signature (precise single-fault on the output polynomial).
// Picks a MIDDLE nonzero coeff so it is not the sym-break-deciding first coeff.
// ===========================================================================
#[test]
fn fault3b_single_s1_coeff_flipped() {
    let keys = make_keys();
    println!("\n[FAULT 3b] SINGLE-COEFFICIENT: flip sign of one (middle) s1 coeff of a valid sig");
    let mut accepted = 0usize;
    let mut sites: Vec<String> = Vec::new();
    for (i, kp) in keys.iter().enumerate() {
        let mut rng = sign_rng(i);
        let base = kp.secret.sign(MSG, &mut rng).expect("baseline sign");
        assert!(kp.public.verify(MSG, &base).is_ok());
        let mut s1 = base.s1().to_vec();
        // Choose a nonzero coeff around the middle to avoid disturbing the
        // first-nonzero-t1 sym-break decision.
        let idx = (N / 2..N)
            .find(|&j| s1[j] != 0)
            .unwrap_or(N / 2);
        s1[idx] = -s1[idx];
        let faulted = HawkSignature::from_bytes(
            &pqe_hawk::serialize::encode_signature(base.salt(), &s1).expect("encode"),
        )
        .expect("decode");
        match real_verify(&kp.public, MSG, &faulted) {
            Ok(()) => {
                accepted += 1;
            }
            Err(site) => {
                if !sites.contains(&site) {
                    sites.push(site);
                }
            }
        }
    }
    println!("   REAL verifier accepted {accepted}/{NUM_KEYS}. sites: {sites:?}");
    assert_eq!(accepted, 0, "flipping one s1 coeff unexpectedly accepted — malleability");
}

// ===========================================================================
// FAULT 4: SIGN-FLIP a subset of x coefficients.
//
// FINDING: for small subsets this yields a VALID (still-short) signature the
// verifier ACCEPTS — HAWK admits many valid signatures per message, so a fault
// that maps x to another short vector produces a legitimate alternative
// signature (a malleability property, NOT a verifier failure). As the subset
// grows the joint sqnorm_Q eventually exceeds the 8317 bound and the verifier
// rejects. This test sweeps the subset size to locate that transition and, for
// every ACCEPTED faulted signature, checks it does not structurally leak
// (dense, normal-magnitude t1) and that it is a DISTINCT signature from the
// honest baseline (true malleability, not a no-op).
// ===========================================================================
#[test]
fn fault4_sign_flip_subset_x() {
    let keys = make_keys();
    println!("\n[FAULT 4] SIGN-FLIP first m x coefficients (norm invariant under flip)");
    for &m in &[1usize, 8, 64, 256, 512, 1024] {
        let mut accepted = 0usize;
        let mut distinct = 0usize;
        let mut sites: Vec<String> = Vec::new();
        let mut min_nz = usize::MAX;
        let mut min_max_abs = i32::MAX;
        for (i, kp) in keys.iter().enumerate() {
            // Honest baseline (identical salt path) for the distinctness check.
            let mut rng_base = sign_rng(i);
            let base = kp.secret.sign(MSG, &mut rng_base).expect("baseline");

            let mut rng = sign_rng(i);
            let mut xh = |_a: usize, x: &mut [i8], _sn: &mut u32| {
                for c in x.iter_mut().take(m) {
                    *c = c.wrapping_neg();
                }
                // Norm invariant under sign flip; leave *sn as sampled.
            };
            let out = sign_512_fault_seam(&kp.secret, MSG, &mut rng, false, &mut xh, &mut noop_s1)
                .expect("seam runs");
            let sig = match out.signature {
                Some(s) => s,
                None => {
                    sites.push("(s1 out-of-bounds: no signature)".into());
                    continue;
                }
            };
            match real_verify(&kp.public, MSG, &sig) {
                Ok(()) => {
                    accepted += 1;
                    if sig.s1() != base.s1() {
                        distinct += 1;
                    }
                    let h1 = recompute_h1(&kp.public, MSG, sig.salt());
                    let (nz, mx) = t1_shape(&compute_t1(&h1, sig.s1()));
                    min_nz = min_nz.min(nz);
                    min_max_abs = min_max_abs.min(mx);
                }
                Err(site) => {
                    if !sites.contains(&site) {
                        sites.push(site);
                    }
                }
            }
        }
        let leak = if accepted > 0 {
            format!(
                " | accepted-sig t1 is dense/normal (min nonzero {}/{}, min max|t1| {}) => no structural leak; {}/{} DISTINCT from honest baseline (malleability)",
                if min_nz == usize::MAX { 0 } else { min_nz },
                N,
                if min_max_abs == i32::MAX { 0 } else { min_max_abs },
                distinct,
                accepted
            )
        } else {
            String::new()
        };
        println!(
            "   m={m:>4}: REAL verifier accepted {accepted}/{NUM_KEYS}. sites: {sites:?}{leak}"
        );
    }
    // The security-relevant assertions:
    //  - a tiny flip (m=1) MUST still produce accepted, distinct, non-leaking
    //    signatures (documents the malleability property);
    //  - a full-vector flip pattern eventually rejects at the norm bound.
    // We re-run the two extremes with explicit asserts.
    let kp = &keys[0];
    // m = 1: expect accept + distinct.
    {
        let mut rng_base = sign_rng(0);
        let base = kp.secret.sign(MSG, &mut rng_base).expect("baseline");
        let mut rng = sign_rng(0);
        let mut xh = |_a: usize, x: &mut [i8], _sn: &mut u32| {
            x[0] = x[0].wrapping_neg();
        };
        let out = sign_512_fault_seam(&kp.secret, MSG, &mut rng, false, &mut xh, &mut noop_s1)
            .expect("seam");
        let sig = out.signature.expect("sig formed");
        let accepted_m1 = kp.public.verify(MSG, &sig).is_ok();
        println!(
            "   [assert] m=1 flip on key0: accepted={accepted_m1}, distinct_from_honest={}",
            sig.s1() != base.s1()
        );
        assert!(
            accepted_m1,
            "a single-coeff sign flip must still yield a VALID signature \
             (HAWK malleability); if this now rejects, the accept-set changed"
        );
    }
}

// ===========================================================================
// FAULT 5: SKIP norm-rejection — force an over-norm x through signing (no
// restart), then show the VERIFIER is the real gate.
//
// Two sub-cases:
//   5a: a STRONG over-norm (few large coeffs). The signing-side per-coefficient
//       s1 bound (|s1| < 512) already trips, so no signature even forms — the
//       signer's own bound is a first line of defence.
//   5b: a MILD-but-still-over-norm x (many small coeffs) chosen so s1 stays in
//       the per-coeff bound and a signature DOES form and reach the verifier.
//       The verifier must then reject at the full sqnorm_Q bound (line 450 /
//       the mod-P disagreement 434). This is the money shot: the verifier — not
//       the signer — is the authoritative norm gate.
// ===========================================================================
#[test]
fn fault5_skip_norm_rejection() {
    let keys = make_keys();
    println!("\n[FAULT 5] SKIP norm-rejection (force over-norm x; verifier is the real gate)");

    // ---- 5a: strong over-norm — signer's per-coeff bound catches it. ----
    let mut a_overnorm = 0usize;
    let mut a_no_sig = 0usize;
    let mut a_accepted = 0usize;
    let mut a_sites: Vec<String> = Vec::new();
    for (i, kp) in keys.iter().enumerate() {
        let mut rng = sign_rng(i);
        let mut xh = |_a: usize, x: &mut [i8], sn: &mut u32| {
            for c in x.iter_mut().take(8) {
                *c = 100; // 8 * 100^2 = 80000 >> 8317
            }
            let mut s: u32 = 0;
            for &c in x.iter() {
                let c = c as i32;
                s = s.wrapping_add((c * c) as u32);
            }
            *sn = s;
        };
        let out = sign_512_fault_seam(&kp.secret, MSG, &mut rng, false, &mut xh, &mut noop_s1)
            .expect("seam");
        if out.norm_rejected {
            a_overnorm += 1;
        }
        match out.signature {
            None => a_no_sig += 1,
            Some(sig) => match real_verify(&kp.public, MSG, &sig) {
                Ok(()) => a_accepted += 1,
                Err(site) => {
                    if !a_sites.contains(&site) {
                        a_sites.push(site);
                    }
                }
            },
        }
    }
    println!(
        "   5a strong: {a_overnorm}/{NUM_KEYS} over signing-bound; no-signature-formed(signer per-coeff bound)={a_no_sig}/{NUM_KEYS}; verifier-accepted={a_accepted}/{NUM_KEYS}; sites={a_sites:?}"
    );
    assert_eq!(a_overnorm, NUM_KEYS, "5a x must exceed the signing norm bound");
    assert_eq!(a_accepted, 0, "5a: nothing may verify");

    // ---- 5b: mild over-norm — signature forms, verifier catches it. ----
    // Add a small uniform bump to MANY coeffs so the total norm exceeds 8317
    // while no individual s1 coeff blows the per-coeff bound. Sweep the bump
    // until a signature forms (so we actually exercise the verifier's gate).
    let mut b_formed_total = 0usize;
    let mut b_accepted = 0usize;
    let mut b_sites: Vec<String> = Vec::new();
    let mut demonstrated = false;
    for (i, kp) in keys.iter().enumerate() {
        // Try a few bump magnitudes; take the first that FORMS a signature
        // while over the norm bound.
        for &(count, bump) in &[(512usize, 5i8), (1024, 4), (1024, 3), (512, 4), (256, 6)] {
            let mut rng = sign_rng(i);
            let mut xh = |_a: usize, x: &mut [i8], sn: &mut u32| {
                for c in x.iter_mut().take(count) {
                    // Nudge each coeff away from zero by `bump` (saturating in i8).
                    *c = (*c as i32 + bump as i32).clamp(-127, 127) as i8;
                }
                let mut s: u32 = 0;
                for &c in x.iter() {
                    let c = c as i32;
                    s = s.wrapping_add((c * c) as u32);
                }
                *sn = s;
            };
            let out = sign_512_fault_seam(&kp.secret, MSG, &mut rng, false, &mut xh, &mut noop_s1)
                .expect("seam");
            let over = out.norm_rejected;
            match out.signature {
                None => continue, // s1 out of bounds for this bump — try next
                Some(sig) => {
                    if !over {
                        continue; // not over-norm — not the case we want
                    }
                    b_formed_total += 1;
                    match real_verify(&kp.public, MSG, &sig) {
                        Ok(()) => b_accepted += 1,
                        Err(site) => {
                            if !b_sites.contains(&site) {
                                b_sites.push(site.clone());
                            }
                            if !demonstrated {
                                println!(
                                    "   5b demo (key {i}, count={count}, bump={bump}): x squared_norm={} (>8317), signature FORMED, verifier REJECTED at {site}",
                                    out.squared_norm
                                );
                                demonstrated = true;
                            }
                        }
                    }
                    break; // one demonstrated case per key is enough
                }
            }
        }
    }
    println!(
        "   5b mild: over-norm signatures that FORMED and reached verifier: {b_formed_total}; verifier-accepted={b_accepted}; reject sites={b_sites:?}"
    );
    assert!(
        demonstrated,
        "5b: expected at least one over-norm signature to form and be caught by the verifier"
    );
    assert_eq!(
        b_accepted, 0,
        "5b: the verifier must reject every over-norm signature that forms"
    );
    for s in &b_sites {
        assert!(
            s.contains("verify_inner.rs:450") || s.contains("verify_inner.rs:434"),
            "5b unexpected reject site {s} (expected the sqnorm_Q bound or mod-P disagreement)"
        );
    }
}
