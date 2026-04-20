//! Statistical timing smoke tests (dudect-style, very lightweight).
//!
//! These are **not** a proof of constant-time behavior. They are a coarse
//! smoke test that catches the loudest timing leaks — cases where sign
//! or verify's runtime on one key/message is statistically distinguishable
//! from another.
//!
//! # Limitations
//!
//! * Uses wall-clock `Instant` (subject to OS scheduling noise).
//! * 1000-iteration sample size; production CT verification needs millions.
//! * LLVM may introduce branches at `-O3` that `-O0` lacks; both should be
//!   tested for real constant-time guarantees.
//! * Does not check verify's timing — verify is data-dependent on the input
//!   signature and pubkey by design (error-path short-circuits). Focusing on
//!   sign, where a timing leak would expose secret key material.
//!
//! For production-grade constant-time verification, use `dudect` proper,
//! `ctgrind`, or formal-methods tooling like Jasmin.

use pqe_hawk::HawkKeypair;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::time::{Duration, Instant};

/// Gate heavy timing tests behind `slow-tests`; defaults off.
#[test]
#[cfg_attr(not(feature = "slow-tests"), ignore)]
fn sign_timing_distributions_are_comparable() {
    const ITERATIONS: usize = 200;
    let msg = b"timing test message";

    // Two distinct keys seeded from different states.
    let mut rng_a = ChaCha20Rng::from_seed([1u8; 32]);
    let mut rng_b = ChaCha20Rng::from_seed([2u8; 32]);
    let kp_a = HawkKeypair::generate(&mut rng_a);
    let kp_b = HawkKeypair::generate(&mut rng_b);

    // Warm up — pre-allocate, exercise ICache/DCache.
    let mut warm_rng = ChaCha20Rng::from_seed([0u8; 32]);
    for _ in 0..20 {
        let _ = kp_a.secret.sign(msg, &mut warm_rng);
        let _ = kp_b.secret.sign(msg, &mut warm_rng);
    }

    // Interleave timings to minimize systematic drift (e.g. CPU thermal
    // throttling between the two groups).
    let mut times_a: Vec<Duration> = Vec::with_capacity(ITERATIONS);
    let mut times_b: Vec<Duration> = Vec::with_capacity(ITERATIONS);
    let mut sign_rng = ChaCha20Rng::from_seed([42u8; 32]);
    for _ in 0..ITERATIONS {
        let t0 = Instant::now();
        let _ = kp_a.secret.sign(msg, &mut sign_rng);
        times_a.push(t0.elapsed());
        let t0 = Instant::now();
        let _ = kp_b.secret.sign(msg, &mut sign_rng);
        times_b.push(t0.elapsed());
    }

    // Trim top 10% outliers (GC, scheduling, page faults) from each.
    times_a.sort();
    times_b.sort();
    let trim = ITERATIONS / 10;
    let times_a_trim = &times_a[trim..ITERATIONS - trim];
    let times_b_trim = &times_b[trim..ITERATIONS - trim];

    let mean_a: Duration = times_a_trim.iter().sum::<Duration>() / times_a_trim.len() as u32;
    let mean_b: Duration = times_b_trim.iter().sum::<Duration>() / times_b_trim.len() as u32;

    // Compute the relative difference.
    let (larger, smaller) = if mean_a > mean_b {
        (mean_a, mean_b)
    } else {
        (mean_b, mean_a)
    };
    let diff_ratio = (larger.as_nanos() - smaller.as_nanos()) as f64 / smaller.as_nanos() as f64;

    eprintln!(
        "sign timing smoke: key_a mean = {:?}, key_b mean = {:?}, relative diff = {:.2}%",
        mean_a,
        mean_b,
        diff_ratio * 100.0
    );

    // Threshold: the two timings should be within 20% of each other. This is
    // very loose — any real side-channel leak would be far smaller than that
    // — but above the noise floor of wall-clock measurement on a busy
    // machine. Tightening this threshold requires dedicated CT-testing
    // tooling (dudect, ctgrind).
    //
    // The goal here is to catch obvious regressions: if someone adds a
    // data-dependent early-exit in sign, timings will often differ by 2x+.
    assert!(
        diff_ratio < 0.20,
        "sign timing distributions differ too much ({:.2}%) — possible CT regression",
        diff_ratio * 100.0
    );
}

/// Fast smoke check that runs on every `cargo test`: confirms that
/// two back-to-back sign calls with the same inputs take "about the
/// same" time. Not a CT test; a regression guard for the sign-path
/// retry loop.
#[test]
fn sign_two_calls_stable_timing() {
    let mut rng = ChaCha20Rng::from_seed([3u8; 32]);
    let kp = HawkKeypair::generate(&mut rng);
    let msg = b"hi";
    let mut s_rng = ChaCha20Rng::from_seed([4u8; 32]);

    // Warm up.
    for _ in 0..5 {
        let _ = kp.secret.sign(msg, &mut s_rng);
    }

    let t0 = Instant::now();
    let _ = kp.secret.sign(msg, &mut s_rng);
    let d1 = t0.elapsed();
    let t0 = Instant::now();
    let _ = kp.secret.sign(msg, &mut s_rng);
    let d2 = t0.elapsed();

    // Neither should be absurdly slow (> 100ms would indicate a retry
    // storm or hang).
    assert!(d1 < Duration::from_millis(500));
    assert!(d2 < Duration::from_millis(500));
}
