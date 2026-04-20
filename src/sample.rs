//! HAWK Gaussian sampler.
//!
//! Port of `sig_gauss_alt` from c-reference/hawk-512/hawk_sign.c:536-650.
//! The HAWK signing process samples `2n` small signed integers from a
//! bimodal Gaussian distribution using a CDT (cumulative distribution
//! table). For HAWK-512 (logn=9), `2n = 1024` and each sample fits in
//! `i8`.
//!
//! This is the "alt" variant that consumes RNG bytes directly (rather than
//! deriving them from internal SHAKE instances). The non-alt variant uses
//! four parallel SHAKE-256 instances; both produce identical distributions.

use crate::params::HAWK_N;
use rand::RngCore;

/// CDT high parts for HAWK-512 (15-bit values stored as u16).
/// Source: hawk_sign.c:339-345 (`sig_gauss_hi_Hawk_512`).
const SIG_GAUSS_HI: [u16; 10] = [
    0x580B, 0x35F9, 0x1D34, 0x0DD7, 0x05B7, 0x020C, 0x00A2, 0x002B, 0x000A, 0x0001,
];

/// CDT low parts for HAWK-512 (full 64-bit words).
/// Source: hawk_sign.c:348-362 (`sig_gauss_lo_Hawk_512`).
const SIG_GAUSS_LO: [u64; 26] = [
    0x0C27920A04F8F267,
    0x3C689D9213449DC9,
    0x1C4FF17C204AA058,
    0x7B908C81FCE3524F,
    0x5E63263BE0098FFD,
    0x4EBEFD8FF4F07378,
    0x56AEDFB0876A3BD8,
    0x4628BC6B23887196,
    0x061E21D588CC61CC,
    0x7F769211F07B326F,
    0x2BA568D92EEC18E7,
    0x0668F461693DFF8F,
    0x00CF0F8687D3B009,
    0x001670DB65964485,
    0x000216A0C344EB45,
    0x00002AB6E11C2552,
    0x000002EDF0B98A84,
    0x0000002C253C7E81,
    0x000000023AF3B2E7,
    0x0000000018C14ABF,
    0x0000000000EBCC6A,
    0x000000000007876E,
    0x00000000000034CF,
    0x000000000000013D,
    0x0000000000000006,
    0x0000000000000000,
];

/// Number of samples produced per call: `2 * HAWK_N`.
pub const NUM_SAMPLES: usize = 2 * HAWK_N;

/// Output of [`sample`]: `2n` signed bytes plus the squared norm.
pub struct GaussianSample {
    pub x: [i8; NUM_SAMPLES],
    pub squared_norm: u32,
}

/// Sample `2n` small signed integers from the bimodal HAWK Gaussian
/// distribution, with per-sample parity selected by `t_parity` (a packed
/// bit vector — the `h0` or `h1` from [`crate::hash::compute_h`]).
///
/// Port of `sig_gauss_alt` (hawk_sign.c:536-650) for HAWK-512 (logn=9).
///
/// Returns the sampled vector and its squared norm. The norm is used by
/// the outer signing loop (Task 22) to decide whether to retry.
pub fn sample<R: RngCore>(rng: &mut R, t_parity: &[u8]) -> GaussianSample {
    assert_eq!(
        t_parity.len(),
        NUM_SAMPLES / 8,
        "t_parity must be NUM_SAMPLES/8 = {} bytes for HAWK-512",
        NUM_SAMPLES / 8
    );

    let mut x = [0i8; NUM_SAMPLES];
    let mut sn: u32 = 0;

    let hi_len = SIG_GAUSS_HI.len();
    let lo_len = SIG_GAUSS_LO.len();

    // Each iteration of the outer loop produces 16 samples and consumes 160 RNG bytes.
    // Buffer layout mirrors the C union (hawk_sign.c:587-599):
    //   lo for sample (j,k): bytes [(j*8 + k*32) .. (j*8 + k*32 + 8)]
    //   hi for sample (j,k): bytes [(j*8 + 128 + k*2) .. (j*8 + 128 + k*2 + 2)]
    // This is an interleaved layout of 4 simulated streams (j=0..3), each
    // contributing 4 samples (k=0..3) per group of 16.
    for u in (0..NUM_SAMPLES).step_by(16) {
        let mut buf = [0u8; 160];
        rng.fill_bytes(&mut buf);

        for j in 0..4usize {
            for k in 0..4usize {
                let v = u + (j << 2) + k;

                // Byte offsets matching C: dec64le(buf.b + (j<<3) + (k<<5))
                //                          dec16le(buf.b + (j<<3) + 128 + (k<<1))
                let lo_off = (j << 3) + (k << 5);
                let hi_off = (j << 3) + 128 + (k << 1);
                let lo = u64::from_le_bytes(buf[lo_off..lo_off + 8].try_into().unwrap());
                let hi = u16::from_le_bytes(buf[hi_off..hi_off + 2].try_into().unwrap()) as u32;

                // Extract sign bit from top of lo; broadcast to u32 mask.
                let neg = (lo >> 63).wrapping_neg() as u32;
                let lo = lo & 0x7FFF_FFFF_FFFF_FFFF;
                let hi = hi & 0x7FFF;

                // Parity bit selects D0 (even column) vs D1 (odd column).
                let pbit = ((t_parity[v >> 3] >> (v & 7)) & 1) as u64;
                let p_odd = pbit.wrapping_neg(); // 0xFFFF...FFFF or 0x0000...0000
                let p_oddw = p_odd as u32;

                let mut r: u32 = 0;

                // Hi-table loop: both hi and lo parts contribute to comparison.
                let mut i = 0;
                while i < hi_len {
                    let tlo0 = SIG_GAUSS_LO[i];
                    let tlo1 = SIG_GAUSS_LO[i + 1];
                    let tlo = tlo0 ^ (p_odd & (tlo0 ^ tlo1));
                    let cc = (lo.wrapping_sub(tlo) >> 63) as u32;
                    let thi0 = SIG_GAUSS_HI[i] as u32;
                    let thi1 = SIG_GAUSS_HI[i + 1] as u32;
                    let thi = thi0 ^ (p_oddw & (thi0 ^ thi1));
                    r = r.wrapping_add(hi.wrapping_sub(thi).wrapping_sub(cc) >> 31);
                    i += 2;
                }

                // Lo-only tail: hi must be zero (hinz gates the contribution).
                let hinz = hi.wrapping_sub(1) >> 31; // 1 iff hi == 0
                while i < lo_len {
                    let tlo0 = SIG_GAUSS_LO[i];
                    let tlo1 = SIG_GAUSS_LO[i + 1];
                    let tlo = tlo0 ^ (p_odd & (tlo0 ^ tlo1));
                    let cc = (lo.wrapping_sub(tlo) >> 63) as u32;
                    r = r.wrapping_add(hinz & cc);
                    i += 2;
                }

                // Multiply by 2 and apply parity.
                let r = (r << 1).wrapping_sub(p_oddw);
                // Apply sign bit.
                let r = (r ^ neg).wrapping_sub(neg);

                // Truncate to i8, matching C: x[v] = (int8_t)*(int32_t *)&r
                x[v] = r as i32 as i8;
                sn = sn.wrapping_add(r.wrapping_mul(r));
            }
        }
    }

    GaussianSample {
        x,
        squared_norm: sn,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn sample_terminates() {
        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
        let t = [0u8; NUM_SAMPLES / 8];
        let result = sample(&mut rng, &t);
        // Just check the function runs to completion without panicking.
        assert_eq!(result.x.len(), NUM_SAMPLES);
    }

    #[test]
    fn sample_deterministic_for_same_seed() {
        let t = [0xAAu8; NUM_SAMPLES / 8];
        let mut rng1 = ChaCha20Rng::from_seed([42u8; 32]);
        let mut rng2 = ChaCha20Rng::from_seed([42u8; 32]);
        let r1 = sample(&mut rng1, &t);
        let r2 = sample(&mut rng2, &t);
        assert_eq!(r1.x, r2.x);
        assert_eq!(r1.squared_norm, r2.squared_norm);
    }

    #[test]
    fn samples_are_small() {
        // The bimodal Gaussian has standard deviation ~2*sigma ≈ 2.6.
        // 99.99% of samples should be within |x| <= 12 or so. We just
        // assert no sample exceeds a generous bound (i8 max 127).
        let mut rng = ChaCha20Rng::from_seed([7u8; 32]);
        let t = [0u8; NUM_SAMPLES / 8];
        let r = sample(&mut rng, &t);
        let max_abs = r.x.iter().map(|&v| v.unsigned_abs() as u32).max().unwrap();
        assert!(max_abs <= 30, "unexpectedly large sample: {}", max_abs);
    }

    #[test]
    fn norm_matches_sum_of_squares() {
        let mut rng = ChaCha20Rng::from_seed([1u8; 32]);
        let t = [0u8; NUM_SAMPLES / 8];
        let r = sample(&mut rng, &t);
        let computed: u32 =
            r.x.iter()
                .map(|&v| {
                    let v = v as i32;
                    v.wrapping_mul(v) as u32
                })
                .sum::<u32>();
        // Note: the C accumulates `r * r` where `r` is the u32 (sign-applied)
        // before being truncated to i8. For samples in [-128, 127], i8→i32→
        // squared matches the C u32 arithmetic.
        assert_eq!(r.squared_norm, computed);
    }
}
