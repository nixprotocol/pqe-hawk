//! Number-Theoretic Transform over Z_q[X]/(X^n+1).
//!
//! Port of `Zq(NTT)` and `Zq(iNTT)` from c-reference/hawk-512/modq.h
//! (specialized for Q=18433). The reference operates on the "Montgomery [1..q]"
//! domain where 0 is represented as q; this module wraps the conversion at the
//! API boundary so that callers see standard `[0, q)` coefficients.
//!
//! Pointwise multiplication uses the formula `tomonty(montymul(a[u], b[u]))`
//! before calling `iNTT`. The formula follows directly from Montgomery-domain
//! semantics: `montymul(a, b)` returns `a*b*R^{-1} mod q` (modq.h:113-116) and
//! `tomonty(·)` multiplies by R (modq.h:127-130), yielding `a*b mod q`
//! in the [1..q] representation that `iNTT` consumes.

use crate::params::{HAWK_LOGN, HAWK_N, HAWK_Q};
use crate::ring::Poly;

// === Constants from modq.h Q=18433 branch (lines 25-32) ===
// Q is HAWK_Q from params.rs; Q0I, R2 are at modq.h:26, 28.
const Q: u32 = HAWK_Q; // 18433
const Q0I: u32 = 3955247103; // -1/q mod 2^32; modq.h:26
const R2: u32 = 806; // 2^64 mod q; modq.h:28

// === Modular helpers ===

/// Modular addition in [1, q]: computes `(x + y) mod q` with result in `[1, q]`.
///
/// Port of `Zq(add)` — modq.h:72–77.
///
/// The C code:
///   x = Q - (x + y);        // may wrap if x+y > Q
///   x += Q & (x >> 16);     // if wrapped (x has high bits), add Q to correct
///   return Q - x;
/// Since Q < 2^15, Q & (x >> 16) is Q or a sub-part — but any non-zero high
/// bits mean x wrapped (x+y > Q), and the mask brings x back in range.
#[inline(always)]
fn mq_add(x: u32, y: u32) -> u32 {
    let t = Q.wrapping_sub(x.wrapping_add(y));
    // t >> 16: non-zero iff x+y > Q (wrapped); Q & that propagates correction
    let t = t.wrapping_add(Q & (t >> 16));
    Q.wrapping_sub(t)
}

/// Modular subtraction: returns `(x - y) mod q`, result in [1, q].
/// Port of `Zq(sub)` from modq.h:80-86.
#[inline(always)]
pub(crate) fn mq_sub(x: u32, y: u32) -> u32 {
    let t = y.wrapping_sub(x);
    let t = t.wrapping_add(Q & (t >> 16));
    Q.wrapping_sub(t)
}

/// Halving in [1, q]: computes `x / 2 mod q` with result in `[1, q]`.
///
/// Port of `Zq(half)` — modq.h:97–100.
///
/// `(Q + 1) >> 1` is the modular inverse of 2 mod q (since q is odd).
#[inline(always)]
fn mq_half(x: u32) -> u32 {
    (x >> 1).wrapping_add(((Q.wrapping_add(1)) >> 1) & (x & 1).wrapping_neg())
}

/// Montgomery reduction: given `x` in `[1, m]` (m = 2^32 - (2^16-1)*(q-1)),
/// returns `x * R^{-1} mod q` in `[1, q]`.
///
/// Port of `Zq(montyred)` — modq.h:104–109.
///
/// Algorithm:
///   x  <- x * Q0I              (mod 2^32, gives the Montgomery correction factor)
///   x  <- (x >> 16) * Q
///   return (x >> 16) + 1
#[inline(always)]
fn mq_montyred(x: u32) -> u32 {
    let x = x.wrapping_mul(Q0I);
    let x = (x >> 16).wrapping_mul(Q);
    (x >> 16).wrapping_add(1)
}

/// Montgomery multiplication: `x * y * R^{-1} mod q` in `[1, q]`.
///
/// Port of `Zq(montymul)` — modq.h:113–116.
#[inline(always)]
pub(crate) fn mq_montymul(x: u32, y: u32) -> u32 {
    mq_montyred(x.wrapping_mul(y))
}

/// Convert `x` (in `[1, q]`) to Montgomery form: `x * R mod q` in `[1, q]`.
///
/// Port of `Zq(tomonty)` — modq.h:127–130.
/// Uses `montyred(x * R2)` = `x * R2 * R^{-1} = x * R` (since R2 = R^2 mod q).
#[inline(always)]
pub(crate) fn mq_tomonty(x: u32) -> u32 {
    mq_montyred(x.wrapping_mul(R2))
}

/// Convert from `[1, q]` representation back to `[0, q-1]`.
///
/// Port of `Zq(unorm)` — modq.h:165–168.
///
/// `x & ((x - q) >> 16)`: when `x == q`, `x - q == 0` so the mask is 0 → returns 0;
/// when `x < q`, `x - q` underflows (bit 16 is set) so the mask is all-ones → returns x.
#[inline(always)]
fn mq_unorm(x: u32) -> u32 {
    x & ((x.wrapping_sub(Q)) >> 16)
}

/// Like `mq_set` but for a small signed input (|x| < q/2).
///
/// Returns the value in `[1, q]` (Montgomery-[1..q] representation).
/// Port of `Zq(set_small)` (modq.h:150-157).
#[inline]
pub fn mq_set_small(x: i32) -> u32 {
    let mut y = (x.wrapping_neg()) as u32;
    y = y.wrapping_add(HAWK_Q & (y >> 16));
    HAWK_Q.wrapping_sub(y)
}

/// Convert a value in `[1, q]` (Montgomery representation) to a signed integer
/// in `[-(q-1)/2, +(q-1)/2]`.
///
/// Port of `Zq(snorm)` (modq.h:174-180).
#[inline]
pub fn mq_snorm(x: u32) -> i32 {
    let x = x.wrapping_sub(HAWK_Q & (((HAWK_Q >> 1).wrapping_sub(x)) >> 16));
    x as i32
}

/// Convert a polynomial of signed 8-bit coefficients to the Montgomery-[1..q]
/// representation used by the q=18433 NTT.
///
/// `d` and `a` must both have length `1 << logn`. Port of
/// `Zq(poly_set_small)` (modq.h:831-839).
pub fn mq_poly_set_small(d: &mut [u16], a: &[i8]) {
    debug_assert_eq!(d.len(), a.len());
    for (di, &ai) in d.iter_mut().zip(a.iter()) {
        *di = mq_set_small(ai as i32) as u16;
    }
}

/// In-place conversion: the first n bytes (i.e. n/2 u16 words) of d contain
/// signed 8-bit coefficients packed little-endian. After this call, d[0..n]
/// holds the corresponding mod-q representation.
///
/// Port of `Zq(poly_set_small_inplace_low)` (modq.h:846-860). The loop
/// iterates from the top down so that earlier writes don't clobber later
/// reads.
pub fn mq_poly_set_small_inplace_low(logn: u32, d: &mut [u16]) {
    let n = 1usize << logn;
    debug_assert!(d.len() >= n);
    let mut u = n;
    while u > 0 {
        u -= 2;
        let x = d[u >> 1];
        let x0 = x as u8 as i8;
        let x1 = (x >> 8) as u8 as i8;
        d[u] = mq_set_small(x0 as i32) as u16;
        d[u + 1] = mq_set_small(x1 as i32) as u16;
    }
}

/// In-place conversion: the last n bytes (i.e. d[n/2..n] interpreted as u8
/// bytes) contain signed 8-bit coefficients; after this call, d[0..n] holds
/// the mod-q representation.
///
/// Port of `Zq(poly_set_small_inplace_high)` (modq.h:867-880). Loop bottom-up.
pub fn mq_poly_set_small_inplace_high(logn: u32, d: &mut [u16]) {
    let n = 1usize << logn;
    let hn = n >> 1;
    debug_assert!(d.len() >= n);
    for u in (0..n).step_by(2) {
        let x = d[hn + (u >> 1)];
        let x0 = x as u8 as i8;
        let x1 = (x >> 8) as u8 as i8;
        d[u] = mq_set_small(x0 as i32) as u16;
        d[u + 1] = mq_set_small(x1 as i32) as u16;
    }
}

/// Convert a polynomial in Montgomery-[1..q] representation to signed
/// coefficients in `[-(q-1)/2, +(q-1)/2]` (stored as u16 via two's complement).
///
/// Port of `Zq(poly_snorm)` (modq.h:882-890).
pub fn mq_poly_snorm(d: &mut [u16]) {
    for x in d.iter_mut() {
        *x = mq_snorm(*x as u32) as u16;
    }
}

// === Domain conversion ===
// The [1..q] representation convention (0 → q) is described in modq.h:50-68.

/// Convert from Poly's [0, q) representation to the [1..q] representation
/// used by mq_ntt/mq_intt internally. The integer 0 is represented as q.
#[inline(always)]
fn lift_to_mq(c: u16) -> u16 {
    if c == 0 {
        Q as u16
    } else {
        c
    }
}

/// Convert from [1..q] back to [0, q) via mq_unorm (modq.h:165-168).
#[inline(always)]
fn lift_from_mq(c: u16) -> u16 {
    mq_unorm(c as u32) as u16
}

// === Twiddle tables ===

/// Forward twiddle factors: `GM[x] = (2^32) * (g^rev(x)) mod q`, where g is a
/// 2048th root of unity mod q (g = 19 for q = 18433).
///
/// Verbatim from c-reference/hawk-512/modq.h Zq(GM) for Q=18433 (lines 379–510).
#[rustfmt::skip]
const MQ_GM: [u16; 1024] = [
     4564, 17110, 12162, 16208, 10701,  9705,  3451,  5078,
    12400, 10202,  8245, 13131,  4631,  3492, 17179,  5622,
     5537,  3399,  2485,  9938,   345, 14064, 10152,   789,
     5092, 15713, 12632,  6516, 16107,  2314, 15385, 17281,
      383,  5515,  5019, 13218,  3293,  4728,  9704, 14263,
     4417,   218, 16011,  2568,  5635,  8516, 18352, 12887,
    13102, 15257, 14316, 12813,  8886, 11051, 13356, 15353,
    12059,  6880, 17926, 11710,  8052,  1737, 16384, 18094,
     3410, 14787, 13788, 14210,  2656, 17550,  7950,  4311,
    18150,  4973, 11548,  7848, 15326, 15517,    97, 11648,
    17990, 17685, 10847, 14695,   282,  1558,  6535, 10743,
     2399,   181, 10165,  8051, 12204, 18401, 13377,  7233,
    10892, 15728, 15002, 11766,  8462, 15245, 12420,  8613,
     5053, 12360, 17415, 12678,  4606,   870,  8429,  9572,
     6542,  1892,  4008, 17045,   371, 10155,   819, 15114,
     1522, 13638, 16576, 17586, 16840,  7671, 13873, 12065,
     7433,  7599,  2497,  5298, 16406,  3443,  9437,  6905,
    14589, 17851,   209, 17496,  1698,  7028,  4444,  8211,
     9159, 16089, 16741,  9085,  2658,  4488,  8650,  3995,
     6532, 11903,   508,   192,  4039, 17347, 12742,  6993,
    14812, 17645,  4527,   695,  8380, 16230,  2153,  3136,
     2133,  4725,  9230, 13213,  6548, 18005,  6108, 16097,
    16952, 13519, 16207, 12802, 16047,  7081, 12818,  8328,
    17091,  8927,  9558,  9273,  9301, 10337, 11142,  5082,
    11846, 15508, 17108,  8498, 16135,  3776,  6752, 12857,
     3590,   486,  3056,  4203, 10364, 17125, 14532,  3025,
    18386, 12029,  1983,  7426, 13553,   623,  6269, 15287,
    16399, 12294,  6987,  8011,  1378, 14019,  3042,  3472,
     4734, 12820, 16363,  7781,  3644, 16472,  3523, 14104,
    11521, 18288, 13956,  4549,  6314, 16320, 16373, 16203,
    15299,  7524,  9080, 15914, 17765, 12520,  5829, 13379,
     9482,  7938,   760, 13350,  9526, 15502, 16160,  6398,
     7067,  1655,  3428,  7827, 10564,  1235, 10800,  8291,
    15614, 14755,  8732,  3010, 12821,  7168,  8131,  1912,
     8093, 10461, 12301, 11616, 13947,  8029, 15138,  8334,
    13345, 13462,  7201, 11285,  8232,  5869,  5652,  8087,
     9232,   151,  5425, 15984,  9061, 10972,   874,  6136,
     9299,  4966, 10442,  5398,  6605, 14398,  7625,  7091,
    10974, 14743,  6836, 17243,   504,  7883, 10503, 12533,
     3111, 13658,  1303,  6153, 12791,   335, 16064,  6652,
    14432, 10970,   558,  5436,   300, 13031, 12835,  7899,
     8435,  7252,  2970, 12879,  2786, 16438, 16584,  2204,
     5974,  6467,  7971, 14624,  9637,  9448, 18144,  7293,
    18121, 10042,  1398, 12430,   157,  6881, 18084, 12060,
     5783,   444, 14853,  7936, 13337, 10411,  4401, 12549,
    17699,  1174,  1162,  5374,  3796,   709,  1424,  8521,
     2238,   991,  9114, 15056,  4900, 16221,   731, 18419,
    14885,  1707, 11644,  7594, 14783,  4281, 12810,  5277,
    10528, 15155, 16633, 13979,  5573,  7912, 15085,  4250,
    13224, 11094,  1717, 11970,  4689, 11787,   613, 14891,
     6892,  1734, 15910, 17044,  1022, 16497,  7473,  4421,
    17061,  2094, 17491, 14013,  1872, 13480, 10045, 17585,
     1562, 10460, 12143, 11266,  2168, 15769,  3047,  7683,
    14260,  9889, 14090, 14179,  4404, 11389, 11461,  4622,
    18349, 14047,  7466, 13272,  5005, 12487,   615,  1829,
    13229, 15305,  3467, 11180,  2855,  8191,  3868,  9735,
    18383, 13189,   933,  7900, 18340, 17527,  4316, 14694,
     5680,  9549, 15669,  5777, 17938,  7070, 11080,  4478,
     1466, 10714, 15409,  8001,  7888,  3707, 14283,  7140,
     3046, 14214, 15419, 16423, 18200, 10217, 10615, 18381,
    13138,  1337,  8483,  7125,  6741, 10966, 18359,  4036,
    11656,  2954,  5907,  1652, 12095, 11393, 12093,  6022,
    11472,  6513, 15239, 12291, 16914,  3635,  2907,   373,
    12897,  8503, 16298,  8337, 10348, 11023,  8932,  5553,
    12984, 11729,  9882, 13024,   556,    65, 10270,  4317,
    14404,  9508,  9191,  9860, 14257, 11049, 13040, 14653,
    13038,  9282, 10349,  4492,  6555,  9154,  8558, 14991,
     4583,  3619,   379, 13206, 11105,  7100, 15820, 14978,
     7277, 12620,  3196, 11513,  7268, 16100,    46, 12935,
    10191,  4142,  9281, 11926, 14900, 14340, 16894,  5224,
     9309, 13388, 13942,  3818,  2937,  7206, 14135, 15212,
     7925,  1689,  8800,  1294,  5524, 14570, 16368, 11992,
     9491,  4458,  3910, 11928, 13598,  1656,  3586,  8177,
    13056,  2322, 16649,  1648, 14699, 18328,  1843,   116,
    10016,  4221,  3330,  2710,  5358, 11169, 13567,  1354,
     8715,  3439,  8805,  5505, 10680, 17825, 14534,  8396,
     4185,  3904,  8543,  2358, 13314, 13160, 14784, 16183,
     3842, 13644, 17524,  1253, 13782, 16530, 12687, 15971,
    13700, 17515,  2420, 10494,  7049,  8615, 15561, 10671,
    10485,  1060,  1583,  2340,  6599, 16718,  5525,  8039,
    12196, 15350, 10577,  8497, 16786, 10118, 13406,  2164,
      696,  7375,  3971,   630, 13829,  4501, 10704,  8545,
     8124, 10763,  4718,  6718, 13636, 11540, 16886,  2173,
    13510,  4961,  9652,  3648,  3009, 16232,  2469,  3836,
     4933,  3461, 12281, 13205, 11756, 13442,  4041,  4285,
     3661, 16043,  9473, 11418, 13814, 10301,  5454, 10915,
     8727, 17232, 13005,  3609,  9965,  5508,  3913, 10768,
    11368,  3716, 15705, 10290, 10822, 12073,  8935,  4393,
     3878, 18157, 11691, 13998, 11637, 16445, 17690,  4654,
    12911,  9234,  2765,  6125, 12586, 12014, 18046,  2176,
    17540,  7355,   811, 12063, 17878, 11837,  8513, 13958,
    16653, 12390,  3722,  4745,  7749,  8299,  2499, 10669,
    16214,  3951, 15969,   375, 13937, 18040, 11638,  9914,
    16136, 15678,  7102, 12699,  9368, 15152, 16159, 12929,
    14186, 13925,  6623,  7438,  5741, 16684,   153, 14572,
    14261,  3358, 14440, 14021, 15097, 18043, 12112, 10964,
     5242, 13012,  9833,  1249, 16386,  5032,  2437, 10065,
     1738,  3850,    11,  1891,  3970,  7161,  7025, 17895,
     6303, 14429, 12523, 17941,  6931,  5087, 11127, 10882,
    13926, 16149,  7788, 11652,  8944,   913, 15223,  6189,
     9511,  2869, 10910,  8768,  6262,  5705, 16606,  5986,
    10784,  2189, 14068, 10397, 14897, 15500, 15844,  5698,
     5743,  3622,   853, 14256,  9576,  2313, 15227, 16931,
     3810,  1440,  6324,  6309,  3400,  6365, 10288, 15790,
    16146,  5667, 10602, 11119,  5700,  7960,  4236,  2617,
    12801,  8757,  1131,  5072, 16068, 17394,  1735,  5010,
     2908, 12275,  3985,  1361, 17206, 13615, 12942,  9536,
    12505,  6468,  8129, 14974,  2983,  1708, 11802,  7944,
    17712,  8436,  5712,  3320, 13774, 13479,  9887, 17235,
     4487,  3873,  3645,  9941, 16825, 13471,  8623, 14435,
     5656,   396,  7269,  9569,   935, 13271, 13889, 18167,
     6320, 14000,    40, 15255,  4382,  7607,  3761,  8098,
    15702, 11450,  2666,  7539, 13722,  2864, 10120,  7018,
    11627,  8023, 14190,  6234, 15359,  2757, 11647,  6434,
     1917, 14513,  7362, 10475,   985,    82, 12956, 10267,
    10798,  2920,   535,  8185, 17135, 16491,  6525,  2321,
    11245, 14410,  9521, 11291,  4326,  4683,  2594, 16946,
    12878,  3561,  9648, 11339,  9944, 13628, 14996, 14086,
    16837,  8831, 12823, 12539,  2930, 16057, 11685, 16318,
    11722, 14300, 10574,  9657, 17379,  8165, 18193,   635,
    17483, 10962, 17727,  2636, 16666,  1219,  8272,  2691,
    15755, 15534,  2783, 17598,  9028,  5299,  7757, 11350,
     9421,   803, 16276,  4555,  2408, 15134, 13315,  6629,
     2575, 12004, 16466, 17109, 14006,  9793, 17355, 17445,
     9993,  6970, 13713,  6344, 17481,  5591, 17027,  2952,
      268,   827,  1635, 12955,  8609, 13704,  8571,  3820,
    15205, 13149, 13046, 12333,  8005, 13766, 18367,  7087,
     5414, 14093, 14734, 10939, 12282,  6674,  3811, 13342,
];

/// Inverse twiddle factors: `iGM[x] = (2^31) * ((1/g)^rev(x)) mod q`.
/// Note: scaled by 2^31 (not 2^32) because the `half()` in `iNTT` absorbs the
/// logn-fold factor, giving an inverse without explicit 1/n scaling.
///
/// Verbatim from c-reference/hawk-512/modq.h Zq(iGM) for Q=18433 (lines 647–779).
#[rustfmt::skip]
const MQ_IGM: [u16; 1024] = [
     2282,  9878, 10329, 12352, 15894,  7491,  4364,  3866,
    15622,   627, 16687,  6901,  2651,  5094, 13332, 12233,
      576,  1524, 17276,  1163, 15175, 12117,  1360, 15887,
     8822, 13357, 11401,  9044, 13464,  7974,  7517,  6448,
     9386, 10241,  8348, 14407, 12578,  9470, 14993,  3187,
     1540, 11755,  3691, 13990,  2810, 11275,  1588, 11882,
     2773,  9257, 14175,  6399, 17149,  1211, 18324,  7008,
     2085, 13581, 16069,  7570, 11824,  6707,  6459,  9025,
     3184,  2280,  5381, 10013,  9640, 10145, 11614, 17672,
    10876,  8807,  4139,  9031,   694, 16429, 17487, 15162,
    13647,  5002, 17998, 16130, 12094,   509, 12253,  6690,
     4910, 12223,  1594, 14202, 12550, 10932, 10569, 12987,
     5600,  2528,    16, 12331,  5191,  4134,  9126,  8017,
     3845,  5949, 17654, 18292,  1869,  3793,   374,  9438,
    12609,  9168,  1458, 10770, 14509, 12659,  6730,  9358,
     7061, 14458,  9658, 17105, 11328, 11539,  1823, 16728,
    15234, 10353, 10682, 13670, 11758, 18053, 14464, 13692,
     2527,  6302, 12173,   334, 10476, 13893, 14671,  1567,
     1115,  1030, 10273, 15276,  6942, 11455,  9289,  3456,
    11381,  7455, 10197, 16611,  5326,  1035, 12023, 16066,
    16697, 16912,  2207, 17744,  5211,  5723, 12286,  1017,
     1573,  6082,  8905,  2440, 14720,  8225,  3202,  9240,
     7704, 11167,   654, 13251,  7115, 16905, 18190, 16638,
     2788, 15057, 16545,  1149, 14184,  9879, 10679, 12510,
    15892, 12862,  4048,  4566,  4580, 13654,  4753,   671,
    14269, 12024,  5676,  1193, 12032,  1113,  2457,  9957,
     1168, 15379,   214, 15159,  2610, 13818,  6854,  8150,
    16865,  8140, 10318, 14243,  8869,  6953,   394, 11027,
     5720, 12062,   543,  7197, 18337, 18179,  3265, 15167,
     7219, 14108, 16189, 17104,  4674,   846,  1172,  4637,
     5111, 16211, 14919, 17584,  9685,  9112,   291,  1922,
     5764,  4498,  7495, 10230, 15784,  7968,  5417,  5500,
     6440, 13967,  3705, 13259,  5048, 10284,  4965,  2768,
     9030,  7763,  7399,  9976,  3071,  1597,  5960, 12697,
    15422,  3170,  3520,  3169, 17607,  6263, 16956, 12605,
    16415,    37, 12950,  5846,  5654,  4975,  8548, 11864,
       26,  3909,  4108,  9333,  1005,  1507, 11326, 16910,
    14863,  2075,  7363, 14489,  5216,  1512, 13076, 17700,
    16194, 12893, 14898,  9464,  6328,  1382,  4442, 15593,
    11086, 16275,   453,  9263, 14483,  8750,  2622,    25,
     4349, 16499,  5121,  7789, 12843,  7483,  1564,  2602,
     8302,  8909,  2973,  6714, 11797, 14700,  2193,    42,
    16122,  3486,  3522, 16231,  2127, 11388,  4272, 11303,
     5375,  7693,  1332, 17349, 12800,  3145, 13203, 17652,
      424,  4194, 11693, 17497,  2210,   471, 17386,   686,
     7006,  5480,   968, 17922,  9911, 10478, 17566, 14987,
     1771,  8910,  3323,  6872, 12448,  8358, 12886, 11821,
    16308,  1674, 14477,  6430,  2227,   900,  1639, 13169,
     6578, 12028,  7076,  1825, 14636, 12611,  8363,  1774,
        7,  8851,  1106, 15983, 10905, 13876,  8721, 17314,
     4956, 17721,  8862, 16535, 15746, 17852, 17846,   367,
     2942,  7016,  4011,  2548, 14465,  1790, 18211,  6325,
    12403,  9391,  5776,  9138, 12218, 17734, 13412,   156,
     5570,  9361, 13709,  4398, 11121,  5231,  5983, 15446,
    17331, 10141, 10214, 17040,  2777, 16948, 14807,  4999,
     5267,  2799,  2701, 18283, 15715, 18154, 12948, 11217,
    15107, 10401,  9049,  2821,  6140,  8565, 11604,  7661,
     2950,  3965,  5275, 18181,   595, 15015,  1845, 12946,
     5671,  5404, 11234,  5914, 15734, 13212, 15950,  4567,
    15365, 17996, 12947,  4686, 10441,  6504,  9141, 13817,
     5173, 15607,  6282, 14317,  3574,  5616, 11702,  2544,
    14266, 10864,  5202,  2243, 12625,  3066,  3986,  5170,
    17477,  5151, 14849,  2806, 16928, 14067,  1839, 10626,
     5071, 13033,  8599, 13151,  5303, 16719,  8389,  5683,
    11762,  7311, 15096, 12292,  3747, 11066,  2170, 15726,
     5673,    33, 11550,  5214,  3050, 11910,  2642,  1614,
    16523,  4931, 11581,  4912,  2739,  8399,  8803, 18299,
    16957,   703,  6421,   476, 15261,  2360, 14948,  4220,
      494,   539,  4320, 11430,   662, 10200, 12431,  7929,
     5902,  2559, 10866, 17229,  6939, 10295,  8815,  4506,
    12758,  5338,  6567, 13919,  9634,  7825, 10666,  1339,
     7871, 14297,  8607, 10100, 17115,   353, 12952,   475,
     8899,   120,  5134,   527,  4388, 13146, 11283, 12572,
    10274,  3374,  1188, 16968,  2947,  2805,  4801,   798,
    11390, 10935, 11619, 13461,  3547, 13609,  7436, 11994,
     9960, 17136,  6875, 16270,  3571,  4456, 11228,  3594,
     8056,  5954,   971,   649,  5124,  8949, 16973, 13034,
     4083, 11955, 18392,  8724,  3979, 14752,  1960,  8258,
    15216,  3393,  7838,  1537, 15316, 11338,  5205,  3403,
    14924, 13373, 17001, 11572,  5447, 17100, 12708, 10582,
    14384,  7336,  5413, 16242,  1589, 18413, 11433, 15273,
      133,  2272,  2581,  8749,  4432,  5582, 18235, 15605,
     1999,  4905,  2481,   804,  4246,  7394,  7280,  6973,
      599,  4273,  2477, 11546, 16773, 15577, 14215,  9577,
    14461, 12532, 17579,  7725, 10946,  5152, 15199,  2964,
    13665, 11962,  2409,  9830,  8536,  7224,  3079, 16979,
    15928,  8349,  9736, 10399, 15897,  8651,  4838,  2816,
     7908, 16315, 14453, 15583,  3657, 13132,  6383, 10360,
    10538, 13289,  6034, 16733,  6062, 15271, 17713, 16528,
      751,  1603,  8060, 13645, 11305,  8790, 16622,  6345,
    15584, 10511, 10683,  1768,  4018, 11399,  8122, 13041,
    15440, 10130,  6364, 15302, 14049, 12978,  7782,  4461,
     6122,  1605,  8760, 13961, 12607, 14539,  1142, 11470,
    12992,  3653,  6673,  5751,   246,  2955,  2002,  6065,
      269,  5704,  5636, 16448,  8271,  9211, 16508, 17564,
     4184,  7998, 15917, 10240,  8592,  4300, 11927, 15812,
    12951, 12377,   195,  1668,  2206, 11213, 16754,  2086,
    11147,  9140, 10091,  6346, 14714,  5905,  2254, 11340,
     2752,  1137, 10857, 13749,  2867, 14882, 10594, 10365,
    13476, 12614,  9413,  2248,  9029,  1232,  7241, 10326,
     3882,  7967,  5067,  5342,  6844, 16572, 12238,   890,
    11454,  4960,  3298,  9494,  3185,  8811,  5539,  9663,
    17345,  9410, 12426, 12140,  6154,  7834, 13816,  2761,
    16106,  9588,   994,  3398, 11434,  3371,   138, 16494,
     7020,  4749,  3180, 13022, 13288,  1364, 16575, 12749,
    13049,  7260, 15679,  4234,  7412,  2714,  9817,  4853,
     3759, 15706,  4066, 11526, 12724,  4480,  1195,  7386,
     7074,  7196, 11712, 12555,  2614,  3076,  7486,  6750,
    16515,  7982, 10317,  7712, 16609, 13607,  6736, 11678,
     8130,  9990, 12663, 11615, 15074, 16074,  3835, 14371,
     4944, 13081,  6966,  2302, 18118,  7231,  5529, 18085,
    17351, 11730, 13374, 10040,  4968,  3928, 10758, 12335,
     5197,  6454, 10074,  5917, 17263,  8425, 17903,  3974,
     3881,  1436,  4909,  5692, 13186, 17223,   459, 11583,
     1231,  2873, 10168, 11542,  8590,  9671, 11611, 16512,
     1125, 11041, 11853, 11776, 17254,  4945, 16481,  7124,
    14235, 11166,   304, 13093,  6464,  4814,  7497,  4859,
    17756,  2433,  3632, 15754, 17078, 16768,  7106, 13425,
    18375,  8295,  9269,  1867, 17609,   892, 17272, 11905,
     5128, 16640, 17605, 11634, 12469, 16478, 16204,  4471,
    12437, 10249, 11148, 15671, 17786, 14033,  8372,  5254,
    10827,  2149, 14830,  7748, 16524, 11462, 11739,  4562,
    15821,  9986, 11263, 10983, 12470,  4576, 16362,  4121,
     2749, 18410, 10383, 14799,  3460, 16835, 12123,  5578,
    10944, 10523, 14883,  3664, 11830,  9027,  7407,  6925,
     1721, 14154, 13856,  5939, 16187,  4042, 13792, 11914,
     1890, 11913,  3692,  2088, 13503,  4621, 13679, 11231,
     7058, 13298,  9184, 18155, 11921, 13492,  3352, 11941,
];

// === Core in-place NTT / iNTT (operate on [1, q] domain) ===

/// Forward NTT in-place on a buffer of n = 2^logn coefficients in [1, q].
///
/// Port of `Zq(NTT)` — modq.h:783–804 (Q=18433 specialization).
pub(crate) fn mq_ntt(logn: usize, a: &mut [u16]) {
    let n = 1usize << logn;
    let mut t = n;
    for lm in 0..logn {
        let m = 1usize << lm;
        let ht = t >> 1;
        let mut v0 = 0usize;
        for u in 0..m {
            let s = MQ_GM[u + m] as u32;
            for v in 0..ht {
                let k1 = v0 + v;
                let k2 = k1 + ht;
                let x1 = a[k1] as u32;
                let x2 = mq_montymul(a[k2] as u32, s);
                a[k1] = mq_add(x1, x2) as u16;
                a[k2] = mq_sub(x1, x2) as u16;
            }
            v0 += t;
        }
        t = ht;
    }
}

/// Inverse NTT in-place on a buffer of n = 2^logn coefficients in [1, q].
///
/// Port of `Zq(iNTT)` — modq.h:808–829 (Q=18433 specialization).
/// The `half()` calls inside absorb the 2^logn factor, so the output is the
/// true inverse without a separate 1/n scaling step.
pub(crate) fn mq_intt(logn: usize, a: &mut [u16]) {
    let mut t = 1usize;
    for lm in 0..logn {
        let hm = 1usize << (logn - 1 - lm);
        let dt = t << 1;
        let mut v0 = 0usize;
        for u in 0..hm {
            let s = MQ_IGM[u + hm] as u32;
            for v in 0..t {
                let k1 = v0 + v;
                let k2 = k1 + t;
                let x1 = a[k1] as u32;
                let x2 = a[k2] as u32;
                a[k1] = mq_half(mq_add(x1, x2)) as u16;
                a[k2] = mq_montymul(s, mq_sub(x1, x2)) as u16;
            }
            v0 += dt;
        }
        t = dt;
    }
}

// === Public API ===

/// Forward NTT: transforms a polynomial from coefficient representation to
/// evaluation representation (NTT domain).
///
/// Input coefficients are in `[0, q)`. Output coefficients are in `[0, q)`.
/// The internal computation uses the `[1, q]` Montgomery domain; conversion
/// happens at the boundary.
pub fn ntt(p: &Poly) -> Poly {
    let mut buf = [0u16; HAWK_N];
    for i in 0..HAWK_N {
        buf[i] = lift_to_mq(p.coeffs[i]);
    }
    mq_ntt(HAWK_LOGN, &mut buf);
    for c in buf.iter_mut() {
        *c = lift_from_mq(*c);
    }
    Poly::new(buf)
}

/// Inverse NTT: transforms a polynomial from evaluation representation back
/// to coefficient representation.
///
/// Input and output coefficients are in `[0, q)`.
pub fn intt(p: &Poly) -> Poly {
    let mut buf = [0u16; HAWK_N];
    for i in 0..HAWK_N {
        buf[i] = lift_to_mq(p.coeffs[i]);
    }
    mq_intt(HAWK_LOGN, &mut buf);
    for c in buf.iter_mut() {
        *c = lift_from_mq(*c);
    }
    Poly::new(buf)
}

/// Pointwise multiplication of two polynomials already in NTT (evaluation) form.
///
/// Uses `tomonty(montymul(a[u], b[u]))` so the result is `a[u] * b[u] mod q`
/// in `[1, q]`, ready for `iNTT`. The formula follows directly from Montgomery-domain
/// semantics: `montymul(a, b)` returns `a*b*R^{-1} mod q` (modq.h:113-116) and
/// `tomonty(·)` multiplies by R (modq.h:127-130), yielding `a*b mod q` in the [1..q]
/// representation that `iNTT` consumes.
///
/// Inputs and output are in `[0, q)` (standard Poly representation).
pub fn pointwise_mul(a: &Poly, b: &Poly) -> Poly {
    let mut out = [0u16; HAWK_N];
    for i in 0..HAWK_N {
        let ai = lift_to_mq(a.coeffs[i]) as u32;
        let bi = lift_to_mq(b.coeffs[i]) as u32;
        let prod = mq_tomonty(mq_montymul(ai, bi));
        out[i] = lift_from_mq(prod as u16);
    }
    Poly::new(out)
}

/// Polynomial multiplication via NTT: computes `a * b mod (X^n + 1, q)`.
///
/// Equivalent to `intt(pointwise_mul(ntt(a), ntt(b)))`.
/// Uses O(n log n) operations instead of the O(n²) schoolbook algorithm.
pub fn mul(a: &Poly, b: &Poly) -> Poly {
    intt(&pointwise_mul(&ntt(a), &ntt(b)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn ntt_intt_is_identity() {
        let mut coeffs = [0u16; HAWK_N];
        for i in 0..HAWK_N {
            coeffs[i] = (i as u16) % HAWK_Q as u16;
        }
        let p = Poly::new(coeffs);
        assert_eq!(intt(&ntt(&p)), p);
    }

    #[test]
    fn ntt_mul_matches_schoolbook() {
        let a_coeffs: [u16; HAWK_N] = std::array::from_fn(|i| ((i as u32 * 3) % HAWK_Q) as u16);
        let b_coeffs: [u16; HAWK_N] = std::array::from_fn(|i| ((i as u32 * 7 + 1) % HAWK_Q) as u16);
        let a = Poly::new(a_coeffs);
        let b = Poly::new(b_coeffs);
        assert_eq!(mul(&a, &b), a.mul_schoolbook(&b));
    }

    #[test]
    fn mq_set_small_roundtrip() {
        // For x in [-3, 3], mq_set_small produces a value in [1, q].
        // Then mq_snorm should recover x (after casting).
        for x in -3i32..=3 {
            let m = mq_set_small(x);
            assert!(m >= 1 && m <= HAWK_Q, "x={} → m={} out of [1, q]", x, m);
            let back = mq_snorm(m);
            assert_eq!(back, x, "roundtrip failed for x={}", x);
        }
        // Boundary: x = (q-1)/2 should work.
        let boundary = ((HAWK_Q - 1) / 2) as i32;
        let m = mq_set_small(boundary);
        let back = mq_snorm(m);
        assert_eq!(back, boundary);
        let m2 = mq_set_small(-boundary);
        let back2 = mq_snorm(m2);
        assert_eq!(back2, -boundary);
    }

    #[test]
    fn mq_poly_set_small_basic() {
        let a: [i8; 8] = [0, 1, -1, 2, -2, 3, -3, 127];
        let mut d = [0u16; 8];
        mq_poly_set_small(&mut d, &a);
        // Verify each coefficient round-trips.
        for i in 0..8 {
            assert_eq!(mq_snorm(d[i] as u32), a[i] as i32, "coef {}", i);
        }
    }

    #[test]
    fn mq_poly_snorm_basic() {
        let a: [i8; 8] = [0, 1, -1, 2, -2, 3, -3, 100];
        let mut d = [0u16; 8];
        mq_poly_set_small(&mut d, &a);
        mq_poly_snorm(&mut d);
        // After snorm, the u16 values are the signed coefficients in two's complement.
        for i in 0..8 {
            assert_eq!(d[i] as i16, a[i] as i16, "coef {}", i);
        }
    }

    proptest! {
        #[test]
        fn proptest_ntt_intt_roundtrip(
            v in prop::collection::vec(any::<u16>(), HAWK_N)
        ) {
            let mut coeffs = [0u16; HAWK_N];
            for i in 0..HAWK_N {
                coeffs[i] = v[i] % HAWK_Q as u16;
            }
            let p = Poly::new(coeffs);
            prop_assert_eq!(intt(&ntt(&p)), p);
        }

        #[test]
        fn proptest_ntt_mul_matches_schoolbook(
            va in prop::collection::vec(any::<u16>(), HAWK_N),
            vb in prop::collection::vec(any::<u16>(), HAWK_N),
        ) {
            let mut a_coeffs = [0u16; HAWK_N];
            let mut b_coeffs = [0u16; HAWK_N];
            for i in 0..HAWK_N {
                a_coeffs[i] = va[i] % HAWK_Q as u16;
                b_coeffs[i] = vb[i] % HAWK_Q as u16;
            }
            let a = Poly::new(a_coeffs);
            let b = Poly::new(b_coeffs);
            prop_assert_eq!(mul(&a, &b), a.mul_schoolbook(&b));
        }
    }
}
