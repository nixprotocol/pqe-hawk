//! Modular arithmetic helpers for 31-bit primes used by the HAWK keygen path.
//!
//! All operations work modulo a prime p with `(4/3)*2^30 < p < 2^31`. Operands
//! are `u32` integers in `[0, p)`. Values "in Montgomery representation" are
//! `x*R mod p` where `R = 2^32 mod p`.
//!
//! Port of the inline modular helpers in c-reference/hawk-512/ng_inner.h
//! (lines 348-488). Keeps all functions constant-time (no branches on
//! secret data).
//!
//! # Visibility note
//!
//! All functions are `pub` (not `pub(crate)`) to allow the FFI cross-check
//! integration tests in `tests/cross_check.rs` to call them directly via
//! `pqe_hawk::keygen::mp31::*`. No API stability guarantees are made beyond
//! this crate's own test suite.

/// Expand the top bit of `x` into a full 32-bit mask (all-ones if
/// `x >= 2^31`, all-zeros otherwise). Constant time.
///
/// Port of `tbmask` from c-reference/hawk-512/ng_inner.h:352-356.
#[inline]
pub fn tbmask(x: u32) -> u32 {
    ((x as i32) >> 31) as u32
}

/// Convert a signed integer in `[-(p-1), p-1]` to the `[0, p)` representation.
///
/// Port of `mp_set` from c-reference/hawk-512/ng_inner.h:362-367.
#[inline]
pub fn mp_set(v: i32, p: u32) -> u32 {
    let w = v as u32;
    w.wrapping_add(p & tbmask(w))
}

/// Convert `[0, p)` to the signed-normalized value in `[-(p-1)/2, +(p-1)/2]`.
///
/// Port of `mp_norm` from c-reference/hawk-512/ng_inner.h:372-377.
#[inline]
pub fn mp_norm(x: u32, p: u32) -> i32 {
    let w = x.wrapping_sub(p & tbmask((p >> 1).wrapping_sub(x)));
    w as i32
}

/// `R = 2^32 mod p`.
///
/// Port of `mp_R` from c-reference/hawk-512/ng_inner.h:398-405.
#[inline]
pub fn mp_r(p: u32) -> u32 {
    (p << 1).wrapping_neg()
}

/// `R/2 = 2^31 mod p`.
///
/// Port of `mp_hR` from c-reference/hawk-512/ng_inner.h:410-417.
#[inline]
pub fn mp_h_r(p: u32) -> u32 {
    (1u32 << 31).wrapping_sub(p)
}

/// Addition mod p.
///
/// Port of `mp_add` from c-reference/hawk-512/ng_inner.h:422-427.
#[inline]
pub fn mp_add(a: u32, b: u32, p: u32) -> u32 {
    let d = a.wrapping_add(b).wrapping_sub(p);
    d.wrapping_add(p & tbmask(d))
}

/// Subtraction mod p.
///
/// Port of `mp_sub` from c-reference/hawk-512/ng_inner.h:432-437.
#[inline]
pub fn mp_sub(a: u32, b: u32, p: u32) -> u32 {
    let d = a.wrapping_sub(b);
    d.wrapping_add(p & tbmask(d))
}

/// Halving mod p (multiplication by the modular inverse of 2).
///
/// Port of `mp_half` from c-reference/hawk-512/ng_inner.h:442-446.
#[inline]
pub fn mp_half(a: u32, p: u32) -> u32 {
    (a.wrapping_add(p & (a & 1).wrapping_neg())) >> 1
}

/// Montgomery multiplication: returns `a*b/R mod p` where `R = 2^32`.
///
/// Port of `mp_montymul` from c-reference/hawk-512/ng_inner.h:460-467.
#[inline]
pub fn mp_montymul(a: u32, b: u32, p: u32, p0i: u32) -> u32 {
    let z = (a as u64).wrapping_mul(b as u64);
    let w = ((z as u32).wrapping_mul(p0i)) as u32;
    let d = (z.wrapping_add((w as u64).wrapping_mul(p as u64)) >> 32) as u32;
    d.wrapping_sub(p)
        .wrapping_add(p & tbmask(d.wrapping_sub(p)))
}

/// Computes `2^(31*e) mod p` using square-and-multiply.
///
/// Port of `mp_Rx31` from c-reference/hawk-512/ng_inner.h:472-488.
pub fn mp_rx31(mut e: u32, p: u32, p0i: u32, r2: u32) -> u32 {
    let mut x = mp_half(r2, p);
    let mut d = 1u32;

    loop {
        if (e & 1) != 0 {
            d = mp_montymul(d, x, p, p0i);
        }
        e >>= 1;
        if e == 0 {
            return d;
        }
        x = mp_montymul(x, x, p, p0i);
    }
}

// === Sub-task 15.4: NTT over a 31-bit prime ===

/// Bit-reversal table over 10 bits (indices 0..1024).
///
/// Port of `REV10` from c-reference/hawk-512/ng_mp31.c:34-121.
#[rustfmt::skip]
pub(crate) const REV10: [u16; 1024] = [
       0,  512,  256,  768,  128,  640,  384,  896,   64,  576,  320,  832,
     192,  704,  448,  960,   32,  544,  288,  800,  160,  672,  416,  928,
      96,  608,  352,  864,  224,  736,  480,  992,   16,  528,  272,  784,
     144,  656,  400,  912,   80,  592,  336,  848,  208,  720,  464,  976,
      48,  560,  304,  816,  176,  688,  432,  944,  112,  624,  368,  880,
     240,  752,  496, 1008,    8,  520,  264,  776,  136,  648,  392,  904,
      72,  584,  328,  840,  200,  712,  456,  968,   40,  552,  296,  808,
     168,  680,  424,  936,  104,  616,  360,  872,  232,  744,  488, 1000,
      24,  536,  280,  792,  152,  664,  408,  920,   88,  600,  344,  856,
     216,  728,  472,  984,   56,  568,  312,  824,  184,  696,  440,  952,
     120,  632,  376,  888,  248,  760,  504, 1016,    4,  516,  260,  772,
     132,  644,  388,  900,   68,  580,  324,  836,  196,  708,  452,  964,
      36,  548,  292,  804,  164,  676,  420,  932,  100,  612,  356,  868,
     228,  740,  484,  996,   20,  532,  276,  788,  148,  660,  404,  916,
      84,  596,  340,  852,  212,  724,  468,  980,   52,  564,  308,  820,
     180,  692,  436,  948,  116,  628,  372,  884,  244,  756,  500, 1012,
      12,  524,  268,  780,  140,  652,  396,  908,   76,  588,  332,  844,
     204,  716,  460,  972,   44,  556,  300,  812,  172,  684,  428,  940,
     108,  620,  364,  876,  236,  748,  492, 1004,   28,  540,  284,  796,
     156,  668,  412,  924,   92,  604,  348,  860,  220,  732,  476,  988,
      60,  572,  316,  828,  188,  700,  444,  956,  124,  636,  380,  892,
     252,  764,  508, 1020,    2,  514,  258,  770,  130,  642,  386,  898,
      66,  578,  322,  834,  194,  706,  450,  962,   34,  546,  290,  802,
     162,  674,  418,  930,   98,  610,  354,  866,  226,  738,  482,  994,
      18,  530,  274,  786,  146,  658,  402,  914,   82,  594,  338,  850,
     210,  722,  466,  978,   50,  562,  306,  818,  178,  690,  434,  946,
     114,  626,  370,  882,  242,  754,  498, 1010,   10,  522,  266,  778,
     138,  650,  394,  906,   74,  586,  330,  842,  202,  714,  458,  970,
      42,  554,  298,  810,  170,  682,  426,  938,  106,  618,  362,  874,
     234,  746,  490, 1002,   26,  538,  282,  794,  154,  666,  410,  922,
      90,  602,  346,  858,  218,  730,  474,  986,   58,  570,  314,  826,
     186,  698,  442,  954,  122,  634,  378,  890,  250,  762,  506, 1018,
       6,  518,  262,  774,  134,  646,  390,  902,   70,  582,  326,  838,
     198,  710,  454,  966,   38,  550,  294,  806,  166,  678,  422,  934,
     102,  614,  358,  870,  230,  742,  486,  998,   22,  534,  278,  790,
     150,  662,  406,  918,   86,  598,  342,  854,  214,  726,  470,  982,
      54,  566,  310,  822,  182,  694,  438,  950,  118,  630,  374,  886,
     246,  758,  502, 1014,   14,  526,  270,  782,  142,  654,  398,  910,
      78,  590,  334,  846,  206,  718,  462,  974,   46,  558,  302,  814,
     174,  686,  430,  942,  110,  622,  366,  878,  238,  750,  494, 1006,
      30,  542,  286,  798,  158,  670,  414,  926,   94,  606,  350,  862,
     222,  734,  478,  990,   62,  574,  318,  830,  190,  702,  446,  958,
     126,  638,  382,  894,  254,  766,  510, 1022,    1,  513,  257,  769,
     129,  641,  385,  897,   65,  577,  321,  833,  193,  705,  449,  961,
      33,  545,  289,  801,  161,  673,  417,  929,   97,  609,  353,  865,
     225,  737,  481,  993,   17,  529,  273,  785,  145,  657,  401,  913,
      81,  593,  337,  849,  209,  721,  465,  977,   49,  561,  305,  817,
     177,  689,  433,  945,  113,  625,  369,  881,  241,  753,  497, 1009,
       9,  521,  265,  777,  137,  649,  393,  905,   73,  585,  329,  841,
     201,  713,  457,  969,   41,  553,  297,  809,  169,  681,  425,  937,
     105,  617,  361,  873,  233,  745,  489, 1001,   25,  537,  281,  793,
     153,  665,  409,  921,   89,  601,  345,  857,  217,  729,  473,  985,
      57,  569,  313,  825,  185,  697,  441,  953,  121,  633,  377,  889,
     249,  761,  505, 1017,    5,  517,  261,  773,  133,  645,  389,  901,
      69,  581,  325,  837,  197,  709,  453,  965,   37,  549,  293,  805,
     165,  677,  421,  933,  101,  613,  357,  869,  229,  741,  485,  997,
      21,  533,  277,  789,  149,  661,  405,  917,   85,  597,  341,  853,
     213,  725,  469,  981,   53,  565,  309,  821,  181,  693,  437,  949,
     117,  629,  373,  885,  245,  757,  501, 1013,   13,  525,  269,  781,
     141,  653,  397,  909,   77,  589,  333,  845,  205,  717,  461,  973,
      45,  557,  301,  813,  173,  685,  429,  941,  109,  621,  365,  877,
     237,  749,  493, 1005,   29,  541,  285,  797,  157,  669,  413,  925,
      93,  605,  349,  861,  221,  733,  477,  989,   61,  573,  317,  829,
     189,  701,  445,  957,  125,  637,  381,  893,  253,  765,  509, 1021,
       3,  515,  259,  771,  131,  643,  387,  899,   67,  579,  323,  835,
     195,  707,  451,  963,   35,  547,  291,  803,  163,  675,  419,  931,
      99,  611,  355,  867,  227,  739,  483,  995,   19,  531,  275,  787,
     147,  659,  403,  915,   83,  595,  339,  851,  211,  723,  467,  979,
      51,  563,  307,  819,  179,  691,  435,  947,  115,  627,  371,  883,
     243,  755,  499, 1011,   11,  523,  267,  779,  139,  651,  395,  907,
      75,  587,  331,  843,  203,  715,  459,  971,   43,  555,  299,  811,
     171,  683,  427,  939,  107,  619,  363,  875,  235,  747,  491, 1003,
      27,  539,  283,  795,  155,  667,  411,  923,   91,  603,  347,  859,
     219,  731,  475,  987,   59,  571,  315,  827,  187,  699,  443,  955,
     123,  635,  379,  891,  251,  763,  507, 1019,    7,  519,  263,  775,
     135,  647,  391,  903,   71,  583,  327,  839,  199,  711,  455,  967,
      39,  551,  295,  807,  167,  679,  423,  935,  103,  615,  359,  871,
     231,  743,  487,  999,   23,  535,  279,  791,  151,  663,  407,  919,
      87,  599,  343,  855,  215,  727,  471,  983,   55,  567,  311,  823,
     183,  695,  439,  951,  119,  631,  375,  887,  247,  759,  503, 1015,
      15,  527,  271,  783,  143,  655,  399,  911,   79,  591,  335,  847,
     207,  719,  463,  975,   47,  559,  303,  815,  175,  687,  431,  943,
     111,  623,  367,  879,  239,  751,  495, 1007,   31,  543,  287,  799,
     159,  671,  415,  927,   95,  607,  351,  863,  223,  735,  479,  991,
      63,  575,  319,  831,  191,  703,  447,  959,  127,  639,  383,  895,
     255,  767,  511, 1023,
];

/// Fill `gm` with powers of `g` (Montgomery form) in bit-reversed order.
///
/// Port of `mp_mkgm` from c-reference/hawk-512/ng_mp31.c:156-179.
/// `gm` must have length at least `1 << logn`.
/// `g` is a primitive 2^10-th root of unity mod p, already precomputed
/// in the `SmallPrime` table.
pub fn mp_mkgm(logn: u32, gm: &mut [u32], g: u32, p: u32, p0i: u32) {
    let n = 1usize << logn;
    debug_assert!(gm.len() >= n, "gm too small: {} < {}", gm.len(), n);

    // Square g enough times so it becomes a primitive 2^logn-th root.
    let mut g = g;
    for _ in logn..10 {
        g = mp_montymul(g, g, p, p0i);
    }
    let k = (10 - logn) as usize;
    let mut x1 = mp_r(p);
    let mut u = 0usize;
    loop {
        let v = REV10[u << k] as usize;
        gm[v] = x1;
        u += 1;
        if u >= n {
            break;
        }
        x1 = mp_montymul(x1, g, p, p0i);
    }
}

/// Fill `igm` with powers of `ig` (Montgomery form) scaled by 1/2, in
/// bit-reversed order.
///
/// Port of `mp_mkigm` from c-reference/hawk-512/ng_mp31.c:195-218.
/// The extra `1/2` factor (via `mp_h_r`) is absorbed into the iNTT
/// normalization so that `mp_intt(mp_ntt(a))` recovers `a`.
pub fn mp_mkigm(logn: u32, igm: &mut [u32], ig: u32, p: u32, p0i: u32) {
    let n = 1usize << logn;
    debug_assert!(igm.len() >= n, "igm too small: {} < {}", igm.len(), n);

    let mut ig = ig;
    for _ in logn..10 {
        ig = mp_montymul(ig, ig, p, p0i);
    }
    let k = (10 - logn) as usize;
    let mut x2 = mp_h_r(p);
    let mut u = 0usize;
    loop {
        let v = REV10[u << k] as usize;
        igm[v] = x2;
        u += 1;
        if u >= n {
            break;
        }
        x2 = mp_montymul(x2, ig, p, p0i);
    }
}

/// Fill both `gm` and `igm` simultaneously (combined mkgm + mkigm).
///
/// Port of `mp_mkgmigm` from c-reference/hawk-512/ng_mp31.c:125-154.
/// `gm` and `igm` must each have length at least `1 << logn`.
/// `g` is a primitive 2^10-th root of unity mod p; `ig` is its inverse.
pub fn mp_mkgmigm(logn: u32, gm: &mut [u32], igm: &mut [u32], g: u32, ig: u32, p: u32, p0i: u32) {
    let n = 1usize << logn;
    debug_assert!(gm.len() >= n, "gm too small: {} < {}", gm.len(), n);
    debug_assert!(igm.len() >= n, "igm too small: {} < {}", igm.len(), n);

    // Square g and ig enough times so they become primitive 2^logn-th roots.
    let mut g = g;
    let mut ig = ig;
    for _ in logn..10 {
        g = mp_montymul(g, g, p, p0i);
        ig = mp_montymul(ig, ig, p, p0i);
    }
    let k = (10 - logn) as usize;
    let mut x1 = mp_r(p);
    let mut x2 = mp_h_r(p);
    let mut u = 0usize;
    loop {
        let v = REV10[u << k] as usize;
        gm[v] = x1;
        igm[v] = x2;
        u += 1;
        if u >= n {
            break;
        }
        x1 = mp_montymul(x1, g, p, p0i);
        x2 = mp_montymul(x2, ig, p, p0i);
    }
}

/// Forward NTT over Z_p\[X\]/(X^n+1). Operates in place on `a`.
///
/// Port of `mp_NTT` from c-reference/hawk-512/ng_mp31.c:220-247.
/// `gm` must have been populated by `mp_mkgm`. `a.len()` must be at
/// least `1 << logn`.
pub fn mp_ntt(logn: u32, a: &mut [u32], gm: &[u32], p: u32, p0i: u32) {
    if logn == 0 {
        return;
    }
    let mut t = 1usize << logn;
    for lm in 0..logn {
        let m = 1usize << lm;
        let ht = t >> 1;
        let mut v0 = 0usize;
        for u in 0..m {
            let s = gm[u + m];
            for v in 0..ht {
                let k1 = v0 + v;
                let k2 = k1 + ht;
                let x1 = a[k1];
                let x2 = mp_montymul(a[k2], s, p, p0i);
                a[k1] = mp_add(x1, x2, p);
                a[k2] = mp_sub(x1, x2, p);
            }
            v0 += t;
        }
        t = ht;
    }
}

/// Inverse NTT over Z_p\[X\]/(X^n+1). Operates in place on `a`.
///
/// Port of `mp_iNTT` from c-reference/hawk-512/ng_mp31.c:249-277.
/// `igm` must have been populated by `mp_mkigm`. After this call,
/// `a` holds the original polynomial (the 1/n factor is folded into
/// the `mp_h_r` seeds of `mp_mkigm`).
pub fn mp_intt(logn: u32, a: &mut [u32], igm: &[u32], p: u32, p0i: u32) {
    if logn == 0 {
        return;
    }
    let mut t = 1usize;
    for lm in 0..logn {
        let hm = 1usize << (logn - 1 - lm);
        let dt = t << 1;
        let mut v0 = 0usize;
        for u in 0..hm {
            let s = igm[u + hm];
            for v in 0..t {
                let k1 = v0 + v;
                let k2 = k1 + t;
                let x1 = a[k1];
                let x2 = a[k2];
                a[k1] = mp_half(mp_add(x1, x2, p), p);
                a[k2] = mp_montymul(mp_sub(x1, x2, p), s, p, p0i);
            }
            v0 += dt;
        }
        t = dt;
    }
}

// === Sub-task 15.20: mp_div ===

/// Constant-time modular division: returns `x / y mod p`. If `y` is not
/// invertible mod `p` (i.e. gcd(y, p) > 1), returns 0.
///
/// Implements a constant-time binary GCD inversion: runs exactly 62 iterations
/// of the half-GCD reduction, accumulating the Bézout coefficient for `y` in
/// `v`. The final mask `tbmask(b - 2)` is all-ones iff b == 1 (gcd = 1).
///
/// Port of `mp_div` (c-reference/hawk-512/ng_mp31.c:4-28).
pub fn mp_div(x: u32, y: u32, p: u32) -> u32 {
    let mut a = y;
    let mut b = p;
    let mut u = x;
    let mut v: u32 = 0;
    for _ in 0..62 {
        let a_odd = (a & 1).wrapping_neg();
        let swap = tbmask(a.wrapping_sub(b)) & a_odd;
        let t1 = swap & (a ^ b);
        a ^= t1;
        b ^= t1;
        let t2 = swap & (u ^ v);
        u ^= t2;
        v ^= t2;
        a = a.wrapping_sub(a_odd & b);
        u = mp_sub(u, a_odd & v, p);
        a >>= 1;
        u = mp_half(u, p);
    }
    v & tbmask(b.wrapping_sub(2))
}

#[cfg(test)]
mod tests {
    use super::*;

    // Use PRIMES[0] as the working prime in all tests.
    // p = 2147473409, p0i = 2042615807, r2 = 419348484
    const P: u32 = 2147473409;
    const P0I: u32 = 2042615807;
    const R2: u32 = 419348484;

    #[test]
    fn tbmask_top_bit_set() {
        assert_eq!(tbmask(0x80000000), 0xFFFFFFFF);
        assert_eq!(tbmask(0xFFFFFFFF), 0xFFFFFFFF);
    }

    #[test]
    fn tbmask_top_bit_clear() {
        assert_eq!(tbmask(0), 0);
        assert_eq!(tbmask(0x7FFFFFFF), 0);
    }

    #[test]
    fn mp_set_positive_and_negative() {
        assert_eq!(mp_set(5, P), 5);
        // mp_set(-1, P) should return P - 1.
        assert_eq!(mp_set(-1, P), P - 1);
        assert_eq!(mp_set(0, P), 0);
    }

    #[test]
    fn mp_norm_roundtrip() {
        // For any v in [-(p-1)/2, (p-1)/2], mp_norm(mp_set(v, p), p) == v
        for v in [-100i32, -1, 0, 1, 100, 1_000_000] {
            assert_eq!(mp_norm(mp_set(v, P), P), v);
        }
        // Values above p/2 should come back negative.
        assert_eq!(mp_norm(P - 1, P), -1);
    }

    #[test]
    fn mp_add_basic() {
        assert_eq!(mp_add(0, 0, P), 0);
        assert_eq!(mp_add(5, 7, P), 12);
        // wrap: (P - 1) + 1 == 0
        assert_eq!(mp_add(P - 1, 1, P), 0);
        // wrap: (P - 1) + (P - 1) == P - 2
        assert_eq!(mp_add(P - 1, P - 1, P), P - 2);
    }

    #[test]
    fn mp_sub_basic() {
        assert_eq!(mp_sub(5, 3, P), 2);
        assert_eq!(mp_sub(0, 0, P), 0);
        // wrap: 0 - 1 == P - 1
        assert_eq!(mp_sub(0, 1, P), P - 1);
    }

    #[test]
    fn mp_half_basic() {
        // For even values, half is straightforward.
        assert_eq!(mp_half(4, P), 2);
        assert_eq!(mp_half(0, P), 0);
        // 1/2 mod P: since P is odd, (P+1)/2 is the inverse of 2.
        // mp_half(1, P) should equal (P + 1) / 2.
        assert_eq!(mp_half(1, P), (P + 1) / 2);
        // Sanity: mp_half(a, P) * 2 mod P == a
        for a in [0u32, 1, 2, 3, 4, 100, P - 1] {
            let h = mp_half(a, P);
            assert_eq!(mp_add(h, h, P), a, "half of {} was {}", a, h);
        }
    }

    #[test]
    fn mp_r_and_h_r() {
        // mp_r(p) == 2^32 mod p (= 2^32 - 2*p for primes in our range)
        let r = mp_r(P);
        // 2*r mod p should equal 0 (since 2*R = 2*(2^32 mod p), and 2^33 mod p needs a further reduction)
        // Easier check: mp_r(p) < p and (p + mp_r(p)) mod 2^32 == 2^32 mod p
        assert!(r < P);
        // mp_h_r: 2^31 mod p == p - 2^31 if 2^31 > p, or 2^31 - p otherwise.
        let h = mp_h_r(P);
        assert!(h < P);
        // 2 * mp_h_r == mp_r (mod p)
        assert_eq!(mp_add(h, h, P), r);
    }

    #[test]
    fn mp_montymul_one_and_r2() {
        // montymul(x, R2) converts x to Montgomery form: x*R mod p
        // In particular montymul(1, R2) == R mod p == mp_r(p)
        assert_eq!(mp_montymul(1, R2, P, P0I), mp_r(P));
    }

    #[test]
    fn mp_montymul_multiplication() {
        // For x = 5, y = 7 in normal form: montymul(x*R, y) == x*y mod p.
        // Easier: montymul(montymul(a, R2), montymul(b, R2)) = montymul(a*b, R2)
        // We can also check: montymul(a, R2) * b normal. Use direct test:
        //   montymul(R2, x) gives x*R mod p. Then montymul that with y normal gives x*y mod p.
        let a_mont = mp_montymul(5, R2, P, P0I);
        let prod = mp_montymul(a_mont, 7, P, P0I);
        assert_eq!(prod, 35);
    }

    #[test]
    fn mp_rx31_zero_e() {
        // mp_rx31(0, p, p0i, r2) should return the multiplicative identity in
        // the Montgomery domain — wait, actually it returns 2^0 mod p = 1.
        // Reading the loop: starts with d = 1, e = 0 never enters the body,
        // returns d = 1.
        assert_eq!(mp_rx31(0, P, P0I, R2), 1);
    }

    #[test]
    fn mp_rx31_one_e() {
        // e = 1: loop runs once, d = montymul(1, mp_half(R2, p), p, p0i).
        // mp_half(R2, p) is 2^63 mod p (Montgomery rep of 2^31).
        // montymul(1, 2^63-mod-p) = 1 * (2^63 mod p) / R = (2^63 / 2^32) mod p = 2^31 mod p.
        // That equals mp_h_r(p).
        assert_eq!(mp_rx31(1, P, P0I, R2), mp_h_r(P));
    }

    // === 15.19a: mp_mkgmigm test ===

    #[test]
    fn mp_mkgmigm_produces_same_as_mkgm_plus_mkigm() {
        use crate::keygen::primes::PRIMES;
        let p = PRIMES[0].p;
        let p0i = PRIMES[0].p0i;
        let g = PRIMES[0].g;
        let ig = PRIMES[0].ig;
        let mut gm_a = [0u32; 512];
        let mut gm_b = [0u32; 512];
        let mut igm_a = [0u32; 512];
        let mut igm_b = [0u32; 512];
        mp_mkgm(9, &mut gm_a, g, p, p0i);
        mp_mkigm(9, &mut igm_a, ig, p, p0i);
        mp_mkgmigm(9, &mut gm_b, &mut igm_b, g, ig, p, p0i);
        assert_eq!(gm_a, gm_b);
        assert_eq!(igm_a, igm_b);
    }

    // === 15.4: NTT correctness tests ===

    #[test]
    fn rev10_size_is_1024() {
        assert_eq!(REV10.len(), 1024);
        // Basic sanity: REV10[1] should be bit-reverse of 0b0000000001 over
        // 10 bits = 0b1000000000 = 512.
        assert_eq!(REV10[1], 512);
        assert_eq!(REV10[0], 0);
        assert_eq!(REV10[1023], 1023);
    }

    #[test]
    fn mp_mkgm_first_entry() {
        use crate::keygen::primes::PRIMES;
        // gm[REV10[0]] = gm[0] should equal mp_r(p) (Montgomery repr of 1).
        let mut gm = [0u32; 512];
        mp_mkgm(9, &mut gm, PRIMES[0].g, PRIMES[0].p, PRIMES[0].p0i);
        assert_eq!(gm[0], mp_r(PRIMES[0].p));
    }

    #[test]
    fn mp_mkigm_first_entry() {
        use crate::keygen::primes::PRIMES;
        // igm[0] should equal mp_h_r(p) (= R/2 mod p; the extra 1/2 factor).
        let mut igm = [0u32; 512];
        mp_mkigm(9, &mut igm, PRIMES[0].ig, PRIMES[0].p, PRIMES[0].p0i);
        assert_eq!(igm[0], mp_h_r(PRIMES[0].p));
    }

    #[test]
    fn mp_ntt_zero_polynomial() {
        use crate::keygen::primes::PRIMES;
        let mut gm = [0u32; 512];
        let mut igm = [0u32; 512];
        mp_mkgm(9, &mut gm, PRIMES[0].g, PRIMES[0].p, PRIMES[0].p0i);
        mp_mkigm(9, &mut igm, PRIMES[0].ig, PRIMES[0].p, PRIMES[0].p0i);
        let mut a = [0u32; 512];
        mp_ntt(9, &mut a, &gm, PRIMES[0].p, PRIMES[0].p0i);
        // NTT of the zero polynomial is zero in evaluation form.
        for i in 0..512 {
            assert_eq!(a[i], 0);
        }
    }

    #[test]
    fn mp_ntt_roundtrip_logn_9() {
        use crate::keygen::primes::PRIMES;
        // Build gm/igm for logn=9 = 512-point NTT over PRIMES[0].
        let mut gm = [0u32; 512];
        let mut igm = [0u32; 512];
        mp_mkgm(9, &mut gm, PRIMES[0].g, PRIMES[0].p, PRIMES[0].p0i);
        mp_mkigm(9, &mut igm, PRIMES[0].ig, PRIMES[0].p, PRIMES[0].p0i);

        // Input polynomial with a simple pattern.
        let mut a = [0u32; 512];
        for i in 0..512 {
            a[i] = (i as u32) % PRIMES[0].p;
        }
        let original = a;

        mp_ntt(9, &mut a, &gm, PRIMES[0].p, PRIMES[0].p0i);
        mp_intt(9, &mut a, &igm, PRIMES[0].p, PRIMES[0].p0i);

        // Verify recovery.
        for i in 0..512 {
            assert_eq!(a[i], original[i], "coeff {} mismatch", i);
        }
    }

    // === Sub-task 15.20: mp_div tests ===

    #[test]
    fn mp_div_basic() {
        use crate::keygen::primes::PRIMES;
        let p = PRIMES[0].p;
        // 7 / 3 mod p: verify that 3 * (7/3) ≡ 7 (mod p).
        let q = mp_div(7, 3, p);
        let prod = ((3u64 * q as u64) % p as u64) as u32;
        assert_eq!(prod, 7);
    }

    #[test]
    fn mp_div_non_invertible_returns_zero() {
        use crate::keygen::primes::PRIMES;
        let p = PRIMES[0].p;
        // y = 0 is certainly non-invertible mod p.
        let q = mp_div(5, 0, p);
        assert_eq!(q, 0);
    }
}
