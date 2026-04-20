//! Feature-gated tests that diff our Rust implementation against the
//! reference C. Requires `--features cross-check-reference-c`.
//!
//! Pattern: each Rust ring/NTT/sign primitive has a corresponding `test_*`
//! C wrapper exported by the patched reference (see patches/). We byte-diff
//! the outputs over many proptest cases.

#![cfg(feature = "cross-check-reference-c")]
#![allow(
    non_upper_case_globals,
    non_camel_case_types,
    non_snake_case,
    dead_code
)]

include!(concat!(env!("OUT_DIR"), "/ffi_bindings.rs"));

use pqe_hawk::keygen::mp31;
use pqe_hawk::keygen::primes::PRIMES;
use pqe_hawk::params::{HAWK_N, HAWK_Q};
use pqe_hawk::ring::Poly;
use pqe_hawk::sample::NUM_SAMPLES;
use proptest::array::uniform32;
use proptest::prelude::*;

/// Deterministic `RngCore` that dispenses bytes from a pre-allocated buffer.
/// Both the Rust sampler and the C wrapper receive an identical byte stream,
/// enabling a byte-exact diff.
struct ByteBufRng<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl rand::RngCore for ByteBufRng<'_> {
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
        let end = self.pos + dst.len();
        dst.copy_from_slice(&self.buf[self.pos..end]);
        self.pos = end;
    }

    fn try_fill_bytes(&mut self, dst: &mut [u8]) -> Result<(), rand::Error> {
        self.fill_bytes(dst);
        Ok(())
    }
}

proptest! {
    #[test]
    fn ffi_poly_add_matches(
        va in prop::collection::vec(any::<u16>(), HAWK_N),
        vb in prop::collection::vec(any::<u16>(), HAWK_N),
    ) {
        let mut a = [0u16; HAWK_N];
        let mut b = [0u16; HAWK_N];
        for i in 0..HAWK_N {
            a[i] = va[i] % HAWK_Q as u16;
            b[i] = vb[i] % HAWK_Q as u16;
        }
        let a_poly = Poly::new(a);
        let b_poly = Poly::new(b);

        let rust_sum = a_poly.add(&b_poly);

        let mut c_sum = [0u16; HAWK_N];
        unsafe {
            test_poly_add(a.as_ptr(), b.as_ptr(), c_sum.as_mut_ptr());
        }

        prop_assert_eq!(rust_sum.coeffs, c_sum);
    }

    #[test]
    fn ffi_ntt_forward_matches(
        v in prop::collection::vec(any::<u16>(), HAWK_N),
    ) {
        let mut coeffs = [0u16; HAWK_N];
        for i in 0..HAWK_N {
            coeffs[i] = v[i] % HAWK_Q as u16;
        }
        let p = Poly::new(coeffs);
        let rust_ntt = pqe_hawk::ntt::ntt(&p);

        let mut c_ntt = [0u16; HAWK_N];
        unsafe {
            test_ntt_forward(coeffs.as_ptr(), c_ntt.as_mut_ptr());
        }
        prop_assert_eq!(rust_ntt.coeffs, c_ntt);
    }

    #[test]
    fn ffi_ntt_inverse_matches(
        v in prop::collection::vec(any::<u16>(), HAWK_N),
    ) {
        let mut coeffs = [0u16; HAWK_N];
        for i in 0..HAWK_N {
            coeffs[i] = v[i] % HAWK_Q as u16;
        }
        let p = Poly::new(coeffs);
        let rust_intt = pqe_hawk::ntt::intt(&p);

        let mut c_intt = [0u16; HAWK_N];
        unsafe {
            test_ntt_inverse(coeffs.as_ptr(), c_intt.as_mut_ptr());
        }
        prop_assert_eq!(rust_intt.coeffs, c_intt);
    }

    #[test]
    fn ffi_poly_mul_ntt_matches(
        va in prop::collection::vec(any::<u16>(), HAWK_N),
        vb in prop::collection::vec(any::<u16>(), HAWK_N),
    ) {
        let mut a = [0u16; HAWK_N];
        let mut b = [0u16; HAWK_N];
        for i in 0..HAWK_N {
            a[i] = va[i] % HAWK_Q as u16;
            b[i] = vb[i] % HAWK_Q as u16;
        }
        let a_poly = Poly::new(a);
        let b_poly = Poly::new(b);
        let rust_mul = pqe_hawk::ntt::mul(&a_poly, &b_poly);

        let mut c_mul = [0u16; HAWK_N];
        unsafe {
            test_poly_mul_ntt(a.as_ptr(), b.as_ptr(), c_mul.as_mut_ptr());
        }
        prop_assert_eq!(rust_mul.coeffs, c_mul);
    }

    #[test]
    fn ffi_compute_hm_matches(
        msg in prop::collection::vec(any::<u8>(), 0..=512usize),
        hpub in prop::collection::vec(any::<u8>(), 0..=64usize),
    ) {
        let rust_hm = pqe_hawk::hash::compute_hm(&msg, &hpub);
        let mut c_hm = [0u8; 64];
        unsafe {
            test_compute_hm(
                msg.as_ptr(), msg.len(),
                hpub.as_ptr(), hpub.len(),
                c_hm.as_mut_ptr(),
            );
        }
        prop_assert_eq!(rust_hm, c_hm);
    }

    #[test]
    fn ffi_compute_h_matches(
        hm in uniform32(any::<u8>()).prop_flat_map(|first32| {
            uniform32(any::<u8>()).prop_map(move |second32| {
                let mut hm = [0u8; 64];
                hm[..32].copy_from_slice(&first32);
                hm[32..].copy_from_slice(&second32);
                hm
            })
        }),
        salt in prop::collection::vec(any::<u8>(), 0..=64usize),
    ) {
        let (rust_h0, rust_h1) = pqe_hawk::hash::compute_h(&hm, &salt);
        let mut c_h0 = [0u8; 64];
        let mut c_h1 = [0u8; 64];
        unsafe {
            test_compute_h(
                hm.as_ptr(),
                salt.as_ptr(), salt.len(),
                c_h0.as_mut_ptr(), c_h1.as_mut_ptr(),
            );
        }
        prop_assert_eq!(rust_h0, c_h0);
        prop_assert_eq!(rust_h1, c_h1);
    }
}

proptest! {
    #[test]
    fn ffi_sample_matches(
        // 10240 random bytes for the RNG stream: (NUM_SAMPLES / 16) * 160 = 64 * 160
        rng_bytes in prop::collection::vec(any::<u8>(), 10240),
        // 128-byte parity vector: NUM_SAMPLES / 8
        t_parity in prop::collection::vec(any::<u8>(), NUM_SAMPLES / 8),
    ) {
        let mut t_arr = [0u8; 128];
        t_arr.copy_from_slice(&t_parity);

        // Rust side: consume bytes via ByteBufRng.
        let mut rng = ByteBufRng { buf: &rng_bytes, pos: 0 };
        let rust_out = pqe_hawk::sample::sample(&mut rng, &t_arr);

        // C side: same byte stream fed via test_sample.
        let mut c_x = [0i8; 1024];
        let c_norm = unsafe {
            test_sample(rng_bytes.as_ptr(), t_arr.as_ptr(), c_x.as_mut_ptr())
        };

        prop_assert_eq!(rust_out.x, c_x);
        prop_assert_eq!(rust_out.squared_norm, c_norm);
    }
}

// === Sub-task 15.3: mp_* helpers byte-exact cross-check against ng_inner.h ===

// Use PRIMES[0] as the primary working prime for proptests:
// p = 2147473409, p0i = 2042615807, r2 = 419348484.
const MP_P: u32 = 2147473409;
const MP_P0I: u32 = 2042615807;
const MP_R2: u32 = 419348484;

// Strategy: produce u32 values in [0, MP_P).
fn mod_p_strategy() -> impl Strategy<Value = u32> {
    any::<u32>().prop_map(|x| x % MP_P)
}

proptest! {
    #[test]
    fn ffi_mp_set_matches(v in -((MP_P - 1) as i32)..=((MP_P - 1) as i32)) {
        let rust = mp31::mp_set(v, MP_P);
        let c = unsafe { test_mp_set(v, MP_P) };
        prop_assert_eq!(rust, c);
    }

    #[test]
    fn ffi_mp_norm_matches(x in mod_p_strategy()) {
        let rust = mp31::mp_norm(x, MP_P);
        let c = unsafe { test_mp_norm(x, MP_P) };
        prop_assert_eq!(rust, c);
    }

    #[test]
    fn ffi_mp_add_matches(a in mod_p_strategy(), b in mod_p_strategy()) {
        let rust = mp31::mp_add(a, b, MP_P);
        let c = unsafe { test_mp_add(a, b, MP_P) };
        prop_assert_eq!(rust, c);
    }

    #[test]
    fn ffi_mp_sub_matches(a in mod_p_strategy(), b in mod_p_strategy()) {
        let rust = mp31::mp_sub(a, b, MP_P);
        let c = unsafe { test_mp_sub(a, b, MP_P) };
        prop_assert_eq!(rust, c);
    }

    #[test]
    fn ffi_mp_half_matches(a in mod_p_strategy()) {
        let rust = mp31::mp_half(a, MP_P);
        let c = unsafe { test_mp_half(a, MP_P) };
        prop_assert_eq!(rust, c);
    }

    #[test]
    fn ffi_mp_montymul_matches(a in mod_p_strategy(), b in mod_p_strategy()) {
        let rust = mp31::mp_montymul(a, b, MP_P, MP_P0I);
        let c = unsafe { test_mp_montymul(a, b, MP_P, MP_P0I) };
        prop_assert_eq!(rust, c);
    }

    #[test]
    fn ffi_mp_rx31_matches(e in 0u32..=100u32) {
        let rust = mp31::mp_rx31(e, MP_P, MP_P0I, MP_R2);
        let c = unsafe { test_mp_rx31(e, MP_P, MP_P0I, MP_R2) };
        prop_assert_eq!(rust, c);
    }
}

// Constant-function tests (no proptest needed):
#[test]
fn ffi_mp_r_matches() {
    let rust = mp31::mp_r(MP_P);
    let c = unsafe { test_mp_r(MP_P) };
    assert_eq!(rust, c);
}

#[test]
fn ffi_mp_h_r_matches() {
    let rust = mp31::mp_h_r(MP_P);
    let c = unsafe { test_mp_h_r(MP_P) };
    assert_eq!(rust, c);
}

// === Sub-task 15.5: mp_NTT / mp_iNTT cross-check ===

proptest! {
    #[test]
    fn ffi_mp_ntt_matches(
        v in prop::collection::vec(any::<u32>(), 512),
    ) {
        // Reduce each coefficient mod P so inputs are in [0, P).
        let mut a = [0u32; 512];
        for i in 0..512 {
            a[i] = v[i] % MP_P;
        }

        // Rust side: build gm and run mp_ntt.
        let mut gm = [0u32; 512];
        pqe_hawk::keygen::mp31::mp_mkgm(
            9, &mut gm, PRIMES[0].g, MP_P, MP_P0I,
        );
        let mut rust_out = a;
        pqe_hawk::keygen::mp31::mp_ntt(
            9, &mut rust_out, &gm, MP_P, MP_P0I,
        );

        // C side.
        let mut c_out = [0u32; 512];
        unsafe {
            test_mp_ntt(
                a.as_ptr(), c_out.as_mut_ptr(),
                PRIMES[0].g, MP_P, MP_P0I,
            );
        }
        prop_assert_eq!(rust_out, c_out);
    }

    #[test]
    fn ffi_mp_intt_matches(
        v in prop::collection::vec(any::<u32>(), 512),
    ) {
        let mut a = [0u32; 512];
        for i in 0..512 {
            a[i] = v[i] % MP_P;
        }

        let mut igm = [0u32; 512];
        pqe_hawk::keygen::mp31::mp_mkigm(
            9, &mut igm, PRIMES[0].ig, MP_P, MP_P0I,
        );
        let mut rust_out = a;
        pqe_hawk::keygen::mp31::mp_intt(
            9, &mut rust_out, &igm, MP_P, MP_P0I,
        );

        let mut c_out = [0u32; 512];
        unsafe {
            test_mp_intt(
                a.as_ptr(), c_out.as_mut_ptr(),
                PRIMES[0].ig, MP_P, MP_P0I,
            );
        }
        prop_assert_eq!(rust_out, c_out);
    }
}

// === Sub-task 15.7: zint primitives cross-check ===

// Strategy for a 31-bit limb (u32 in [0, 2^31)).
fn limb31_strategy() -> impl Strategy<Value = u32> {
    any::<u32>().prop_map(|x| x & 0x7FFF_FFFF)
}

// Strategy for a big-int of N limbs.
fn bigint_strategy(n: usize) -> impl Strategy<Value = Vec<u32>> {
    prop::collection::vec(limb31_strategy(), n..=n)
}

proptest! {
    #[test]
    fn ffi_zint_mul_small_matches(
        m in bigint_strategy(8),
        x in limb31_strategy(),
    ) {
        let mut rust_m = m.clone();
        let rust_cc = pqe_hawk::keygen::zint31::zint_mul_small(&mut rust_m, x);

        let mut c_m = m.clone();
        let c_cc = unsafe {
            test_zint_mul_small(c_m.as_mut_ptr(), c_m.len(), x)
        };

        prop_assert_eq!(rust_cc, c_cc);
        prop_assert_eq!(rust_m, c_m);
    }

    #[test]
    fn ffi_zint_mod_small_unsigned_matches(
        d in bigint_strategy(8),
    ) {
        let p = PRIMES[0].p;
        let p0i = PRIMES[0].p0i;
        let r2 = PRIMES[0].r2;

        let rust_r = pqe_hawk::keygen::zint31::zint_mod_small_unsigned(
            &d, d.len(), 1, p, p0i, r2,
        );
        let c_r = unsafe {
            test_zint_mod_small_unsigned(
                d.as_ptr(), d.len(), 1, p, p0i, r2,
            )
        };
        prop_assert_eq!(rust_r, c_r);
    }

    #[test]
    fn ffi_zint_add_mul_small_matches(
        x_initial in bigint_strategy(9),  // len+1 = 8+1
        y in bigint_strategy(8),
        s in limb31_strategy(),
    ) {
        let len = 8usize;
        let mut rust_x = x_initial.clone();
        pqe_hawk::keygen::zint31::zint_add_mul_small(
            &mut rust_x, len, 1, &y, s,
        );

        let mut c_x = x_initial.clone();
        unsafe {
            test_zint_add_mul_small(
                c_x.as_mut_ptr(), len, 1, y.as_ptr(), s,
            );
        }
        prop_assert_eq!(rust_x, c_x);
    }

    #[test]
    fn ffi_zint_norm_zero_matches(
        x_initial in bigint_strategy(8),
    ) {
        let p = {
            let mut pbuf = [0u32; 8];
            // Use PRIMES[0].p in limb 0 and zero elsewhere — positive odd integer,
            // valid for zint_norm_zero's comparison logic.
            pbuf[0] = PRIMES[0].p;
            pbuf
        };

        let mut rust_x = x_initial.clone();
        let xlen = rust_x.len();
        pqe_hawk::keygen::zint31::zint_norm_zero(&mut rust_x, xlen, 1, &p);

        let mut c_x = x_initial.clone();
        unsafe {
            test_zint_norm_zero(
                c_x.as_mut_ptr(), c_x.len(), 1, p.as_ptr(),
            );
        }
        prop_assert_eq!(rust_x, c_x);
    }
}

// === Sub-task 15.10b: zint_bezout + scaled add/sub cross-check ===

/// Strategy for a positive odd 8-limb big-int. We force the low bit of the
/// lowest limb to 1 to ensure oddness.
fn odd_bigint_strategy(n: usize) -> impl Strategy<Value = Vec<u32>> {
    bigint_strategy(n).prop_map(|mut v| {
        if !v.is_empty() {
            v[0] |= 1; // ensure odd
        }
        // Force top limb's sign bit (bit 30) to 0 to keep positive.
        if let Some(last) = v.last_mut() {
            *last &= 0x3FFF_FFFF;
            if *last == 0 {
                *last = 1; // ensure non-zero (positive)
            }
        }
        v
    })
}

proptest! {
    #[test]
    fn ffi_zint_bezout_matches(
        x in odd_bigint_strategy(8),
        y in odd_bigint_strategy(8),
    ) {
        let len = 8usize;
        let mut rust_u = vec![0u32; len];
        let mut rust_v = vec![0u32; len];
        let mut rust_tmp = vec![0u32; 4 * len];
        let rust_r = pqe_hawk::keygen::zint31::zint_bezout(
            &mut rust_u, &mut rust_v, &x, &y, len, &mut rust_tmp,
        );

        let mut c_u = vec![0u32; len];
        let mut c_v = vec![0u32; len];
        let mut c_tmp = vec![0u32; 4 * len];
        let c_r = unsafe {
            test_zint_bezout(
                c_u.as_mut_ptr(), c_v.as_mut_ptr(),
                x.as_ptr(), y.as_ptr(), len, c_tmp.as_mut_ptr(),
            )
        } as u32;
        prop_assert_eq!(rust_r, c_r);
        prop_assert_eq!(rust_u, c_u);
        prop_assert_eq!(rust_v, c_v);
    }

    #[test]
    fn ffi_zint_add_scaled_mul_small_matches(
        x_init in bigint_strategy(8),
        y in bigint_strategy(4),
        k in any::<i32>(),
        sch in 0u32..=3u32,
        scl in 0u32..=30u32,
    ) {
        let xlen = 8usize;
        let ylen = 4usize;
        let mut rust_x = x_init.clone();
        pqe_hawk::keygen::zint31::zint_add_scaled_mul_small(
            &mut rust_x, xlen, &y, ylen, 1, k, sch, scl,
        );

        let mut c_x = x_init.clone();
        unsafe {
            test_zint_add_scaled_mul_small(
                c_x.as_mut_ptr(), xlen, y.as_ptr(), ylen, 1, k, sch, scl,
            );
        }
        prop_assert_eq!(rust_x, c_x);
    }

    #[test]
    fn ffi_zint_sub_scaled_matches(
        x_init in bigint_strategy(8),
        y in bigint_strategy(4),
        sch in 0u32..=3u32,
        scl in 0u32..=30u32,
    ) {
        let xlen = 8usize;
        let ylen = 4usize;
        let mut rust_x = x_init.clone();
        pqe_hawk::keygen::zint31::zint_sub_scaled(
            &mut rust_x, xlen, &y, ylen, 1, sch, scl,
        );

        let mut c_x = x_init.clone();
        unsafe {
            test_zint_sub_scaled(
                c_x.as_mut_ptr(), xlen, y.as_ptr(), ylen, 1, sch, scl,
            );
        }
        prop_assert_eq!(rust_x, c_x);
    }
}

// === Sub-task 15.13: Fxr (Q32.32 fixed-point) cross-check ===

proptest! {
    #[test]
    fn ffi_fxr_mul_matches(x in any::<u64>(), y in any::<u64>()) {
        use pqe_hawk::keygen::fxr::Fxr;
        let rust = Fxr::of_scaled32(x).mul(Fxr::of_scaled32(y)).0;
        let c = unsafe { test_fxr_mul(x, y) };
        prop_assert_eq!(rust, c);
    }

    #[test]
    fn ffi_fxr_div_matches(
        x in any::<u64>(),
        // y must be non-zero to avoid division by zero.
        y in any::<u64>().prop_filter("y != 0", |&y| y != 0),
    ) {
        use pqe_hawk::keygen::fxr::Fxr;
        let rust = Fxr::of_scaled32(x).div(Fxr::of_scaled32(y)).0;
        let c = unsafe { test_fxr_div(x, y) };
        prop_assert_eq!(rust, c);
    }

    #[test]
    fn ffi_fxr_inv_matches(
        x in any::<u64>().prop_filter("x != 0", |&x| x != 0),
    ) {
        use pqe_hawk::keygen::fxr::Fxr;
        let rust = Fxr::of_scaled32(x).inv().0;
        let c = unsafe { test_fxr_inv(x) };
        prop_assert_eq!(rust, c);
    }

    #[test]
    fn ffi_fxr_round_matches(x in any::<u64>()) {
        use pqe_hawk::keygen::fxr::Fxr;
        let rust = Fxr::of_scaled32(x).round();
        // C returns int64_t; Rust returns i32. Cast C result to i32 for comparison.
        let c = unsafe { test_fxr_round(x) } as i32;
        prop_assert_eq!(rust, c);
    }
}

// === Sub-task 15.16: vect_FFT / vect_iFFT / vect_mul_fft cross-check ===

/// Strategy for a length-N vector of Fxr bit patterns.
fn fxr_vec_strategy(n: usize) -> impl Strategy<Value = Vec<u64>> {
    prop::collection::vec(any::<u64>(), n..=n)
}

proptest! {
    #[test]
    fn ffi_vect_fft_matches(
        f_bits in fxr_vec_strategy(512),
    ) {
        use pqe_hawk::keygen::fxr::Fxr;

        // Rust side: convert u64 bits to Fxr, run vect_fft, extract u64 bits.
        let mut rust_f: Vec<Fxr> = f_bits.iter().map(|&b| Fxr(b)).collect();
        pqe_hawk::keygen::fxp::vect_fft(9, &mut rust_f);
        let rust_out: Vec<u64> = rust_f.iter().map(|x| x.0).collect();

        // C side: pass the raw u64 buffer.
        let mut c_out = f_bits.clone();
        unsafe {
            test_vect_fft(9, c_out.as_mut_ptr());
        }

        prop_assert_eq!(rust_out, c_out);
    }

    #[test]
    fn ffi_vect_ifft_matches(
        f_bits in fxr_vec_strategy(512),
    ) {
        use pqe_hawk::keygen::fxr::Fxr;

        let mut rust_f: Vec<Fxr> = f_bits.iter().map(|&b| Fxr(b)).collect();
        pqe_hawk::keygen::fxp::vect_ifft(9, &mut rust_f);
        let rust_out: Vec<u64> = rust_f.iter().map(|x| x.0).collect();

        let mut c_out = f_bits.clone();
        unsafe {
            test_vect_ifft(9, c_out.as_mut_ptr());
        }

        prop_assert_eq!(rust_out, c_out);
    }

    #[test]
    fn ffi_vect_mul_fft_matches(
        a_bits in fxr_vec_strategy(512),
        b_bits in fxr_vec_strategy(512),
    ) {
        use pqe_hawk::keygen::fxr::Fxr;

        let mut rust_a: Vec<Fxr> = a_bits.iter().map(|&b| Fxr(b)).collect();
        let rust_b: Vec<Fxr> = b_bits.iter().map(|&b| Fxr(b)).collect();
        pqe_hawk::keygen::fxp::vect_mul_fft(9, &mut rust_a, &rust_b);
        let rust_out: Vec<u64> = rust_a.iter().map(|x| x.0).collect();

        let mut c_a = a_bits.clone();
        unsafe {
            test_vect_mul_fft(9, c_a.as_mut_ptr(), b_bits.as_ptr());
        }

        prop_assert_eq!(rust_out, c_a);
    }
}

// === Sub-task 15.18b: solve_NTRU_deepest cross-check ===

use pqe_hawk::keygen::ntru_profile::{NtruProfile, SOLVE_HAWK_512};

/// Strategy for a random i8 polynomial of exactly `n` coefficients, with
/// coefficients in [-8, 8] matching HAWK's binary Gaussian distribution.
fn small_i8_poly_strategy(n: usize) -> impl Strategy<Value = Vec<i8>> {
    prop::collection::vec(-8i8..=8i8, n..=n)
}

proptest! {
    #[test]
    fn ffi_solve_ntru_deepest_matches(
        f in small_i8_poly_strategy(512),
        g in small_i8_poly_strategy(512),
    ) {
        let prof = SOLVE_HAWK_512;
        let logn_top: u32 = 9;
        let tmp_size = 32 * 512usize;

        // Rust side.
        let mut rust_tmp = vec![0u32; tmp_size];
        let rust_r = pqe_hawk::keygen::ntru::solve_ntru_deepest(
            &prof, logn_top, &f, &g, &mut rust_tmp,
        );

        // C side.
        let mut c_tmp = vec![0u32; tmp_size];
        let c_r = unsafe {
            test_solve_ntru_deepest(
                &prof as *const NtruProfile as *const std::ffi::c_void,
                logn_top,
                f.as_ptr(),
                g.as_ptr(),
                c_tmp.as_mut_ptr(),
            )
        };

        prop_assert_eq!(rust_r, c_r);
        // Only compare the output region on a successful solve — error paths
        // may leave scratch values in tmp that differ.
        if rust_r == 0 {
            // F and G are each len = max_bl_small[logn_top] u32 words.
            let len = prof.max_bl_small[logn_top as usize] as usize;
            prop_assert_eq!(&rust_tmp[..2 * len], &c_tmp[..2 * len]);
        }
    }
}

// === Sub-task 15.19d: solve_NTRU_intermediate (via full pipeline from deepest) ===

proptest! {
    #[test]
    fn ffi_solve_ntru_intermediate_depth_8_matches(
        f in small_i8_poly_strategy(512),
        g in small_i8_poly_strategy(512),
    ) {
        let prof = SOLVE_HAWK_512;
        let logn_top: u32 = 9;
        let depth: u32 = 8;
        let n = 512usize;
        // Generous tmp size for the pipeline.
        let tmp_size = 64 * n;

        // Rust side.
        let mut rust_tmp = vec![0u32; tmp_size];
        let rust_r_deepest = pqe_hawk::keygen::ntru::solve_ntru_deepest(
            &prof, logn_top, &f, &g, &mut rust_tmp,
        );
        let rust_r = if rust_r_deepest == 0 {
            pqe_hawk::keygen::ntru::solve_ntru_intermediate(
                &prof, logn_top, &f, &g, depth, &mut rust_tmp,
            )
        } else {
            rust_r_deepest
        };

        // C side.
        let mut c_tmp = vec![0u32; tmp_size];
        let c_r = unsafe {
            test_solve_ntru_intermediate_after_deepest(
                &prof as *const NtruProfile as *const std::ffi::c_void,
                logn_top, f.as_ptr(), g.as_ptr(), depth,
                c_tmp.as_mut_ptr(),
            )
        };

        prop_assert_eq!(rust_r, c_r);
        // On success, compare the output F at the start of tmp.
        // At depth 8, n at this level = 1 shifted, but the F output is a
        // polynomial of degree 2^(logn_top - depth) = 2^1 = 2, with llen
        // words per coefficient. Compare the first llen*2 words.
        if rust_r == 0 {
            // At depth=8, logn=1, n=2, and llen = max_bl_large[8] = 121.
            // So output is 2 * 121 = 242 words.
            let llen = prof.max_bl_large[depth as usize] as usize;
            let output_len = (1usize << (logn_top - depth)) * llen;
            prop_assert_eq!(&rust_tmp[..output_len], &c_tmp[..output_len]);
        }
    }
}

// === Sub-task 15.21: solve_NTRU full pipeline cross-check ===

proptest! {
    #![proptest_config(proptest::test_runner::Config {
        cases: 20,
        ..Default::default()
    })]
    #[test]
    fn ffi_solve_ntru_matches(
        f in small_i8_poly_strategy(512),
        g in small_i8_poly_strategy(512),
    ) {
        let prof = SOLVE_HAWK_512;
        let logn: u32 = 9;
        let n = 512usize;
        let tmp_size = 64 * n;

        let mut rust_tmp = vec![0u32; tmp_size];
        let rust_r = pqe_hawk::keygen::ntru::solve_ntru(
            &prof, logn, &f, &g, &mut rust_tmp,
        );

        let mut c_tmp = vec![0u32; tmp_size];
        let c_r = unsafe {
            test_solve_ntru(
                &prof as *const NtruProfile as *const std::ffi::c_void,
                logn, f.as_ptr(), g.as_ptr(),
                c_tmp.as_mut_ptr(),
            )
        };

        prop_assert_eq!(rust_r, c_r);
        if rust_r == 0 {
            // F and G packed as i8 into tmp[0..2*n/4] words (2*n bytes total).
            let packed_words = 2 * n / 4; // n=512 -> 256 words
            prop_assert_eq!(&rust_tmp[..packed_words], &c_tmp[..packed_words]);
        }
    }
}

// === Sub-task 15.22b: hawk_regen_fg cross-check ===

proptest! {
    #[test]
    fn ffi_hawk_regen_fg_matches(
        seed in prop::collection::vec(any::<u8>(), 24..=24),
    ) {
        let mut rust_f = vec![0i8; 512];
        let mut rust_g = vec![0i8; 512];
        pqe_hawk::keygen::regen_fg::hawk_regen_fg(&mut rust_f, &mut rust_g, &seed);

        let mut c_f = vec![0i8; 512];
        let mut c_g = vec![0i8; 512];
        unsafe {
            test_hawk_regen_fg(
                9,
                c_f.as_mut_ptr(),
                c_g.as_mut_ptr(),
                seed.as_ptr() as *const std::ffi::c_void,
            );
        }
        prop_assert_eq!(rust_f, c_f);
        prop_assert_eq!(rust_g, c_g);
    }
}

// === Sub-task 15.9: zint_rebuild_crt and zint_negate cross-check ===

proptest! {
    #[test]
    fn ffi_zint_negate_matches(
        a in bigint_strategy(8),
        ctl in 0u32..=1u32,
    ) {
        let mut rust_a = a.clone();
        pqe_hawk::keygen::zint31::zint_negate(&mut rust_a, ctl);

        let mut c_a = a.clone();
        unsafe {
            test_zint_negate(c_a.as_mut_ptr(), c_a.len(), ctl);
        }
        prop_assert_eq!(rust_a, c_a);
    }

    #[test]
    fn ffi_zint_rebuild_crt_matches(
        // 4 interleaved big-ints (n=4), each of xlen=4 limbs, single set.
        // Total storage: 1 * 4 * 4 = 16 u32 limbs.
        xx in bigint_strategy(16),
        normalize in 0u32..=1u32,
    ) {
        let xlen = 4usize;
        let n = 4usize;
        let num_sets = 1usize;

        // We need each limb to be a valid residue mod the corresponding prime.
        // The C relies on limb u being in [0, PRIMES[u].p). Reduce accordingly.
        let mut normalized_xx = xx.clone();
        for u in 0..xlen {
            let p = PRIMES[u].p;
            for v in 0..n {
                // Slot is xx[u * n + v] (stride n across limbs, offset v for the v-th bigint).
                let idx = u * n + v;
                normalized_xx[idx] = normalized_xx[idx] % p;
            }
        }

        let mut rust_xx = normalized_xx.clone();
        let mut rust_tmp = vec![0u32; xlen];
        pqe_hawk::keygen::zint31::zint_rebuild_crt(
            &mut rust_xx, xlen, n, num_sets, normalize != 0, &mut rust_tmp,
        );

        let mut c_xx = normalized_xx.clone();
        let mut c_tmp = vec![0u32; xlen];
        unsafe {
            test_zint_rebuild_crt(
                c_xx.as_mut_ptr(), xlen, n, num_sets, normalize as i32,
                c_tmp.as_mut_ptr(),
            );
        }
        prop_assert_eq!(rust_xx, c_xx);
        prop_assert_eq!(rust_tmp, c_tmp);
    }
}

// === Sub-task 15.24: Hawk_keygen end-to-end cross-check ===

proptest! {
    #![proptest_config(proptest::test_runner::Config {
        cases: 5,
        ..Default::default()
    })]
    #[test]
    fn ffi_hawk_keygen_matches(
        seed_bytes in prop::collection::vec(any::<u8>(), 65_536..=65_536),
    ) {
        let logn = 9u32;
        let n = 512usize;

        // Rust side: consume bytes from seed_bytes sequentially via closure.
        let rust_rng_bytes = seed_bytes.clone();
        let mut rust_pos = 0usize;
        let rust_kp = pqe_hawk::keygen::hawk_keygen::hawk_keygen_512(|buf| {
            buf.copy_from_slice(&rust_rng_bytes[rust_pos..rust_pos + buf.len()]);
            rust_pos += buf.len();
        });

        // C side: same byte stream fed via test_hawk_keygen.
        let mut c_f = vec![0i8; n];
        let mut c_g = vec![0i8; n];
        let mut c_f_cap = vec![0i8; n];
        let mut c_g_cap = vec![0i8; n];
        let mut c_q00 = vec![0i16; n];
        let mut c_q01 = vec![0i16; n];
        let mut c_q11 = vec![0i32; n];
        let mut c_seed = vec![0u8; 24];
        let mut c_tmp = vec![0u32; 32 * n];
        let c_r = unsafe {
            test_hawk_keygen(
                logn,
                seed_bytes.as_ptr(), seed_bytes.len(),
                c_f.as_mut_ptr(), c_g.as_mut_ptr(),
                c_f_cap.as_mut_ptr(), c_g_cap.as_mut_ptr(),
                c_q00.as_mut_ptr(), c_q01.as_mut_ptr(), c_q11.as_mut_ptr(),
                c_seed.as_mut_ptr(),
                c_tmp.as_mut_ptr(), c_tmp.len() * 4,
            )
        };
        // Emit both byte positions for diagnosing retry-count divergence.
        prop_assert_eq!(c_r, 0, "C keygen should succeed (rust_bytes_consumed={})", rust_pos);
        // Check that both consumed the same number of bytes (same retry count).
        prop_assert_eq!(rust_kp.seed, c_seed.clone(),
            "seed mismatch — likely retry-count divergence (rust_bytes={} c_seed={:?})",
            rust_pos, &c_seed[..4]);
        prop_assert_eq!(rust_kp.f, c_f, "f mismatch");
        prop_assert_eq!(rust_kp.g, c_g, "g mismatch");
        prop_assert_eq!(rust_kp.f_cap, c_f_cap, "F mismatch");
        prop_assert_eq!(rust_kp.g_cap, c_g_cap, "G mismatch");
        prop_assert_eq!(rust_kp.q00, c_q00, "q00 mismatch");
        prop_assert_eq!(rust_kp.q01, c_q01, "q01 mismatch");
        prop_assert_eq!(rust_kp.q11, c_q11, "q11 mismatch");
    }
}

// === Task 19: HawkKeypair::generate + to_bytes() cross-check ===
//
// Validates that the Rust public-API wrapper (HawkKeypair::generate) plus the
// Task 17 serialization (encode_public / encode_private) produce byte-identical
// output to the reference C's hawk_keygen (Zh(keygen)) + encode_public +
// encode_private pipeline.

proptest! {
    #![proptest_config(proptest::test_runner::Config {
        cases: 3,
        ..Default::default()
    })]
    #[test]
    fn ffi_hawk_keypair_encoded_matches(
        seed_bytes in prop::collection::vec(any::<u8>(), 65_536..=65_536),
    ) {
        // Rust side: HawkKeypair::generate via ByteBufRng, then to_bytes().
        let mut rng = ByteBufRng { buf: &seed_bytes, pos: 0 };
        let kp = pqe_hawk::HawkKeypair::generate(&mut rng);
        let rust_pk = kp.public.to_bytes().expect("encode_public");
        let rust_sk = kp.secret.to_bytes();

        // C side: same byte stream fed to hawk_keygen (the public API that
        // wraps Hawk_keygen + encode_public + encode_private in a retry loop).
        let mut c_sk = vec![0u8; 184];
        let mut c_pk = vec![0u8; 1024];
        // hawk_keygen needs HAWK_TMPSIZE_KEYGEN(9) = 22*512+7 = 11271 bytes,
        // but the public-api wrapper also stores f,g,F,G,q00,q01,q11,seed,tpriv,tpub
        // inside tmp before calling the inner keygen. 64*512 = 32768 bytes is safe.
        let mut c_tmp = vec![0u8; 64 * 512];
        let c_r = unsafe {
            test_hawk_keygen_encoded(
                seed_bytes.as_ptr(), seed_bytes.len(),
                c_sk.as_mut_ptr(), c_pk.as_mut_ptr(),
                c_tmp.as_mut_ptr() as *mut std::ffi::c_void, c_tmp.len(),
            )
        };
        prop_assert_eq!(c_r, 1, "C hawk_keygen should succeed (returns 1)");
        prop_assert_eq!(rust_pk.as_slice(), c_pk.as_slice(), "public key mismatch");
        prop_assert_eq!(rust_sk.as_slice(), c_sk.as_slice(), "secret key mismatch");
    }
}

// === Sub-task 22a: bp_mulmod_512 and basis_m2_mul cross-check ===

fn u8_vec_64() -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(any::<u8>(), 64..=64)
}

proptest! {
    #[test]
    fn ffi_bp_mulmod_512_matches(
        a in u8_vec_64(),
        b in u8_vec_64(),
    ) {
        let mut rust_d = [0u8; 64];
        let mut rust_tmp = [0u8; 160];
        pqe_hawk::sign::bp::bp_mulmod_512(&mut rust_d, &a, &b, &mut rust_tmp);

        let mut c_d = [0u8; 64];
        let mut c_tmp = [0u8; 160];
        unsafe {
            test_bp_mulmod_512(
                c_d.as_mut_ptr(),
                a.as_ptr(), b.as_ptr(),
                c_tmp.as_mut_ptr(),
            );
        }
        prop_assert_eq!(rust_d, c_d);
    }

    #[test]
    fn ffi_basis_m2_mul_matches(
        h0 in u8_vec_64(),
        h1 in u8_vec_64(),
        f2 in u8_vec_64(),
        g2 in u8_vec_64(),
        f_cap2 in u8_vec_64(),
        g_cap2 in u8_vec_64(),
    ) {
        let mut rust_t0 = [0u8; 64];
        let mut rust_t1 = [0u8; 64];
        let mut rust_tmp = [0u8; 224];
        pqe_hawk::sign::bp::basis_m2_mul(
            &mut rust_t0, &mut rust_t1,
            &h0, &h1, &f2, &g2, &f_cap2, &g_cap2,
            &mut rust_tmp,
        );

        let mut c_t0 = [0u8; 64];
        let mut c_t1 = [0u8; 64];
        let mut c_tmp = [0u8; 224];
        unsafe {
            test_basis_m2_mul(
                c_t0.as_mut_ptr(), c_t1.as_mut_ptr(),
                h0.as_ptr(), h1.as_ptr(),
                f2.as_ptr(), g2.as_ptr(),
                f_cap2.as_ptr(), g_cap2.as_ptr(),
                c_tmp.as_mut_ptr(),
            );
        }
        prop_assert_eq!(rust_t0, c_t0);
        prop_assert_eq!(rust_t1, c_t1);
    }
}

// === Sub-task 24.5: fx32_FFT / fx32_iFFT cross-check ===

proptest! {
    #[test]
    fn ffi_fx32_fft_matches(
        v in prop::collection::vec(any::<u32>(), 512..=512),
    ) {
        // Rust side.
        let mut rust_a = v.clone();
        pqe_hawk::verify::fx32::fx32_fft(9, &mut rust_a);

        // C side.
        let mut c_a = v.clone();
        unsafe {
            test_fx32_fft(9, c_a.as_mut_ptr());
        }

        prop_assert_eq!(rust_a, c_a);
    }

    #[test]
    fn ffi_fx32_ifft_matches(
        v in prop::collection::vec(any::<u32>(), 512..=512),
    ) {
        let mut rust_a = v.clone();
        pqe_hawk::verify::fx32::fx32_ifft(9, &mut rust_a);

        let mut c_a = v.clone();
        unsafe {
            test_fx32_ifft(9, c_a.as_mut_ptr());
        }

        prop_assert_eq!(rust_a, c_a);
    }
}

// === Sub-task 24.9: verify's mp_div / mp_ntt_autoadj / mp_poly_to_ntt[_autoadj] ===

proptest! {
    #[test]
    fn ffi_verify_mp_div_matches(
        x in 0u32..pqe_hawk::verify::consts::P1,
        y in 0u32..pqe_hawk::verify::consts::P1,
    ) {
        let p = pqe_hawk::verify::consts::P1;
        let p0i = pqe_hawk::verify::consts::P1_0I;
        let m16 = pqe_hawk::verify::consts::P1_M16;

        let rust = pqe_hawk::verify::mp::mp_div(x, y, p, p0i, m16);
        let c = unsafe { test_verify_mp_div(x, y, p, p0i, m16) };
        prop_assert_eq!(rust, c);
    }

    #[test]
    fn ffi_mp_poly_to_ntt_matches_p1(
        a in prop::collection::vec(any::<i16>(), 512..=512),
    ) {
        use pqe_hawk::verify::consts::{P1, P1_0I};
        use pqe_hawk::verify::mp::{GM_P1, mp_poly_to_ntt};

        let mut rust_d = vec![0u32; 512];
        mp_poly_to_ntt(9, &mut rust_d, &a, P1, P1_0I, &GM_P1);

        let mut c_d = vec![0u32; 512];
        unsafe {
            test_mp_poly_to_ntt(
                9, c_d.as_mut_ptr(), a.as_ptr(),
                P1, P1_0I, GM_P1.as_ptr(),
            );
        }
        prop_assert_eq!(rust_d, c_d);
    }

    #[test]
    fn ffi_mp_poly_to_ntt_autoadj_matches_p1(
        a_half in prop::collection::vec(any::<i16>(), 256..=256),
    ) {
        use pqe_hawk::verify::consts::{P1, P1_0I};
        use pqe_hawk::verify::mp::{GM_P1, mp_poly_to_ntt_autoadj};

        let mut rust_d = vec![0u32; 256];
        mp_poly_to_ntt_autoadj(9, &mut rust_d, &a_half, P1, P1_0I, &GM_P1);

        let mut c_d = vec![0u32; 256];
        unsafe {
            test_mp_poly_to_ntt_autoadj(
                9, c_d.as_mut_ptr(), a_half.as_ptr(),
                P1, P1_0I, GM_P1.as_ptr(),
            );
        }
        prop_assert_eq!(rust_d, c_d);
    }
}

// === Sub-task 15.23b: make_q001 cross-check ===

proptest! {
    #[test]
    fn ffi_make_q001_matches(
        f in small_i8_poly_strategy(512),
        g in small_i8_poly_strategy(512),
        f_cap in prop::collection::vec(-127i8..=127i8, 512..=512),
        g_cap in prop::collection::vec(-127i8..=127i8, 512..=512),
    ) {
        let logn = 9u32;
        let lim00 = 1i32 << 9;
        let lim01 = 1i32 << 11;
        let lim11 = 1i32 << 13;

        // Rust side.
        let mut rust_tmp = vec![0u32; 16 * 512];
        let rust_r = pqe_hawk::keygen::q_derive::make_q001(
            logn, lim00, lim01, lim11,
            &f, &g, &f_cap, &g_cap,
            &mut rust_tmp,
        );

        // C side.
        let mut c_tmp = vec![0u32; 16 * 512];
        let mut c_q00 = vec![0i16; 512];
        let mut c_q01 = vec![0i16; 512];
        let mut c_q11 = vec![0i32; 512];
        let c_r = unsafe {
            test_make_q001(
                logn,
                lim00,
                lim01,
                lim11,
                f.as_ptr(),
                g.as_ptr(),
                f_cap.as_ptr(),
                g_cap.as_ptr(),
                c_q00.as_mut_ptr(),
                c_q01.as_mut_ptr(),
                c_q11.as_mut_ptr(),
                c_tmp.as_mut_ptr(),
            )
        };

        let rust_ok = rust_r.is_ok();
        match (rust_r, c_r) {
            (Ok(rust_out), 1) => {
                prop_assert_eq!(rust_out.q00, c_q00, "q00 mismatch");
                prop_assert_eq!(rust_out.q01, c_q01, "q01 mismatch");
                prop_assert_eq!(rust_out.q11, c_q11, "q11 mismatch");
            }
            (Err(()), 0) => {
                // Both reported failure; agree.
            }
            _ => {
                prop_assert!(
                    false,
                    "Rust {:?} vs C return {}",
                    rust_ok,
                    c_r
                );
            }
        }
    }
}

// === Sub-task 24.13: HawkPublicKey::verify + FFI cross-check ===

proptest! {
    #![proptest_config(proptest::test_runner::Config {
        cases: 5,
        ..Default::default()
    })]

    /// Rust signs a message, C verifies via test_hawk_verify.
    /// Expect all C returns to be 1 (valid).
    #[test]
    fn ffi_verify_accepts_rust_signature(
        key_seed in any::<[u8; 32]>(),
        msg in prop::collection::vec(any::<u8>(), 0..=256usize),
    ) {
        use rand::SeedableRng;
        use rand_chacha::ChaCha20Rng;

        let mut rng = ChaCha20Rng::from_seed(key_seed);
        let kp = pqe_hawk::HawkKeypair::generate(&mut rng);

        let mut sign_rng = ChaCha20Rng::from_seed([99u8; 32]);
        let sig = kp.secret.sign(&msg, &mut sign_rng).expect("sign ok");

        let pk_bytes = kp.public.to_bytes().unwrap();
        let sig_bytes = sig.to_bytes().unwrap();

        // tmp: 10*512*4 = 20480 bytes + 8 alignment headroom.
        let mut tmp = vec![0u8; 10 * 512 * 4 + 8];
        let c_result = unsafe {
            test_hawk_verify(
                pk_bytes.as_ptr(),
                pk_bytes.len(),
                msg.as_ptr(),
                msg.len(),
                sig_bytes.as_ptr(),
                sig_bytes.len(),
                tmp.as_mut_ptr() as *mut std::ffi::c_void,
                tmp.len(),
            )
        };
        prop_assert_eq!(c_result, 1, "C rejected Rust signature");
    }

    /// Rust signs a message, Rust verifies with HawkPublicKey::verify.
    /// Self-consistency check: verify must accept every signature sign produces.
    #[test]
    fn rust_verify_accepts_rust_signature(
        key_seed in any::<[u8; 32]>(),
        msg in prop::collection::vec(any::<u8>(), 0..=256usize),
    ) {
        use rand::SeedableRng;
        use rand_chacha::ChaCha20Rng;

        let mut rng = ChaCha20Rng::from_seed(key_seed);
        let kp = pqe_hawk::HawkKeypair::generate(&mut rng);

        let mut sign_rng = ChaCha20Rng::from_seed([99u8; 32]);
        let sig = kp.secret.sign(&msg, &mut sign_rng).expect("sign ok");

        prop_assert!(
            kp.public.verify(&msg, &sig).is_ok(),
            "Rust verify rejected Rust signature"
        );
    }
}
