//! Wire-format encoding and decoding for HAWK-512 keys.
//!
//! Byte layouts match the reference C (`hawk_kgen.c` + `hawk_vrfy.c`) so that
//! serialized keys interchange byte-for-byte with the NIST KAT vectors.
//!
//! Key encodings (HAWK-512):
//!   - Public key (1024 bytes): Golomb-Rice encoded q00 (half-size) +
//!     Golomb-Rice encoded q01 (full-size), zero-padded to 1024.
//!   - Secret key (184 bytes): kgseed(24) || F_mod2(64) || G_mod2(64) || hpub(32).

use crate::error::HawkError;
use crate::params::{HAWK_N, HAWK_PUBLIC_KEY_BYTES, HAWK_SECRET_KEY_BYTES};

const HAWK_512_LOWBITS_Q00: u32 = 5;
const HAWK_512_LOWBITS_Q01: u32 = 9;
/// Extra bits for q00[0]: 16 - (low00 + 4) = 16 - (5 + 4) = 7.
const HAWK_512_EB00_LEN: u32 = 7;
const HAWK_512_HPUB_LEN: usize = 32;
const HAWK_512_KGSEED_LEN: usize = 24;

/// Pack n int8 low bits into n/8 bytes (little-endian within each byte).
/// Port of `extract_lowbit` (hawk_kgen.c:7-21).
pub fn extract_lowbit(f: &[i8], out: &mut [u8]) {
    let n = f.len();
    assert_eq!(out.len(), n / 8);
    for u in (0..n).step_by(8) {
        let fu = f[u..u + 8].iter().map(|&x| x as u8);
        let mut b = 0u8;
        for (i, byte) in fu.enumerate() {
            b |= (byte & 1) << i;
        }
        out[u / 8] = b;
    }
}

/// Top-bit mask: all-ones if x has bit 31 set, else 0.
/// Matches the C `tbmask` macro used in hawk_kgen.c and hawk_vrfy.c.
#[inline]
fn tbmask(x: u32) -> u32 {
    ((x as i32) >> 31) as u32
}

/// Golomb-Rice encoder for a polynomial of length `1 << logn` int16 values.
/// `low` = fixed-size prefix bits per coefficient.
///
/// On success writes bytes into `dst` starting at offset 0, and returns
/// `(written_bytes, num_ignored_bits_in_last_byte)`.
/// On buffer overflow returns `Err(())`.
///
/// Layout written: sign_bits (n/8 bytes) || fixed_parts (low*n/8 bytes) ||
///                 variable_parts (unary tails, variable length).
///
/// Port of `encode_gr` (hawk_kgen.c:86-178).
pub(crate) fn encode_gr(
    logn: u32,
    dst: &mut [u8],
    a: &[i16],
    low: u32,
) -> Result<(usize, u32), ()> {
    let n = 1usize << logn;
    let dst_len = dst.len();
    let min_len = ((low + 1) as usize) << (logn - 3);
    if dst_len < min_len {
        return Err(());
    }
    let mut out_off: usize = 0;

    // --- Sign bits: n/8 bytes ---
    for u in (0..n).step_by(8) {
        let mut x = 0u8;
        for v in 0..8 {
            // Extract the sign bit (MSB of the i16, treated as u16).
            x |= ((a[u + v] as u16 >> 15) as u8) << v;
        }
        dst[u >> 3] = x;
    }
    out_off += n >> 3;

    // --- Fixed-size parts: low*n/8 bytes total ---
    let low_mask: u32 = (1u32 << low) - 1;
    if low <= 8 {
        // Each group of 8 coefficients contributes low*8 bits = low bytes.
        for u in (0..n).step_by(8) {
            let mut x: u64 = 0;
            let mut vv: u32 = 0;
            for v in 0..8 {
                let mut w = a[u + v] as u32; // sign-extended to 32 bits
                w ^= tbmask(w); // if negative: w = |w|-1 (two's complement abs - 1)
                x |= ((w & low_mask) as u64) << vv;
                vv += low;
            }
            // Write low bytes (low*8 bits, one per coefficient group).
            for i in 0..(low as usize) {
                dst[out_off] = (x >> (i * 8)) as u8;
                out_off += 1;
            }
        }
    } else {
        // low > 8: split into two halves of 4 coefficients each.
        for u in (0..n).step_by(8) {
            let mut x0: u64 = 0;
            let mut vv: u32 = 0;
            for v in 0..4 {
                let mut w = a[u + v] as u32;
                w ^= tbmask(w);
                x0 |= ((w & low_mask) as u64) << vv;
                vv += low;
            }
            let mut x1: u64 = 0;
            vv = 0;
            for v in 4..8 {
                let mut w = a[u + v] as u32;
                w ^= tbmask(w);
                x1 |= ((w & low_mask) as u64) << vv;
                vv += low;
            }
            // Merge: x0 holds 4*low bits; shift x1 up to pack after x0.
            // low << 2 = 4*low (bits used in x0).
            x0 |= x1 << (low << 2);
            x1 >>= 64 - (low << 2);
            // Write first 8 bytes from x0.
            for i in 0..8usize {
                dst[out_off] = (x0 >> (i * 8)) as u8;
                out_off += 1;
            }
            // Write remaining bytes from x1: (low*8 - 64) / 8 = low - 8 bytes.
            let extra_bytes = (low as i32 - 8) as usize;
            for i in 0..extra_bytes {
                dst[out_off] = (x1 >> (i * 8)) as u8;
                out_off += 1;
            }
        }
    }

    // --- Variable-size parts: unary prefix (k zeros then a 1) per coefficient ---
    let mut acc: u32 = 0;
    let mut acc_len: i32 = 0;
    for u in 0..n {
        let w = a[u] as u32;
        // k = magnitude >> low (unary count).
        let k = ((w ^ tbmask(w)) >> low) as i32;
        // Write unary code: k zeros then a 1 (bit at position acc_len + k).
        acc |= 1u32 << (acc_len + k);
        acc_len += 1 + k;
        while acc_len >= 8 {
            if out_off >= dst_len {
                return Err(());
            }
            dst[out_off] = acc as u8;
            out_off += 1;
            acc >>= 8;
            acc_len -= 8;
        }
    }
    if acc_len > 0 {
        if out_off >= dst_len {
            return Err(());
        }
        dst[out_off] = acc as u8;
        out_off += 1;
    }
    // Number of unused (padding) bits in the last byte.
    let num_ignored = ((-acc_len) & 7) as u32;
    Ok((out_off, num_ignored))
}

/// Encode a HAWK-512 public key.
/// `q00` and `q01` must each have length HAWK_N (512).
/// Writes exactly `HAWK_PUBLIC_KEY_BYTES` (1024) bytes.
///
/// Port of `encode_public` (hawk_kgen.c:189-247) for HAWK-512 (logn=9).
pub fn encode_public(q00: &[i16], q01: &[i16]) -> Result<[u8; HAWK_PUBLIC_KEY_BYTES], HawkError> {
    assert_eq!(q00.len(), HAWK_N);
    assert_eq!(q01.len(), HAWK_N);

    let logn = 9u32;
    let low00 = HAWK_512_LOWBITS_Q00;
    let low01 = HAWK_512_LOWBITS_Q01;
    let eb00_len = HAWK_512_EB00_LEN;

    let mut out = [0u8; HAWK_PUBLIC_KEY_BYTES];

    // --- Encode q00 (half-size: only q00[0..n/2]) ---
    // q00 is auto-adjoint; only the first n/2 coefficients are encoded.
    // q00[0] is specially downscaled: encode q00[0] >> eb00_len, then
    // append the low eb00_len bits separately.
    let sav_q00_0 = q00[0];
    let mut q00_half: Vec<i16> = q00[..HAWK_N / 2].to_vec();
    q00_half[0] = sav_q00_0 >> (eb00_len as i16); // arithmetic right-shift

    let (mut len00, ni) = encode_gr(logn - 1, &mut out[..], &q00_half, low00)
        .map_err(|()| HawkError::MalformedPublicKey("encode_gr q00 overflow".into()))?;

    // Append extra-bits of q00[0] (the low eb00_len bits).
    let eb00 = (sav_q00_0 as u32) & ((1u32 << eb00_len) - 1);
    if eb00_len <= ni {
        // The extra bits fit in the unused bits of the last byte.
        out[len00 - 1] |= (eb00 << (8 - ni)) as u8;
    } else {
        // Need one more byte for the overflow bits.
        if len00 >= HAWK_PUBLIC_KEY_BYTES {
            return Err(HawkError::MalformedPublicKey(
                "q00 extra bits overflow buffer".into(),
            ));
        }
        out[len00 - 1] |= (eb00 << (8 - ni)) as u8;
        out[len00] = (eb00 >> ni) as u8;
        len00 += 1;
    }

    // --- Encode q01 (full-size: q01[0..n]) ---
    let (_len01, _ni01) = encode_gr(logn, &mut out[len00..], q01, low01)
        .map_err(|()| HawkError::MalformedPublicKey("encode_gr q01 overflow".into()))?;

    // The remainder of `out` is already zeroed (array initialised to 0).
    Ok(out)
}

/// Golomb-Rice decoder.
///
/// Reads `1 << logn` int16 values from `buf`, writing into `d`.
/// `low` = fixed-size bits per coefficient; `lim_bits` = upper bound on
/// unary count (lim_hi = 1 << (lim_bits - low)).
///
/// The layout in `buf` is: sign_bits || fixed_parts || variable_parts.
/// The decoder reads variable_parts first (from offset (low+1)*n/8 onward),
/// then goes back to read sign_bits and fixed_parts from the start.
///
/// Returns `(consumed_bytes, num_ignored_bits_in_last_byte)` on success,
/// or `Err(())` on malformed input.
///
/// Port of `decode_gr` (hawk_vrfy.c:1215-1323).
pub(crate) fn decode_gr(
    logn: u32,
    d: &mut [i16],
    buf: &[u8],
    low: u32,
    lim_bits: u32,
) -> Result<(usize, u32), ()> {
    let n = 1usize << logn;
    let buf_len = buf.len();
    let min_len = ((low + 1) as usize) << (logn - 3);
    if buf_len < min_len {
        return Err(());
    }
    // Variable-part bytes start right after sign_bits + fixed_parts.
    let mut voff: usize = min_len;

    // --- NTZ (number of trailing zeros) lookup table for 8-bit values ---
    // ntz[0] = 8 (sentinel: all-zeros byte has 8 trailing zeros).
    static NTZ: [u8; 256] = {
        let mut t = [0u8; 256];
        t[0] = 8;
        let mut i = 1usize;
        while i < 256 {
            let mut c = 0u8;
            let mut x = i as u8;
            while x & 1 == 0 {
                x >>= 1;
                c += 1;
            }
            t[i] = c;
            i += 1;
        }
        t
    };

    // --- Variable part: unary prefixes ---
    // Each coefficient u has unary code: k trailing zeros then a 1-bit.
    // k gives the high bits: magnitude >> low = k, so magnitude = k<<low | lp.
    let mut acc: u32 = 0;
    let mut acc_off: i32 = 0;
    let lim_hi: i32 = 1i32 << (lim_bits - low);
    for u in 0..n {
        // Refill accumulator while it's empty.
        while acc == 0 {
            if acc_off >= lim_hi {
                return Err(());
            }
            if voff >= buf_len {
                return Err(());
            }
            acc |= (buf[voff] as u32) << acc_off;
            voff += 1;
            acc_off += 8;
        }
        // Count trailing zeros in the low 8 bits of acc.
        //
        // Soundness note: the first-byte NTZ returns k in [0, 8]. For k < 8
        // there is no explicit k >= lim_hi rejection here — that's safe for
        // every HAWK parameter set because lim_hi ∈ {8, 16}:
        //
        //   decode_gr_5_9  (q00, s1):  lim_hi = 16, so any k < 8 is valid.
        //   decode_gr_9_12 (q01):      lim_hi = 8,  and k < 8 is the only path
        //                              that can succeed (k == 8 falls into the
        //                              secondary check, k > 8 via second byte
        //                              triggers rejection there).
        //
        // This matches the reference C exactly (hawk_vrfy.c:1265-1271); the
        // invariant `lim_hi >= 8` is what makes the single-byte branch
        // trivially safe. Any future parameter set with lim_hi < 8 would need
        // an additional guard here.
        let mut k = NTZ[(acc & 0xFF) as usize] as i32;
        if k == 8 {
            // Low byte was all zero; look at next byte.
            k += NTZ[((acc >> 8) & 0xFF) as usize] as i32;
            if k >= lim_hi {
                return Err(());
            }
        }
        // Store the high part of the decoded magnitude.
        d[u] = ((k as u16) << low) as i16;
        acc >>= (k + 1) as u32;
        acc_off -= k + 1;
    }
    let num_ignored = acc_off as u32;

    // --- Sign bits and fixed-width parts ---
    // Sign bits:  buf[0 .. n/8]
    // Fixed parts: buf[n/8 .. (low+1)*n/8]
    let mut loff: usize = 1 << (logn - 3); // skip over sign-bit bytes
    let lmask: u32 = (1u32 << low) - 1;

    if low <= 8 {
        for u in (0..n).step_by(8) {
            let sbb = buf[u >> 3] as u32;
            // Read low*8 bits = low bytes for this group of 8 coefficients.
            let mut lpp: u64 = 0;
            for j in 0..(low as usize) {
                lpp |= (buf[loff] as u64) << (j * 8);
                loff += 1;
            }
            // Apply sign and fixed-part bits to each coefficient.
            for i in 0..8usize {
                let lp = ((lpp >> (i as u32 * low)) as u32) & lmask;
                let sm = ((sbb >> i) & 1).wrapping_neg(); // 0 or 0xFFFF_FFFF
                                                          // d[u+i] starts as k<<low (from variable part).
                                                          // XOR trick: if positive, d ^= lp; if negative, d ^= 0xFFFF ^ lp.
                                                          // This gives: positive → (k<<low)|lp, negative → ~((k<<low)|lp).
                let cur = d[u + i] as u16 as u32;
                d[u + i] = (cur ^ (sm & 0xFFFF) ^ lp) as u16 as i16;
            }
        }
    } else {
        // low > 8 (e.g., low=9 for q01 in HAWK-512).
        for u in (0..n).step_by(8) {
            let sbb = buf[u >> 3] as u32;
            // Read first 4 bytes for lpp0 (covers 3 full coefficients + partial 4th).
            let lpp0 = u32::from_le_bytes([buf[loff], buf[loff + 1], buf[loff + 2], buf[loff + 3]]);
            loff += 4;
            // Read remaining (low - 4) bytes for lpp1.
            let mut lpp1: u64 = 0;
            for k in 0..(low as usize - 4) {
                lpp1 |= (buf[loff] as u64) << (k * 8);
                loff += 1;
            }
            // First 3 coefficients: fit entirely within lpp0.
            for i in 0..3usize {
                let lp = (lpp0 >> (i as u32 * low)) & lmask;
                let sm = ((sbb >> i) & 1).wrapping_neg();
                let cur = d[u + i] as u16 as u32;
                d[u + i] = (cur ^ (sm & 0xFFFF) ^ lp) as u16 as i16;
            }
            // Coefficients 3..7: combine remaining bits of lpp0 with lpp1.
            // lpp0 >> (3*low) gives the leftover bits; shift lpp1 up to join.
            let combined: u64 = (lpp1 << (32 - 3 * low)) | (lpp0 >> (3 * low)) as u64;
            for i in 3..8usize {
                let lp = ((combined >> ((i - 3) as u32 * low)) as u32) & lmask;
                let sm = ((sbb >> i) & 1).wrapping_neg();
                let cur = d[u + i] as u16 as u32;
                d[u + i] = (cur ^ (sm & 0xFFFF) ^ lp) as u16 as i16;
            }
        }
    }

    Ok((voff, num_ignored))
}

/// Decode a HAWK-512 public key from bytes.
///
/// Returns `(q00, q01)` where `q00` has length HAWK_N (with
/// `q00[HAWK_N/2..]` zeroed — the auto-adjoint half is not encoded) and
/// `q01` has length HAWK_N.
///
/// Port of `decode_q00 + decode_q01` (hawk_vrfy.c:1345-1423).
pub fn decode_public(buf: &[u8; HAWK_PUBLIC_KEY_BYTES]) -> Result<(Vec<i16>, Vec<i16>), HawkError> {
    let logn = 9u32;
    let low00 = HAWK_512_LOWBITS_Q00;
    let low01 = HAWK_512_LOWBITS_Q01;
    let eb00_len = HAWK_512_EB00_LEN;
    // lim_bits are specific to each Golomb-Rice instantiation in the C macros:
    //   - q00 at logn=9: `decode_gr_5_9`  → lim_bits = 9
    //   - q01 at logn=9: `decode_gr_9_12` → lim_bits = 12
    // The `low + 4` relation coincidentally holds for q00 but not q01, so we
    // hardcode both directly.
    let lim_00: u32 = 9;
    let lim_01: u32 = 12;

    let mut q00 = vec![0i16; HAWK_N];
    let mut q01 = vec![0i16; HAWK_N];

    // Decode q00: half-size (n/2 = 256 coefficients) at logn-1=8.
    let (mut len00, ni) = decode_gr(logn - 1, &mut q00[..HAWK_N / 2], buf, low00, lim_00)
        .map_err(|()| HawkError::MalformedPublicKey("decode_gr q00 failed".into()))?;

    // Recover q00[0] extra bits.
    // The low eb00_len bits of the original q00[0] were packed into the
    // unused bits of the last byte (and possibly one overflow byte).
    let eb00: i32;
    let last: u32;
    if eb00_len <= ni {
        let top = (buf[len00 - 1] as u32) >> (8 - ni);
        eb00 = (top & ((1u32 << eb00_len) - 1)) as i32;
        last = top >> eb00_len;
    } else {
        if len00 >= HAWK_PUBLIC_KEY_BYTES {
            return Err(HawkError::MalformedPublicKey(
                "q00 extra bits buffer overflow".into(),
            ));
        }
        let lo_bits = (buf[len00 - 1] as u32) >> (8 - ni);
        let next_byte = buf[len00] as u32;
        eb00 = (lo_bits | (next_byte << ni)) as i32 & ((1i32 << eb00_len) - 1);
        last = next_byte >> (eb00_len - ni);
        len00 += 1;
    }
    if last != 0 {
        return Err(HawkError::MalformedPublicKey(
            "q00 extra bits nonzero padding".into(),
        ));
    }
    // Reconstruct q00[0]: the encoded value was (q00[0] >> eb00_len); combine.
    q00[0] = (((q00[0] as i32) << eb00_len) | eb00) as i16;

    // Decode q01: full-size (n=512 coefficients) at logn=9.
    let pub_remaining = &buf[len00..];
    let (len01, ni01) = decode_gr(logn, &mut q01, pub_remaining, low01, lim_01)
        .map_err(|()| HawkError::MalformedPublicKey("decode_gr q01 failed".into()))?;

    // Check that unused bits in the last q01 byte are zero.
    if len01 > 0 {
        let last_byte_top = (pub_remaining[len01 - 1] as u32) >> (8 - ni01);
        if last_byte_top != 0 {
            return Err(HawkError::MalformedPublicKey(
                "q01 trailing bits nonzero".into(),
            ));
        }
    }

    // Check zero padding.
    for &b in &pub_remaining[len01..] {
        if b != 0 {
            return Err(HawkError::MalformedPublicKey("padding nonzero".into()));
        }
    }

    Ok((q00, q01))
}

/// Encode a HAWK-512 secret key.
///
/// Layout: kgseed(24) || F_mod2(64) || G_mod2(64) || hpub(32) = 184 bytes.
/// `hpub` is `SHAKE256(encoded_public_key)[..32]` and must be pre-computed
/// by the caller.
///
/// Port of `encode_private` (hawk_kgen.c:26-62).
pub fn encode_private(
    seed: &[u8],
    f_cap: &[i8],
    g_cap: &[i8],
    hpub: &[u8],
) -> [u8; HAWK_SECRET_KEY_BYTES] {
    assert_eq!(seed.len(), HAWK_512_KGSEED_LEN);
    assert_eq!(f_cap.len(), HAWK_N);
    assert_eq!(g_cap.len(), HAWK_N);
    assert_eq!(hpub.len(), HAWK_512_HPUB_LEN);

    let mut out = [0u8; HAWK_SECRET_KEY_BYTES];
    let mut off = 0usize;

    out[off..off + HAWK_512_KGSEED_LEN].copy_from_slice(seed);
    off += HAWK_512_KGSEED_LEN;

    extract_lowbit(f_cap, &mut out[off..off + HAWK_N / 8]);
    off += HAWK_N / 8;

    extract_lowbit(g_cap, &mut out[off..off + HAWK_N / 8]);
    off += HAWK_N / 8;

    out[off..off + HAWK_512_HPUB_LEN].copy_from_slice(hpub);

    out
}

/// Decode a HAWK-512 secret key.
///
/// Returns `(seed, F_mod2, G_mod2, hpub)` as byte vectors.
/// The caller reconstructs f, g via `Hawk_regen_fg(seed)` and lifts
/// F, G from the mod-2 bitmasks.
pub fn decode_private(buf: &[u8; HAWK_SECRET_KEY_BYTES]) -> (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>) {
    let mut off = 0usize;
    let seed = buf[off..off + HAWK_512_KGSEED_LEN].to_vec();
    off += HAWK_512_KGSEED_LEN;
    let f_mod2 = buf[off..off + HAWK_N / 8].to_vec();
    off += HAWK_N / 8;
    let g_mod2 = buf[off..off + HAWK_N / 8].to_vec();
    off += HAWK_N / 8;
    let hpub = buf[off..off + HAWK_512_HPUB_LEN].to_vec();
    (seed, f_mod2, g_mod2, hpub)
}

// ---------------------------------------------------------------------------
// Signature encode / decode
// ---------------------------------------------------------------------------

/// HAWK-512 signature s1 Golomb-Rice low-bits parameter.
/// Source: `hawk_sign.c:684` — `low = (logn == 10) ? 6 : 5`.
const HAWK_512_LOWBITS_S1: u32 = 5;
/// lim_bits for decode_gr_5_9.
const HAWK_512_LIMBITS_S1: u32 = 9;

/// Encode a HAWK-512 signature: salt || GR(s1) || zero-padding, 555 bytes.
///
/// Port of `encode_sig` (hawk_sign.c:680-755) for HAWK-512 (logn=9).
pub fn encode_signature(
    salt: &[u8; crate::params::HAWK_SALT_BYTES],
    s1: &[i16],
) -> Result<[u8; crate::params::HAWK_SIGNATURE_BYTES], HawkError> {
    assert_eq!(s1.len(), HAWK_N);
    let logn = 9u32;
    let low = HAWK_512_LOWBITS_S1;

    let mut out = [0u8; crate::params::HAWK_SIGNATURE_BYTES];
    let salt_len = crate::params::HAWK_SALT_BYTES;

    // Salt.
    out[..salt_len].copy_from_slice(salt);

    // Golomb-Rice of s1. encode_gr writes (signs || fixed-parts || variable)
    // and returns (bytes_written, num_ignored_bits).
    let (enc_len, _ni) = encode_gr(logn, &mut out[salt_len..], s1, low)
        .map_err(|()| HawkError::MalformedSignature("encode_gr s1 overflow".into()))?;
    let written = salt_len + enc_len;
    // Zero-pad remainder.
    for b in out[written..].iter_mut() {
        *b = 0;
    }
    Ok(out)
}

/// Decode a HAWK-512 signature: extract salt and s1 from 555 bytes.
/// Rejects malformed encodings (trailing bits nonzero, padding nonzero).
///
/// Port of `decode_sig_inner` (hawk_vrfy.c:1475-1491) for HAWK-512.
pub fn decode_signature(
    bytes: &[u8; crate::params::HAWK_SIGNATURE_BYTES],
) -> Result<crate::sign::HawkSignature, HawkError> {
    let logn = 9u32;
    let low = HAWK_512_LOWBITS_S1;
    let lim = HAWK_512_LIMBITS_S1;
    let salt_len = crate::params::HAWK_SALT_BYTES;

    let mut salt = [0u8; crate::params::HAWK_SALT_BYTES];
    salt.copy_from_slice(&bytes[..salt_len]);

    let mut s1 = vec![0i16; HAWK_N];
    let buf = &bytes[salt_len..];
    let (enc_len, ni) = decode_gr(logn, &mut s1, buf, low, lim)
        .map_err(|()| HawkError::MalformedSignature("decode_gr s1 failed".into()))?;
    // Check that the unused bits in the last used byte are zero.
    if enc_len > 0 {
        let last_byte = (buf[enc_len - 1] as u32) >> (8 - ni);
        if last_byte != 0 {
            return Err(HawkError::MalformedSignature(
                "trailing bits nonzero".into(),
            ));
        }
    }
    // Check padding beyond encoded content is zero.
    for &b in buf[enc_len..].iter() {
        if b != 0 {
            return Err(HawkError::MalformedSignature("padding nonzero".into()));
        }
    }

    Ok(crate::sign::HawkSignature { salt, s1 })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extract_lowbit_basic() {
        let f: [i8; 16] = [1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1];
        let mut out = [0u8; 2];
        extract_lowbit(&f, &mut out);
        // byte 0: bits f[0]..f[7] LSB-first = 1,0,1,1,0,0,1,0 = 0b01001101
        assert_eq!(out[0], 0b01001101);
        // byte 1: bits f[8]..f[15] LSB-first = 1,1,1,0,0,1,0,1 = 0b10100111
        assert_eq!(out[1], 0b10100111);
    }

    #[test]
    fn encode_private_then_decode_roundtrips() {
        let seed = vec![7u8; 24];
        let f_cap: Vec<i8> = (0i16..512).map(|i| i as i8).collect();
        let g_cap: Vec<i8> = (0i16..512).map(|i| (i as i8).wrapping_mul(3)).collect();
        let hpub = vec![0xABu8; 32];

        let enc = encode_private(&seed, &f_cap, &g_cap, &hpub);
        let (s2, f_mod2, g_mod2, h2) = decode_private(&enc);

        assert_eq!(s2, seed);
        assert_eq!(h2, hpub);
        // Each bit of f_mod2 must equal the LSB of the corresponding f_cap element.
        for i in 0..512usize {
            let bit = (f_mod2[i / 8] >> (i % 8)) & 1;
            assert_eq!(bit as i8, f_cap[i] & 1, "f_mod2 bit {} mismatch", i);
            let bit = (g_mod2[i / 8] >> (i % 8)) & 1;
            assert_eq!(bit as i8, g_cap[i] & 1, "g_mod2 bit {} mismatch", i);
        }
    }
}
