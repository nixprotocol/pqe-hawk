//! Constants for HAWK-512 verification.
//!
//! Port of c-reference/hawk-512/hawk_vrfy.c:257-266, 1194-1202.

/// First prime modulus used by verify's dual-prime NTT: 2^31 - 10239.
pub const P1: u32 = 2147473409;
/// -1/P1 mod 2^32.
pub const P1_0I: u32 = 2042615807;
/// 2^64 mod P1.
pub const P1_R2: u32 = 419348484;
/// 2^96 mod P1.
pub const P1_R3: u32 = 1819566170;
/// 16*2^32 mod P1.
pub const P1_M16: u32 = 327648;

/// Second prime modulus used by verify's dual-prime NTT: 2^31 - 94207.
pub const P2: u32 = 2147389441;
/// -1/P2 mod 2^32.
pub const P2_0I: u32 = 1862176767;
/// 2^64 mod P2.
pub const P2_R2: u32 = 1141604340;
/// 2^96 mod P2.
pub const P2_R3: u32 = 976758995;
/// 16*2^32 mod P2.
pub const P2_M16: u32 = 3014624;

/// `bits_lim01[logn]` — max bit size of q01 coefficients (excluding sign bit).
pub const BITS_LIM01: [i8; 11] = [0, 0, 0, 0, 0, 0, 0, 0, 11, 12, 14];
/// `bits_lims0[logn]` — limit for s0 reconstruction.
pub const BITS_LIMS0: [i8; 11] = [0, 0, 0, 0, 0, 0, 0, 0, 12, 13, 14];
/// `bits_lims1[logn]` — limit for s1 (same as q00 encoding lim + 4).
pub const BITS_LIMS1: [i8; 11] = [0, 0, 0, 0, 0, 0, 0, 0, 9, 9, 10];

/// Maximum tnorm for HAWK-512 verify. (sigma_ver^2)*8*n^2.
/// Source: hawk_vrfy.c:1618 — `max_tnorm = 8317` for logn=9.
pub const HAWK_512_MAX_TNORM: u32 = 8317;

/// KeyNormCheck floor for `q00[0]` (Dao, eprint 2026/1298).
///
/// For an auto-adjoint public key, `q00[0] = ||(f,g)||^2`. Honest keygen
/// resamples until `||(f,g)||^2 >= 2080` (see `keygen::hawk_keygen`,
/// `HAWK_512_L2LOW`, from reference `ng_hawk.c:462`), so every legitimately
/// generated HAWK-512 key has `q00[0] >= 2080`. The reference verifier — and
/// this port — never re-checked this floor, so a maliciously-formed public key
/// with a tiny `q00[0]` (e.g. the constant `q00=1` weak keys) is accepted and
/// breaks the BUFF properties (message-bound signatures, exclusive ownership).
/// Rejecting `q00[0] < 2080` at verify time closes that gap and is a no-op for
/// honest keys. This is the correct HAWK-512 threshold; do not use the
/// HAWK-256 value (556) — it would leave a `[556, 2080)` gap.
///
/// The floor is exactly tight, not conservative: measured over 48 honest
/// HAWK-512 keys, `q00[0]` ranged `[2080, 2290]`, hitting the floor exactly.
/// So it costs zero honest-key margin, and lowering it to 556 would open a
/// band that no honest key ever occupies. Regression test: `tests/buff_poc.rs`.
pub const HAWK_512_Q00_FLOOR: i32 = 2080;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primes_have_correct_relations() {
        // P0 (from PRIMES[0]) should equal P1 here (they're the same prime).
        assert_eq!(P1, 2147473409);
        assert!(P1 < (1u32 << 31));
        assert!(P2 < (1u32 << 31));
        assert!(P1 > P2);
        // P1 - 1 multiple of 2048.
        assert_eq!((P1 - 1) % 2048, 0);
        assert_eq!((P2 - 1) % 2048, 0);
    }
}
