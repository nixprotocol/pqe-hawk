//! Regenerate the (f, g) polynomials of a HAWK secret key from a seed.
//!
//! Port of c-reference/hawk-512/ng_hawk.c:95-182. Only the HAWK-512 (logn=9)
//! variant is ported; logn=8 and logn=10 are unused at this parameter set.

use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;

/// Seed length in bytes for HAWK-512 regen_fg. Port of `seed_len = 24`
/// from ng_hawk.c:98.
pub const HAWK_REGEN_FG_SEED_LEN: usize = 24;

/// Regenerate f and g polynomials from a 24-byte seed (HAWK-512).
/// Both f and g must be preallocated as `[i8; 512]`.
///
/// Port of `regen_fg_9` (ng_hawk.c:95-128).
pub fn hawk_regen_fg(f: &mut [i8], g: &mut [i8], seed: &[u8]) {
    assert_eq!(f.len(), 512, "f must be 512 i8 for HAWK-512");
    assert_eq!(g.len(), 512, "g must be 512 i8 for HAWK-512");
    assert_eq!(
        seed.len(),
        HAWK_REGEN_FG_SEED_LEN,
        "seed must be {} bytes",
        HAWK_REGEN_FG_SEED_LEN
    );

    for j in 0..4u8 {
        // Four parallel SHAKE instances, each initialized with seed || j_byte.
        let mut shake = Shake256::default();
        shake.update(seed);
        shake.update(&[j]);
        let mut reader = shake.finalize_xof();

        // For each of 32 chunks (1024 / 32 = 32 chunks per instance):
        //   extract 8 bytes
        //   popcount each byte (via the 3-step mask-add trick)
        //   store 8 signed coefficients into f or g with stride 8 and offset j*8
        //
        // u ranges 0, 32, 64, ..., 992 — 32 iterations.
        let mut u = 0usize;
        while u < 1024 {
            let mut qb = [0u8; 8];
            reader.read(&mut qb);
            let mut q = u64::from_le_bytes(qb);

            // Popcount each byte: reduce bit-pairs, then nibbles, then byte-wise.
            q = (q & 0x5555_5555_5555_5555) + ((q >> 1) & 0x5555_5555_5555_5555);
            q = (q & 0x3333_3333_3333_3333) + ((q >> 2) & 0x3333_3333_3333_3333);
            q = (q & 0x0F0F_0F0F_0F0F_0F0F) + ((q >> 4) & 0x0F0F_0F0F_0F0F_0F0F);

            // q now holds 8 byte-popcounts in its 8 bytes. Each popcount is in
            // [0, 8] (since a byte has 8 bits). Subtract 4 to center, giving
            // coefficient values in [-4, +4].
            let mut vv = [0i8; 8];
            let mut qq = q;
            for i in 0..8 {
                vv[i] = ((qq & 0xFF) as i32 - 4) as i8;
                qq >>= 8;
            }

            // Target buffer and offset.
            // Source offset (within the "f" half or "g" half): u_local = u if u<512 else u-512.
            // Target stride: 8 per j; within a chunk, 8 consecutive bytes.
            let offset_within = if u < 512 { u } else { u - 512 };
            let base_offset = offset_within + (j as usize) * 8;
            if u < 512 {
                f[base_offset..base_offset + 8].copy_from_slice(&vv);
            } else {
                g[base_offset..base_offset + 8].copy_from_slice(&vv);
            }

            u += 32;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regen_fg_produces_coefficients_in_range() {
        let seed = [0u8; 24];
        let mut f = [0i8; 512];
        let mut g = [0i8; 512];
        hawk_regen_fg(&mut f, &mut g, &seed);
        // All coefficients should be in [-4, +4].
        for (i, &c) in f.iter().enumerate() {
            assert!(c >= -4 && c <= 4, "f[{}] = {} out of [-4, 4]", i, c);
        }
        for (i, &c) in g.iter().enumerate() {
            assert!(c >= -4 && c <= 4, "g[{}] = {} out of [-4, 4]", i, c);
        }
    }

    #[test]
    fn regen_fg_is_deterministic() {
        let seed = [42u8; 24];
        let mut f1 = [0i8; 512];
        let mut g1 = [0i8; 512];
        let mut f2 = [0i8; 512];
        let mut g2 = [0i8; 512];
        hawk_regen_fg(&mut f1, &mut g1, &seed);
        hawk_regen_fg(&mut f2, &mut g2, &seed);
        assert_eq!(f1, f2);
        assert_eq!(g1, g2);
    }

    #[test]
    fn different_seeds_yield_different_outputs() {
        let seed_a = [1u8; 24];
        let seed_b = [2u8; 24];
        let mut fa = [0i8; 512];
        let mut ga = [0i8; 512];
        let mut fb = [0i8; 512];
        let mut gb = [0i8; 512];
        hawk_regen_fg(&mut fa, &mut ga, &seed_a);
        hawk_regen_fg(&mut fb, &mut gb, &seed_b);
        assert_ne!(fa, fb);
        assert_ne!(ga, gb);
    }
}
