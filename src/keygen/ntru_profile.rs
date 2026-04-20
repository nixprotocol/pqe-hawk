//! NTRU solver profile: per-degree tuning parameters controlling the
//! multi-depth NTRU equation solver.
//!
//! Port of `ntru_profile` from ng_inner.h:1395-1404 and the per-degree
//! profile constants from ng_hawk.c:3-34.

/// Per-degree solver tuning parameters.
///
/// Port of `ntru_profile` from c-reference/hawk-512/ng_inner.h:1395-1404.
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NtruProfile {
    pub q: u32,
    pub min_logn: u32,
    pub max_logn: u32,
    pub max_bl_small: [u16; 11],
    pub max_bl_large: [u16; 10],
    pub word_win: [u16; 10],
    pub reduce_bits: u32,
    pub coeff_fg_limit: [u8; 11],
    pub min_save_fg: [u16; 11],
}

/// Solver result: no error.
pub const SOLVE_OK: i32 = 0;
/// Solver result: gcd(Res(f,X^n+1), Res(g,X^n+1)) ≠ 1.
pub const SOLVE_ERR_GCD: i32 = -1;
/// Solver result: reduction failed (NTRU equation no longer satisfied).
pub const SOLVE_ERR_REDUCE: i32 = -2;
/// Solver result: (F, G) coefficients outside the allowed range.
pub const SOLVE_ERR_LIMIT: i32 = -3;

/// HAWK-256 solver profile. Port of `SOLVE_Hawk_256` from ng_hawk.c:3-12.
pub const SOLVE_HAWK_256: NtruProfile = NtruProfile {
    q: 1,
    min_logn: 8,
    max_logn: 8,
    max_bl_small: [1, 1, 1, 2, 3, 5, 9, 17, 34, 0, 0],
    max_bl_large: [1, 1, 2, 4, 7, 13, 26, 50, 0, 0],
    word_win: [1, 1, 1, 2, 3, 3, 3, 4, 0, 0],
    reduce_bits: 14,
    coeff_fg_limit: [0, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127],
    min_save_fg: [0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3],
};

/// HAWK-512 solver profile. Port of `SOLVE_Hawk_512` from ng_hawk.c:14-23.
pub const SOLVE_HAWK_512: NtruProfile = NtruProfile {
    q: 1,
    min_logn: 9,
    max_logn: 9,
    max_bl_small: [1, 1, 1, 2, 3, 6, 11, 21, 41, 82, 0],
    max_bl_large: [1, 2, 3, 5, 8, 16, 31, 61, 121, 0],
    word_win: [1, 1, 1, 2, 2, 3, 3, 4, 6, 0],
    reduce_bits: 11,
    coeff_fg_limit: [0, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127],
    min_save_fg: [0, 0, 1, 2, 2, 2, 2, 2, 2, 2, 3],
};

/// HAWK-1024 solver profile. Port of `SOLVE_Hawk_1024` from ng_hawk.c:25-34.
pub const SOLVE_HAWK_1024: NtruProfile = NtruProfile {
    q: 1,
    min_logn: 10,
    max_logn: 10,
    max_bl_small: [1, 1, 2, 2, 4, 7, 13, 25, 48, 96, 191],
    max_bl_large: [1, 2, 3, 5, 10, 19, 37, 72, 143, 284],
    word_win: [1, 1, 2, 2, 3, 3, 3, 4, 4, 7],
    reduce_bits: 9,
    coeff_fg_limit: [0, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127],
    min_save_fg: [0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3],
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hawk_512_profile_fields() {
        let p = SOLVE_HAWK_512;
        assert_eq!(p.q, 1);
        assert_eq!(p.min_logn, 9);
        assert_eq!(p.max_logn, 9);
        assert_eq!(p.reduce_bits, 11);
        assert_eq!(p.max_bl_small[8], 41);
        assert_eq!(p.max_bl_large[8], 121);
        assert_eq!(p.coeff_fg_limit[9], 127);
    }

    #[test]
    fn profile_sizes_are_fixed() {
        assert_eq!(SOLVE_HAWK_512.max_bl_small.len(), 11);
        assert_eq!(SOLVE_HAWK_512.max_bl_large.len(), 10);
        assert_eq!(SOLVE_HAWK_512.word_win.len(), 10);
        assert_eq!(SOLVE_HAWK_512.coeff_fg_limit.len(), 11);
        assert_eq!(SOLVE_HAWK_512.min_save_fg.len(), 11);
    }

    #[test]
    fn error_codes() {
        assert_eq!(SOLVE_OK, 0);
        assert_eq!(SOLVE_ERR_GCD, -1);
        assert_eq!(SOLVE_ERR_REDUCE, -2);
        assert_eq!(SOLVE_ERR_LIMIT, -3);
    }
}
