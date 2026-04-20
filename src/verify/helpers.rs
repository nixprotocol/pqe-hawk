//! Small verify-time helpers.

/// In-place compute `t1 = h1 - 2*s1` where `d` starts as `s1` and `h1` is
/// packed bits (n bits in n/8 bytes). After this call `d` holds the i16
/// coefficients of `t1`.
///
/// Port of `make_t1` (hawk_vrfy.c:1565-1576).
pub fn make_t1(logn: u32, d: &mut [i16], h1: &[u8]) {
    let n = 1usize << logn;
    debug_assert!(d.len() >= n);
    debug_assert!(h1.len() >= n / 8);
    for u in (0..n).step_by(8) {
        let mut h1b = h1[u >> 3] as u32;
        for v in 0..8 {
            let x = d[u + v] as i32 as u32;
            let new_x = (h1b & 1).wrapping_sub(x << 1);
            d[u + v] = new_x as i32 as i16;
            h1b >>= 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn make_t1_zero_s1_zero_h1() {
        let mut d = vec![0i16; 512];
        let h1 = vec![0u8; 64];
        make_t1(9, &mut d, &h1);
        assert_eq!(d, vec![0i16; 512]);
    }

    #[test]
    fn make_t1_simple_case() {
        // s1 = [1, 2, 3, ...], h1 bit = 1 everywhere.
        // After: d[u] = 1 - 2*s1[u] = 1 - 2u.
        let mut d: Vec<i16> = (0..512).map(|i| i as i16).collect();
        let h1 = vec![0xFFu8; 64]; // all bits set
        make_t1(9, &mut d, &h1);
        for u in 0..512 {
            let expected = 1i32 - 2 * (u as i32);
            assert_eq!(d[u] as i32, expected as i16 as i32, "u={}", u);
        }
    }
}
