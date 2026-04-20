//! NTRU equation solver for HAWK key generation.
//!
//! Port of c-reference/hawk-512/ng_ntru.c. The solver has a multi-depth
//! recursive structure:
//!   - `make_fg_step` reduces (f,g) from depth d to depth d+1 in RNS form.
//!   - `make_fg_deepest` runs all reductions down to degree 1 (resultants).
//!   - `solve_NTRU_deepest` applies Bezout to the resultants to get the
//!     base-case (F,G).
//! Higher-depth lift operations come in later sub-tasks.

#[allow(unused_imports)]
use crate::keygen::mp31::{mp_add, mp_div, mp_norm, mp_set, mp_sub};
use crate::keygen::mp31::{mp_intt, mp_mkgm, mp_mkigm, mp_montymul, mp_ntt, mp_rx31};
use crate::keygen::ntru_profile::{NtruProfile, SOLVE_ERR_GCD, SOLVE_ERR_REDUCE, SOLVE_OK};
use crate::keygen::poly::poly_mp_set_small;
use crate::keygen::primes::PRIMES;
use crate::keygen::zint31::{zint_bezout, zint_mod_small_signed, zint_mul_small, zint_rebuild_crt};

/// Offset in tmp[] for saving intermediate (f,g) values — expressed in u32
/// words, keyed by `logn_top`. Port of `FG_SAVE_OFFSET` (ng_ntru.c:220).
pub(crate) fn fg_save_offset(logn_top: u32) -> usize {
    5usize << logn_top
}

/// Convert (f,g) to RNS+NTT form at the top depth, using PRIMES[0].
/// Output written at the start of tmp as [ft|gt|gm] (each n u32 limbs).
///
/// Port of `make_fg_zero` (ng_ntru.c:9-25).
pub(crate) fn make_fg_zero(logn: u32, f: &[i8], g: &[i8], tmp: &mut [u32]) {
    let n = 1usize << logn;
    let p = PRIMES[0].p;
    let p0i = PRIMES[0].p0i;
    // Split tmp into ft | gt | gm regions.
    let (ft, rest) = tmp.split_at_mut(n);
    let (gt, rest) = rest.split_at_mut(n);
    let gm = &mut rest[..n];
    poly_mp_set_small(logn, ft, f, p);
    poly_mp_set_small(logn, gt, g, p);
    mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);
    mp_ntt(logn, ft, gm, p, p0i);
    mp_ntt(logn, gt, gm, p, p0i);
}

/// One depth step of (f,g) reduction. See ng_ntru.c:36-130.
///
/// Input: tmp[0..2*n*slen] holds (f,g) at depth `depth` in RNS+NTT form.
/// Output: tmp[0..2*hn*tlen] holds (f',g') at depth `depth+1` in RNS+NTT form.
///
/// where:
///   n    = 2^(logn_top - depth)
///   hn   = n / 2
///   slen = prof.max_bl_small[depth]
///   tlen = prof.max_bl_small[depth + 1]
///
/// Port of `make_fg_step` (ng_ntru.c:36-130).
pub(crate) fn make_fg_step(prof: &NtruProfile, logn_top: u32, depth: u32, tmp: &mut [u32]) {
    let logn = logn_top - depth;
    let n = 1usize << logn;
    let hn = n >> 1;
    let slen = prof.max_bl_small[depth as usize] as usize;
    let tlen = prof.max_bl_small[depth as usize + 1] as usize;

    // Layout after memmove:
    //   fd[0..hn*tlen]   output f'
    //   gd[0..hn*tlen]   output g'
    //   fs[0..n*slen]    source f (moved here from tmp[0])
    //   gs[0..n*slen]    source g (moved here from tmp[n*slen])
    //   t1[0..n]         NTT support
    //   t2[0..]          extra
    let fd_off = 0usize;
    let gd_off = fd_off + hn * tlen;
    let fs_off = gd_off + hn * tlen;
    let gs_off = fs_off + n * slen;
    let t1_off = gs_off + n * slen;
    let t2_off = t1_off + n;

    // memmove(fs, tmp, 2 * n * slen): copy the existing (f,g) to the fs/gs region.
    tmp.copy_within(0..2 * n * slen, fs_off);

    // === First slen primes: process directly from source (RNS+NTT → plain RNS). ===
    // For each prime u in 0..slen:
    //   - Compute yf[v] = montymul(montymul(xf[2v], xf[2v+1]), R2) for v in 0..hn.
    //   - Inverse-NTT xf and xg to get them in plain RNS.
    //   - Advance pointers: xf += n, xg += n, yf += hn, yg += hn.
    //
    // We work via index arithmetic, since split_at_mut chains would be
    // unwieldy here due to the interleaved reads and writes.
    for u in 0..slen {
        let p = PRIMES[u].p;
        let p0i = PRIMES[u].p0i;
        let r2 = PRIMES[u].r2;

        let xf_base = fs_off + u * n; // tmp[xf_base..xf_base+n] = row u of fs
        let xg_base = gs_off + u * n; // tmp[xg_base..xg_base+n] = row u of gs
        let yf_base = fd_off + u * hn; // tmp[yf_base..yf_base+hn] = row u of fd
        let yg_base = gd_off + u * hn; // tmp[yg_base..yg_base+hn] = row u of gd

        for v in 0..hn {
            let a = tmp[xf_base + 2 * v];
            let b = tmp[xf_base + 2 * v + 1];
            tmp[yf_base + v] = mp_montymul(mp_montymul(a, b, p, p0i), r2, p, p0i);
            let a = tmp[xg_base + 2 * v];
            let b = tmp[xg_base + 2 * v + 1];
            tmp[yg_base + v] = mp_montymul(mp_montymul(a, b, p, p0i), r2, p, p0i);
        }

        // Build iGM table in t1, then iNTT the source rows in place.
        {
            let (left, right) = tmp.split_at_mut(t1_off);
            let t1 = &mut right[..n];
            mp_mkigm(logn, t1, PRIMES[u].ig, p, p0i);
            mp_intt(logn, &mut left[xf_base..xf_base + n], t1, p, p0i);
            mp_intt(logn, &mut left[xg_base..xg_base + n], t1, p, p0i);
        }
    }

    // === Rebuild fs and gs from RNS into plain integer form. ===
    // zint_rebuild_crt(fs, slen, n, 2, true, t1)
    // fs spans tmp[fs_off..fs_off + 2*n*slen] with num_sets=2, n cols, slen rows.
    // t1 starts right after (t1_off), length at least slen.
    {
        let (_, right) = tmp.split_at_mut(fs_off);
        let (fs_gs, rest) = right.split_at_mut(2 * n * slen);
        // rest starts at t1_off (relative to fs_off) — pass it in full so
        // zint_rebuild_crt can index up to rest[slen-1].
        zint_rebuild_crt(fs_gs, slen, n, 2, true, rest);
    }

    // === Remaining primes slen..tlen: compute new NTT rows from plain integers. ===
    for u in slen..tlen {
        let p = PRIMES[u].p;
        let p0i = PRIMES[u].p0i;
        let r2 = PRIMES[u].r2;
        let rx = mp_rx31(slen as u32, p, p0i, r2);

        let yf_base = fd_off + u * hn;
        let yg_base = gd_off + u * hn;

        // Build GM in t1.
        {
            let (_, right) = tmp.split_at_mut(t1_off);
            let t1 = &mut right[..n];
            mp_mkgm(logn, t1, PRIMES[u].g, p, p0i);
        }

        // t2[v] = zint_mod_small_signed(fs + v, slen, n, p, p0i, r2, rx) for v in 0..n.
        // fs is at tmp[fs_off..], stride n.
        for v in 0..n {
            let val = zint_mod_small_signed(&tmp[fs_off + v..], slen, n, p, p0i, r2, rx);
            tmp[t2_off + v] = val;
        }

        // NTT t2 using gm from t1.
        // Split at t2_off so left = tmp[0..t2_off] (contains t2) and
        // right = tmp[t2_off..] (t2 is at left[t2_off..t2_off+n] which is
        // left[t2_off - 0 .. ] relative to left).
        // Actually: split at t1_off, then sub-split right into t1 | t2.
        {
            let (left_t1, right_t1) = tmp.split_at_mut(t1_off);
            let (t1, t2_slice) = right_t1.split_at_mut(n);
            let _ = left_t1;
            mp_ntt(logn, &mut t2_slice[..n], t1, p, p0i);
        }

        // yf[v] = montymul(montymul(t2[2v], t2[2v+1]), r2).
        for v in 0..hn {
            let a = tmp[t2_off + 2 * v];
            let b = tmp[t2_off + 2 * v + 1];
            tmp[yf_base + v] = mp_montymul(mp_montymul(a, b, p, p0i), r2, p, p0i);
        }

        // t2[v] = zint_mod_small_signed(gs + v, slen, n, ...) for v in 0..n.
        for v in 0..n {
            let val = zint_mod_small_signed(&tmp[gs_off + v..], slen, n, p, p0i, r2, rx);
            tmp[t2_off + v] = val;
        }

        // NTT t2 using gm from t1.
        {
            let (left_t1, right_t1) = tmp.split_at_mut(t1_off);
            let (t1, t2_slice) = right_t1.split_at_mut(n);
            let _ = left_t1;
            mp_ntt(logn, &mut t2_slice[..n], t1, p, p0i);
        }

        // yg[v] = montymul(montymul(t2[2v], t2[2v+1]), r2).
        for v in 0..hn {
            let a = tmp[t2_off + 2 * v];
            let b = tmp[t2_off + 2 * v + 1];
            tmp[yg_base + v] = mp_montymul(mp_montymul(a, b, p, p0i), r2, p, p0i);
        }
    }
}

/// Compute (f,g) at the deepest level (Res(f,X^n+1) and Res(g,X^n+1)).
/// Saves intermediate (f,g) at depths >= prof.min_save_fg[logn_top] to the
/// tail of tmp (at positions counted back from sav_off) for later lift use.
///
/// Returns `true` if f is invertible mod PRIMES[0].p (precondition for
/// the overall keygen to succeed).
///
/// Port of `make_fg_deepest` (ng_ntru.c:163-201).
pub(crate) fn make_fg_deepest(
    prof: &NtruProfile,
    logn_top: u32,
    f: &[i8],
    g: &[i8],
    tmp: &mut [u32],
    mut sav_off: usize,
) -> bool {
    make_fg_zero(logn_top, f, g, tmp);

    // Check invertibility of f: all NTT coefficients of f (stored at
    // tmp[0..n]) must be non-zero.
    let n = 1usize << logn_top;
    let mut b: u32 = 0;
    for u in 0..n {
        b |= tmp[u].wrapping_sub(1);
    }
    let r = (b >> 31) == 0; // r = 1 - (b >> 31) in C, true if invertible

    for d in 0..logn_top {
        make_fg_step(prof, logn_top, d, tmp);

        // make_fg_step computed (f,g) for depth d+1. Save it if:
        //   d+1 < logn_top  (not the deepest level itself)
        //   d+1 >= prof.min_save_fg[logn_top]
        let d2 = d + 1;
        if d2 < logn_top && d2 >= prof.min_save_fg[logn_top as usize] as u32 {
            let slen = prof.max_bl_small[d2 as usize] as usize;
            let logn_d2 = logn_top - d2;
            let nd2 = 1usize << logn_d2;
            let fglen = slen * nd2 * 2; // 2 * n_d2 * slen words (f and g)
            sav_off -= fglen;
            // memmove(tmp + sav_off, tmp, fglen)
            tmp.copy_within(0..fglen, sav_off);
        }
    }

    r
}

/// Solve the NTRU equation at the deepest level.
/// On success, writes (F, G) to the start of tmp (len u32 words each) and
/// returns `SOLVE_OK`. Returns negative error code on failure.
///
/// Port of `solve_NTRU_deepest` (ng_ntru.c:231-286).
pub fn solve_ntru_deepest(
    prof: &NtruProfile,
    logn_top: u32,
    f: &[i8],
    g: &[i8],
    tmp: &mut [u32],
) -> i32 {
    let sav_off = fg_save_offset(logn_top);
    if !make_fg_deepest(prof, logn_top, f, g, tmp, sav_off) {
        return SOLVE_ERR_GCD;
    }

    let len = prof.max_bl_small[logn_top as usize] as usize;

    // Memory layout:
    //   Fp[0..len]     output F
    //   Gp[0..len]     output G
    //   fp[0..len]     Res(f,X^n+1) — moved here from tmp[0..len]
    //   gp[0..len]     Res(g,X^n+1) — moved here from tmp[len..2*len]
    //   t1[...]        scratch
    //
    // memmove(fp, tmp, 2*len): copy current tmp[0..2*len] to tmp[2*len..4*len].
    tmp.copy_within(0..2 * len, 2 * len);

    // Rebuild the resultants from RNS into plain integers.
    // fp and gp are adjacent big-ints of len limbs each, num_sets=2, n=1.
    {
        let (_, rest) = tmp.split_at_mut(2 * len);
        let (fp_gp, t1) = rest.split_at_mut(2 * len);
        zint_rebuild_crt(fp_gp, len, 1, 2, false, t1);
    }

    // Apply binary GCD to get (F, G) with f*G - g*F = 1.
    // C layout: Fp = tmp[0..len], Gp = tmp[len..2*len].
    // zint_bezout(Gp, Fp, fp, gp, ...) writes Gp to first arg, Fp to second.
    // fp = tmp[2*len..3*len], gp = tmp[3*len..4*len], t1 = tmp[4*len..].
    {
        let (fp_gp_out, rest) = tmp.split_at_mut(2 * len);
        // Slice t1 to exactly 4*len so that zint_bezout's internal b slice
        // ends up exactly len elements (the debug_assert in zint_co_reduce
        // checks b.len() == len).
        let (fp_gp_in, t1_full) = rest.split_at_mut(2 * len);
        let t1 = &mut t1_full[..4 * len];
        let (fp_out, gp_out) = fp_gp_out.split_at_mut(len); // Fp=tmp[0..len], Gp=tmp[len..2*len]
        let (fp_in, gp_in) = fp_gp_in.split_at_mut(len);
        // zint_bezout(u=Gp, v=Fp, x=fp, y=gp, ...): fp*Gp - gp*Fp = 1
        if zint_bezout(gp_out, fp_out, fp_in, gp_in, len, t1) == 0 {
            return SOLVE_ERR_GCD;
        }
    }

    // If q != 1, multiply Fp by q.
    // Fp = tmp[0..len].
    if prof.q != 1 {
        let fp_out = &mut tmp[..len];
        if zint_mul_small(fp_out, prof.q) != 0 {
            return SOLVE_ERR_REDUCE;
        }
    }

    SOLVE_OK
}

/// Compute (f,g) at intermediate depth `depth` in RNS+NTT form, writing the
/// result at the start of `tmp`.
///
/// Port of `make_fg_intermediate` (ng_ntru.c:142-152).
pub(crate) fn make_fg_intermediate(
    prof: &NtruProfile,
    logn_top: u32,
    f: &[i8],
    g: &[i8],
    depth: u32,
    tmp: &mut [u32],
) {
    make_fg_zero(logn_top, f, g, tmp);
    for d in 0..depth {
        make_fg_step(prof, logn_top, d, tmp);
    }
}

/// Minimum logn at which the NTT variant of poly_sub_scaled is faster
/// than the schoolbook variant. See `MIN_LOGN_FGNTT` (ng_ntru.c:293).
const MIN_LOGN_FGNTT: u32 = 4;

/// Solving the NTRU equation at an intermediate recursion level.
/// Input: F from the deeper level (half-degree), plain representation,
/// at the start of tmp. Output: F at this level, written at the start of tmp.
///
/// Port of `solve_NTRU_intermediate` (ng_ntru.c:303-707).
///
/// Variables named `scale_FG`, `scale_FG_init` deliberately mirror the C's
/// identifiers to aid byte-diff review against hawk_sign.c; the
/// `non_snake_case` lint is silenced for that reason.
#[allow(non_snake_case)]
pub fn solve_ntru_intermediate(
    prof: &NtruProfile,
    logn_top: u32,
    f: &[i8],
    g: &[i8],
    depth: u32,
    tmp: &mut [u32],
) -> i32 {
    use crate::keygen::fxp::{vect_fft, vect_ifft, vect_inv_mul2e_fft, vect_mul_fft};
    use crate::keygen::fxr::Fxr;
    use crate::keygen::mp31::{mp_intt, mp_mkgm, mp_mkgmigm, mp_montymul, mp_ntt, mp_rx31, tbmask};
    use crate::keygen::poly::{
        divrem31, poly_big_to_fixed, poly_max_bitlength, poly_sub_kf_scaled_depth1,
        poly_sub_scaled, poly_sub_scaled_ntt,
    };
    use crate::keygen::zint31::{zint_mod_small_signed, zint_rebuild_crt};

    let logn = logn_top - depth;
    let n = 1usize << logn;
    let hn = n >> 1;

    // slen   size for (f,g) at this level (also size of output F)
    // llen   size for unreduced F at this level
    // dlen   size for input F from deeper level
    let slen = prof.max_bl_small[depth as usize] as usize;
    let llen = prof.max_bl_large[depth as usize] as usize;
    let dlen = prof.max_bl_small[(depth + 1) as usize] as usize;

    // =========================================================================
    // Step 1: Get (f,g) for this level in RNS+NTT form.
    //
    // Initial layout: [Fd(dlen*hn) | fgt_placeholder ...]
    // We write (f,g) RNS+NTT into tmp starting at fgt_start = dlen*hn.
    // =========================================================================
    let fgt_start = dlen * hn;

    if depth < prof.min_save_fg[logn_top as usize] as u32 {
        // Compute from scratch via make_fg_intermediate, writing to tmp[fgt_start..].
        let (_fd, fgt_buf) = tmp.split_at_mut(fgt_start);
        make_fg_intermediate(prof, logn_top, f, g, depth, fgt_buf);
    } else {
        // Restore from the saved region at the end of tmp.
        // The saves are packed from fg_save_offset() backwards, deepest first.
        let mut sav_off = fg_save_offset(logn_top);
        for d in prof.min_save_fg[logn_top as usize] as u32..=depth {
            sav_off -= (prof.max_bl_small[d as usize] as usize) << (logn_top + 1 - d);
        }
        // Move tmp[sav_off .. sav_off + 2*slen*n] to tmp[fgt_start ..].
        let copy_len = 2 * slen * n;
        tmp.copy_within(sav_off..sav_off + copy_len, fgt_start);
    }

    // =========================================================================
    // Step 2: Rearrange buffers to make room for the unreduced F.
    //
    // New layout:
    //   Ft  = tmp[0          .. llen*n]          (unreduced F placeholder)
    //   ft  = tmp[ft_start   .. ft_start+slen*n]
    //   gt  = tmp[gt_start   .. gt_start+slen*n]
    //   Fd  = tmp[fd_start   .. fd_start+dlen*hn]
    //   t1  = tmp[t1_start   ..]
    //
    // Old layout:
    //   tmp[0        .. fgt_start]        Fd  (dlen*hn)
    //   tmp[fgt_start.. fgt_start+2*n*slen] fgt = [ft|gt]
    //
    // Steps:
    //   1. Copy fgt (at tmp[fgt_start..]) to new ft (at tmp[llen*n..]).
    //      Since llen >= dlen, llen*n >= dlen*hn = fgt_start, so destination
    //      >= source: use copy_within (safe for forward overlap).
    //   2. Copy old Fd (at tmp[0..dlen*hn]) to new Fd (at tmp[llen*n + 2*slen*n..]).
    //      This must happen AFTER step 1 because the source could be overwritten.
    // =========================================================================
    let ft_start = llen * n;
    let gt_start = ft_start + slen * n;
    let fd_start = gt_start + slen * n;
    let t1_start = fd_start + dlen * hn;

    // Step 2a: move fgt -> ft (copy_within handles overlapping forward copy).
    tmp.copy_within(fgt_start..fgt_start + 2 * n * slen, ft_start);

    // Step 2b: move old Fd (tmp[0..dlen*hn]) -> new position at fd_start.
    // The destination (fd_start = llen*n + 2*slen*n) is always > dlen*hn
    // so no overlap with source when llen >= 1.
    tmp.copy_within(0..hn * dlen, fd_start);

    // =========================================================================
    // Step 3: Convert Fd to RNS, storing values in the *last hn* slots of each
    //         n-word row of Ft (i.e. Ft[u*n + hn .. u*n + n]).
    //
    // After this, Fd is no longer needed; t1 is moved to fd_start.
    // =========================================================================
    for u in 0..llen {
        let p = PRIMES[u].p;
        let p0i = PRIMES[u].p0i;
        let r2 = PRIMES[u].r2;
        let rx = mp_rx31(dlen as u32, p, p0i, r2);
        for v in 0..hn {
            // Read from Fd (now at fd_start), write to Ft[u*n + hn + v].
            // These regions are disjoint: u*n + hn + v < llen*n <= fd_start.
            let z = zint_mod_small_signed(&tmp[fd_start + v..], dlen, hn, p, p0i, r2, rx);
            tmp[u * n + hn + v] = z;
        }
    }
    // Fd no longer needed; t1 now starts at fd_start.
    let t1_start = fd_start;

    // =========================================================================
    // Step 4: Compute unreduced F modulo llen small primes.
    //         Also inverse-NTT (f,g) as we go; after slen primes, rebuild CRT.
    //
    // For each prime u in 0..llen:
    //   a. build gm/igm in t1
    //   b. get g mod p (either from gt[u*n] directly, or computed from plain gt)
    //   c. NTT the hn-entry slice of Ft[u*n+hn..]  (logn-1 NTT)
    //   d. multiply to get unreduced F row
    //   e. iNTT the F row (logn iNTT)
    // =========================================================================
    for u in 0..llen {
        if u == slen {
            // ft and gt are now in plain RNS; rebuild them via CRT.
            // ft and gt together occupy tmp[ft_start .. ft_start + 2*slen*n].
            // zint_rebuild_crt needs adjacent scratch starting at fd_start.
            let (head, tail) = tmp.split_at_mut(ft_start);
            let (ft_gt_region, rest) = tail.split_at_mut(2 * slen * n);
            zint_rebuild_crt(ft_gt_region, slen, n, 2, true, rest);
            let _ = head;
        }

        let p = PRIMES[u].p;
        let p0i = PRIMES[u].p0i;
        let r2 = PRIMES[u].r2;

        // Step 4a: build gm, igm at t1_start and t1_start+n.
        {
            let (head, rest) = tmp.split_at_mut(t1_start);
            let (gm, rest2) = rest.split_at_mut(n);
            let (igm, _) = rest2.split_at_mut(n);
            mp_mkgmigm(logn, gm, igm, PRIMES[u].g, PRIMES[u].ig, p, p0i);
            let _ = head;
        }

        // Step 4b: get gx (g mod p in NTT form) at t1_start + 2*n.
        if u < slen {
            // Copy gt[u*n .. (u+1)*n] into gx area (t1_start + 2*n).
            let gx_start = t1_start + 2 * n;
            tmp.copy_within(gt_start + u * n..gt_start + (u + 1) * n, gx_start);

            // iNTT ft[u*n..(u+1)*n] in place.
            {
                let (head, rest) = tmp.split_at_mut(t1_start);
                let (_gm, rest2) = rest.split_at_mut(n);
                let (igm, _) = rest2.split_at_mut(n);
                let ft_u = &mut head[ft_start + u * n..ft_start + (u + 1) * n];
                mp_intt(logn, ft_u, igm, p, p0i);
            }
            // iNTT gt[u*n..(u+1)*n] in place.
            {
                let (head, rest) = tmp.split_at_mut(t1_start);
                let (_gm, rest2) = rest.split_at_mut(n);
                let (igm, _) = rest2.split_at_mut(n);
                let gt_u = &mut head[gt_start + u * n..gt_start + (u + 1) * n];
                mp_intt(logn, gt_u, igm, p, p0i);
            }
        } else {
            // gt is in plain representation; compute g mod p via zint_mod_small_signed + NTT.
            let rx = mp_rx31(slen as u32, p, p0i, r2);
            let gx_start = t1_start + 2 * n;
            for v in 0..n {
                let z = zint_mod_small_signed(&tmp[gt_start + v..], slen, n, p, p0i, r2, rx);
                tmp[gx_start + v] = z;
            }
            // NTT gx using gm.
            let (head, rest) = tmp.split_at_mut(t1_start);
            let (gm, rest2) = rest.split_at_mut(n);
            let (_igm, rest3) = rest2.split_at_mut(n);
            let gx = &mut rest3[..n];
            mp_ntt(logn, gx, gm, p, p0i);
            let _ = head;
        }

        // Step 4c: NTT the hn entries Fe+hn = Ft[u*n + hn .. u*n + n] with logn-1.
        {
            let (fe_region, rest) = tmp.split_at_mut(t1_start);
            let (gm, _) = rest.split_at_mut(n);
            let fe_hn = &mut fe_region[u * n + hn..u * n + n];
            mp_ntt(logn - 1, fe_hn, gm, p, p0i);
        }

        // Step 4d: compute F (unreduced) mod p.
        // For v in 0..hn:
        //   ga = gx[2v], gb = gx[2v+1]
        //   mFp = mp_montymul(Fe[v + hn], r2, p, p0i)
        //   Fe[2v]   = mp_montymul(gb, mFp, p, p0i)
        //   Fe[2v+1] = mp_montymul(ga, mFp, p, p0i)
        {
            // Read gx from t1+2*n; write Fe = tmp[u*n..(u+1)*n].
            // Split at t1_start to get mutable Fe and shared gx.
            let (fe_region, rest) = tmp.split_at_mut(t1_start);
            let (_gm, rest2) = rest.split_at_mut(n);
            let (_igm, rest3) = rest2.split_at_mut(n);
            let gx = &rest3[..n];
            let fe = &mut fe_region[u * n..(u + 1) * n];
            for v in 0..hn {
                let ga = gx[(v << 1)];
                let gb = gx[(v << 1) + 1];
                let mfp = mp_montymul(fe[v + hn], r2, p, p0i);
                fe[v << 1] = mp_montymul(gb, mfp, p, p0i);
                fe[(v << 1) + 1] = mp_montymul(ga, mfp, p, p0i);
            }
        }

        // Step 4e: iNTT the new F row (full logn).
        {
            let (ft_region, rest) = tmp.split_at_mut(t1_start);
            let (_gm, rest2) = rest.split_at_mut(n);
            let (igm, _) = rest2.split_at_mut(n);
            let fe = &mut ft_region[u * n..(u + 1) * n];
            mp_intt(logn, fe, igm, p, p0i);
        }
    }

    // g is no longer needed; in the C, t1 is moved to gt. We leave t1_start at
    // fd_start (same region) since we already have t1_start = fd_start.
    // The C sets t1 = gt, which is ft_start + slen*n = gt_start.  We do the same.
    let t1_start = gt_start;

    // =========================================================================
    // Step 5: Edge case — if slen == llen, f hasn't been CRT-rebuilt yet.
    // =========================================================================
    if slen == llen {
        let (head, tail) = tmp.split_at_mut(ft_start);
        let (ft_region, rest) = tail.split_at_mut(slen * n);
        zint_rebuild_crt(ft_region, slen, n, 1, true, rest);
        let _ = head;
    }

    // =========================================================================
    // Step 6: Rebuild unreduced F (plain representation) from RNS.
    //
    // CRT scratch must start at gt_start (after ft) to avoid clobbering ft.
    // In C: t1 = gt = tmp + ft_start + slen*n = gt_start.
    // We split at gt_start so that `right` starts after ft.
    // =========================================================================
    {
        let (left, right) = tmp.split_at_mut(gt_start);
        let ft_only = &mut left[..llen * n];
        zint_rebuild_crt(ft_only, llen, n, 1, true, right);
    }

    // =========================================================================
    // Step 7: Babai round-off reduction loop.
    //
    // If depth > 1 and logn >= MIN_LOGN_FGNTT, use poly_sub_scaled_ntt and
    // convert f to NTT over slen+1 primes (needs extra n words after ft).
    // =========================================================================
    let use_sub_ntt = depth > 1 && logn >= MIN_LOGN_FGNTT;
    // If use_sub_ntt, t1 advances by n more (to leave room for the extra NTT word).
    let t1_start = if use_sub_ntt { t1_start + n } else { t1_start };

    // Sizes for fixed-point approximation window.
    let rlen = {
        let r = prof.word_win[depth as usize] as usize;
        if r > slen {
            slen
        } else {
            r
        }
    };
    let blen = slen - rlen;
    let ftb_off = ft_start + blen * n; // top rlen words of ft
    let scale_fg: u32 = 31 * blen as u32;
    let scale_FG_init: u32 = 31 * llen as u32;

    // Compute scale_x from the maximum bit-length of the top rlen words.
    let scale_x = poly_max_bitlength(logn, &tmp[ftb_off..ftb_off + rlen * n], rlen);
    // scale_t: target scale so that coefficients don't overflow in FFT.
    // scale_t <= 15 - logn, but clamp to scale_x.
    let scale_t_target: u32 = 15u32.wrapping_sub(logn);
    let scale_t =
        scale_t_target ^ (scale_t_target ^ scale_x) & tbmask(scale_x.wrapping_sub(scale_t_target));
    let scdiff = scale_x.wrapping_sub(scale_t);

    // Step 7a: Compute rt3 = poly_big_to_fixed(f top rlen words, scdiff).
    // Then: rt3 <- vect_fft(rt3); rt3 <- vect_inv_mul2e_fft(rt3, scale_t).
    // rt3 lives in a local Vec<Fxr> to avoid unsafe pointer casts.
    let mut rt3 = vec![Fxr(0); n];
    {
        let ftb = &tmp[ftb_off..ftb_off + rlen * n];
        poly_big_to_fixed(logn, &mut rt3, ftb, rlen, scdiff);
    }
    vect_fft(logn, &mut rt3);
    vect_inv_mul2e_fft(logn, &mut rt3, scale_t);

    // Step 7b: If use_sub_ntt, convert f to NTT over slen+1 primes and store
    //          the result in tmp[ft_start..ft_start + (slen+1)*n].
    //
    // t2 region starts at t1_start (after the potential +n shift).
    // We use a temporary local Vec for gm and the slen+1 NTT rows.
    if use_sub_ntt {
        let mut gm_buf = vec![0u32; n];
        let mut tn_buf = vec![0u32; (slen + 1) * n];
        for u in 0..=slen {
            let p = PRIMES[u].p;
            let p0i = PRIMES[u].p0i;
            let r2 = PRIMES[u].r2;
            let rx = mp_rx31(slen as u32, p, p0i, r2);
            mp_mkgm(logn, &mut gm_buf, PRIMES[u].g, p, p0i);
            let base = u * n;
            for v in 0..n {
                tn_buf[base + v] =
                    zint_mod_small_signed(&tmp[ft_start + v..], slen, n, p, p0i, r2, rx);
            }
            mp_ntt(logn, &mut tn_buf[base..base + n], &gm_buf, p, p0i);
        }
        tmp[ft_start..ft_start + (slen + 1) * n].copy_from_slice(&tn_buf);
    }

    // Step 7c: Reduction loop.
    let mut fg_len = llen;
    let mut scale_FG = scale_FG_init;
    let mut k_buf = vec![0i32; n];
    let mut rt1 = vec![Fxr(0); n];

    loop {
        // 7c.i: Convert F to fixed-point at scale_FG + scale_x.
        let (tlen, toff) = divrem31(scale_FG);
        let tlen = tlen as usize;
        {
            let src = &tmp[tlen * n..(tlen + (fg_len - tlen)) * n];
            poly_big_to_fixed(
                logn,
                &mut rt1,
                src,
                fg_len - tlen,
                scale_x.wrapping_add(toff),
            );
        }

        // 7c.ii: rt1 <- (F * adj(f)) / (f * adj(f))
        vect_fft(logn, &mut rt1);
        vect_mul_fft(logn, &mut rt1, &rt3);
        vect_ifft(logn, &mut rt1);

        // 7c.iii: k <- round(rt1)
        for u in 0..n {
            k_buf[u] = rt1[u].round();
        }

        // 7c.iv: Apply poly_sub variant.
        let scale_k = scale_FG.wrapping_sub(scale_fg);

        if depth == 1 {
            // poly_sub_kf_scaled_depth1's `k` parameter is interpreted as
            // signed on input and holds Montgomery-domain u32 products on
            // output; declaring it `&mut [u32]` matches the C API shape.
            //
            // Safe conversion: copy into a Vec<u32> (bit-preserving cast
            // from i32), pass, discard — avoids unsafe slice aliasing.
            let mut k_as_u32: Vec<u32> = k_buf.iter().map(|&v| v as u32).collect();
            let t2_start = t1_start;
            let (head, tail) = tmp.split_at_mut(t2_start);
            let f_big = &mut head[..fg_len * n];
            poly_sub_kf_scaled_depth1(logn_top, f_big, fg_len, &mut k_as_u32, scale_k, f, tail);
            // `k_buf` no longer read after this branch in the iteration; the
            // write-back into `k_as_u32` is intentionally discarded.
        } else if use_sub_ntt {
            // ft NTT rows are at tmp[ft_start .. ft_start + (slen+1)*n].
            // Ft (read-write) is at tmp[0 .. fg_len*n].
            // t2 scratch is at tmp[ft_start + (slen+1)*n ..].
            //
            // Split at ft_start: head = Ft region, tail = ft rows + t2 scratch.
            // Copy ft rows to a local vec so we can pass Ft as &mut and t2 as &mut
            // simultaneously (the borrow checker can't see their disjointness within tail).
            let ft_ntt_copy = tmp[ft_start..ft_start + (slen + 1) * n].to_vec();
            let (head, tail) = tmp.split_at_mut(ft_start);
            let f_big = &mut head[..fg_len * n];
            // t2 scratch starts after the ft rows (which we've already copied out).
            let t2_scratch = &mut tail[(slen + 1) * n..];
            poly_sub_scaled_ntt(
                logn,
                f_big,
                fg_len,
                &ft_ntt_copy,
                slen,
                &k_buf,
                scale_k,
                t2_scratch,
            );
        } else {
            // poly_sub_scaled: Ft at tmp[0..fg_len*n], ft at tmp[ft_start..].
            let (head, tail) = tmp.split_at_mut(ft_start);
            let f_big = &mut head[..fg_len * n];
            let ft_slice = &tail[..slen * n];
            poly_sub_scaled(logn, f_big, fg_len, ft_slice, slen, &k_buf, scale_k);
        }

        // 7c.v: Adjust scale_FG; break when fully reduced.
        if scale_FG <= scale_fg {
            break;
        }
        if scale_FG <= scale_fg + prof.reduce_bits {
            scale_FG = scale_fg;
        } else {
            scale_FG -= prof.reduce_bits;
        }

        // 7c.vi: Possibly shrink fg_len and verify redundant high words.
        while fg_len > slen && 31u32 * (fg_len - slen) as u32 > scale_FG.wrapping_sub(scale_fg) + 30
        {
            fg_len -= 1;
            // xp = Ft + (fg_len - 1) * n
            let xp_start = (fg_len - 1) * n;
            for u in 0..n {
                let sw = (tmp[xp_start + u] >> 30).wrapping_neg() >> 1;
                if tmp[xp_start + u + n] != sw {
                    return SOLVE_ERR_REDUCE;
                }
            }
        }
    }

    SOLVE_OK
}

// === Sub-task 15.20: solve_ntru_depth0 ===

/// Solve the NTRU equation at the top recursion level (depth == 0).
///
/// On entry, `tmp[0..hn]` holds `Fd` (the half-degree F from depth 1),
/// in plain integer representation (each word is a signed value stored in
/// u32 two's-complement).
///
/// On success returns `SOLVE_OK` with `tmp[0..n]` = F and `tmp[n..2n]` = G,
/// both in plain signed-integer representation.
///
/// Port of `solve_NTRU_depth0` (c-reference/hawk-512/ng_ntru.c:718-977).
pub fn solve_ntru_depth0(
    prof: &NtruProfile,
    logn: u32,
    f: &[i8],
    g: &[i8],
    tmp: &mut [u32],
) -> i32 {
    use crate::keygen::fxp::{vect_div_autoadj_fft, vect_fft, vect_ifft};
    use crate::keygen::fxr::Fxr;
    use crate::keygen::mp31::{mp_add, mp_div, mp_norm, mp_set, mp_sub};
    use crate::keygen::poly::{poly_mp_norm, poly_mp_set};

    let n = 1usize << logn;
    let hn = n >> 1;
    let p = PRIMES[0].p;
    let p0i = PRIMES[0].p0i;
    let r2 = PRIMES[0].r2;

    assert!(
        tmp.len() >= 5 * n,
        "tmp too small: {} < {}",
        tmp.len(),
        5 * n
    );

    // ----------------------------------------------------------------
    // Step 1: Convert Fd (at tmp[0..hn]) to RNS+NTT, into t3[0..hn].
    //   gm lives in t4 for this step.
    // ----------------------------------------------------------------
    {
        let (t1, rest) = tmp.split_at_mut(n);
        let (_t2, rest) = rest.split_at_mut(n);
        let (t3, rest) = rest.split_at_mut(n);
        let (t4, _) = rest.split_at_mut(n);
        mp_mkgm(logn, t4, PRIMES[0].g, p, p0i);
        // poly_mp_set converts t1[0..hn] (plain int) -> RNS in [0, p).
        poly_mp_set(logn - 1, &mut t1[..hn], p);
        mp_ntt(logn - 1, &mut t1[..hn], t4, p, p0i);
        t3[..hn].copy_from_slice(&t1[..hn]);
    }

    // ----------------------------------------------------------------
    // Step 2: Compute unreduced F (RNS+NTT) into t1.
    //   t2 <- g (RNS+NTT); use g(X^2) trick to extend from hn -> n.
    //   t1[2u]   = g[2u+1] * Fd[u]   (mod p, Montgomery)
    //   t1[2u+1] = g[2u]   * Fd[u]   (mod p, Montgomery)
    // ----------------------------------------------------------------
    {
        let (t1, rest) = tmp.split_at_mut(n);
        let (t2, rest) = rest.split_at_mut(n);
        let (t3, rest) = rest.split_at_mut(n);
        let (t4, _) = rest.split_at_mut(n);
        poly_mp_set_small(logn, t2, g, p);
        mp_ntt(logn, t2, t4, p, p0i);
        for u in 0..hn {
            let ga = t2[(u << 1)];
            let gb = t2[(u << 1) + 1];
            let mf = mp_montymul(t3[u], r2, p, p0i);
            t1[(u << 1)] = mp_montymul(gb, mf, p, p0i);
            t1[(u << 1) + 1] = mp_montymul(ga, mf, p, p0i);
        }
    }

    // ----------------------------------------------------------------
    // Step 3: Convert f to RNS+NTT into t3; check all evaluations != 0.
    // ----------------------------------------------------------------
    {
        let (_t1, rest) = tmp.split_at_mut(n);
        let (_t2, rest) = rest.split_at_mut(n);
        let (t3, rest) = rest.split_at_mut(n);
        let (t4, _) = rest.split_at_mut(n);
        poly_mp_set_small(logn, t3, f, p);
        mp_ntt(logn, t3, t4, p, p0i);
        for u in 0..n {
            if t3[u] == 0 {
                return SOLVE_ERR_REDUCE;
            }
        }
    }

    // ----------------------------------------------------------------
    // Step 4: Compute kn (numerator) and kd (denominator) for k.
    //   G = (q + g*F) / f  (exact mod p)
    //   kn[u] = f*[n-1-u] * F[u] + g*[n-1-u] * G[u]  (in RNS+NTT)
    //   kd[u] = f[u] * f*[n-1-u] + g[u] * g*[n-1-u]  (denominator)
    // The auto-adjoint denominator satisfies kd[u] == kd[n-1-u].
    // Results stored into t2 (kn) and t3 (kd).
    // ----------------------------------------------------------------
    let q = prof.q;
    {
        let (t1, rest) = tmp.split_at_mut(n);
        let (t2, rest) = rest.split_at_mut(n);
        let (t3, _) = rest.split_at_mut(n);
        for u in 0..hn {
            let tf0 = t3[u];
            let tf1 = t3[n - 1 - u];
            let tg0 = t2[u];
            let tg1 = t2[n - 1 - u];
            let tf0_large = t1[u];
            let tf1_large = t1[n - 1 - u];
            let mf0 = mp_montymul(tf0, r2, p, p0i);
            let mf1 = mp_montymul(tf1, r2, p, p0i);
            let mg0 = mp_montymul(tg0, r2, p, p0i);
            let mg1 = mp_montymul(tg1, r2, p, p0i);
            let tg0_val = mp_div(mp_add(q, mp_montymul(mg0, tf0_large, p, p0i), p), tf0, p);
            let tg1_val = mp_div(mp_add(q, mp_montymul(mg1, tf1_large, p, p0i), p), tf1, p);
            let kn0 = mp_add(
                mp_montymul(mf1, tf0_large, p, p0i),
                mp_montymul(mg1, tg0_val, p, p0i),
                p,
            );
            let kn1 = mp_add(
                mp_montymul(mf0, tf1_large, p, p0i),
                mp_montymul(mg0, tg1_val, p, p0i),
                p,
            );
            let kd = mp_add(
                mp_montymul(mf0, tf1, p, p0i),
                mp_montymul(mg0, tg1, p, p0i),
                p,
            );
            t2[u] = kn0;
            t2[n - 1 - u] = kn1;
            t3[u] = kd;
            t3[n - 1 - u] = kd;
        }
    }

    // ----------------------------------------------------------------
    // Step 5: iNTT t2, t3; convert to signed via mp_norm.
    // ----------------------------------------------------------------
    {
        let (_t1, rest) = tmp.split_at_mut(n);
        let (t2, rest) = rest.split_at_mut(n);
        let (t3, rest) = rest.split_at_mut(n);
        let (t4, _) = rest.split_at_mut(n);
        mp_mkigm(logn, t4, PRIMES[0].ig, p, p0i);
        mp_intt(logn, t2, t4, p, p0i);
        mp_intt(logn, t3, t4, p, p0i);
        for u in 0..n {
            t2[u] = mp_norm(t2[u], p) as u32;
            t3[u] = mp_norm(t3[u], p) as u32;
        }
    }

    // ----------------------------------------------------------------
    // Step 6: Convert numerator (t2) and denominator (t3) to FFT
    //   representation, downscaled by 2^DOWNSCALE.
    //   rt4 = FFT of denominator (auto-adjoint: copy [0..hn] to [hn..n]).
    //   rt3 = FFT of numerator.
    // ----------------------------------------------------------------
    const DOWNSCALE: u32 = 10;
    // Read t3 (denominator, at tmp[2n..3n]) and t2 (numerator, at tmp[n..2n])
    // into local Vecs before mutating tmp further.
    let denom_plain: Vec<i32> = tmp[2 * n..3 * n].iter().map(|&w| w as i32).collect();
    let numer_plain: Vec<i32> = tmp[n..2 * n].iter().map(|&w| w as i32).collect();

    let mut rt4: Vec<Fxr> = denom_plain
        .iter()
        .map(|&s| {
            let x = ((s as i64) << (32 - DOWNSCALE)) as u64;
            Fxr::of_scaled32(x)
        })
        .collect();
    vect_fft(logn, &mut rt4);
    // Denominator is auto-adjoint: FFT has half-size symmetry.
    for u in 0..hn {
        rt4[hn + u] = rt4[u];
    }

    let mut rt3: Vec<Fxr> = numer_plain
        .iter()
        .map(|&s| {
            let x = ((s as i64) << (32 - DOWNSCALE)) as u64;
            Fxr::of_scaled32(x)
        })
        .collect();
    vect_fft(logn, &mut rt3);

    // ----------------------------------------------------------------
    // Step 7: Divide numerator by denominator in FFT domain, round to k.
    //   rt5 = rt4[hn..n] (second half of denominator FFT, half-size).
    // ----------------------------------------------------------------
    {
        let rt5: &[Fxr] = &rt4[hn..n];
        vect_div_autoadj_fft(logn, &mut rt3, rt5);
    }
    vect_ifft(logn, &mut rt3);

    // ----------------------------------------------------------------
    // Step 8: Round rt3 into t2 (RNS, using gm still in t4 below).
    // ----------------------------------------------------------------
    for u in 0..n {
        tmp[n + u] = mp_set(rt3[u].round(), p);
    }

    // ----------------------------------------------------------------
    // Step 9: Reload gm into t5; reload f, g into t3, t4.
    // ----------------------------------------------------------------
    {
        let (_t1_t2, rest) = tmp.split_at_mut(2 * n);
        let (t3, rest) = rest.split_at_mut(n);
        let (t4, t5) = rest.split_at_mut(n);
        mp_mkgm(logn, t5, PRIMES[0].g, p, p0i);
        poly_mp_set_small(logn, t3, f, p);
        poly_mp_set_small(logn, t4, g, p);
    }

    // ----------------------------------------------------------------
    // Step 10: NTT of t2 (k), t3 (f), t4 (g) using gm in t5.
    // ----------------------------------------------------------------
    {
        let (_t1, rest) = tmp.split_at_mut(n);
        let (t2, rest) = rest.split_at_mut(n);
        let (t3, rest) = rest.split_at_mut(n);
        let (t4, t5) = rest.split_at_mut(n);
        mp_ntt(logn, t2, t5, p, p0i);
        mp_ntt(logn, t3, t5, p, p0i);
        mp_ntt(logn, t4, t5, p, p0i);
    }

    // ----------------------------------------------------------------
    // Step 11: Reduce F by F -= k*f; recompute G = (q + g*F) / f.
    // ----------------------------------------------------------------
    {
        let (t1, rest) = tmp.split_at_mut(n);
        let (t2, rest) = rest.split_at_mut(n);
        let (t3, rest) = rest.split_at_mut(n);
        let (t4, _) = rest.split_at_mut(n);
        for u in 0..n {
            let tf = t3[u];
            let tg = t4[u];
            let tk = t2[u];
            let mf = mp_montymul(tf, r2, p, p0i);
            let mg = mp_montymul(tg, r2, p, p0i);
            let tf_new = mp_sub(t1[u], mp_montymul(mf, tk, p, p0i), p);
            let tg_new = mp_div(mp_add(q, mp_montymul(mg, tf_new, p, p0i), p), tf, p);
            t1[u] = tf_new;
            t2[u] = tg_new;
        }
    }

    // ----------------------------------------------------------------
    // Step 12: iNTT t1 (F), t2 (G); normalize to signed integers.
    // ----------------------------------------------------------------
    {
        let (t1, rest) = tmp.split_at_mut(n);
        let (t2, rest) = rest.split_at_mut(n);
        let (t3, _) = rest.split_at_mut(n);
        mp_mkigm(logn, t3, PRIMES[0].ig, p, p0i);
        mp_intt(logn, t1, t3, p, p0i);
        mp_intt(logn, t2, t3, p, p0i);
        poly_mp_norm(logn, t1, p);
        poly_mp_norm(logn, t2, p);
    }

    SOLVE_OK
}

// === Sub-task 15.21: solve_ntru orchestrator ===

/// Solve the full NTRU equation end-to-end.
///
/// Calls `solve_ntru_deepest`, then `solve_ntru_intermediate` for each
/// depth from `logn-1` down to 1, then `solve_ntru_depth0`.
/// On success, converts (F, G) from signed-31-bit to int8 with a range
/// check, and writes them packed in `tmp[0..2*n]` (one i8 per u32, low byte).
///
/// Returns `SOLVE_OK` on success, or a negative error code on failure.
///
/// Port of `solve_NTRU` (ng_ntru.c:980-1019).
pub fn solve_ntru(prof: &NtruProfile, logn: u32, f: &[i8], g: &[i8], tmp: &mut [u32]) -> i32 {
    use crate::keygen::ntru_profile::SOLVE_ERR_LIMIT;
    use crate::keygen::poly::poly_big_to_small;

    let n = 1usize << logn;

    let err = solve_ntru_deepest(prof, logn, f, g, tmp);
    if err != SOLVE_OK {
        return err;
    }

    let mut depth = logn;
    while depth > 1 {
        depth -= 1;
        let err = solve_ntru_intermediate(prof, logn, f, g, depth, tmp);
        if err != SOLVE_OK {
            return err;
        }
    }

    let err = solve_ntru_depth0(prof, logn, f, g, tmp);
    if err != SOLVE_OK {
        return err;
    }

    // Post-processing: convert F and G from signed-31-bit to int8 with range
    // check, then pack back as raw bytes (matching C's memmove).
    //
    // C:
    //   int8_t *F = (int8_t *)(tmp + 2*n);  // byte ptr at word-offset 2*n
    //   int8_t *G = F + n;
    //   poly_big_to_small(logn, F, tmp,   lim) -> n bytes at F
    //   poly_big_to_small(logn, G, tmp+n, lim) -> n bytes at G
    //   memmove(tmp, F, 2*n)  // copy 2*n bytes back to tmp[0]
    //
    // The memmove packs i8 values 4-per-u32 (little-endian on x86).
    let lim = prof.coeff_fg_limit[logn as usize] as i32;

    let mut f_out = vec![0i8; n];
    let mut g_out = vec![0i8; n];
    if !poly_big_to_small(logn, &mut f_out, &tmp[..n], lim) {
        return SOLVE_ERR_LIMIT;
    }
    if !poly_big_to_small(logn, &mut g_out, &tmp[n..2 * n], lim) {
        return SOLVE_ERR_LIMIT;
    }

    // Pack f_out[0..n] || g_out[0..n] as 2*n bytes into tmp[0..n/2].
    // n i8 values occupy n bytes; packed into n/4 u32 words (little-endian).
    let words_per = n / 4; // words needed for n i8 values
    for w in 0..words_per {
        let i = w * 4;
        tmp[w] = u32::from_le_bytes([
            f_out[i] as u8,
            f_out[i + 1] as u8,
            f_out[i + 2] as u8,
            f_out[i + 3] as u8,
        ]);
    }
    for w in 0..words_per {
        let i = w * 4;
        tmp[words_per + w] = u32::from_le_bytes([
            g_out[i] as u8,
            g_out[i + 1] as u8,
            g_out[i + 2] as u8,
            g_out[i + 3] as u8,
        ]);
    }

    SOLVE_OK
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::keygen::ntru_profile::SOLVE_HAWK_512;

    #[test]
    fn solve_ntru_deepest_on_trivial_input_does_not_panic() {
        // Use f = [1, 0, 0, ...] (polynomial 1) and g = [1, 0, 0, ...].
        // This is a degenerate case that may fail gcd, but shouldn't crash.
        let n = 512usize;
        let mut f = vec![0i8; n];
        let mut g = vec![0i8; n];
        f[0] = 1;
        g[0] = 1;
        // tmp needs room for lots — use a generous size.
        let mut tmp = vec![0u32; 32 * n];
        let r = solve_ntru_deepest(&SOLVE_HAWK_512, 9, &f, &g, &mut tmp);
        // Accept any outcome (gcd failure is expected; just verify no panic).
        let _ = r;
    }

    #[test]
    fn solve_ntru_intermediate_at_depth_0_no_panic() {
        // Placeholder test that ensures the function doesn't panic on zero input.
        // Real correctness validation comes from FFI cross-check (15.19d).
        let n = 512;
        let mut f = vec![0i8; n];
        let mut g = vec![0i8; n];
        f[0] = 1;
        g[0] = 1;
        // tmp needs generous size.
        let mut tmp = vec![0u32; 64 * n];
        let r = solve_ntru_intermediate(&SOLVE_HAWK_512, 9, &f, &g, 0, &mut tmp);
        let _ = r;
    }

    // === Sub-task 15.20: solve_ntru_depth0 tests ===

    #[test]
    fn solve_ntru_depth0_no_panic() {
        // Dummy input; ensures the code path runs without panic/OOB.
        // Real correctness is validated in 15.20b (FFI cross-check).
        let n = 512;
        let mut f = vec![0i8; n];
        let mut g = vec![0i8; n];
        f[0] = 1;
        g[0] = 1;
        let mut tmp = vec![0u32; 64 * n];
        // tmp[0..hn] = Fd from depth 1 (all zeros is a valid no-panic input).
        let _ = solve_ntru_depth0(&SOLVE_HAWK_512, 9, &f, &g, &mut tmp);
    }
}
