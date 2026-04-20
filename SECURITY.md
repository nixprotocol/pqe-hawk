# Security policy for pqe-hawk

## Current status: testnet-ready, not audited

This crate is a faithful Rust port of the upstream HAWK reference C implementation (github.com/hawk-sign/dev@1b9fef5, MIT licensed). It is **testnet-ready** — suitable for non-custodial testnet use where a failure has no monetary impact — but has **not** been audited for production deployment.

The upstream HAWK scheme itself has not yet received full security review (see upstream notes: "no security review yet — use at own risk").

**Do not use this crate to secure live funds, production authentication, or any setting where a compromise has material consequences, without performing your own security analysis first.**

## Hardening posture

Testnet-readiness has been reached via the following measures:

- **Secret key zeroization.** `HawkSecretKey` derives `zeroize::ZeroizeOnDrop`, so its `f, g, F, G, seed` fields are zeroed when dropped. Verified by an explicit unit test.
- **Constant-time public-key equality.** `HawkPublicKey::ct_eq` uses `subtle::ConstantTimeEq` to avoid early-exit byte-by-byte comparison.
- **FFI cross-check.** 51 proptests (256 cases each) validate byte-exact equivalence with the vendored C reference at every primitive: mp31 arithmetic, NTT, big-int (`zint31`), fixed-point FFT, Gaussian sampler, NTRU solver, keygen, sign, verify.
- **NIST KAT pin.** The generated PQCsignKAT_HAWK-512.rsp from the reference harness is committed; case 0's 1024-byte pubkey and 184-byte secret key are asserted byte-for-byte against the Rust implementation.
- **Fuzz targets.** Three `cargo-fuzz` harnesses (`fuzz_decode_pubkey`, `fuzz_decode_signature`, `fuzz_verify`) assert that decoders and the verify pipeline never panic on arbitrary inputs. Stable-Rust proptest mirrors run in CI.
- **Malformed-input proptests.** Corrupting the trailing zero-padding region of a valid pubkey or signature is asserted to return `Err(HawkError::Malformed*)` with no panic. Pathological all-zeros / all-ones buffers exercise the unary tail decoder.
- **Timing smoke tests.** A dudect-style coarse test asserts that `sign` timings for two different keys stay within 20% of each other; catches gross timing regressions (e.g., an accidental data-dependent early-exit).

## Constant-time posture

The port preserves the branchless / mask-based idioms of the upstream C reference. In particular:

- Modular arithmetic (`mp_add`, `mp_sub`, `mp_half`, `mp_montymul`, `mp_div`) uses the upstream's sign-bit-mask pattern for conditional operations.
- The Gaussian sampler (`sig_gauss_alt`) uses CDT comparisons without data-dependent branches.
- The sym-break check (`poly_symbreak`) is a linear scan with carry-forward, not an early-exit.
- Public-key equality is constant-time (`ct_eq`).
- Secret material is zeroed on drop (`ZeroizeOnDrop`).

However:

- Rust's default bounds checks on slice indexing panic on out-of-range access. Since our code uses `wrapping_*` arithmetic and explicit index math, OOB accesses shouldn't occur with well-formed inputs — and the fuzz harness + malformed-input proptests exercise adversarial inputs against the decode/verify boundary.
- Allocator behavior (via `Vec<u16>`, `Vec<u32>`, etc.) is implementation-defined. Timing side-channels involving heap allocation are theoretically possible.
- No formal constant-time audit (e.g., dudect, ctgrind, Jasmin) has been performed. The built-in timing smoke test is a regression guard, not a proof.

## Upstream caveats

The upstream HAWK team explicitly notes (per the C reference `README`):

> WARNING: This code has not been audited. HAWK itself is a relatively recent scheme; the security reduction is not yet as well-studied as, e.g., for Falcon. Use at your own risk.

## Reporting vulnerabilities

If you find a security issue in this port, please open a GitHub issue or contact the maintainer at the address in `Cargo.toml`.
