# Security policy for pqe-hawk

## Current status: SCHEME BROKEN — research/reference only

### 🚫 DO NOT USE IN PRODUCTION — the HAWK scheme is cryptographically broken.

> **HAWK was withdrawn from NIST standardization on 2026-07-29.** The HAWK design
> team (Léo Ducas et al.) withdrew HAWK from NIST's additional-signature
> standardization process following the key-recovery attack of Strážnickas & Weis
> (Anthropic), *"HAWK-n Key Recovery Reduces to SVP in Dimension n/2+1."* The team
> confirmed the attack **approximately halves the lattice-reduction block size**
> required to recover an equivalent secret key, and stated that naïve
> countermeasures (doubling parameters, higher-rank modules) make HAWK
> uncompetitive.

This crate is a faithful Rust port of the upstream HAWK reference C implementation (github.com/hawk-sign/dev@1b9fef5, MIT licensed). The port is byte-exact against that reference, but **the scheme it implements is cryptographically broken.** Earlier releases described this crate as "testnet-ready"; that framing is retracted — the crate is now suitable only for research, reproducibility, and reference.

**Do not use this crate to secure live funds, production authentication, or any setting where a compromise has material consequences.**

### The break, concretely

- **What it is.** A public, deterministic, polynomial-time reduction: HAWK-*n* key recovery reduces to one exact-SVP call in dimension *n*/2+1 (vs. the designers' assumed ≈dimension *n*). It reads only the public key and exploits the scheme's structure — `det B = 1` and public key = Gram matrix `B*B`.
- **Impact on HAWK-512** (this crate's only parameter set): key-recovery cost drops from the designers' claimed ≈2¹⁵⁰ gates to **≤2¹⁰⁸ gates** (AGPS20 model). **HAWK-512 no longer meets its claimed NIST Level I security.** HAWK-256 has been recovered end-to-end in practice.
- **Why the port cannot fix it.** The attack targets the scheme's mathematics, not any implementation detail. Nothing in this port (sampler, encoding, zeroization, constant-time work below) affects it; a byte-exact port faithfully reproduces the vulnerability. It does not transfer to Falcon, and is evaded only by odd-prime-power conductors — not a change available within HAWK's parameter sets.

## Implementation hardening posture

The measures below concern the *implementation* only. They do not — and cannot — address the scheme-level break described above; they are documented because the port's fidelity and implementation hygiene remain useful for research and reference:

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

As of 2026-07-29 the upstream HAWK team has **withdrawn HAWK from NIST
standardization** (see the status section above). Their pre-existing warning, per
the C reference `README`, was already:

> WARNING: This code has not been audited. HAWK itself is a relatively recent scheme; the security reduction is not yet as well-studied as, e.g., for Falcon. Use at your own risk.

That warning is now superseded by the confirmed key-recovery break, not merely an
absence of review.

## Fixed in 0.1.1

Two issues found and fixed in 0.1.1 (see `CHANGELOG.md` for detail). Both were
reproduced against this crate before being fixed, and both carry regression
tests. If you are on 0.1.0, upgrade — it is a drop-in change with no API break.

- **Weak-key BUFF break at verify** (Dao, eprint 2026/1298). The verifier
  accepted maliciously formed public keys with a tiny `q00[0]`, under which a
  single trivial signature verified for essentially every message (measured:
  256 of 256). Fixed by a `KeyNormCheck` floor rejecting `q00[0] < 2080`.
  Honestly generated keys were never affected, and EUF/SUF-CMA security under
  honest keys was not affected — only the BUFF add-on properties.
- **Signing could emit an unserializable signature.** An in-bounds `s1` can
  still overflow the fixed 555-byte Golomb-Rice buffer in aggregate; signing
  now resigns on overflow instead of returning a signature that fails to
  encode. Observed in practice stalling a live testnet on block commit.

## Reporting vulnerabilities

If you find a security issue in this port, please report it privately via
[GitHub Security Advisories](https://github.com/nixprotocol/pqe-hawk/security/advisories/new)
rather than opening a public issue, so a fix can be prepared before disclosure.

If the issue is in HAWK itself rather than this port, please also contact the
upstream HAWK team (https://github.com/hawk-sign/dev).
