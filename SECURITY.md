# Security policy for pqe-hawk

## Current status: SCHEME BROKEN — research/reference only

### 🚫 DO NOT USE IN PRODUCTION — the HAWK scheme is cryptographically broken.

> **HAWK was withdrawn from NIST standardization on 2026-07-29.** The HAWK design
> team (Léo Ducas et al.) withdrew HAWK from NIST's additional-signature
> standardization process following the key-recovery attack of Straznickas & Weis
> (Anthropic), *"HAWK-n Key Recovery Reduces to SVP in Dimension n/2+1."* The team
> confirmed the attack **approximately halves the lattice-reduction block size**
> required to recover an equivalent secret key, and stated that naïve
> countermeasures (doubling parameters, higher-rank modules) make HAWK
> uncompetitive.

This crate is a faithful Rust port of the upstream HAWK reference C implementation (github.com/hawk-sign/dev@1b9fef5, MIT licensed). The port is byte-exact against that reference, but **the scheme it implements is cryptographically broken.** Earlier releases described this crate as "testnet-ready"; that framing is retracted — the crate is now suitable only for research, reproducibility, and reference.

**Do not use this crate to secure live funds, production authentication, or any setting where a compromise has material consequences.**

### The break, concretely

- **What it is.** Straznickas & Weis (Anthropic), [eprint 2026/1593](https://eprint.iacr.org/2026/1593): a public, unconditional, deterministic polynomial-time reduction from HAWK-*n* key recovery to poly(*n*) calls to an exact-SVP oracle in dimension *n*/2+1 (vs. the designers' assumed ≈dimension *n*). The lever is a **nontrivial automorphism of the key lattice** — the Galois involution τ: ζ ↦ −ζ — which is recoverable as a *shortest vector of a public rank-n lattice* (Ducas's block reduction on that near-hypercubic class); the van Gent–Pulles descent then extracts an equivalent secret key. It reads only the public key.
- **Impact on HAWK-512** (this crate's only parameter set): key-recovery cost drops from the designers' claimed ≈2¹⁵⁰ gates to **≤2¹⁰⁸ gates** (AGPS20 model). **HAWK-512 no longer meets its claimed NIST Level I security.** HAWK-256 has been recovered end-to-end in practice, in a few hours on a single server.
- **Why the port cannot fix it.** The attack targets the scheme's mathematics, not any implementation detail. Nothing in this port (sampler, encoding, zeroization, constant-time work below) affects it; a byte-exact port faithfully reproduces the vulnerability. It does not transfer to Falcon, and is evaded only by conductors *m* ∈ {*p*ᵏ, 2*p*ᵏ} for odd prime *p* (equivalently, cyclotomic conductors *m* > 4 with cyclic (ℤ/*m*)ˣ) — not a change available within HAWK's parameter sets.

### Separately: key recovery from signing-side leakage

The break above needs only the public key. Independent work shows HAWK is *also* fragile when signing leaks, which is a different threat model and is **not** addressed by anything in this port:

- **Sign leakage.** Brinkmann, Kraus & May, *"Halfspace Learning for Lattice Signature Key Recovery from Signs"* ([eprint 2026/1366](https://eprint.iacr.org/2026/1366)): leaking the sign — or merely the Hamming weight — of a randomness coordinate recovers a HAWK key from **30 signatures in about 10 minutes** at the 128-bit level in the noise-free model, tolerating up to 35% sign error. The method is not HAWK-specific (Falcon falls at 100 signatures, ML-DSA at 190,000), but HAWK is the scheme this crate ships.
- **Sampler leakage.** Chi, Lee & Lee, *"HAWK with Hint: Algebraic Key Recovery from Side-Channel Leakage"* ([eprint 2026/699](https://eprint.iacr.org/2026/699)): only *partial* sampler leakage is needed, because HAWK's public algebraic structure amplifies it. Full-coefficient leakage recovers a HAWK-1024 key from a single signature; the exact-sign hint model reached 100% success over 100 independently generated keys from 14 signatures on the reference implementation.

Both require observing the signer, so they do not change the conclusion for a signer in isolation. They do mean the ≤2¹⁰⁸-gate figure is not what an attacker who can watch signing actually faces.

## Implementation hardening posture

The measures below concern the *implementation* only. They do not — and cannot — address the scheme-level break described above; they are documented because the port's fidelity and implementation hygiene remain useful for research and reference:

- **Secret key zeroization.** `HawkSecretKey` derives `zeroize::ZeroizeOnDrop`, so its `f`, `g`, `f_cap` (F), `g_cap` (G), and `seed` fields are zeroed when dropped. Verified by the explicit unit test `zeroize_on_drop_clears_secret_material`.
- **Constant-time public-key equality.** `HawkPublicKey::ct_eq` uses `subtle::ConstantTimeEq` to avoid early-exit byte-by-byte comparison.
- **Opaque verify errors.** `verify_inner` collapses every internal rejection into `HawkError::InvalidSignature`. The `VerifyReject` variants that name *which* check failed are `#[doc(hidden)]` and test-only, so the public API never reveals the failing check (which would aid malleability attacks).
- **FFI cross-check.** 50 tests in `tests/cross_check.rs` cover every primitive — mp31 arithmetic, NTT, big-int (`zint31`), fixed-point FFT, Gaussian sampler, NTRU solver, keygen, sign, verify. 49 byte-diff against the vendored C reference (47 proptest-driven at proptest's 256-case default, a few blocks lower; plus 2 direct diffs); the 50th is a Rust-only round-trip. All are gated behind `--features cross-check-reference-c`, which builds the C reference and runs bindgen, so they do **not** run under a default `cargo test` or in CI.
- **Fault-injection survey** (`tests/fault_survey.rs`). Drives the real signing inner loop through a doc-hidden seam and injects plausible faults — truncated noise, forced and skipped sym-break, single-coefficient flips, sign flips, skipped norm rejection — then feeds each faulted signature to the real verifier. Every model runs over 8 deterministic keys alongside an un-faulted control that must accept, proving the harness reproduces production. The suite also asserts that the opaque `verify()` and the labelled internal path never disagree on accept/reject.
- **Residual BUFF probing** (`tests/buff_residual.rs`). The 0.1.1 `KeyNormCheck` floor closes the constant-`q00` weak-key class but does not prove full BUFF restoration. This suite probes the explicitly-open residual classes against the real verifier: S-CEO / S-DEO key substitution, non-constant weak keys that pass the floor, and PolyQnorm soundness on `(q00, q01)` pairs that are not a genuine `||(f,g)||²` form.
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
- The signing-side leakage attacks described above are out of scope for this section: they exploit information the *signer* emits while signing, not the branch structure of this port. A constant-time implementation does not defend against them.

## Upstream caveats

As of 2026-07-29 the upstream HAWK team has **withdrawn HAWK from NIST
standardization** (see the status section above). Their pre-existing warning, as
recorded verbatim in `c-reference/PROVENANCE.md` (the upstream `README` itself is
not part of what we vendored), was already:

> This code was written for the on-ramp call for post-quantum signature schemes,
> organized by NIST and currently serves research purposes. There was no
> security review yet so use at own risk.

That warning is now superseded by the confirmed key-recovery break, not merely an
absence of review.

## Release notes: what to run

**Current release: 0.1.2.** It is documentation-only — no code behaviour
changed — so everything in this document applies to it unchanged. The latest
*security-relevant* changes shipped in **0.1.1**: two issues found and fixed
(see `CHANGELOG.md` for detail). Both were reproduced against this crate before
being fixed, and both carry regression tests. If you are on 0.1.0, upgrade — it
is a drop-in change with no API break.

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
