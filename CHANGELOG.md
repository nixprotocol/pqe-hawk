# Changelog

All notable changes to this crate are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Security

- **HAWK withdrawn from NIST standardization (2026-07-29); the scheme is
  broken.** The HAWK design team (Léo Ducas et al.) withdrew HAWK from NIST's
  additional-signature standardization process following the key-recovery attack
  of Strážnickas & Weis (Anthropic), *"HAWK-n Key Recovery Reduces to SVP in
  Dimension n/2+1."* The attack gives a public, deterministic, polynomial-time
  reduction from HAWK-*n* key recovery to one exact-SVP call in dimension
  *n*/2+1, approximately halving the lattice-reduction block size versus the
  designers' assumption. For HAWK-512 it lowers key recovery from a claimed
  ≈2^150 gates to ≤2^108 gates (AGPS20 model), so **HAWK-512 no longer meets its
  claimed NIST Level I security**; HAWK-256 has been recovered end-to-end. The
  attack reads only the public key and exploits the scheme's structure
  (`det B = 1`, public key = Gram matrix `B*B`), so it **cannot** be mitigated in
  this port — a byte-exact port faithfully reproduces the vulnerability.

### Changed

- Documentation retracts the prior "testnet-ready" / "NIST Level I" framing
  across `README.md`, `SECURITY.md`, `Cargo.toml`, and the `src/lib.rs` crate
  header: this crate is now **research/reference only**. No code behavior
  changed — keygen, sign, and verify remain byte-exact against the reference C.

## 0.1.1 - 2026-07-21

Two correctness fixes. No API changes, so this is a drop-in upgrade from 0.1.0.
Both issues were confirmed by reproduction against this crate before fixing, and
each fix is covered by a regression test proven to fail without it.

### Fixed

- **Weak public keys broke the BUFF properties at verify time** (Dao, [eprint
  2026/1298](https://eprint.iacr.org/2026/1298)). The verifier rejected only a
  negative `q00[0]`, never a small positive one, so a maliciously formed public
  key with a tiny `q00[0]` (the constant `q00 = 1` weak-key class, preimage
  `(1,0,a,1)`) was accepted. A single trivial signature (`salt = 0`, `s1 = 0`)
  then verified under such a key for essentially every message, breaking
  message-bound signatures (MBS) and M-S-UEO. Measured on this crate before the
  fix: the witness verified for **256 of 256** distinct messages; after, 0 of
  256.

  `verify_inner` now enforces a `KeyNormCheck` floor, rejecting `q00[0] < 2080`.
  For an auto-adjoint key `q00[0] = ||(f,g)||^2`, and keygen resamples until
  that reaches `HAWK_512_L2LOW = 2080`, so the check is exact for HAWK-512 and a
  no-op for honest keys (measured over 48 honest keys, `q00[0]` spans
  `[2080, 2290]`). Note this is the HAWK-512 floor, not the HAWK-256 value 556,
  which would leave a `[556, 2080)` band no honest key occupies.

  This does **not** affect EUF/SUF-CMA security under honestly generated keys;
  only the BUFF add-on properties over verifier-accepted-but-malformed keys were
  affected.

- **Signing could return a signature that fails to serialize.** `sign_512`
  encoded outside its rejection loop, assuming an in-bounds `s1` always fits the
  fixed 555-byte wire buffer. The Golomb-Rice variable section grows with
  `sum(|s1[i]| >> 5)`, so a signature can be in bounds per coefficient
  (`|s1[i]| < 512`) yet overflow in aggregate; `encode_gr` then returned an
  error only later, when a caller tried to serialize. This was observed in
  practice, stalling a live single-validator testnet on block commit.

  Encoding now happens inside the retry loop, and an overflow triggers a resign
  with a fresh salt rather than returning an unserializable signature. This
  restores the invariant that every returned `HawkSignature` encodes to 555
  bytes, and matches the reference C, which restarts on `encode_sig` failure
  (`hawk_sign.c:1142-1158`). Overflow shares the existing
  `HAWK_SAMPLER_RETRY_BUDGET`; measured upstream at 0 overflow-driven resigns
  over 20,000 signings, so the shared budget has ample headroom.

### Security posture

Unchanged: this crate is **not audited**. See `SECURITY.md`. The BUFF fix closes
the constant-`q00` weak-key class; Dao does not claim full BUFF restoration, and
the S-CEO/S-DEO games remain open in the literature.

## 0.1.0 - 2026-04-20

Initial public release: pure-Rust HAWK-512 keygen, sign, and verify, validated
byte-for-byte against the reference C via an FFI cross-check harness.
