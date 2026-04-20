# pqe-hawk

Pure-Rust port of the [HAWK](https://hawk-sign.info/) post-quantum signature scheme at the HAWK-512 parameter set (NIST Level I).

## Status

**Testnet-ready, not audited.** This crate is a faithful Rust port of the upstream reference C implementation (github.com/hawk-sign/dev@1b9fef5, MIT licensed). The entire keygen, sign, and verify pipeline has been validated byte-for-byte against the C reference via an FFI cross-check harness.

- **51 FFI cross-check proptests** (256 cases each) validate primitive-level equivalence with C: modular arithmetic, NTT, big-integer operations, fixed-point FFT, Gaussian sampler, NTRU solver, keygen, sign, verify.
- **Self-consistency proptests** verify keygen determinism, key/sig serialization round-trips, and sign/verify round-trips.
- **NIST KAT pin** against the committed `PQCsignKAT_HAWK-512.rsp` from the reference harness: case 0's pk (1024 B) and sk (184 B) match byte-for-byte.
- **Secret-key zeroization** via `zeroize::ZeroizeOnDrop`; constant-time public-key equality via `subtle::ConstantTimeEq`.
- **Fuzz targets** (`cargo-fuzz`): decoder and verify pipeline never panic on arbitrary inputs. Mirrored by stable-Rust proptests and malformed-input tests that corrupt the trailing padding region of valid pubkey/signature encodings.
- **Timing smoke tests** assert that `sign` timings between two keys stay within 20% of each other.

**Not audited.** Do not use to secure live funds, production authentication, or any setting where a compromise has material consequences. See SECURITY.md for the full posture.

## Usage

```rust
use pqe_hawk::HawkKeypair;
use rand::rngs::OsRng;

let kp = HawkKeypair::generate(&mut OsRng);

let msg = b"hello world";
let sig = kp.secret.sign(msg, &mut OsRng).unwrap();

assert!(kp.public.verify(msg, &sig).is_ok());
```

## Parameter set

HAWK-512 only (NIST Level I):

| Parameter          | Value      |
|--------------------|-----------:|
| Ring dimension `n` | 512        |
| Modulus `q`        | 18433      |
| Public key size    | 1024 bytes |
| Secret key size    | 184 bytes  |
| Signature size     | 555 bytes  |
| Salt size          | 24 bytes   |

HAWK-256 and HAWK-1024 are not currently implemented.

## Feature flags

- `default`: pure Rust, no build-time C dependency.
- `cross-check-reference-c`: builds the vendored C reference and runs bindgen to enable FFI byte-diff tests. Requires a C toolchain (cc) at build time. Development-only.

## Deterministic signing

In addition to randomized `sign(msg, rng)`, the crate exposes `sign_deterministic(msg, nonce_seed)` for crash-recovery scenarios:

```rust
let sig = kp.secret.sign_deterministic(b"hello", &[42u8; 32]).unwrap();
// Same inputs always produce the same signature bytes.
```

The "randomness" for deterministic signing is derived from SHAKE-256 over `sk_bytes || nonce_seed || msg`.

## Running tests

```bash
# Self-consistency tests (no C toolchain required)
cargo test

# Full FFI cross-check against reference C (requires cc)
cargo test --features cross-check-reference-c

# Slow / expensive proptests (keygen-heavy fuzz cases, padding corruption,
# statistical timing smoke)
cargo test --features slow-tests
```

## License

Dual-licensed under Apache-2.0 or MIT, at your option. See LICENSE-APACHE and LICENSE-MIT in the repo root.

The vendored C reference in `c-reference/hawk-512/` is MIT-licensed (copy of github.com/hawk-sign/dev@1b9fef5).
