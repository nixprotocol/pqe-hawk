# pqe-hawk fuzz targets

Coverage-guided fuzz targets for the decode + verify boundary. Run with
[`cargo-fuzz`](https://github.com/rust-fuzz/cargo-fuzz):

```bash
cargo install cargo-fuzz
cd fuzz

# Decode tests — must never panic on arbitrary bytes.
cargo fuzz run fuzz_decode_pubkey
cargo fuzz run fuzz_decode_signature

# Full verify pipeline — must Err gracefully for any input.
cargo fuzz run fuzz_verify
```

Each run will fuzz until you stop it (Ctrl+C). For CI integration, add
a time budget:

```bash
cargo fuzz run fuzz_verify -- -max_total_time=60   # 60 seconds
```

Any panic, crash, or out-of-memory indicates a vulnerability: verify
must always return `Ok(())` or `Err(HawkError::*)` for any input,
without aborting the process. Found crashes are saved under
`fuzz/artifacts/<target>/`.

**Note:** fuzz targets require nightly Rust. The parent crate builds on
stable; these targets are gated in a sub-workspace and only compile
under `cargo fuzz`. They are excluded from the main workspace via the
`[package.metadata] cargo-fuzz = true` marker.
