# Reference C provenance

**Upstream:** https://github.com/hawk-sign/dev
**Commit SHA:** 1b9fef52559273fe7b40fe3e22968eaedd3a4c2a
**Date retrieved:** 2026-04-16
**License:** MIT (see `LICENSE` file in this directory)
**Authors:** Thomas Pornin <thomas.pornin@nccgroup.com>, Ludo Pulles <ludo.pulles@cwi.nl>
**Paper:** [eprint 2022/1155](https://eprint.iacr.org/2022/1155) — "Hawk: Module LIP makes Lattice Signatures Fast, Compact and Simple" (Ducas, Postlethwaite, Pulles, van Woerden)

## What we vendored

The `hawk-512/` directory contains the generated Reference_Implementation for
HAWK-512. The upstream `src/` directory contains macro-configurable templates;
running `python3 build.py` (or `make`) generates per-degree NIST packages under
`NIST/Reference_Implementation/hawk512/`. We vendored that generated output —
NOT the upstream `src/` generator.

Vendored files:
- `hawk-512/*.c` and `hawk-512/*.h` — verbatim copy from upstream
  `NIST/Reference_Implementation/hawk512/` (includes NIST API wrappers
  `api.c`, `api.h`, `rng.c`, `rng.h`, and KAT generator `PQCgenKAT_sign.c`).
- `LICENSE` — upstream MIT license, unchanged.

Note on directory naming: the upstream output path for the NIST package is
`NIST/Reference_Implementation/hawk512/` (no hyphen), not
`Reference_Implementation/hawk-512/` as the standard NIST submission layout
would suggest. We normalised to `hawk-512/` in our repo for consistency with
naming conventions in this workspace.

## KAT files

KAT files were not shipped in the upstream repo at this commit. The repo ships
`PQCgenKAT_sign.c` (a C program that generates KAT `.rsp` files when compiled
and run), but the pre-generated `.rsp` files are absent. Task 27 of the
implementation plan will either:

- (a) Regenerate KATs by compiling and running `PQCgenKAT_sign.c` via the
      reference C test harness, or
- (b) Obtain KATs from the original NIST PQC on-ramp submission package.

## Patches applied

See `patches/` directory. Patches are applied at build time by `build.rs`
when the `cross-check-reference-c` feature is enabled. Patches are numbered
and applied in order.

Currently none (Task 3 adds the first HAWK_DEBUG_TRACE patch).

## Upstream disclaimer

From the upstream README, reproduced verbatim:

> This code was written for the on-ramp call for post-quantum signature schemes,
> organized by NIST and currently serves research purposes. There was no
> security review yet so use at own risk.

This disclaimer applies to our port as well.

## How to regenerate (for future maintainers)

```bash
cd /tmp
git clone https://github.com/hawk-sign/dev hawk-dev
cd hawk-dev && git checkout <desired-commit-sha>
python3 build.py   # or: make
# copy NIST/Reference_Implementation/hawk512/ into c-reference/hawk-512/
# copy LICENSE.txt into c-reference/LICENSE
# update this PROVENANCE.md with new SHA and date
# if parameters or encoding changed upstream, regenerate KAT vectors and update tests/kat_vectors/
```
