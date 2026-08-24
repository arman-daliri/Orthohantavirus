# Orthohantavirus pipeline — working notes for Claude sessions

## RESUME HERE (last session: 2026-08-25)

**State: everything green and committed.** HEAD `76f68ed`. Suite **169 tests, all pass**
(~14 min). Validation harness **PASS** (sensitivity 18/18, specificity 27/27, blind spot
5/5, adversarial 4/4). `main` and `03-pipeline-phylogeny-reassortment` are identical;
the feature branch can be deleted with `-D` (`-d` refuses only because its upstream is
stale — everything is reachable from `main` and from `origin/main`).

**FIRST ACTION TOMORROW: `git push`** — 1 commit (`76f68ed`) is unpushed.

**Just resolved (do not redo):**
- ViennaRNA **2.7.2 installed from PyPI** — matches the `environment.yml` pin exactly.
  That pin is no longer unverified.
- `_parse_rnacofold_dg`'s regex is **validated** — parses `( -5.30)`, `(-12.30)` and
  FASTA-header forms correctly. Remove the "ASSUMED FORMAT, NOT VALIDATED" caveat in its
  docstring when convenient.

**Open, in priority order:**
1. **M_pack is not usable as scaled** (new finding, blocks enabling it). Real cofold
   energies: cognate panhandle dG ≈ -22.4, non-cognate ≈ -3.6 to 0.0, so realistic
   ddG ≈ 20 kcal/mol. With RT ≈ 0.616, `exp(-ddG/RT)` gives M_pack ≈ **1e-14**. Because
   `rpi_generation` is multiplicative, switching M_pack on as written would annihilate
   the index for every non-cognate pair. Needs a rescaling decision before it can ship.
   Harmless today: M_pack is hardcoded to 1.0.
2. **Gate 4 still unvalidated.** `mafft` and `iqtree` are NOT on PyPI (`pymafft` is an
   unrelated package — do not install it). Needs conda or manual Windows binaries, which
   requires the user's go-ahead. Gate 4 is the only layer that can catch intra-species
   reassortment — the class Rodriguez 1998, Klempa 2018 and Park 2025 all describe, and
   the one Ortiz 2025 shows is dense (nine named genotypes inside *O. andesense*, which
   has 45 complete triplets here).
3. **`dendropy` reconciliation.** pip 5.0.13 in system Python, not the conda env.
   Options discussed: add a `requirements-pip.txt` documenting the system-Python path so
   `environment.yml` is not silently contradicted.
4. Unfixed review findings N4/N5/N6 (all LOW) — see the review notes: nan-ddG behaviour
   flip, `_json_safe` np.float32 hole, and the thin ordering test.

**User asked at sign-off:** pick between (a) `requirements-pip.txt` reconciliation or
(b) the M_pack scaling question. That choice is still open.


## Repository policy

**This repo stays PRIVATE.** It is the single working repo: code, data, docs, notes and
session backups all live here together. Do not propose splitting them out.

**Never rewrite history. Never force-push.** `git filter-repo`, `filter-branch`,
`rebase` onto pushed commits, and `push --force` are all off the table. The remote is
`arman-daliri/Orthohantavirus` — not the primary user's account — so a rewrite would
break other people's clones. If something sensitive lands in history, the answer is
"the repo is private," not a rewrite.

**Public release is a SEPARATE, CLEAN repository**, built by copying an explicit file
list into a fresh repo. Never by flipping this repo's visibility. The list lives in
`PUBLIC_RELEASE.md` — read it before any release work, and keep it current as the
project grows.

**Saving session backups and working notes under `docs/` is fine** while the repo is
private. They are excluded from public release by `PUBLIC_RELEASE.md`, not by deletion.

**If unsure whether a file belongs in the public release, ask.** Do not guess.

## Environment

- Active interpreter is **system Python 3.14**, not the `hanta_phylo` conda env in
  `environment.yml`. `conda`/`mamba`/`micromamba` are not on PATH here.
- **`dendropy` was installed with pip** into that system Python to get the suite green.
  This is a deviation from `environment.yml` and is not yet reconciled.
- **`mafft` and `iqtree` are NOT installed.** Gate 4 (per-segment phylogenetics)
  therefore cannot be exercised for real. Tests substitute fake binaries; on Windows
  those need a `.cmd` shim (`_write_fake_binary`) because `CreateProcess` cannot exec a
  `.py` directly.
- **ViennaRNA/RNAcofold is NOT installed.** `environment.yml` pins `viennarna=2.7.2`,
  never reconciled against a real build. `M_pack` is hardcoded to 1.0 as a result.

## Running things

```bash
python -m unittest reassortment_detector          # 169 tests, ~14 min (gate 4 is slow)
python validate_reassortment_detector.py --json output/validation/validation_report.json
python reassortment_detector.py --rpi-report ...  # writes rpi_validation_report.json
```

## Known limitations to keep in mind

- **Gate 1 cannot see intra-species reassortment.** It assigns each segment to its
  nearest reference *species*; same-species parents all assign to that species, so it
  reports `concordant`. This is the class Rodriguez 1998, Klempa 2018 and Park 2025 all
  describe. Gate 4 is the layer that would catch it, and gate 4 has never been validated
  against real binaries. Do not read "0 false positives" as "detects all reassortment."
- **`rpi_generation` is currently `E × K`**, not the four-factor product the schema
  advertises — `M_pack` and `M_pol` are both scaffolded at 1.0.
- **`rpi_k_stratum` is derived from a protein k-mer proxy**, not parental nucleotide
  divergence. Every row carries `rpi_k_stratum_provisional = true`.
- See the README's "Reassortment Plausibility Index — limitations" section (items 1–10)
  before quoting any number.

## A recurring defect class in this codebase

Three separate reviews have caught the same shape: **an invariant asserted on an
intermediate value while the value that actually ships quietly violates it.**

- `rpi_m_pol_status` reported `"scored"` for a constant 1.0
- the sensitivity sweep reported rank correlations that were arithmetic identities
- `rpi_e_total` rounded to a literal `0.0` while the test asserted on the raw dict

**Assert on the row that ships, not on the function that feeds it.** When adding a test
for a guarantee, check it against the CSV/JSON artifact, not the intermediate return.
Also add a guard-the-guard assertion where a fixture could drift into making the test
trivially true.
