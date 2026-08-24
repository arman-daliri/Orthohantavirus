# Orthohantavirus pipeline — working notes for Claude sessions

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
