# Public release policy

**Status of this repository: PRIVATE, and it stays private.** It is the main working
repo — code, data, notes, session backups all live here together.

**Public release is a SEPARATE, CLEAN repository.** It is built by *copying* the files in
the COPY list below into a fresh repo. It is never made by flipping this repo's
visibility, and never by rewriting this repo's history.

## The decision, and why

A push included `docs/claude-sessions-backup-20260812-183855/` (66 session transcript
files, ~16 MB) and two PDFs. Rather than rewrite history and force-push a shared repo,
we keep this repo private and control what becomes public at copy time.

Two things were checked before settling on this:

- **The transcripts contain no credentials.** Scanned for `sk-ant-`, `ghp_`,
  `github_pat_`, `AKIA`, `AIza`, `xox[baprs]-`, private-key headers and JWTs — zero hits.
  No `api_key`, no `Bearer`, no `.env`. A `Brave Passwords.csv` appears only as a
  filename inside an `ls -l` listing; its contents were never read and no credential rows
  exist anywhere in the tree. What the transcripts do expose is one contributor's home
  directory layout (`/home/victorfrankestein/...`) — a reason to keep them private, not
  an emergency.
- **Neither PDF was a licensing problem.** Both are Open Access CC BY (see below). They
  were removed to keep the repo small, not for copyright.

**Consequence to keep in mind:** because history is not rewritten, anything ever
committed here stays in this repo's history. That is acceptable *only* while the repo is
private. It is the reason public release must be a copy, not a visibility change.

---

## COPY to the future public repo

- `reassortment_detector.py`, `all_hanta.py`, `phylogeny.py` — the pipeline
- `validate_reassortment_detector.py` — validation harness
- `rpi_calibration.py` and `data/rpi_calibration.json` — the pinned ΔΔG calibration panel
  and the τ derived from it. **Both, not just the script.** The spec's defence of τ is
  that it comes from a versioned panel shipped in the repo rather than being picked by
  eye, and a reviewer cannot check that without the artifact. It records the fold
  mechanism and ViennaRNA version, the panel's composition, and the pre-registered D1
  decision rule with every input to it. Contains no sequence data — only public
  GenBank-derived accessions and isolate keys.
- `README.md` — including the RPI limitations and validation sections
- `environment.yml`
- `Mydocs/*.md` — reader-facing docs (`MyReport.md`, `Orthohantavirus_01.md`,
  `Pipeline_Report.md`, `Reassortment_Report.md`,
  `Reassortment_Candidate_Scoring_Pipeline_v2.md`,
  `reference-notes-goodfellow-ortiz.md`)
- `Mydocs/TEACHME.html`
- `docs/superpowers/specs/` and `docs/superpowers/plans/` — design specs and plans
  (**review each before copying**; they are internal working documents)
- **Small example data only** — a handful of candidate FASTAs sufficient to run the
  README's example invocation end to end. Not the full package.
- `LICENSE` — *not yet present, must be added before release*
- `CITATION.cff` — *not yet present, must be added before release*
- `PUBLIC_RELEASE.md` — this file, trimmed to the COPY list

## NEVER copy

- `docs/claude-sessions-backup-*/` — session transcripts; exposes a contributor's local
  paths, and is working material, not a research artifact
- `Mydocs/*.pdf` — publisher PDFs. Both current ones are CC BY so redistribution would be
  *permitted*, but cite the DOI instead of shipping bytes; the rule stays simple so a
  future non-CC-BY PDF is not waved through by habit.
- `all_hanta.zip` — 14 MB, duplicates the extracted `all_hanta/` tree
- `all_hanta/` — the full NCBI package, ~100 MB and regenerable with
  `datasets download virus genome taxon 1980442`
- `output/` — every pipeline product (~40 MB+), all regenerable by running the pipeline
- `output/*.log` — carries absolute paths from whichever machine ran it
- `.superpowers/` — SDD ledger; already gitignored
- `suggesting.md`, `docs/My Library.csv`, `docs/*.xlsx` — scratch/working notes
- `__pycache__/`, `*.pyc`
- Anything the pipeline regenerates. **If it can be rebuilt, it does not ship.**

## Before release day

- [ ] Add `LICENSE` (pick one — code and derived data may want different terms)
- [ ] Add `CITATION.cff`
- [ ] Confirm the public README's example invocation runs against the small example data
- [ ] Re-read `docs/superpowers/specs/` and `plans/` for anything internal
- [ ] Confirm no absolute local paths in any copied file
- [ ] Confirm `dendropy`, `mafft`, `iqtree`, `viennarna` pins in `environment.yml` match a
      real, built environment

## When unsure

**Ask.** Do not guess a file into the public list.

---

## Reference PDFs removed from this repo

Both Open Access; re-downloadable from their DOIs. Content extracted to
`Mydocs/reference-notes-goodfellow-ortiz.md`.

| File | Paper | Licence |
|---|---|---|
| `GoodFellow.pdf` | Goodfellow et al. (2024) npj Viruses 2:11, doi:10.1038/s44298-024-00016-6 | CC BY 4.0 |
| `Ortiz2025.pdf` | Ortiz et al. (2025) PLOS Negl Trop Dis 19(9):e0013489, doi:10.1371/journal.pntd.0013489 | CC BY |

---

*Keep this file current as the project grows, so publication day is just copying the
first list.*
