# Reassortment Candidate Detector — Design Spec

**Date:** 2026-08-12
**Branch:** `05-feature/reassortment-detector`
**Status:** approved design, pending implementation plan

## 1. Purpose

`Reassortment_Report.md` (§6) already defines a four-gate framework for judging whether a
*future* candidate reassortant genome is a plausible natural event or an implausible/junk
combination that could never arise or survive. That framework was written for humans applying
it by hand to the 12 untestable isolates in this dataset. This spec turns it into a runnable
tool: given a brand-new S+M+L genome (real new isolate, or ML-generated), score it against the
same four gates — plus a new Gate 0 the report doesn't need because it always assumed
well-formed input — and return a verdict with full supporting evidence.

Plain-language version: you hand the tool three new genome pieces. It checks — is this even a
real-looking hantavirus genome (Gate 0)? Does the pattern of similarity match how nature
actually swaps segments (Gate 1)? Would the pieces still work together physically (Gate 2,
informational only)? Could the "parent" viruses that contributed each piece have ever actually
met (Gate 3)? And, built as real evolutionary trees, does the statistical evidence hold up
(Gate 4)? Only evidence that survives all of this earns the label "plausible natural
reassortant" — anything else gets a specific, honest reason why not.

## 2. Architecture

New standalone script, `reassortment_detector.py`, with its own CLI entry point. It **imports**
`all_hanta.py` as a module and reuses its existing machinery rather than duplicating it:
`resolve_inputs`, `parse_fasta`, `parse_cds_fasta`, `build_isolate_records`, `kmer_profile`,
`cosine_distance_matrix`, `run_mafft`, `run_iqtree`, `classify_reassortment_call`,
`_analyze_group`, `PhyloCandidate`, `IsolateRecord`, `SEGMENT_LENGTH_RANGE`,
`length_anomaly`, `species_risk_profile`, `PHYLOGROUP`. **No changes to `all_hanta.py` itself.**

Rationale for a separate file: this is a forward-looking, exploratory tool judging genomes that
aren't part of the reproducible NCBI-ingest pipeline's tracked dataset. Keeping it separate
avoids growing the already-4000+-line core pipeline file and keeps the reproducibility
guarantees of `all_hanta.py` (checksums, atomic writes, stable schemas) scoped to what they were
built for.

## 3. Input contract

A candidate is three **full genomic** FASTA sequences (S, M, L) — the same representation
`ortho_ml_sequence_table.csv` stores (`s_sequence`/`m_sequence`/`l_sequence`), and what a real
new NCBI isolate download looks like. No parent labels required; the tool infers them.

CDS is derived internally, not supplied: a longest-ORF extractor per segment, sanity-checked
against the expected protein length (derived empirically from the reference set's
`s/m/l_protein_aa` lengths — not hardcoded), because Gate 4's alignment must match the CDS
boundaries the reference tree-building already uses (`select_primary_cds` reads
NCBI-precomputed CDS; it cannot produce CDS from a genome that isn't in the NCBI package).

**Batch mode:** a manifest CSV (`candidate_id,s_fasta,m_fasta,l_fasta`), one row per candidate.
Output is one CSV, schema-stable like `ortho_reassortment_table.csv`, written next to wherever
`--output` points (default: `<output>/candidate_reassortment_table.csv`), never inside
`output/phylogenetics/`.

## 4. Reference loading

Runs the existing ingest once per invocation (`resolve_inputs` → the same accession-building
loop `run_pipeline` uses → `build_isolate_records` → `parse_cds_fasta`) to get the real
9,852-isolate reference set with both genomic and CDS sequences available. Cached in-process so
batch mode doesn't re-ingest per candidate (~16s one-time cost, not per-candidate).

For any Gate 3 (ecological) lookup, reference isolates flagged `is_lab_derived` or
`is_vaccine_strain` are excluded from consideration as a "wild" parent — a putative parent's
host must never resolve to a cell line.

## 5. Gate 0 — Sanity / junk screen

New relative to the report (which assumes well-formed input). Per segment:

- Length within `SEGMENT_LENGTH_RANGE` (existing constant, genomic-length ranges).
- A plausible ORF found via longest-ORF extraction, at roughly the expected protein length
  (empirical, from reference data).
- Terminal panhandle motif present, checked against the **consensus region at positions ~2–12**
  (not the literal first/last base) — archived "complete" sequences can carry 1–2 artificially
  added terminal bases (Liu et al. 2024; matches `Orthohantavirus_01.md` point 7's note about
  GenBank submission adding missing terminal bases). A single terminal-base mismatch is a
  warning, not a failure. The motif itself is derived empirically as the consensus first/last
  ~11nt across reference isolates for that segment — not asserted from memory. **Checked and
  measured against the real dataset during implementation (Task 11's integration test):** only
  ~22% of this project's own reference genomic S sequences start with the recognizable motif in
  forward orientation — many real NCBI submissions are deposited on the opposite strand (e.g.
  isolate OR915064.1/"vindeln" starts `GGGAGTCTACTACTA`, whose reverse complement,
  `TAGTAGTAGACTCCC`, is a near-exact panhandle match) or simply omit the untranslated terminal
  region entirely. So the check now also tries the reverse complement (mirroring how the ORF
  search already checks both strands), and even after that, a panhandle mismatch is
  **informational only, not a Gate 0 failure** — measured directly against this dataset, treating
  it as a hard failure would have called real reference isolates `junk_sequence`. ORF length and
  overall genomic length remain the load-bearing structural checks.
- Minimum k-mer-identity floor to anything in the reference set (existing `kmer_profile` /
  cosine distance machinery).

**Split outcome:**

- Fails length / ORF structure → verdict `junk_sequence`, with the specific failing check named
  (terminal panhandle mismatch alone no longer triggers this — see above). Short-circuits — no
  tree-building wasted.
- Passes structure but fails the k-mer identity floor → verdict-track `divergent_novel_candidate`,
  **not junk**. A genuinely novel lineage (new shrew/mole-borne species are still being
  described) looks exactly like this. Routes into Gate 4 with the comparator set built from
  whatever partial similarity exists, rather than being discarded.

## 6. Gate 1 — Fast genetic-compatibility pre-screen

Per segment, k-mer nearest-neighbor against the reference set (existing `kmer_profile` /
`cosine_distance_matrix`) identifies the best-matching species. Reports the **top-3 species
hits with distances**, not just the winner — k-mer composition distance tracks reference sample
size and diversity (the exact mechanism that produced the withdrawn classifier's inverted
ranking), so a close top-3 is itself informative and all three feed Gate 4's comparator set, not
just the top hit.

- All three segments' best hit is the same species → fast "no apparent reassortment" path
  (still confirmed/refuted by Gate 4, not decided here).
- One segment's best hit differs from the other two → preliminary M-only-vs-S/L-involving
  pattern. Also computes parent-divergence magnitude:
  - Too divergent, no phylogenetic bridge between putative parents → flagged as more consistent
    with an artifact or independent evolution than a real reassortment event.
  - **Too similar → labeled `undetectable_by_design`, not suspicious.** Reassortment between
    near-identical parents is invisible to any tree-based method; that's a detection limit, not
    a finding either way.

**This gate can never by itself produce a `supported_reassortment_candidate` verdict.** It only
decides which species/comparators Gate 4 pulls in — exactly the discipline that separates this
tool from the withdrawn classifier the report describes in §2.

## 7. Gate 2 — Structural / functional compatibility (annotation only, no veto)

Per `Reassortment_Report.md` §7's own stated limitation, only sequence-derivable checks are in
scope — no real structural/docking modeling. **Gate 2 never fails a candidate outright**; it
annotates the row with informational flags a reader can weigh, consistent with the report's
framing that a definitive non-viability call needs real structural modeling this tool doesn't
do.

### Track A — Universal structural integrity

A small per-species **structural reference table**, not a single fixed reference species.
Researched directly (fetched and cross-checked primary sources) during this design session:

| Species | Segment/region | Accession | Basis |
|---|---|---|---|
| PUUV (Sotkamo) | M (Gc) | MN832783.1 | `solved` — Tscherne et al. 2025 §4.2/§4.4 |
| ANDV (Chile-9717869) | GPC precursor, Gc aa 652–1107 | NP_604472.1 | `solved` — Guardado-Calvo et al., *"The Hantavirus Surface Glycoprotein Lattice and Its Fusion Control Mechanism,"* PMC7572791 |
| HTNV (76-118) | Gn only | P08668.1 (PMC7572791); AIL25321.1 residues 18–371 (an older, separate construct) | `partial` — Gn structure solved, WAASA confirmed, but no citable TMIS/pH-switch residue numbers found for HTNV Gc; those positions are `homology_inferred` from the PUUV/ANDV anchors if used at all, explicitly flagged as unconfirmed |

Key finding from cross-checking two independent papers: PUUV's TMIS/pH-switch numbers
(Tscherne et al. 2025) turned out to be **identical** to ANDV's own independently solved
numbers (Guardado-Calvo et al., PMC7572791) — Trp766/Tyr745/Phe900 (TMIS) and Glu757/Asp759
(pH-switch) appear in both papers at the same positions. This isn't coincidence: Guardado-Calvo
co-authors both papers and the group reuses its own ANDV-established Gc numbering convention for
PUUV, because the two species' Gc ectodomains align with essentially no indels in this region.
So the table has **two independently `solved` anchors that agree**, not one PUUV-only reference
being stretched across the genus — meaningfully stronger grounding than the original v2 draft
assumed, and exactly the reason the check below is alignment-based rather than a hardcoded
integer, since agreement between two solved structures is not a guarantee every candidate's
indel pattern will match.

- **WAASA cleavage motif** — confirmed independently for both PUUV (Tscherne et al. 2025 §4.4)
  and HTNV (2017 *J Virol* PMC5640846: "WAASA recognition sequence during protein
  biosynthesis") — genus-general, not anchor-specific.
- **TMIS fusion residues** — Trp766, Tyr745, Phe900, confirmed independently in both the PUUV
  (Tscherne et al. 2025 §4.2) and ANDV (PMC7572791) papers.
- **pH-switch residues** — Glu757/Asp759, same dual confirmation. PMC7572791 additionally notes
  His953 makes contacts with Asp679/Lys833/Asn955 around the post-fusion trimer's 3-fold axis in
  ANDV specifically (`H953F` is a known trimer-destabilizing mutation) — recorded as an optional
  ANDV-specific extra check, not required for the core Track A pass/fail.
- Sotkamo accessions per Tscherne et al. 2025 Figure 1: MN832782.1 (L), **MN832783.1 (M)**,
  MN832784.1 (S) — corrected from v2's draft, which cited a single accession with a `.1–.4`
  range; there are three accessions, one per segment, and only the M-segment one
  (MN832783.1) is relevant to these Gn/Gc features.

Checked by **alignment, not fixed integer position**: candidate M-segment protein is aligned to
whichever anchor Gate 1's implicated species is closer to (PUUV or ANDV row; HTNV candidates
align to whichever of the two gives the higher-identity alignment, since HTNV's own Gc numbering
for these positions isn't independently confirmed); the motif/residues are checked at the
aligned position, so the check is correct regardless of indels between the candidate and the
reference row. Output per feature: `structural_feature_conserved` (bool) for WAASA / TMIS /
pH-switch, plus an N-linked glycosylation site count delta (loss of a conserved glycan is a
known trafficking-impairment signal; gain of a novel site is rare and separately noted) and ORF
continuity at the segment junctions (no premature stop introduced by the swap).

**No confidence penalty for using the ANDV vs. PUUV anchor** — both are independently `solved`
and agree, so the check runs the same way regardless of which one a given candidate aligns
better to. The only "unavailable" case is when neither anchor aligns usably (candidate's
implicated species/clade is too distant from both PUUV and ANDV for a meaningful alignment):
output `structural_check_not_available`, Track A skipped for that candidate, not treated as a
weak fail.

### Track B — Lineage-specific signature check (species-matched)

Only runs if Gate 1 implicated a specific parental species with a known marker set. Aligns
candidate NP/GP to *that* species' own reference (e.g., HTNV 76-118 if HTNV implicated) and
checks the documented lineage-defining residues for that species specifically — e.g., L179 NP,
V71/G700/D967 GP per `Orthohantavirus_01.md` — never applied cross-species. Output:
`lineage_marker_match` (bool), `matched_species`, `marker_set_used`. No marker set for the
implicated species → `not_informative`, never a fail.

**Implementation note:** `Orthohantavirus_01.md` names these residues but doesn't cite which
strain's numbering they're counted against, and a targeted search during this design session
didn't turn up a citable HTNV NP reference accession for L179 specifically (unlike Track A's
residues, which were run down to primary sources). Track B ships as a working, tested mechanism
with its marker table **empty** rather than guessing a reference accession — every candidate
gets `not_informative` until a sourced entry is added, which is exactly the designed fallback
for "no marker set for the implicated species." Filling in the real HTNV entry is future work,
not silently assumed.

### Track C — Interaction compatibility (informational)

- `sl_mismatch_high_risk` (bool): flagged when Gate 1 implicates different parents for S and L
  — experimental systems mostly show instability here, but documented natural exceptions exist
  (Park et al. 2025 HTNV L-segment reassortants), so this is a flag, never a fail.
- `receptor_class_note`: if M moved, checked against literature more current than
  `Orthohantavirus_01.md`'s "β3-integrins vs. Protocadherin-1" framing, which turns out to be an
  oversimplification — Jangra et al. 2018 (*Nature*, doi:10.1038/s41586-018-0702-1) established
  PCDH1 as required for New World hantavirus (ANDV, SNV) entry specifically, and Gerrard et al.
  2021's genetic-depletion study (eLife, PMID 34232859) found β3/β1 integrin loss had "little or
  no effect" on entry for a broad panel of hantaviruses — challenging the integrin-as-receptor
  assumption rather than confirming an Old-World counterpart to PCDH1. So the check is scoped to
  what's actually established: `new_world_pcdh1_dependent` vs. `old_world_pcdh1_independent`
  (by genus clade, using the existing `PHYLOGROUP`/geography data — no per-species integrin claim
  asserted), noting only whether an M-segment donor/recipient pair crosses that New World/Old
  World line, not a specific receptor-molecule match.

### Gate 2 output columns

`panhandle_match_moved_segment`, `waasa_intact`, `tmis_conserved`, `ph_switch_conserved`,
`glycan_count_delta`, `orf_intact`, `structural_check_not_available`, `lineage_marker_match`,
`matched_species`, `marker_set_used`, `sl_mismatch_high_risk`, `receptor_class_note`,
`structural_red_flags` (count of the above that came back negative), `viability_concern` (bool,
`structural_red_flags > 0`) — informational, never gates the verdict directly.

## 8. Gate 3 — Ecological opportunity

Pulls host_order / country / region / year for the reference isolates behind each implicated
segment's putative parent (lab/vaccine-excluded per §4), checks host/geographic/temporal
overlap the same way the report's own worked examples do (Bruges/Nova coinfection; New Mexico
*Peromyscus* sympatry).

**Missing metadata is not failure.** If either putative parent's metadata is absent or too
sparse to judge → `ecology_unknown`. `no_ecological_route` is reserved for **positive evidence
of non-overlap** — different continents, non-sympatric host orders, or implausibly separated in
time. A Gate 3 pass never upgrades a verdict by itself — ecological opportunity is necessary but
not sufficient (the report's own Bruges/Nova case: real, frequent coinfection, still no
reassortment between those two viruses).

## 9. Gate 4 — Real phylogenetic evidence

The methodologically load-bearing gate, matching how the report itself validates known isolates
(§2) rather than repeating the withdrawn classifier's mistake.

- **Comparator set** = union of reference isolates from every species Gate 1's top-3 implicated
  (not a single species tree) — a pure M-exchange candidate needs both the S/L-matching
  species' isolates *and* the M-matching species' isolates in the same tree to be visible at
  all; a single-species tree structurally cannot show an inter-species donor.
- Candidate is injected as one more `PhyloCandidate`, built with a synthetic minimal
  `IsolateRecord` (only the fields `_analyze_group`/`group_candidates`/`_phylo_row` actually
  read are populated: `isolate_key="CANDIDATE:<id>"`, `species_name` = Gate 1's best-guess
  label, `host_order=""`, `is_lab_derived=False`, `is_vaccine_strain=False`,
  `species_status="classified"`, `has_complete_triplet=True`; everything else empty/default —
  confirmed safe by reading `collect_phylo_candidates`/`group_candidates`/`_analyze_group`
  directly, which only touch those fields).
- Calls `_analyze_group()` directly (bypassing `collect_phylo_candidates`/`group_candidates`,
  since the comparator set here is Gate-1-driven, not the pipeline's own selection criteria) —
  real MAFFT alignment, real IQ-TREE ML tree, real bootstrap-gated `classify_reassortment_call`,
  fixed `--phylo-seed 42` so a candidate scores identically on rerun.
- **Codon-aware alignment for divergent pairs:** within-species nucleotide MAFFT is used when
  the comparator set is single-species; when it spans two-or-more Gate-1-implicated species,
  align at the amino-acid level and back-translate to codon-aligned nucleotides before IQ-TREE,
  so synonymous-site saturation across the deeper divergence doesn't erase the incongruence
  signal.
- **Min-taxa rule carried over:** an implicated species with fewer than `--phylo-min-taxa`
  dataset comparators → `insufficient_comparators` (mirrors the existing
  `species_below_min_taxa`), never a forced tree of two or three.
- **Recombination control:** a lightweight RDP5-style screen (≥2 methods agreeing) on the
  segment Gate 1 flagged as discordant, so intra-segment recombination isn't misread as
  whole-segment reassortment (Liu et al. 2024; Park et al. 2025's recommended protocol, also
  referenced in `Orthohantavirus_01.md` §31 goal).
- Output lands in a scratch directory (`<output>/candidate_phylogenetics/<candidate_id>/`),
  never touching `output/phylogenetics/`.

## 10. Verdict rubric

| Verdict | Trigger |
|---|---|
| `junk_sequence` | Gate 0 structural fail (length / ORF / terminal motif) |
| `divergent_novel_candidate` | Gate 0 structural pass, k-mer identity floor fail |
| `no_conflict` | Gates 1 and 4 find no discordant signal |
| `supported_reassortment_candidate` | Gates 1–4 pass; Gate 4 bootstrap ≥ `--min-bootstrap-support` (default 95); recombination screen negative |
| `supported_reassortment_candidate_with_structural_warning` | As above, but Gate 2 Track A/C flags a red flag (`viability_concern=True`) |
| `no_ecological_route` | Gate 3 finds positive evidence of host/geographic/temporal non-overlap |
| `ecology_unknown` | Gate 3 metadata missing for one or both putative parents (not a failure) |
| `insufficient_evidence` | Gate 4 ran but bootstrap support is below threshold |
| `insufficient_comparators` | An implicated species has too few dataset comparators for a meaningful tree |
| `undetectable_by_design` | Gate 1's implicated parents are too similar for any tree-based method to resolve |

`dead_end_non_viable` (v1's verdict) is retired — Gate 2 has no veto power, per the report's
own §7 limitation; its findings ride as annotation columns (`structural_red_flags`,
`viability_concern`, `lineage_marker_match`) instead of a hard verdict. Every candidate's output
row carries every gate's raw values plus the verdict — CSV-joinable like the rest of `output/`.

## 11. CLI

```
python reassortment_detector.py \
  --reference-input all_hanta \
  --output output \
  --candidate CAND01:cand01_S.fasta,cand01_M.fasta,cand01_L.fasta
# or batch:
python reassortment_detector.py --reference-input all_hanta --output output \
  --candidates-manifest manifest.csv
```

Reuses `all_hanta.py`'s existing flags where they mean the same thing:
`--mafft-bin`/`--iqtree-bin`, `--bootstrap`, `--min-bootstrap-support`, `--phylo-min-taxa`,
`--phylo-seed` (fixed default 42), `--phylo-threads`, `--verbose`.

## 12. Testing strategy

- Unit tests for the new pieces: ORF extraction against known reference protein lengths,
  empirical panhandle-consensus derivation, Gate 0 length/floor checks, Gate 1 top-3 k-mer
  logic, Gate 2 alignment-mapped residue lookup (against the verified Tscherne et al. 2025
  positions), Gate 3 overlap logic (using the report's own Bruges/Nova and New Mexico
  *Peromyscus* cases as fixtures with known expected outcomes), verdict-rubric combination
  logic (table-driven, one test per row in §10).
- Integration test building synthetic test candidates from **real held-out reference
  sequences**: e.g., a true no-conflict candidate (all three segments from one real isolate),
  and an artificial M-swap candidate (S+L from one species' isolate, M from another's) as a
  positive control that Gate 1/4 should flag — mirrors how the report validates the pipeline's
  own real isolates in §3.
- No live MAFFT/IQ-TREE run in unit tests (too slow); integration tests that do exercise
  `_analyze_group` are marked slow/opt-in, matching how `--phylogenetics` is already excluded
  from the default fast pipeline run.

## 13. Out of scope

- No synthetic-candidate generator — you bring the FASTA.
- No real structural/docking modeling for Gate 2 Track A — alignment-mapped motif/residue
  checks only, explicitly informational.
- No changes to `all_hanta.py` itself.
- SEOV and any species beyond PUUV/ANDV/HTNV have no solved-structure anchor found; Track A
  is `structural_check_not_available` for candidates that don't align usably to one of the
  three table rows.

## 14. Open items for implementation — resolved

HTNV's own Gc TMIS/pH-switch numbers were not confirmed by any paper found during research,
so this was verified empirically instead: fetched the three actual GPC proteins (PUUV
QOC68549.1, ANDV NP_604472.1, HTNV P08668.1), anchored each on its own WAASA motif position,
and checked whether the WAASA-relative offset that works for ANDV's raw citation (Y745/E757/
D759/W766/F900) reproduces the same residue identities in the other two:

| Species | Accession | WAASA at (1-based) | Offset vs. ANDV | Y745/E757/D759/W766/F900 match |
|---|---|---|---|---|
| ANDV | NP_604472.1 | 647 | 0 (reference) | 5/5 — exact at cited integers |
| PUUV | QOC68549.1 | 654 | −7 | 5/5 — but only after the −7 shift; **not** at Tscherne et al. 2025's literal integers, confirming that citation reuses ANDV's own numbering convention rather than counting from QOC68549.1 itself |
| HTNV | P08668.1 | 644 | +3 | 4/5 (Y745/E757/D759/W766 match; F900's equivalent position is N, likely a local indel near that specific loop, not a broken anchor) |

This directly validates the alignment-based (not fixed-integer) implementation this spec already
called for, with a concrete counterexample (PUUV) showing why: a naive fixed-position lookup
using Tscherne's literal numbers against QOC68549.1 would have been silently wrong. HTNV's Gc
row is upgraded from "unconfirmed" to `homology_inferred_verified` — this project's own
alignment check, not a cited paper, confirmed 4 of 5 positions directly against primary sequence
data; F900's HTNV equivalent should be treated as a known exception, not extended confidently to
novel candidates without its own check.

The three reference sequences (fetched via NCBI eutils, embedded verbatim in the implementation
plan so no network access is needed to build Track A) are the actual data files Track A aligns
candidates against — see the implementation plan for the full FASTA content.

Remaining open item:

- `select_primary_cds`/ORF-extractor edge cases (ambiguous bases, frameshift-looking candidates)
  need explicit handling: Gate 0 should catch these as an ORF-length failure, not crash.

## 15. References for Gate 2 Track A (verified this session)

- Tscherne, A., Guardado-Calvo, P., Clark, J., Krause, R., & Krammer, F. (2025). Puumala
  orthohantavirus: prevalence, biology, disease, animal models and recent advances in
  therapeutics development and structural biology. *Frontiers in Immunology*, 16.
  doi:10.3389/fimmu.2025.1575112 — WAASA (§4.4), TMIS W766/Y745/F900 and pH-switch E757/D759
  (§4.2), PUUV Sotkamo accessions MN832782–4.1 (L/M/S) (Figure 1).
- Guardado-Calvo, P. et al. *The Hantavirus Surface Glycoprotein Lattice and Its Fusion Control
  Mechanism.* PMC7572791 — ANDV strain Chile-9717869 (NP_604472.1, Gn aa 21–484, Gc aa
  652–1107), same TMIS/pH-switch residue numbers as Tscherne et al. 2025, plus ANDV-specific
  His953/Asp679/Lys833/Asn955 post-fusion trimer contacts; HTNV GnH accession P08668.1, strain
  76-118.
- 2017 *J Virol*, *Structural Transitions of the Conserved and Metastable Hantaviral
  Glycoprotein Envelope*, PMC5640846 — HTNV Gn construct (AIL25321.1, residues 18–371), WAASA
  motif confirmed for HTNV independently of PUUV/ANDV.
- Fusion-loop mutagenesis (Glu677/Asp679/Arg951/His953 abolish ANDV fusion; His953 protonation
  drives pH-triggered Gc homodimer dissociation) — consistent with the PMC7572791 numbering
  above.
