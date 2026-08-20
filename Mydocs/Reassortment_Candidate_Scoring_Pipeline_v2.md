# Candidate Reassortment Scoring Pipeline — Revised Plan (v2)

Revised after biological review against Reassortment_Report.md, Orthohantavirus_01.md, Liu et al. 2024, and Tscherne et al. 2025. Changes from v1: Gate 0 junk/novel split, Gate 1 database-bias mitigation, Gate 2 full redesign (species-matched, no veto), Gate 3 missing-data handling, Gate 4 codon alignment + recombination control, corrected verdict rubric.

## 1. Input

- A candidate is three full genomic FASTA sequences (S, M, L) — matching what `ortho_ml_sequence_table.csv` stores and what a real new isolate download looks like. No parent labels required; the tool infers them.
- Batch mode: a manifest CSV (`candidate_id, s_fasta, m_fasta, l_fasta`) so multiple future samples run in one pass, output as one CSV sibling to `ortho_reassortment_table.csv`.

## 2. Reference loading

- Runs the existing ingest once (`resolve_inputs` → `build_isolate_records` → `parse_cds_fasta`) to get the real 9,852-isolate reference set with both genomic and CDS sequences available — the same objects `all_hanta.py` already builds. Cached so batch mode doesn't re-ingest per candidate.
- **New:** for any metadata-driven lookup (Gate 3), exclude reference records whose host/isolation-source fields indicate laboratory material (Vero, cell culture, passage, vaccine) so a putative parent's "host" never reads as a primate cell line.

## 3. Gate 0 — Sanity / junk screen

Per segment:

- Length within `SEGMENT_LENGTH_RANGE` (already in `all_hanta.py`).
- A plausible ORF at roughly the expected protein length (derived empirically from reference s/m/l_protein_aa lengths, not hardcoded).
- Terminal panhandle motif present — derived empirically as the consensus terminal ~11 nt across reference isolates for that segment. **Check the consensus region (positions ~2–12), not the literal terminal base**, because archived "complete" sequences can carry 1–2 artificially added terminal bases (Liu et al. 2024, Table 1 footnote). A single terminal-base mismatch logs a warning, not a failure.
- Minimum k-mer identity floor to anything in the reference set.

**Split outcome (new):**

- Fails length / ORF / terminals → `junk_sequence` with the specific reason — no tree-building wasted.
- Passes structure but fails the k-mer floor → `divergent_novel_candidate`, not junk. Genuinely novel species (new shrew/mole-borne lineages are described regularly) look exactly like this; route to Gate 4 with an outgroup-aware comparator set instead of discarding.

## 4. Gate 1 — Fast genetic-compatibility pre-screen

- Per segment, k-mer nearest-neighbor against the reference set identifies the best-matching species per segment (existing `kmer_profile` / `cosine_distance_matrix`).
- **Database-bias mitigation (new):** report the **top-3 species hits with distances**, not just the winner. K-mer nearest neighbor tracks reference sample size and diversity — the artifact that killed the original composition classifier. When top hits are close, all of them feed Gate 4.
- If S/M/L all point to the same species → "no apparent reassortment" fast path.
- If one segment's nearest species differs from the other two → preliminary M-only vs. S/L-involving pattern, plus a parent-divergence-magnitude check:
  - Too divergent (no phylogenetic bridge) → flag as more consistent with artifact/independent evolution.
  - **Too similar → label `undetectable_by_design`, not suspicious.** Reassortment between near-identical parents is invisible to phylogenetics; that is a detection limit, not evidence.
- This gate can never produce a "supported" verdict. It only decides which species get pulled in as Gate 4 comparators.

## 5. Gate 2 — Structural / functional compatibility (redesigned)

No single-species reference. Two tracks plus an informational layer. **Gate 2 has no veto power** — it annotates; it never kills a candidate (per Reassortment_Report §7: definitive non-viability calls need real structural modeling, which is out of scope).

### Track A — Universal structural integrity (all candidates)

Mapped to the candidate's own sequence via alignment, using PUUV Sotkamo (MN832782.1–.4) only as the **structural numbering donor** for Gn/Gc features, because its solved structures make it the numbering standard:

- **Panhandle complementarity** of the moved segment vs. the recipient backbone's own termini (NP recognition is sequence-specific; a mismatch predicts failed encapsidation).
- **WAASA cleavage motif** present and correctly positioned in M.
- **TMIS fusion residues** (W766, Y745, F900 in PUUV numbering) mapped through alignment — never asserted as fixed integer positions.
- **pH-switch residues** (E757/D759 equivalents).
- **N-linked glycosylation** count/position vs. alignment (loss of conserved glycans impairs spike trafficking; novel sites are rare).
- **ORF continuity** of the moved segment (no premature stops at reassortment junctions).

### Track B — Lineage-specific signature check (species-matched)

- Only if Gate 1 implicated a specific parental species: align candidate NP/GP to the **reference of that species** (e.g., HTNV 76-118 if HTNV implicated, PUUV Sotkamo if PUUV implicated).
- Check the documented lineage-defining residues **for that species** (e.g., L179 NP, V71/G700/D967 GP are HTNV Lineage-J markers from Liu et al. 2024 — they apply to HTNV candidates only, never as a cross-species test).
- Output: `lineage_marker_match` (bool), `matched_species`, `marker_set_used`.
- If no signature set exists for the implicated species → `not_informative`. Never a fail.

### Track C — Interaction compatibility (informational)

- **NP–RdRp:** if S and L come from different parents, flag `sl_mismatch_high_risk` (experimental systems show this is usually unstable — Rodriguez SNV×BCCV, McElroy SAS-11, Handke PHPUV; but natural exceptions exist, e.g., Park et al. 2025 HTNV L-segment reassortants, so flag, don't fail).
- **Receptor class:** if M moved, note whether donor and recipient species share a receptor class (β3-integrins vs. PCDH1) per literature. A mismatch flags possible host-range change, not non-viability.

### Gate 2 outputs (annotation columns only)

`panhandle_match`, `waasa_intact`, `tmis_conserved`, `ph_switch_conserved`, `glycan_count_delta`, `orf_intact`, `lineage_marker_match`, `sl_mismatch_high_risk`, `receptor_class_note`, `structural_red_flags` (count), `viability_concern` (bool).

## 6. Gate 3 — Ecological opportunity

- Pulls host_order / country / region / year for the nearest-neighbor reference isolates behind each implicated segment's putative parent, from `ortho_isolate_table.csv`-equivalent fields (lab-host records excluded per §2).
- Checks host / geographic / temporal overlap as in the report's worked examples (Bruges/Nova; New Mexico Peromyscus).
- **Missing data is not failure (new):** if metadata for a parent is absent → `ecology_unknown`. Reserve `no_ecological_route` for **positive evidence of non-overlap** (different continents, non-sympatric hosts, implausibly separated in time).
- Ecology remains necessary-but-not-sufficient: a Gate 3 pass never upgrades a verdict by itself (Bruges and Nova coinfect the same moles and still don't reassort).

## 7. Gate 4 — Real phylogenetic evidence

- Comparator set = union of reference isolates from **every species Gate 1 implicated** (top-3 hits included when distances are close) — a pure M-exchange candidate needs both the S/L-matching species' and the M-matching species' isolates in the same tree to be visible.
- Inject the candidate as one more `PhyloCandidate` with a synthetic minimal `IsolateRecord`; call `_analyze_group()` directly — real MAFFT alignment, real IQ-TREE ML trees, real bootstrap-gated `classify_reassortment_call`, fixed `--phylo-seed` so a candidate scores identically on rerun.
- **Codon-aware alignment for divergent pairs (new):** within-species nucleotide MAFFT is fine; across two divergent species, use amino-acid-guided codon alignment before IQ-TREE so saturation doesn't destroy the incongruence signal.
- **Min-taxa rule carried over (new):** if an implicated species has too few dataset comparators → `insufficient_comparators` (mirrors `species_below_min_taxa`), never a forced tree of two.
- **Recombination control (new):** lightweight RDP5-style screen (≥2 methods) on the moved segment, so intra-segment recombination isn't misread as whole-segment reassortment (Liu et al. 2024; Park et al. 2025 protocol).
- Output lands in a scratch directory, never touching `output/phylogenetics/`.

## 8. Verdict rubric (corrected)

| Verdict | Trigger |
|---|---|
| `junk_sequence` | Gate 0 structural fail (length / ORF / terminals) |
| `divergent_novel_candidate` | Gate 0 structural pass, k-mer floor fail |
| `no_conflict` | Gates 1/4 find no discordant signal |
| `supported_reassortment_candidate` | Gates 1–4 pass; Gate 4 bootstrap ≥ 95; recombination screen negative |
| `supported_reassortment_candidate_with_structural_warning` | As above, but Gate 2 Track A flags panhandle / WAASA / TMIS disruption |
| `no_ecological_route` | Positive evidence of host / geographic / temporal non-overlap |
| `ecology_unknown` | Gate 3 metadata missing (not a failure) |
| `insufficient_evidence` | Gate 4 below bootstrap threshold |
| `insufficient_comparators` | Too few same-species references for a meaningful tree |

- `dead_end_non_viable` is **retired as a verdict**. Gate 2 findings ride as annotation columns (`structural_red_flags`, `viability_concern`, `lineage_marker_match`), which is the honest reflection of what sequence-only analysis can claim.
- Each candidate's output row carries every gate's raw values plus the verdict — CSV-joinable like the rest of `output/`.

## 9. Out of scope (unchanged)

- No synthetic-candidate generator (you bring the FASTA).
- No real structural/docking modeling for Gate 2 Track A (alignment-mapped motif checks only).
- No changes to `all_hanta.py` itself.
