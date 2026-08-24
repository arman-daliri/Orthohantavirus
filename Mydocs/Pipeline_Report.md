# A Three-Stage Pipeline for Orthohantavirus Data Curation, Phylogenetic Screening, and Reassortment Candidate Scoring

**Repository state this report describes:** branch `03-pipeline-phylogeny-reassortment`, commit `a1cae12` (2026-08-15 19:52:06 +04:00). Stage 1 figures below were independently re-run and re-verified on 2026-08-20 against this exact commit. Stage 2 figures come from the run committed at `7851ffc` (2026-08-14 23:46:55 – 2026-08-15 13:22:04, wall time ~13h35m); the only commit between that run and `a1cae12` (`a1cae12` itself) is documented as a mechanical relocation of test code with no behavioral change ("Same 99 tests, all passing"), so the committed Stage 2 output is treated here as current. Stage 2's own wall-clock cost (mostly IQ-TREE's ModelFinder and ultrafast-bootstrap steps on trees up to 755 taxa) made an independent re-run impractical within the scope of this report; every Stage 2/3 number below is cited to its source file rather than re-derived.

---

## Abstract

This report describes a three-stage bioinformatics pipeline (`all_hanta.py`, `phylogeny.py`, `reassortment_detector.py`) that processes a 16,028-record NCBI Datasets package for the genus *Orthohantavirus* into curated isolate tables, real per-segment phylogenetic trees, and a five-gate scoring framework for candidate reassortant genomes. Stage 1 verifies package integrity, reconstructs 9,867 physical isolates from 16,028 accessions via a five-level metadata cascade, and yields 206 coding-complete, non-laboratory, non-vaccine isolates. Stage 2 replaces an earlier k-mer-based reassortment classifier — withdrawn after it was shown to invert the literature's own species-level reassortment ranking — with real MAFFT alignments and IQ-TREE maximum-likelihood trees, screening 2,995 (isolate, comparator-scope) pairs. None of the 206 curated isolates receive a `supported_reassortment_candidate` call. Outside that strict curated set, however, 7 unique isolates (previously unreported in this project's own findings document) do receive a `supported_reassortment_candidate` call and 6 receive `weak_signal_low_bootstrap_support`, all newly attributable to a phylogroup-level comparator scope added since that document was last written. Stage 3 provides a five-gate scoring rubric (genetic compatibility, structural plausibility, ecological opportunity, and phylogenetic evidence) for scoring a new candidate genome one at a time; it was not run in bulk over the 13 Stage-2-flagged isolates for this report, and doing so is the paper's principal recommended follow-up.

---

## 1. Introduction

Orthohantaviruses carry a tripartite, negative-sense RNA genome: S (nucleocapsid protein), M (glycoprotein precursor, cleaved to Gn/Gc), and L (RNA-dependent RNA polymerase). Reassortment — the exchange of one or more whole genome segments between two co-infecting parental lineages during replication in a shared host cell — is mechanistically and evolutionarily distinct from intra-segment recombination (template switching within one segment), which the pipeline's source literature explicitly treats as requiring separate detection methods (e.g., RDP5-style screening) and which this pipeline does not screen for.

The literature this project draws on establishes several constraints that motivate the pipeline's design (full annotated bibliography in `Orthohantavirus_01.md`):

- **S and L tend to co-segregate.** The nucleocapsid protein must physically interact with the polymerase for replication and transcription, and must recognize the conserved terminal "panhandle" structure of the genomic RNA to encapsidate it; reassortants that break the native S–L pairing are frequently unstable experimentally (Rizvanov et al., 2004; Tscherne et al., 2025).
- **M is comparatively exchange-tolerant.** Documented natural and experimental reassortants disproportionately swap M while S and L stay matched to one parent (Tscherne et al., 2025).
- **The rule is strong but not absolute** — S- or L-involving reassortment is documented in HTNV, SEOV, and PUUV (Park et al., 2025).
- **Ecology gates plausibility independently of genetics** — niche overlap and co-infection susceptibility between the two parental lineages are a separate, additional requirement (Rizvanov et al., 2004).
- **Most reassortants are less fit than either parent**, since reassortment can disrupt co-adapted packaging signals and RNA–protein interactions (McDonald et al., 2016; Lowen, 2018).

This report describes, and verifies against the code and the last recorded pipeline runs, how the three-stage pipeline operationalizes these constraints — from raw file integrity through a per-candidate scoring rubric.

---

## 2. Data and Methods

### 2.1 Data source and integrity verification

The input is an NCBI Datasets virus data package for taxon ID 1980442 (*Orthohantavirus*), containing seven files (`annotation_report.jsonl`, `biosample_report.jsonl`, `cds.fna`, `data_report.jsonl`, `dataset_catalog.json`, `genomic.fna`, `protein.faa`) covering 16,028 accessioned records. Package integrity is verified against NCBI's own published `md5sum.txt` before any record is parsed or transformed — file-level MD5 for each of the seven files, plus content-level verification of the genomic FASTA against NCBI's `nucleotide.sequenceHash`/CRC32 field. Both the 2026-08-20 re-run and the committed run report all 7/7 files verified; no run used the `--allow-checksum-mismatch` override, so a failed check would have aborted the pipeline rather than being suppressed.

### 2.2 Isolate reconstruction: a five-level metadata cascade

Individual sequence accessions (one per segment) are grouped into physical isolates by `build_group_key`, a five-level cascade evaluated in order, each level attached to a numeric confidence:

| Tier | Key | Confidence |
|---|---|---|
| 1 | species + BioSample accession + isolate/strain label | 1.00 |
| 2 | species + normalized isolate name | 0.95 |
| 3 | species + normalized isolate-or-strain label (isolate name absent) | 0.90 |
| 4 | species + strain + collection date + location | 0.75 |
| 5 | species + strain + location | 0.65 |
| — (fallback) | accession alone | 0.10 |

BioSample linkage (tier 1) is rare in this package — only 99 of 16,028 records (0.6%) carry a BioSample accession — so tier 1 alone would shatter the dataset into near-singleton groups; tier 2, keyed on the submitter-supplied isolate name (present on 15,528 records, 96.9%), does the great majority of the grouping work. Empirically, in the current 16,028-record package, only three of these six levels are ever reached: tier 1 resolves 99 records (0.6%), tier 2 resolves 15,429 records (96.3%), and the remaining 500 records (3.1%) fall all the way through to the accession-only fallback. Tiers 3–5 exist in the code to handle records with a strain-and-location combination but no isolate name, and are not triggered by this particular package's metadata completeness pattern. This grouping reconstructs 9,867 isolates from the 16,028 input accessions.

### 2.3 Completeness classification and the curated wild-isolate set

Each accession's segment is independently assessed for coding-completeness (NCBI's own `complete` flag, qualified — not overridden — by a minimum plausible length per segment). An isolate has a `coding_complete_triplet` when its reconstructed S, M, and L members are all independently coding-complete; 223 of the 9,867 reconstructed isolates meet this bar.

`build_ml_sequence_rows` (Stage 1) applies three further isolate-level filters to produce the sequence-level ML table (`ortho_ml_sequence_table.csv`): the isolate must have a coding-complete triplet; it must not be laboratory-host-derived (`isLabHost` true on any of its three member accessions, or an explicit lab-host record); and it must not be a vaccine strain (`isVaccineStrain` true on any member). Applying this restated, corrected version of the pipeline's own filter chain to the 223 complete triplets yields **206 complete, wild-type isolates** — 17 excluded for lab/vaccine status; the "missing genomic FASTA for one accession despite a complete-triplet call" edge case the code also guards against did not fire in either the 2026-08-14 or the 2026-08-20 run (0 isolates dropped for that reason).

### 2.4 Sequence divergence and the withdrawn k-mer reassortment classifier

Stage 1 additionally computes a cosine-distance k-mer profile per segment protein and performs single-linkage clustering within each species (threshold 0.02) to assign a `kmer_cluster_id`, and — for species with at least 8 assessable isolates (`--min-species-size`, default 8) — a per-segment divergence diagnostic (1,713 isolates assessable project-wide in the current run, 34 flagged with cross-segment taxonomic incongruence).

An earlier version of this pipeline used this same k-mer metric to directly *call* reassortment (`m_segment_reassortment_candidate`, `sl_cosegregation_violation`). That classifier was withdrawn after validation against the literature showed it produced the **inverse** of the published species-level ranking — flagging 72.1% of *Orthohantavirus andesense* (ANDV) isolates, a species the literature reports as reassortment-absent-or-sparse, while flagging only 7.7% of *O. puumalaense* (PUUV) isolates, where the literature reports 19.1–32%. The cause is intrinsic to the metric: cosine distance on amino-acid k-mer composition tracks within-species sample size and diversity, not tree topology, so it cannot distinguish "many closely related isolates" from "a topologically discordant segment." The diagnostic is retained only as an honest input to phylogenetic candidate selection (a `requires_phylogenetic_confirmation` flag, §2.5) — it no longer produces a reassortment call of its own, and `ortho_risk_table.csv` is unconditionally phylo-blind regardless of whether Stage 2 has been run.

### 2.5 Stage 2 — real phylogenetic screening (`phylogeny.py`)

Stage 2 is a standalone script sharing Stage 1's ingestion and isolate-reconstruction code but performing no output-table generation of its own beyond `ortho_reassortment_table.csv` and the tree/alignment files under `output/phylogenetics/` (species scope) and `output/phylogenetics_phylogroup/` (phylogroup scope).

**Candidate selection.** An isolate enters the phylogenetics candidate pool if it is human-affecting-or-likely, or was flagged `requires_phylogenetic_confirmation` by the (retained) k-mer diagnostic, or has a coding-complete S+M+L triplet — intersected with a hard requirement that CDS sequence is available for all three segments. This is deliberately broader than the strict §2.3 curated set: an isolate can enter Stage 2's candidate pool without being coding-complete or without being wild-type, provided it is human-affecting-or-likely.

**Two independent comparator scopes**, controlled by `--phylo-scope` (default `both`): a **species** scope, comparing an isolate only against other members of the same ICTV species; and a **phylogroup** scope, comparing isolates that could not be assigned an ICTV species against others sharing a host-order-derived phylogroup classification. An isolate can receive up to two independent rows in the output table, one per scope, and the two are not required to agree (they answer different questions).

**Alignment and tree-building.** For each qualifying group (≥4 candidate taxa, `--phylo-min-taxa`), each segment's CDS nucleotide sequences are aligned with `mafft --auto` and a maximum-likelihood tree is built with `iqtree -m MFP -B <replicates> --seqtype DNA` (automatic substitution-model selection plus ultrafast bootstrap, 1,000 replicates by default, fixed random seed 42 for reproducibility). Trees are unrooted unless `--phylo-outgroup` is supplied — the committed run did not supply one, so the resulting calls can show *that* segments from different lineages cluster together, but not *which* lineage was donor and which recipient.

**Per-isolate call.** For each candidate, patristic (branch-length) nearest-neighbour and rank-separation diagnostics are computed across the three segment trees. The discordant segment is whichever tree disagrees with the other two on nearest neighbour. A call of `supported_reassortment_candidate` requires both (a) a discordant-segment pattern and (b) the relevant clade's bootstrap support at or above `--min-bootstrap-support` (95% by default); below that threshold the call is `weak_signal_low_bootstrap_support` instead — explicitly not treated as a finding. A per-species/per-group Robinson–Foulds distance between each pair of segment trees (via `dendropy`) additionally summarizes overall tree disagreement independent of any single isolate's call. Every candidate isolate receives a row regardless of outcome (`phylo_status` of `analyzed`, `failed`, or `species_below_min_taxa`) — nothing is silently dropped.

### 2.6 Stage 3 — candidate reassortment scoring (`reassortment_detector.py`)

Stage 3 is not a bulk-dataset stage; it scores one externally supplied candidate genome (or a manifest of several) against the reconstructed reference dataset through five gates, combined into a single verdict from a closed, ten-entry rubric (`docs/superpowers/specs/2026-08-12-reassortment-detector-design.md`, §10):

| Verdict | Trigger |
|---|---|
| `junk_sequence` | Gate 0 structural check fails (length / ORF / terminal-motif) |
| `divergent_novel_candidate` | Gate 0 structure passes, k-mer identity floor fails |
| `no_conflict` | Gates 1 and 4 find no discordant signal |
| `supported_reassortment_candidate` | Gates 1–4 pass; Gate 4 bootstrap ≥ threshold |
| `supported_reassortment_candidate_with_structural_warning` | As above, with a Gate 2 structural red flag |
| `no_ecological_route` | Gate 3 finds positive evidence of host/geographic/temporal non-overlap |
| `ecology_unknown` | Gate 3 metadata missing for one or both putative parents |
| `insufficient_evidence` | Gate 4 ran but bootstrap support is below threshold |
| `insufficient_comparators` | An implicated species has too few dataset comparators for a tree |
| `undetectable_by_design` | Gate 1's implicated parent species are too similar for any tree-based method to resolve |

Briefly: **Gate 0** checks genomic length, longest-ORF plausibility, and (as an informational, non-failing signal — only ~22% of this dataset's own real S-segment records carry a recognizable terminal panhandle motif) terminal panhandle conservation, plus a k-mer identity floor against the reference set. **Gate 1** finds each segment's nearest reference species by k-mer distance and flags the discordant segment, additionally distinguishing "too divergent to place" from "undetectable by design" (top-2 candidate parent species too close to distinguish). **Gate 2** has three tracks — Track A (universal structural markers: TM/IS motifs, pH-switch residues, WAASA site), Track B (species-specific lineage marker residues), Track C (informational S/L-mismatch and receptor-class flags) — and, per the project's own design decision, none of Gate 2's findings carry independent veto power over the verdict; they surface only as a `_with_structural_warning` qualifier. **Gate 3** checks ecological co-occurrence (host, geography, time) of the implicated parent species in the reference isolate table. **Gate 4** reuses Stage 2's exact tree-building and call logic (`phylogeny.py`'s `_analyze_group`) to place the candidate against real reference comparators.

### 2.7 Reproducibility

Stage 1 writes every output table atomically and emits a SHA-256 manifest (`ortho_output_checksums.txt`, verifiable with `sha256sum -c`) alongside the run. IQ-TREE's random seed is fixed (`--phylo-seed`, default 42) across Stages 2 and 3 for tree reproducibility. `output/` is git-tracked by project convention; the numbers cited in §3 are tied to the specific commits named at the top of this report, not to "the pipeline" in the abstract.

---

## 3. Results

### 3.1 Stage 1 attrition

Independently re-run on 2026-08-20 against commit `a1cae12` (identical to the committed 2026-08-14 run in every figure below):

| Step | Count |
|---|---|
| Input accessions | 16,028 |
| ICTV species resolved | 39 |
| Reconstructed isolates | 9,867 |
| Coding-complete S+M+L triplets | 223 |
| **Curated wild-type, complete-triplet isolates (`ortho_ml_sequence_table.csv`)** | **206** |
| Divergence-assessable isolates (≥8 per species) | 1,713 |
| — with cross-segment taxonomic incongruence | 34 |

### 3.2 Stage 2 phylogenetic screening

The committed run (§2.5) considered 28 species-level groups (14 met the ≥4-taxa threshold and were analyzed; 14 excluded, each with 1–3 candidates) and 5 phylogroup-level groups (4 met the threshold; `chiroptera_borne` excluded with 3 candidates), for a combined 2,995 (isolate, scope) rows — 1,470 from species scope, 1,525 from phylogroup scope:

| `reassortment_call` | Rows | Share |
|---|---|---|
| `no_conflict` | 2,961 | 98.9% |
| `not_analyzed` (species/group below min-taxa) | 20 | 0.7% |
| `supported_reassortment_candidate` | 8 | 0.27% |
| `weak_signal_low_bootstrap_support` | 6 | 0.20% |

Joining `isolate_key` to `ortho_isolate_table.csv` reduces the 8 `supported_reassortment_candidate` rows to **7 unique isolates** (one isolate, LN06/*O. seoulense*, is flagged in both scopes) and the 6 `weak_signal_low_bootstrap_support` rows to 6 unique isolates, none overlapping the supported set — 13 unique flagged isolates project-wide:

| Isolate | Species | Call | Discordant segment | Matches canonical M-exchange pattern | Bootstrap |
|---|---|---|---|---|---|
| ROKA16-10 | *O. hantanense* | supported | L | no | 100 |
| GAN22/2021 | *O. seoulense* | supported | L | no | 100 |
| Hunan2020C007-45 | *O. seoulense* | supported | L | no | 100 |
| LN06 | *O. seoulense* | supported (both scopes) | L | no | 100 |
| YY27 | *O. seoulense* | supported | L | no | 100 |
| 202003197 | *O. sinnombreense* | supported | L | no | 100 |
| MSB89863 | unclassified *Orthohantavirus* | supported | M | **yes** | 100 |
| Chile-9717869 | *O. andesense* | weak signal | S | no | 73 |
| Xianning200 | *O. dabieshanense* | weak signal | M | yes | 61 |
| 1 gen14: Z56527 | *O. puumalaense* | weak signal | L | no | 92 |
| WFS_ChangWei | *O. wufangense* | weak signal | M | yes | 72 |
| 3877 | unclassified *Orthohantavirus* | weak signal | S | no | 88 |
| KBM15 | unclassified *Orthohantavirus* | weak signal | S | no | 74 |

Cross-referencing against `ortho_isolate_table.csv`'s own `has_complete_triplet` / `is_lab_derived` / `is_vaccine_strain` flags: only 2 of these 13 isolates (Chile-9717869 and WFS_ChangWei, both `weak_signal_low_bootstrap_support`) have a coding-complete triplet and wild-type status — and neither is a `supported_reassortment_candidate`. **None of the 206 curated wild-type isolates (§2.3, §3.1) receive a `supported_reassortment_candidate` call.** All 7 supported-candidate isolates entered Stage 2's broader candidate pool via the human-affecting-or-likely or `requires_phylogenetic_confirmation` routes (§2.5), not via triplet completeness.

This project's own prior findings document (`Reassortment_Report.md`) states that a project-wide `--phylogenetics` run "analyzed 1,447 candidates and returned zero `supported_reassortment_candidate` calls anywhere." That figure predates the phylogroup-scope analysis (1,447 ≈ this report's 1,470-candidate species scope alone) and is now **superseded** by the committed run described here: 6 of the 7 currently-supported isolates are *O. seoulense* / *O. hantanense* / *O. sinnombreense* members of the species-scope analysis that were not flagged in the earlier run, and the seventh (MSB89863) is only reachable through the phylogroup scope that did not exist when that document was written. §4 discusses what these calls do, and do not, establish.

---

## 4. Discussion

**The absence of a reassortment call is the pipeline's central, hard-won result, not a footnote.** The project's own history includes a classifier that *did* produce confident-looking reassortment calls and was withdrawn specifically because those calls were shown to be an artifact of sample size rather than tree topology (§2.4). The pipeline's design since that withdrawal treats "no signal" as the default, correct answer unless a real tree, with real bootstrap support, says otherwise — and `ortho_risk_table.csv` is deliberately phylo-blind so that a user who never runs Stage 2 cannot be shown a reassortment-adjacent score that was never actually tested against a tree.

**A `supported_reassortment_candidate` call from Stage 2 is a topology-and-bootstrap statement, not a completed biological claim.** It means: across real per-segment ML trees with ≥95% ultrafast bootstrap support on the discordant clade, one segment's nearest neighbour disagrees with the other two. It does not, on its own, mean the pipeline's own structural (Gate 2), ecological (Gate 3), or recombination-screen considerations have been checked — those live in Stage 3, and Stage 3 was not run in bulk over the 7 newly-identified supported candidates for this report. Doing so — feeding each of ROKA16-10, GAN22/2021, Hunan2020C007-45, LN06, YY27, 202003197, and MSB89863 through `reassortment_detector.py` as a candidate against its implicated parent species — is this report's principal recommended follow-up, and is the only way to know whether any of the 7 also clears Gate 3's ecological-opportunity check or trips Gate 2's structural-warning qualifier.

**One candidate is more literature-consistent than the other six.** MSB89863's call is M-segment discordant with `matches_canonical_m_exchange = true` — the dominant, comparatively fitness-tolerant natural pattern the literature describes (§1). The other six supported calls are all L-segment discordant, the rarer and (per Rizvanov et al., 2004 and Tscherne et al., 2025) mechanistically harder-to-sustain pattern, which — per this project's own framework document — warrants stronger corroborating evidence before acceptance, not automatic skepticism; L-involving reassortment is documented in real HTNV/SEOV/PUUV surveillance (Park et al., 2025), and 6 of these 7 calls are themselves SEOV/HTNV/SNV-lineage isolates.

**Unrooted trees limit what any of this can say about direction.** Because the committed run supplied no `--phylo-outgroup`, both the species- and phylogroup-scope trees can show that a segment clusters with an unexpected lineage, but cannot say which lineage donated and which received it.

---

## 5. Limitations

- Stage 2's topology-and-bootstrap test is not a formal statistical topology test (e.g. no SH-test/AU-test on constrained topologies); a `supported_reassortment_candidate` call is consistent with reassortment, not a certainty.
- Species represented by only one dataset isolate cannot be tested by any tree-based method regardless of literature depth — 20 species/group rows in §3.2 were excluded on exactly this basis, and this report's earlier statistics (§3.1) show 223 complete triplets against only 39 resolved species, so single-representative species are common in this package.
- Pathogenicity and species-risk reference tables consumed throughout the pipeline are explicitly hand-curated literature lookups (`all_hanta.py`, in-code documentation), not values computed from the sequence data — they are inputs to the pipeline, not outputs of it, and a reviewer should treat them accordingly.
- Gate 2 (Stage 3) is a sequence-marker checklist against curated reference residues; it does not compute receptor-binding or packaging-signal compatibility from structural modeling.
- The 7 newly-identified `supported_reassortment_candidate` isolates (§3.2) have not been run through Stage 3's full gate framework or checked against the species-specific literature the way `Reassortment_Report.md` did for the 12 single-representative species — that cross-check, and a corresponding update to that document's now-superseded "zero project-wide" claim, is outstanding work.
- This report's Stage 2/3 figures were not independently re-run (§ header note); they are cited to the committed output of a run whose wall-clock cost (~13.5 hours, dominated by IQ-TREE on trees up to 755 taxa) was outside this report's scope to reproduce.

---

## References

1. ICTV. Genus: *Orthohantavirus*. International Committee on Taxonomy of Viruses. <https://ictv.global/report/chapter/hantaviridae/hantaviridae/mammantavirinae/orthohantavirus>
2. Chen R-X, Gong H-Y, Wang X, Sun M-H, Ji Y-F, Tan S-M, Chen J-M, Shao J-W, Liao M (2023). Zoonotic *Hantaviridae* with Global Public Health Significance. *Viruses* 15(8):1705. DOI: [10.3390/v15081705](https://doi.org/10.3390/v15081705)
3. Klempa B (2018). Reassortment events in the evolution of hantaviruses. *Virus Genes* 54(5):638–646. DOI: [10.1007/s11262-018-1590-z](https://doi.org/10.1007/s11262-018-1590-z) / [PMC6153690](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6153690/)
4. Laenen L, Vergote V, Calisher CH, Klempa B, Klingström J, Kuhn JH, Maes P (2019). *Hantaviridae*: Current Classification and Future Perspectives. *Viruses* 11(9):788. DOI: [10.3390/v11090788](https://doi.org/10.3390/v11090788)
5. Rivero R, Simons D, Damodaran L, Karegi I, Gurev S, Becker DJ, Warren DL, Müller NF, Rasmussen DA, Seifert SN (2026). Molecular and ecological determinants of effective reassortment in orthohantaviruses. *bioRxiv* preprint. DOI: [10.64898/2026.06.10.731004](https://doi.org/10.64898/2026.06.10.731004)
6. European Centre for Disease Prevention and Control. Hantavirus infections factsheet. <https://www.ecdc.europa.eu/en/infectious-disease-topics/hantavirus-infection/factsheet-orthohantavirus-infections>
7. NCBI. Datasets virus data report schema. <https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/virus/>
8. Rizvanov AA et al. (2004). *Virology*. (Full citation: `Orthohantavirus_01.md`)
9. Tscherne DM et al. (2025). *Frontiers in Immunology*. (Full citation: `Orthohantavirus_01.md`)
10. Park S et al. (2025). *Scientific Reports*. (Full citation: `Orthohantavirus_01.md`)
11. Liu Y et al. (2024). *PLOS Neglected Tropical Diseases*. (Full citation: `Orthohantavirus_01.md`)
12. McDonald SM et al. (2016). Reassortment in segmented RNA viruses: mechanisms and outcomes. *Nature Reviews Microbiology*. (Full citation: `Orthohantavirus_01.md`)
13. Lowen AC (2018). It's in the mix: Reassortment of segmented viral genomes. *PLOS Pathogens*. (Full citation: `Orthohantavirus_01.md`)
14. Ortiz N et al. (2025). Orthohantavirus rodent hosts and genotypes in Southern South America: A narrative review. *PLOS Neglected Tropical Diseases* 19(9): e0013489. DOI: [10.1371/journal.pntd.0013489](https://doi.org/10.1371/journal.pntd.0013489)

*Species-specific literature for single-representative isolates (rockportense, brugesense, fugongense, luxiense, montanoense, carrizalense, delgaditoense, maporalense, tatenalense, boweense, kenkemeense, wufangense) is compiled in `Reassortment_Report.md` §4 and is not repeated here; it was not re-verified for this report and speaks to a different (species-level literature) question than this report's Stage 1/2/3 code-verification scope.*
