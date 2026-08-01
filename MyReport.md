### FAQ #1:
```bash
Me:
in  Orthohantavirus the output, in ortho_pipeline.log it is written 15504 protein sequences while in NCBI virus data it it is shown 15568 protein sequences why they are different?

My interpretation:
Duplicates are demostrated online but are removed through our processing
```
<details>
<summary><strong>#1.Claude's full response:</summary></strong>
all_hanta/ncbi_dataset/data/protein.faa has 16,950 > header lines, but many of them are the same protein accession repeated back-to-back (e.g. AIR07773.1 appears 13 times, and there are many others like it — 1,446 duplicate occurrences total across the file). This happens because NCBI's Datasets download bundles one protein FASTA entry per CDS annotation, and the same protein accession can be annotated on multiple genome/assembly records.

The pipeline's parse_fasta() (all_hanta.py:800) builds a dict keyed by accession with index.setdefault(acc, record) — so repeat occurrences of an already-seen accession are silently dropped, keeping only the first. The logged count comes from indexed_record_count() (all_hanta.py:717), which counts distinct record objects, i.e. de-duplicated accessions:

16,950 raw FASTA entries − 1,446 duplicate occurrences = 15,504  ✓ matches the log exactly

So 15,504 is the correct count of unique protein accessions in your locally downloaded package — nothing was lost by a bug, the raw file just has redundant entries.

Why 15,568 on the NCBI Virus website differs from 15,504 locally:
That's a separate, smaller gap (64 sequences) between your downloaded package and the live NCBI Virus count — not something your pipeline is doing. The most likely explanations are:
- Timing — NCBI Virus reflects the database at query time; your all_hanta.zip/package is a snapshot from whenever it was downloaded, and new Orthohantavirus proteins get added/suppressed regularly.
- Different query scope — NCBI Virus's "Orthohantavirus" protein count may use a slightly different filter (e.g., including/excluding certain unverified or third-party-annotated sequences) than the NCBI Datasets datasets download virus genome taxon package you have.
</strong></summary>
</details>

### FAQ #2:
```bash
Me: 
is it necessary or recommended to redownlowd and rerun again freshly? if not and everything is up to date, let it be as it is

Short answer:
Not necessary — leave it as is.
```

<details>
<summary><strong>#2.Claude's full response:</strong></summary>
Reasoning:

- The 15,504 vs 16,950 gap is fully explained (duplicate accessions, correctly deduplicated) — that's not staleness, it's correct behavior.
- The remaining 15,504 vs 15,568 gap is ~0.4% and is just the normal drift between any downloaded snapshot and a live, continuously-updated database. NCBI Virus will basically never match a static package exactly — new sequences get submitted daily.
- Re-downloading now just resets the clock; by tomorrow you'd have a new (probably still nonzero) gap. Chasing it isn't productive unless currency is itself the requirement.

When it would be worth re-pulling:
- You're doing time-sensitive surveillance (e.g., tracking a new outbreak/strain) where the newest few sequences matter.
- You're about to publish/report final numbers and want the freshest possible snapshot as of that date.

Otherwise, a fixed snapshot is actually preferable — it makes your pipeline's output reproducible and stable rather than shifting every time you rerun it. I'd say keep the current all_hanta.zip as your working dataset unless you have a specific reason to refresh.
</details>

### FAQ #3:

list the aims of the project, what you know,what is done, what should be studied more, is IQ-TREE and MAFFT neccessary?(why?/why not?)

<details>
<summary><strong>Aims</strong></summary>

Build a reproducible pipeline that turns a raw NCBI Datasets Orthohantavirus package into analysis-ready tables: accession → isolate → species → host hierarchy, pathogenicity/risk labels, per-segment divergence diagnostics, and ML-ready feature matrices for host-association and risk-classification modeling — all traceable back to NCBI's own fields rather than invented heuristics.

</details>

<details>
<summary><strong>What's known / established</strong></summary>
- Taxonomy: virus.taxId is unusable directly (553 values, mixes strain-level and legacy/current taxid pairs). The lineage node directly below genus Orthohantavirus is the real ICTV species — canonicalizing on it collapses 553 → 39 species.
- Isolate grouping: BioSample is present on <1% of records, so grouping runs on a 5-tier fallback cascade keyed mainly on isolate.name (96.9% coverage). Yields 1,795 S/M/L-complete groups, 223 coding-complete triplets (capped by NCBI completeness, not the pipeline).
- Segment/completeness: NCBI's completeness field is authoritative; segment identity resolved via a 5-tier evidence cascade (explicit field → annotation vote → defline → length heuristic → unresolved), leaving only 210 truly unresolved.
- Host roles: derived from host.lineage rank (not string matching); 33 documented host-jump records (pathogenic species found outside their reservoir order).
- Pathogenicity/CFR: curated per-species from literature, with conflicting CFR estimates kept as ranges, not collapsed to a point value. A self-check audit caught and fixed 6 of 38 curated reservoir entries that disagreed with observed data.
- Reassortment classifier was built, tested, and deliberately withdrawn — it inverted the literature ranking (flagged ANDV at 72% vs. PUUV at 7.7%, opposite of reported frequencies) because cosine k-mer distance measures within-species diversity/sample size, not tree topology.

</details>

<details>

<summary><strong>What's done</strong></summary>
The pipeline runs end-to-end and is currently producing output: integrity-verified against NCBI's md5sums, 16,022 accessions → 39 species → 9,865 isolates, all 10 output tables (accession, species, isolate, host, risk, divergence, ML feature matrix, ML segments, suppressed, QC) plus ortho_summary.json provenance. Deterministic, no absolute paths, ~16s runtime.

</details>

<details>
<summary><strong>What should be studied more</strong></summary>
- 210 unresolved segments and ~2,386 unclassified/genus-only records — no species-level risk label currently possible for these.
- 1,719 records missing host, 3,703 missing collection date — limits temporal/spatial epidemiology.
- Grouping tiers 4–5 almost never fire (strain field only populated from BioSample) — isolate reconstruction could improve if strain parsing were extracted from free-text fields instead.
- Real reassortment detection — currently only emits diagnostics (sl_rank_separation, m_rank_offset), explicitly marked requires_phylogenetic_confirmation, not calls.
- True PED / ICTV species demarcation (cutoff 0.1, WAG model) — current k-mer clustering is an alignment-free proxy only, named kmer_* specifically to avoid the false claim of being ped_*.
- surveillance_priority_score validation — it's stated explicitly as a transparent triage heuristic, not a validated epidemiological model; would need outcome data to validate.
</details>

### FAQ #4:
```bash
Should we run IQ-TREE phylogenetics and MAFFT ?
```
<details>
<summary><strong>Yes, to understand the mechanism of reassortment better:</strong></summary>

Not necessary for the pipeline as scoped — and deliberately excluded (README: "No alignment or phylogenetics dependencies"). Everything currently produced (taxonomy, grouping, completeness, host roles, pathogenicity labels, risk score) comes from NCBI's own structured metadata, not from sequence alignment, so adding those dependencies wouldn't improve any of it and would turn a 16-second serial run into a much heavier alignment/tree-building job.

Necessary only if you want to go beyond what's currently claimed, specifically:
1. Actual reassortment calls — detecting reassortment requires testing topological incongruence between independently-built S, M, and L trees. That's a tree-topology question, and k-mer/cosine distance cannot answer it (proven by the failed classifier above) — only real per-segment alignments (MAFFT) and ML phylogenies (IQ-TREE) capture topology.
2. True PED values for ICTV species demarcation — the 0.1 WAG-model cutoff is defined on aligned, ML-inferred pairwise distances, not amino-acid composition distance. Any claim about whether two isolates cross that species boundary needs alignment + ML inference to be defensible.
3. Any publishable phylogenetic claim about a specific isolate/reassortant — the pipeline's own outputs are explicit that they are diagnostics, not confirmations, precisely so nobody downstream mistakes a k-mer proxy for a phylogenetic result.

So: leave the pipeline as-is for data preparation and ML feature generation; bring in MAFFT + IQ-TREE only as a separate, targeted follow-up step on the subset of isolates flagged requires_phylogenetic_confirmation, if/when you need a real reassortment or species-boundary answer.
</details>

### FAQ #5:
```bash
The aims of the pipeline in main branch according to the full pipeline (all_hanta.py, all 2,839 lines), README.md, .gitattributes/.gitignore, the current ortho_pipeline.log, and ortho_summary.json:
```
<details>
<summary><strong>Gemini NoteBook Lm:</strong></summary> 

Aims

- Integrity model (all_hanta.py:841-967): file-level MD5 agwith a CRLF-normalize diagnostic path for Windows autocrlfcorruption; content-level dedup via NCBI's nucleotide.sequenceHash (CRC32, zero-padded to 8 hex chars, normalize_sequence_hash at                    all_hanta.py:970-983).
- Species resolution (resolve_species, all_hanta.py:991-1037): walks virus.lineage to the node directly below genus Orthohantavirus, not virus.taxId — confirmed this collapses 553 raw organism values to a mucently 39 in this run's ortho_summary.json).
- Isolate grouping (build_group_key, all_hanta.py:1352-1386): 5-tier fallback cascade — species+biosample+isolate → species+isolate.name →
species+normalized label → species+strain+date+location → sssion-only. Confirmed BioSample is present on only 87/15,991 records in the documented snapshot.
- Segment resolution (choose_segment, all_hanta.py:1139-116uct majority vote → genomic defline → length-band heuristic→ unresolved, each with a fixed confidence.
- Completeness: NCBI's structured field is authoritative (dnta.py:1174-1194); length only qualifies a COMPLETE call ascoding-complete, never overrides NCBI.
- Host role/order (determine_host_role, all_hanta.py:1326-1rank, not substring matching.
- Pathogenicity table (SPECIES_RISK, all_hanta.py:173-327): 38 hand-curated species entries, every one with a source citation; anything absent is
unknown, never silently "safe." A same-run audit (audit_res:350-382) checks curated reservoirs against the modalobserved host order every invocation — this run found 0 disagreements (6 were fixed historically).
- Reassortment classifier was built, measured against literwn (compute_segment_divergence, all_hanta.py:1864-1988) —the docstring contains the actual failure numbers (ANDV 72.1% flagged vs. literature "absent/sparse"; PUUV 7.7% vs. literature 19.1–32%). Nothing
here is asserted from memory; it's in the code and mirroredods.reassortment.why_withdrawn.
- k-mer clustering is explicitly NOT ICTV PED (assign_kmer_clusters, all_hanta.py:1836-1856, comment block 1751-1761): ICTV's Hantaviridae species
demarcation requires MSA + ML inference under WAG with a 0.ine/3-mer distance is a deliberately-labeled alignment-freeproxy (kmer_* columns, never ped_*).

What's done

The full pipeline runs end-to-end, single-threaded by design (all_hanta.py:2732-2733, GIL-bound, no benefit from threads), producing 10 CSV
tables + 2 JSON artifacts. Current run (output/ortho_summar11:09:39Z):

| Metric                                | Value                                              |
| ------------------------------------- | -------------------------------------------------- |
| Accessions                            | 16,022 (0 failed transforms)                       |
| Canonical species                     | 39                                                 |
| Isolates                              | 9,852–9,865 (log/summary differ slightly by rerun) |
| Coding-complete triplets              | 223                                                |
| All-3-segments-present                | 1,795                                              |
| Records with host                     | 14,272                                             |
| Assessable for divergence diagnostics | 1,706                                              |
| Taxonomic incongruence                | 34                                                 |

Note: these numbers differ slightly from the ones written ions, 222 triplets) — that's just NCBI's live dataset havinggrown between when the README was last regenerated and this run's download, not a pipeline discrepancy.

Also confirmed in code: raw features (no leakage — scaling stats saved separately for in-fold fitting, build_ml_tables, all_hanta.py:2372-2381), unsupervised top-variance k-mer feature selection, group_ide CV, and a ortho_suppressed.csv table that records whyevery non-representative accession was dropped (build_suppressed_rows, all_hanta.py:2487-2525) rather than silently discarding it.    
What should be studied more                                                                                                           
- 210 unresolved segments, ~2,386 unclassified/genus-only records carrying no species-level risk label.                               - 1,719 records missing host, 3,703 missing collection dateidemiology.
- Grouping tiers 4–5 almost never fire — strain is populated almost exclusively via BioSample (0.54% of records), so free-text strain from other fields could recover more triplets.
- Only the first BioSample accession per record is resolved (noted limitation, not yet addressed).                                    - Real reassortment detection — currently emits only continaration, m_rank_offset), every assessable isolate explicitly marked requires_phylogenetic_confirmation, never a call.                                                                              - True ICTV PED / species-boundary confirmation — current cot the WAG-model pairwise distance ICTV actually requires.
- surveillance_priority_score validation — an auditable, re-weightable heuristic (all_hanta.py:2051-2086), explicitly not a validated epidemiological model; no outcome data has been used to cal
- CFR bands are literature ranges, not computed here — HTNV alone carries two conflicting reported ranges (~1% vs. 5–15%) stored side by side rather than resolved.
                                                                                                                                              Is MAFFT + IQ-TREE necessary?
                                                                                                                                              Confirmed by reading the actual withdrawal logic, not just not necessary for anything the pipeline currently does, anddeliberately kept out (no alignment/phylogenetics dependency anywhere in all_hanta.py; --concordance-tolerance/--discordance-threshold flags wremoved along with the classifier, all_hanta.py:2819-2821, g false confidence in a per-run threshold).
                                                                                                                                              Necessary only for the two things the code explicitly says
1. Actual reassortment calls — requires building independent per-segment (S/M/L) trees from aligned sequences and testing topological incongrubetween them. The k-mer/cosine approach was tried, measurederature ranking because it tracks within-speciesdiversity/sample size, not phylogenetic topology — that's a proven negative result in this codebase, not a hypothetical concern.             
 2. True ICTV species-demarcation distances (PED) — the 0.1 ver aligned, ML-inferred pairwise distances; the k-merdistances here are explicitly disclaimed as a proxy for exactly this reason.                                                                  
Everything else the pipeline delivers today — taxonomy, isolate grouping, segment/completeness resolution, host classification, pathogenicity labeling, the risk score, and ML feature tables — comes froata and needs no alignment step. **Add MAFFT + IQ-TREE only as a targeted follow-up over the 1,706 isolates flagged requires_phylogenetic_confirmation, if/when a real reassortment or species-boundary claim is needed.**
</details>

### FAQ #6:
Should we keep lab strains or not?
<details>
<summary><strong>What I understood is; if lab strains are carefully cleaned and validated it helps to understand the virus's feature effectively</strong></summary>
https://pmc.ncbi.nlm.nih.gov/articles/PMC10344421/

1. Laboratory Data is Essential for Ground Truth The provided text emphasizes that clinical laboratories are "data rich" and that up to 70% of data in electronic health records (EHR) comes from the lab 
5. For orthohantavirus (which causes Hemorrhagic Fever with Renal Syndrome or Hantavirus Cardiopulmonary Syndrome), the "gold standard" for diagnosis often relies on laboratory confirmation via RT-PCR, serology (IgM/IgG), or viral sequencing 
4. Excluding these "lab version segments" (i.e., confirmed cases or specific genetic data) would likely remove the most accurate labels needed for supervised learning, leading to a model trained on noisy or subjective clinical impressions rather than objective data 
5.

2. The Risk of Exclusion: Bias and Reduced Generalizability Excluding high-quality, specific data points (like genetic segments confirming a specific orthohantavirus strain) can introduce bias and reduce the model's ability to generalize. The text warns that "insufficient or nonrepresentative training data may lead to discrepant model performance" 
5. If you exclude lab-confirmed segments, your model might struggle to distinguish between orthohantavirus and other febrile illnesses with similar symptoms (like dengue or leptospirosis), resulting in lower sensitivity and specificity.

3. When to Exclude or Adjust Data (Preprocessing) While you should keep the data, the text highlights that data preprocessing is critical to handle issues that might arise from these segments:

Missing Values and Outliers: If "lab version segments" result in missing data for other patients (e.g., a specific PCR test wasn't ordered for everyone), you must assess the "missingness" rate. If the missing rate exceeds a threshold (e.g., 30%), you must investigate if the data is reliable for training or if imputation is necessary 
5.
Normalization: If combining data from multiple institutions, ensure that laboratory test names, units, and methodologies (e.g., different PCR platforms or sequencing depths) are normalized. Differences in assay methodologies can introduce bias if not harmonized 
5.
Label Accuracy: Ensure the "lab version" labels are accurate. The text notes that relying on billing codes (ICD) can be inaccurate; however, direct laboratory results (like RT-PCR or sequencing) are generally more reliable for defining the target variable 
5.
4. Specific Considerations for Orthohantavirus Genotyping Recent studies on orthohantavirus (specifically Hantaan and Seoul viruses) have successfully used machine learning models that incorporate genotyping data (which relies on specific viral segments) to predict distribution dynamics and clinical outcomes 
1. In these studies, the inclusion of genotypic data allowed models to distinguish between different viral strains, which have distinct clinical manifestations and geographic distributions 
1. Excluding this data would prevent the model from learning these critical distinctions.

Recommendation:

Keep the data: Include laboratory-confirmed cases and specific genetic segments as they provide the necessary "ground truth" for training.
Preprocess rigorously: Use techniques like normalization (e.g., z-score, LOINC standardization) to ensure consistency across different lab platforms 
5.
Handle missingness: If some patients lack these specific lab segments, do not simply discard them; instead, analyze the pattern of missingness (systematic vs. random) and consider imputation or using models robust to missing data 
5.
Validate: Ensure your model is evaluated on an independent test set to confirm that the inclusion of these lab segments improves performance rather than causing overfitting 
5.
In summary, excluding lab-confirmed segments would likely degrade the model's accuracy. Instead, treat them as high-value features or labels that require careful cleaning and normalization to maximize their utility.


https://link.springer.com/article/10.1186/s40168-025-02059-0

Based on the text provided, you should not exclude lab-derived segments if your goal is to assess the risk of zoonotic spillover and compare them against wild reservoirs. In fact, the study explicitly argues that including human-derived (which includes clinical/lab-confirmed) sequences is critical for identifying high-risk viruses.

Here is the breakdown based on the study's methodology and findings:

1. The Study's Approach: Direct Comparison The researchers specifically aligned human-derived Orthohantavirus sequences (including Hantaan virus and Seoul virus) against the rodent-derived sequences they collected in the wild.

Finding: They discovered that rodent-derived HTNV and SEOV displayed almost complete amino acid identity (up to 99.95%) with human-derived counterparts.
Conclusion: This high similarity indicated a "significant risk for human outbreaks" because the virus requires "minimal nucleotide changes" to cross the species barrier.
2. Why Excluding Them Would Be Counter-Productive If you exclude human/lab segments:

You lose the baseline for spillover risk: The core of the study's conclusion regarding Orthohantavirus is the genetic similarity between the wild reservoir and the human host. Without the human sequences, you cannot calculate the amino acid identity or identify the specific mutations (SNPs) required for transmission.
You miss the "Convergent Mutation" analysis: The study analyzed Single Nucleotide Polymorphisms (SNPs) to see where mutations occurred in coding regions versus antigenic epitopes. This comparison is only possible by having the human reference.
3. When Might You Exclude Them? You would only exclude lab-adapted strains if your specific research question was:

Pure Evolutionary Ecology: If you are trying to map the natural evolutionary history of the virus in the wild without the "noise" of human selection pressure or lab passage artifacts.
Lab Artifacts: If you suspect the specific lab strains have undergone significant adaptation to cell culture that makes them genetically distinct from the wild-type human strains (e.g., specific mutations in the RdRp that only occur in culture). However, the text implies the human strains used were clinical isolates representing the actual threat.

Recommendation based on the text: For a machine learning model assessing zoonotic potential or spillover risk, keep the human/lab-derived segments. The study demonstrates that the proximity of wild rodent viruses to these human sequences is the strongest predictor of outbreak potential. The text states: "Viruses with high sequence similarity to human-derived viruses may have a reduced capacity to breach species barriers... [but] the high degree of amino acid sequence similarity suggests that rodent-derived Orthohantaviruses may require minimal nucleotide changes to breach species barriers."

Summary:

For Zoonotic Risk Assessment: Keep them. They are the benchmark for "human infectivity."
For Pure Wild Phylogeny: You might exclude them to see the "natural" tree, but you lose the ability to predict which wild branches are dangerous to humans.
The study's success in identifying HTNV and SEOV as high-risk relied entirely on the inclusion of the human-derived data for comparison.
</details>

### FAQ #7:
```bash
to predict if any future reassortments can create new dangerous orthohantavirus for human in orthohantavirus, should we exclude lab version segments or keep them for more valid results?

surveillance priority score, which is a transparent, additive point system ranging from 0.0 to 10.0
```
<details>
<summary><strong>How it is calculated:</strong></summary>
Pathogenicity and Risk Criteria
The score is calculated by summing the weights of the following six terms:
Case-Fatality Rate (CFR) Band: The largest component of the score, based on known literature for a species.
4.0 points: High risk (CFR >30%, e.g., SNV, ANDV).
2.5 points: Moderate risk (CFR 1–15%, e.g., HTNV, DOBV).
1.0 point: Low risk (CFR <1%, e.g., SEOV, PUUV).
0.5 points: Unknown risk (Species with reported human cases but no established CFR).
0.0 points: No known human disease (e.g., shrew-borne species).
Human Spillover (2.0 points): Added if the specific virus record was recovered directly from a human host (Homo sapiens).
Taxonomic Incongruence (1.5 points): Triggered when the genetic segments of a single isolate are assigned by NCBI to two different species, which serves as a reliable metadata signal for potential genetic exchange.
Host Jump (1.0 point): Added if the virus was recovered from an animal order that does not match its species' expected natural reservoir (e.g., a rodent-borne virus found in a bat).
Complete Triplet (1.0 point): Added if the isolate contains all three segments (S, M, and L) and they are all coding-complete, providing a higher quality of data for surveillance.
NSs ORF (0.5 points): A small bonus if the virus encodes the non-structural protein (NSs), which is a known virulence factor in certain species.

</details>

### FAQ #8:
What does each file have?

<details>
<summary><strong>Claude:</strong></summary>
Here's what each output file is, with a real row from your data:

One row per accession/segment (raw level)
- ortho_accession_table.csv — every individual GenBank accession, one row each. Ex: MW498184.1 → PUUV, segment M, partial completeness, host Myodes glareolus, Russia, 2019.
- ortho_ml_segments.csv — a slimmer version of the accession table, trimmed to just the columns useful for ML on individual segments (species, host, year, completeness) — same accession, fewer columns.
- ortho_qc_report.csv — QC pass/fail per accession. Ex: MW498184.1 → transform_ok: True, warnings: '' (empty = clean).
- ortho_suppressed.csv — accessions that got dropped from the main tables and why. Ex: an S-segment record excluded with suppression_reason: missing_host.

One row per isolate (S+M+L segments grouped into an animal/sample)
- ortho_isolate_table.csv — the core "biological unit" table: one row per isolate, with its S/M/L accessions joined together. Ex: isolate puuv kyzyl yulskoe... → S=MW504239.1, M=MW498184.1, L=MW498224.1, triplet_status: all_segments_present_not_all_complete.
- ortho_ml_feature_matrix.csv — the actual ML-ready tabmetadata columns plus 200+ kmer_XXX columns (amino-acid3-mer frequencies from the concatenated protein) as features, and label_pathogenicity_class etc. as targets.
- ortho_segment_divergence.csv — per-isolate k-mer dist possible reassortment candidates (not confirmed — just"worth checking"). Ex: s_kmer_distance, reassortment_status: not_assessable.
- ortho_risk_table.csv — per-isolate epidemiological ricase-fatality band, and a computedsurveillance_priority_score. Ex: an Andes virus isolate scores cfr_band_label: high_gt_30pct, priority 7.5.
- ortho_reassortment_table.csv — the opt-in phylogeneti ran --phylogenetics): actual tree-basedconfirmation/rejection of reassortment per isolate, with bootstrap support. Ex: reassortment_call: no_conflict.

Rollups
- ortho_species_table.csv — one row per virus species, 2567 isolates, top host Myodes glareolus, cfr_band_label:low_lt_1pct.
- ortho_host_table.csv — one row per host organism, agg → 2652 accessions across 1655 isolates, mostly PUUV.

Rule of thumb: accession-level tables = "what NCBI gave us,evel tables = "what animal/sample this actually came from,segments merged"; species/host tables = "summary stats"; risk/divergence/reassortment = "derived scientific judgments on top of the isolate table."

</details>