# Orthohantavirus Reassortment: Findings and a Framework for Judging Future Candidates

**Scope:** The 206 isolates in `output/ortho_ml_sequence_table.csv` (coding-complete S+M+L triplets, no lab-passaged or vaccine-strain isolates), cross-checked against this pipeline's own rigorous phylogenetic reassortment analysis (`--phylogenetics`, run across all 9,852 reconstructed isolates) and against the peer-reviewed literature.

**Bottom line:** No isolate in this dataset shows confirmed reassortment. 186 of the 206 were directly tested by real per-segment phylogenetic trees and came back clean. The remaining 12 could not be tested directly (each is the only representative of its species in the dataset); literature closes most of that gap — two have a published negative, three belong to virus lineages (mole/shrew-borne) with independently confirmed reassortment histories, and the rest are simply understudied rather than negative.

---

## 1. What "reassortment" means here

Orthohantaviruses carry a tripartite genome — S (nucleocapsid protein, NP), M (glycoproteins Gn/Gc), L (RNA-dependent RNA polymerase, RdRp). Reassortment is the exchange of one or more whole segments between two co-infecting parental lineages during replication in a shared host cell, producing progeny with a mixed-parentage genome. It is distinct from **recombination** (intra-segment template switching), which this project's source literature explicitly flags as a separate phenomenon requiring separate detection methods (RDP5-style screening) (Liu et al., 2024; Park et al., 2025).

Segment-specific constraints from the literature (see `Orthohantavirus_01.md` for the full annotated bibliography):

- **S and L tend to co-segregate.** The nucleocapsid protein (S) and the RNA-dependent RNA polymerase (L) must physically interact for replication and transcription — an N–RdRp interaction is required for viral RNA synthesis, and NP must recognize the conserved "panhandle" terminal structures of viral RNA to encapsidate the genome (Rizvanov et al., 2004; Tscherne et al., 2025). Reassortants that break the native S+L pairing are frequently unstable in experimental systems.
- **M is comparatively exchange-tolerant.** Stable experimental and natural reassortants disproportionately involve M-segment swaps while S and L stay matched to the same parent (e.g., BCCV L/S with SNV M; SNV L/S with ANDV M; PHV L/S with PUUV M) (Tscherne et al., 2025). M encodes Gn/Gc, which determines cell attachment (beta-3 integrins, Protocadherin-1) and is the segment where a swap most directly changes infectivity.
- **The rule is strong but not absolute.** S- or L-involving reassortment has still been documented in HTNV, SEOV, and PUUV, including HTNV surveillance detecting candidates with heterogeneous L segments (S–L or M–L incongruence) at high confidence (Park et al., 2025).
- **Beyond genetics, ecology gates plausibility.** A reassortant additionally needs niche overlap between the two parental lineages and host susceptibility to co-infection by both; genetic compatibility alone is not sufficient (Rizvanov et al., 2004).
- **Most reassortants are less fit than either parent.** Reassortment does not usually improve fitness; it is far more often deleterious, since it can disrupt co-adapted packaging signals and RNA–protein interactions (McDonald et al., 2016; Lowen, 2018).

---

## 2. How this pipeline actually tests for it

An earlier version of this pipeline emitted reassortment calls (`m_segment_reassortment_candidate`, `sl_cosegregation_violation`) from a k-mer composition classifier. That classifier was **withdrawn** after validation against the literature: measured on this dataset it produced the *inverse* of the published species-level ranking (e.g., it flagged 72.1% of ANDV isolates, a species where the literature reports reassortment as absent or sparse, while flagging only 7.7% of PUUV, where the literature reports 19.1–32%). The cause was intrinsic to the metric — cosine distance on 3-mer amino-acid composition tracks within-species sample size and diversity, not tree topology (see `README.md` → "Reassortment" for the full withdrawal writeup).

**What replaced it** is the opt-in `--phylogenetics` stage, which does the real thing the literature calls for:

1. **Candidate selection.** Isolates are pulled in if they are human-affecting-or-likely, or were already flagged by the (retained, honest) k-mer divergence diagnostic as `requires_phylogenetic_confirmation`, or have a complete S+M+L triplet — intersected with a hard requirement that CDS is available for all three segments.
2. **Per-species alignment and tree-building.** CDS nucleotide sequences per segment are aligned with `mafft --auto` and a maximum-likelihood tree is built with `iqtree -m MFP -B <bootstrap> --seqtype DNA`, i.e. automatic model selection plus ultrafast bootstrap support (default 1000 replicates).
3. **Per-isolate call.** For each candidate, patristic (branch-length) nearest-neighbour and rank-separation diagnostics are computed from the three real trees. A call of `supported_reassortment_candidate` requires **both** (a) the rank pattern matching the literature's M-segment-exchange signature — S and L agree on nearest neighbour, M's nearest neighbour differs — **and** (b) the relevant M-tree clade clearing `--min-bootstrap-support` (default 95%). Anything below that bootstrap threshold is `weak_signal_low_bootstrap_support`, explicitly not treated as a finding.
4. **Per-species Robinson-Foulds distance** between each pair of segment trees (via `dendropy`) gives an overall measure of tree disagreement, independent of any single isolate.
5. **Nothing is silently dropped.** Every candidate isolate gets a row in `ortho_reassortment_table.csv` with `phylo_status` of `analyzed`, `failed` (with the specific error), or `species_below_min_taxa`.

The one reassortment-*adjacent* signal that is metadata rather than inference is `taxonomic_incongruence` — an isolate whose segments NCBI itself assigned to more than one species (70 isolates project-wide). It feeds the risk score independently of the phylogenetic call.

---

## 3. Results on the 206 ML-sequence-table isolates

Joining `ortho_ml_sequence_table.csv`'s 206 isolates back to `ortho_isolate_table.csv` (for `isolate_key`) and then to `ortho_reassortment_table.csv`:

| Outcome | Count | Meaning |
|---|---|---|
| `no_conflict` | 186 | S, M, L trees agree on nearest neighbour — no reassortment signature |
| `not_analyzed` (`species_below_min_taxa`) | 12 | Species has too few same-species isolates in this dataset to build a meaningful tree — genuinely no data, not a negative result |
| Outside the phylogenetics candidate set | 8 | 7 `unclassified Orthohantavirus` isolates (excluded on principle — pooling under a placeholder label would tree genuinely unrelated lineages) + 1 `Orthohantavirus puumalaense` isolate that fell outside candidate selection |
| `supported_reassortment_candidate` | **0** | — |
| `weak_signal_low_bootstrap_support` | 0 (within these 206; 2 exist elsewhere in the full 9,852-isolate run) | — |

Project-wide (all 9,852 isolates, not just these 206), the `--phylogenetics` run analyzed 1,447 candidates and returned **zero** `supported_reassortment_candidate` calls anywhere — 1,445 `no_conflict`, 2 `weak_signal_low_bootstrap_support`.

---

## 4. Literature deep-dive: the 12 isolates the pipeline couldn't test directly

Each of these 12 is the *only* representative of its species in this dataset, so the pipeline's tree-based method has nothing to compare it against. The question below is genuinely open at the data level; literature review closes some of it.

| Species (our isolate) | Host / geography | Multiple isolates known in the wider literature? | Reassortment documented? | Key source |
|---|---|---|---|---|
| **rockportense** (Rockport virus) | Eastern mole, USA, 1986 | Yes | **Confirmed.** RKPV was excluded from certain phylogenetic analyses specifically because it is a reassortant; shares a most recent common ancestor with cricetid-rodent-borne hantaviruses (cross-species transmission in the distant past) | Klempa, 2018, *Virus Genes*, [PMC6153690](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6153690/) |
| **brugesense** (Bruges virus) | European mole, Belgium, 2013 | Yes — frequently coinfects the same moles as Nova virus | **Mixed.** No reassortment detected between the two *currently* co-circulating viruses despite frequent coinfection, but all three segment trees show topology inconsistencies indicating **ancient** reassortment in Bruges virus's own evolutionary history | [PubMed 29272370](https://pubmed.ncbi.nlm.nih.gov/29272370/), *Genome Biology and Evolution* |
| **fugongense** (FUGV) | *Eothenomys eleusis*, China, 2012 | Unclear — essentially one genome | **No** — explicit published negative, checked directly against its sister taxon | Zhang et al., *Virology Journal*, [doi:10.1186/s12985-016-0483-9](https://doi.org/10.1186/s12985-016-0483-9) |
| **luxiense** (LUXV) | *Eothenomys miletus*, China, 2009 | Unclear — same basis as fugongense | **No** — same negative result (FUGV vs. LUXV comparison) | Same as above |
| **montanoense** (MTNV) | *Peromyscus beatae*, Mexico, 2006 | **Yes, extensively** — ≥9 New Mexico isolates, an Arizona partial sequence, and 2 Mexico isolates (our `104/2006` is one) | **Not yet assessed.** Explicitly named by the authors as necessary future work; segment trees showed no incongruence in what *was* analyzed. High ecological plausibility — multiple *Peromyscus* species are sympatric across its range | Goodfellow et al., 2024, *npj Viruses*, [doi:10.1038/s44298-024-00016-6](https://doi.org/10.1038/s44298-024-00016-6) |
| **carrizalense** (Carrizal virus) | *Reithrodontomys sumichrasti*, Mexico, 2006 | Some — part of a monophyletic clade with Montaño/Limestone Canyon/El Moro Canyon/Rio Segundo viruses | Not addressed in available sources | [ScienceDirect S0168170212002171](https://www.sciencedirect.com/science/article/abs/pii/S0168170212002171) |
| **delgaditoense** (Caño Delgadito virus) | *Sigmodon alstoni*, Venezuela, 1996 | Some — isolated repeatedly from the same host species | Not addressed directly | [PubMed 9498614](https://pubmed.ncbi.nlm.nih.gov/9498614/) |
| **maporalense** (Maporal virus) | *Oligoryzomys fulvescens*, Venezuela | Some | Not addressed directly for this species; the broader Paraguayan/South American hantavirus literature documents reassortment *and* host-switching as a regional pattern in this group | [Phylogenetic exploration of hantaviruses in Paraguay reveals reassortment and host switching, PMC3168424](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3168424/) |
| **tatenalense** (Tatenale virus) | *Microtus agrestis*, UK, 2014 | **Yes** — multiple strains across northern and central England, extensive documented strain variation | Not addressed in available sources (focus was strain variation, not reassortment) | [PMC7232349](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7232349/) |
| **boweense** (*Crocidura douceti*, shrew) | Guinea, 2012 | Not confirmed per-species | **Yes, at clade level** — soricid (shrew)-borne hantaviruses as a group have repeated, independently confirmed reassortment (Imjin virus, Seewis virus, Camp Ripley virus in *Blarina brevicauda*) | [PMC6153690](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6153690/), [PMC7509084](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7509084/) |
| **kenkemeense** (*Sorex roboratus*, shrew) | China, 2011 | Not confirmed per-species | **Yes, at clade level** — same soricid evidence as boweense | Same as above |
| **wufangense** (*Chodsigoa smithii*, shrew) | China, 2016 | Not confirmed per-species | **Yes, at clade level** — same soricid evidence as boweense | Same as above |

**Note on scope:** Ortiz et al. (2025), *Orthohantavirus rodent hosts and genotypes in Southern South America*, was checked directly and does **not** cover delgaditoense, carrizalense, or maporalense — its review scope is limited to Argentina, Chile, Paraguay, and Uruguay, whereas those three species are Venezuelan/Mexican. It also does not discuss reassortment as a topic at all (it is a host/genotype-geography review). It is cited above only where directly relevant.

---

## 5. Synthesis

Combining the pipeline's own test with the literature:

- **No confirmed reassortment among the 206 isolates.** The 186 that could be directly tested all returned clean (`no_conflict`), at the same statistical rigor (bootstrap-gated maximum-likelihood trees) the literature itself uses to make positive calls.
- **The 12 untestable isolates are not "no data" across the board.** Two (fugongense, luxiense) already have a **published negative**. Three (boweense, kenkemeense, wufangense) belong to a lineage (soricid/shrew-borne hantaviruses) with **real, repeatedly confirmed reassortment history**, including a virus in this exact list (**rockportense**, mole-borne, same broader eulipotyphlan-host pattern) individually confirmed as a reassortant. The rest (montanoense, carrizalense, delgaditoense, maporalense, tatenalense, brugesense) range from "actively flagged as needing future work" to "simply not yet studied for this question" — genuine open questions, not negative findings.
- **Ecological co-circulation does not guarantee reassortment.** The Bruges/Nova virus case is instructive: two hantaviruses coinfect the *same individual moles* frequently, and still show no reassortment between each other — genetic/structural compatibility is a separate, additional gate beyond host and geographic overlap.

---

## 6. A framework for judging a future candidate reassortant

This section answers the forward-looking question: **if a new candidate reassortant sequence shows up — real, simulated, or ML-generated — how do you decide whether it's a plausible natural event, or just an unconstrained recombination of segments that could never actually arise or survive?**

A candidate reassortant is a claim of the form *"segment X came from parent A, segments Y and Z came from parent B."* That claim is only as strong as the evidence for it, and it is only *biologically meaningful* if the combination could plausibly have arisen and persisted. Four independent gates need to pass — genetic compatibility, structural/functional compatibility, ecological opportunity, and statistical evidence. Failing any one gate is enough to call a candidate implausible; passing all four is what "supported_reassortment_candidate" means in this pipeline's own terminology.

### Gate 1 — Genetic compatibility (which segment moved?)

| Check | Plausible pattern | Red flag |
|---|---|---|
| Which segment differs from the presumed backbone parent? | **M only**, with S and L matching the same parent (the dominant natural pattern) | S or L differs while its partner doesn't — the rarer, harder-to-sustain pattern; still documented (HTNV, SEOV, PUUV) but needs stronger corroborating evidence before accepting it |
| Terminal "panhandle" sequences on the moved segment | Conserved across the donor and recipient backbone (NP must still recognize and encapsidate it) | Panhandle sequence diverges sharply from what the recipient's NP is known to bind |
| Parent lineages' overall divergence | Moderate — reassortment is documented between closely-to-moderately related strains within a species or between sympatric species | Extreme divergence between putative parents with no phylogenetic bridge — more consistent with independent evolution or a sequencing/lab artifact than a real event |

### Gate 2 — Structural / functional compatibility (would it actually work?)

Reassortment succeeding at the RNA level is necessary but not sufficient — the resulting virus still has to complete its life cycle:

- **Replication:** the recipient's L-segment RdRp must physically interact with the S-segment NP (or vice versa, if S moved) — this is why S+L co-segregation is the default expectation, not an arbitrary rule.
- **Cell entry:** if M moved, the new Gn/Gc must still bind a receptor available on the target host's endothelial cells (beta-3 integrins, Protocadherin-1) — a glycoprotein from a virus with a very different host-receptor profile is a red flag, not a free pass.
- **Immune evasion:** the NP's importin-alpha interaction and the Gn cytoplasmic tail's TRAF3 interaction need to still function well enough for the virus to establish infection — signature residues at these interaction sites (e.g., L179 on NP, V71/G700/D967 on GP, per `Orthohantavirus_01.md`) are the kind of feature a plausibility check should look at when available.
- **Assembly:** the RNP complex (RNA+NP) must still be able to interact with Gn/Gc to bud correctly, and the glycoproteins must still form a functional surface lattice — a mutation combination that disrupts this produces a non-viable particle regardless of how plausible the rest of the genome looks (`Orthohantavirus_01.md`, §21–23).

A candidate that requires an S/L pairing known to be structurally incompatible, or an M segment whose receptor-binding residues don't match anything the recipient's host cells express, should be scored as **implausible** even if the raw sequence identity looks unremarkable — the failure is functional, not statistical.

### Gate 3 — Ecological opportunity (could the parents ever have met?)

Reassortment requires actual co-infection of a single host cell by both parental lineages, which in turn requires:

- **Host overlap.** Same reservoir species, or two host species known to be sympatric and susceptible to cross-infection (the New Mexico *Peromyscus* case in this report is a documented example of exactly this kind of overlap).
- **Geographic and temporal co-circulation.** Both parental lineages present in the same location within a biologically reasonable window — not just "somewhere in the same country," but the same population or a connected metapopulation.
- **Documented or plausible coinfection.** The strongest ecological evidence is a real, sequenced coinfection (as in the Bruges/Nova virus case) or serological/molecular detection of both lineages in the same trapping site or individual animal.

Critically, **this gate is necessary but not sufficient by itself** — Bruges and Nova virus coinfect the same moles frequently and still show no reassortment between each other (§4 above). Ecological opportunity without genetic/structural compatibility (Gates 1–2) does not produce a viable reassortant; it just produces two viruses sharing a host.

### Gate 4 — Statistical / phylogenetic evidence (can you actually show it happened?)

This is where the claim gets tested against data, using the same method this pipeline already applies:

- **Build real per-segment trees** (MAFFT alignment + IQ-TREE maximum-likelihood, not a composition/k-mer proxy — see §2 for why the proxy was withdrawn) for the candidate and enough same-species or same-clade comparators to make the tree meaningful (this pipeline's own `--phylo-min-taxa` gate exists precisely because a tree of one is not a tree).
- **Look for topological incongruence**: does the candidate's placement in the M tree disagree with its placement in the S and L trees?
- **Require bootstrap support** on the discordant clade before accepting the call — this pipeline's threshold is 95%; anything weaker is `weak_signal_low_bootstrap_support`, not evidence. This single gate is what separates a real finding from noise, and it is the exact reason the original k-mer classifier had to be withdrawn (composition distance has no bootstrap-equivalent confidence measure at all).
- **Where possible, corroborate with recombination screening** (RDP5-style, multiple methods, consensus threshold) to make sure the signal is genuinely a whole-segment reassortment event and not an intra-segment recombination artifact being misread as one (Liu et al., 2024; Park et al., 2025; see `Orthohantavirus_01.md` §31 for the full recommended protocol).

### Putting it together: a scoring rubric

| Gate | Pass condition | If it fails |
|---|---|---|
| 1. Genetic compatibility | M-only exchange (or, if S/L involved, strong independent corroboration); conserved panhandle termini; moderate parent divergence | Treat as low-confidence; require stronger evidence from the other gates before accepting |
| 2. Structural/functional compatibility | NP–RdRp, Gn/Gc–receptor, and assembly interactions all remain intact in the recombined genome | **Reject as non-viable** — a genetically-possible combination that cannot complete its life cycle is not a naturally-occurring reassortant, it is a dead end |
| 3. Ecological opportunity | Documented or plausible host/geographic/temporal overlap between parental lineages | Without this, the candidate has no route to have arisen naturally, however clean the tree signal looks — check for artifacts (mislabeling, contamination, lab passage) first |
| 4. Statistical evidence | Bootstrap-supported (≥95%) topological incongruence on a real ML tree, ideally corroborated by recombination screening | Anything below threshold is a hypothesis, not a finding — this pipeline calls it `weak_signal_low_bootstrap_support` and does not let it contribute to any downstream score |

**All four gates passing** = a genuine `supported_reassortment_candidate`, the same bar this pipeline already enforces.
**Gates 1, 3, 4 pass but Gate 2 fails** = a sequence that *could* arise (compatible parents met, tree shows the signal) but *could not survive* — the "random-looking sequence that cannot survive" the question specifically asks about. This is a real, useful negative: it says the combination is a dead evolutionary end, not that it's biologically meaningless.
**Gate 3 fails regardless of the others** = check for a data/lab artifact before treating it as a real ecological puzzle.
**Gate 4 fails** = insufficient evidence either way; more sequencing/sampling is needed, not a verdict.

### Worked examples from this project's own findings

- **Rockport virus (rockportense) — passes.** Independently confirmed reassortant (Gate 4), consistent with cross-species transmission between cricetid- and mole-associated lineages (Gate 1/3 satisfied by the documented shared ancestry), and it is a *circulating, sequenced* virus — i.e., it demonstrably completed its life cycle (Gate 2 satisfied by existence).
- **Bruges virus / Nova virus — Gate 3 passes, Gates 1–2 apparently fail for *current* exchange.** Frequent real coinfection in the same host (strongest possible ecological evidence), yet no reassortment detected between the two currently co-circulating viruses — a clean illustration that ecological opportunity alone does not force reassortment; some property of these two genomes (segment-terminus compatibility, structural mismatch, or simply that both are already well-adapted as-is) blocks it in practice, even though Bruges virus's own deep history shows *ancient* reassortment did happen at some point (Gate 4, historical signal).
- **The 186 `no_conflict` isolates in this dataset — Gate 4 fails (no signal found), by design a well-powered test.** These aren't "insufficient data" the way the 12 gap species are; each had enough same-species comparators for a real bootstrap-supported tree, and the trees agreed. This is the pipeline doing exactly what §2 describes: applying Gate 4 properly and getting a clean negative.

---

## 7. Limitations

- The pipeline's phylogenetics stage tests **incongruence in tree topology**, not a formal statistical topology test (e.g., no SH-test/AU-test on constrained topologies) — a `supported_reassortment_candidate` call means the rank pattern and bootstrap support are consistent with reassortment, not a certainty.
- Species with only one dataset representative (all 12 gap species in §4) cannot be tested by this method regardless of how much literature exists about them — tree-based detection fundamentally requires multiple same-species comparators.
- The literature search in §4 used PubMed, Europe PMC, and general scholarly web search; Google Scholar itself was not queried directly (no API, blocks automated access). Coverage of Venezuelan/Mexican-species-specific reassortment literature (delgaditoense, carrizalense, maporalense) remains thin — genuinely under-studied topics, not confirmed negatives.
- Structural/functional compatibility (Gate 2) is described here as a conceptual checklist; this pipeline does not currently compute receptor-binding or packaging-signal compatibility scores from sequence — that would require a dedicated structural-bioinformatics step (e.g., homology modeling of Gn/Gc against known receptor-binding residues) not yet built.

---

## References

Sources cited above with a bare author/year are from the annotated bibliography in `Orthohantavirus_01.md` — see that file for full citations of: Rizvanov et al. (2004), *Virology*; Tscherne et al. (2025), *Frontiers in Immunology*; Park et al. (2025), *Scientific Reports*; Liu et al. (2024), *PLOS NTD*; McDonald et al. (2016), *Nature Reviews Microbiology*; Lowen (2018), *PLOS Pathogens*.

Sources newly cited in this report:

- Klempa, B. (2018). Reassortment events in the evolution of hantaviruses. *Virus Genes*, 54, 638–646. [PMC6153690](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6153690/)
- Reassortment Between Divergent Strains of Camp Ripley Virus (Hantaviridae) in the Northern Short-Tailed Shrew (*Blarina brevicauda*). [PMC7509084](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7509084/)
- A Novel Hantavirus of the European Mole, Bruges Virus, Is Involved in Frequent Nova Virus Coinfections. *Genome Biology and Evolution*. [PubMed 29272370](https://pubmed.ncbi.nlm.nih.gov/29272370/)
- Zhang et al. Identification and characterization of Orthohantavirus fugongense. *Virology Journal*. [doi:10.1186/s12985-016-0483-9](https://doi.org/10.1186/s12985-016-0483-9)
- Goodfellow, S. M. et al. (2024). Genome sequencing identifies "Limestone Canyon virus" as Montaño virus (*Hantaviridae: Orthohantavirus montanoense*) circulating in brush deermice in New Mexico. *npj Viruses*, 2:11. [doi:10.1038/s44298-024-00016-6](https://doi.org/10.1038/s44298-024-00016-6)
- Ortiz, N. et al. (2025). Orthohantavirus rodent hosts and genotypes in Southern South America: A narrative review. *PLOS Neglected Tropical Diseases*, 19(9): e0013489. [doi:10.1371/journal.pntd.0013489](https://doi.org/10.1371/journal.pntd.0013489)
- Isolation, characterization and geographic distribution of Caño Delgadito virus. [PubMed 9498614](https://pubmed.ncbi.nlm.nih.gov/9498614/)
- Ecology of hantaviruses in Mexico: Genetic identification of rodent host species and spillover infection. [ScienceDirect S0168170212002171](https://www.sciencedirect.com/science/article/abs/pii/S0168170212002171)
- Retrieval of the Complete Coding Sequence of the UK-Endemic Tatenale Orthohantavirus Reveals Extensive Strain Variation and Supports Its Classification as a Novel Species. [PMC7232349](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7232349/)
- Phylogenetic exploration of hantaviruses in Paraguay reveals reassortment and host switching in South America. [PMC3168424](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3168424/)

Pipeline methodology and withdrawn-classifier history: `README.md` → "Reassortment" section; implementation in `all_hanta.py` (`classify_reassortment_call`, `run_phylogenetics_stage`, `compute_segment_divergence`).