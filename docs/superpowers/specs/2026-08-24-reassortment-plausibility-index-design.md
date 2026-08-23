# Reassortment Plausibility Index (RPI) — Design Spec

**Date:** 2026-08-24
**Branch:** `03-pipeline-phylogeny-reassortment`
**Status:** Design approved in brainstorming; not yet implemented
**Supersedes:** nothing. **Amends:** nothing — see §10 (Integration).

---

## 0. Limitations, stated up front

This section is deliberately first. Every number this index produces inherits these
constraints, and any Methods paragraph describing it must reproduce them.

1. **The index is relative, not an absolute probability.** `rpi_generation = 0.4` does not
   mean a 40% chance of anything. It means this constellation scores 0.4 on a bounded
   index whose factors are anchored to published measurements. Ranking is meaningful;
   the absolute value is not.
2. **Every in vitro anchor is Vero E6** — an African green monkey kidney line, not
   reservoir rodent cells. Strandin et al. (2020), already in the project library, exists
   precisely because host-cell context changes hantavirus behaviour. Reassortment
   frequencies measured in Vero E6 may not transfer to *Myodes*, *Peromyscus* or
   *Apodemus* cells.
3. **Total observed reassortants across all direct measurements ≈ 35.** Rodriguez et al.
   (1998) contributed 25 + 1; Rizvanov et al. (2004) contributed 10 monoploid. The
   constellation prior (§7) rests on this, and no amount of smoothing manufactures
   statistical power that is not there.
4. **MOI is confounded across studies.** Rizvanov used MOI 3; Rodriguez used MOI 1.
   Rizvanov names this explicitly as the likely reason their cross-species frequency
   (8.9%) exceeded Rodriguez's within-species frequency (8.5%). **Cross-study reassortment
   frequencies are therefore not comparable**, and this spec never compares them. Only
   Rodriguez's internal contrast (both arms MOI 1, shared NMR11 parent) is controlled.
5. **Compatibility is directional; our metrics are largely symmetric.** McDonald et al.
   (2016) report that φ13 packages φ6's M segment while φ6 cannot package φ13's M unless
   φ6 packaging sequences are appended. Sequence- and ΔG-based metrics computed
   symmetrically will lose this asymmetry. Documented, not solved.
6. **No gold-standard label set exists.** The pipeline's own output contains 8
   `supported_reassortment_candidate` calls out of 2,995 rows. Nothing can be fitted to 8
   positives. This is the reason the index is literature-anchored rather than learned,
   and it is also why §11 orders the validation evidence the way it does.
7. **`rpi_retention` models immediate post-reassortment fitness only.** McDonald's Fig. 4
   mode (c) — compensatory mutations restoring interaction interfaces after reassortment —
   is **out of scope** and is a named limitation, not a term in the model. A reassortant
   scored low may still emerge after adaptation.

---

## 1. Problem statement

`reassortment_detector.py` currently answers a categorical question: *does this candidate
show a supported reassortment signal?* — returning one of a closed 10-verdict list from
`combine_verdict` (`reassortment_detector.py:894`).

It does not answer: *how plausible is it that this constellation could arise and persist
in the real world?* That requires a continuous, mechanistically-grounded score.

**Deliverable:** two continuous columns per candidate, in [0, 1]:

- `rpi_generation` — relative plausibility that this segment constellation is physically
  formed in a co-infected cell and replicates.
- `rpi_retention` — relative plausibility that, once formed, it persists well enough to be
  sampled.

The split is not decorative. McDonald et al. (2016) state that failure to recover
reassortants after co-infection "might reflect the poor fitness of hybrid progeny caused
by mismatched alleles, rather than restrictions on the actual generation of the reassortant
during co-infection." The two failure modes are confounded in any single number. The
empirical signature is visible in the data: Rodriguez's tissue culture produced M-odd and
S-odd reassortants in equal numbers (12 vs 12), while Henderson's field data gave 14 M-odd
against 2 S-odd. Same virus, opposite constellation profile — because one measures
generation and the other measures what survived.

## 2. What this is not

- **Not a detector.** Detection remains Stage 2's bootstrap-gated ML topology comparison
  plus Gate 4. RPI scores candidates; it never flags them.
- **Not a virulence predictor.** Klempa (2018) records that the SNV/ANDV reassortant
  carrying ANDV's M segment did *not* confer lethal HPS (McElroy et al. 2004), and that
  Handke's PUUV/PHV reassortant retained the innate-response phenotype of whichever parent
  supplied S and L. A high `rpi_generation` for an M-swap says nothing about pathogenicity.
  Output columns must keep these separate.
- **Not an import of any published model's coefficients.** See §3.

## 3. Evidence base and its handling

Sources are tiered by how directly this spec relies on them.

**Tier 1 — read in full during design, directly load-bearing:**

| Source | Contribution |
|---|---|
| Rodriguez et al. 1998, *Virology* 242:99–106 | Controlled within-study frequencies; constellation counts; divergence table |
| Rizvanov et al. 2004, *Virology* 327:225–232, doi:10.1016/j.virol.2004.07.012 | Distant-pair constellation collapse; MOI confound; divergence figures |
| Klempa 2018, *Virus Genes* 54:638–646, doi:10.1007/s11262-018-1590-z | Hantavirus in vitro synthesis (all recovered reassortants M-exchanged); sympatry ecology |
| Lowen 2018, *PLoS Pathog* 14(9):e1007200, doi:10.1371/journal.ppat.1007200 | Multi-scale sequential filter (Fig. 1); fitness filter (Fig. 2); bunyavirus polymerase–nucleoprotein incompatibility |
| McDonald et al. 2016, *Nat Rev Microbiol* 14:448–460, doi:10.1038/nrmicro.2016.46 | Two-mechanism restriction framework (Fig. 3); rotavirus terminal-motif precedent; fitness modes (Fig. 4); cognate-competition |

**Tier 2 — preprint, used for structural corroboration only, never for coefficients:**

- Rivero et al. 2026, bioRxiv, doi:10.64898/2026.06.10.731004, posted 12 June 2026.
  **Not peer reviewed.** Its sequential-filter architecture is independently established by
  Lowen (2018) Fig. 1, and its S–L coupling finding by Klempa, Rodriguez and Lowen. Its
  *fitted coefficients are not used anywhere in this spec* — partly on preprint grounds,
  and partly because Rivero themselves report that their models transfer poorly across
  species, which disqualifies transplanting them onto a different dataset.

**Tier 3 — cited within Tier 1 sources, not read directly. Must be marked secondhand:**

- Henderson et al. 1995, *Virology* 214:602–610 — natural SNV constellation counts (via Rodriguez 1998).
- Kirsanovs et al. 2010, *Virus Genes* 41:319–328 — DOBV in vitro 65/207 (via Klempa 2018).
- Handke et al. 2010, *J Gen Virol* 91:2351–2359 — PUUV/PHV reassortant phenotype (via Klempa 2018). **PDF is in the project library and should be read before the spec is implemented.**
- Zhang, S. et al. 2014, *J Virol* 88:9842–9863 — cross-segment protein co-adaptation via mutual information (McDonald 2016, ref. 112). Verified against McDonald's reference list; not independently read.
- Kang et al. 2002 — reported 68.19% HTNV×SEOV progeny reassortant (via Rizvanov 2004). Chinese-language source, unverified, and a substantial outlier. **Do not use as an anchor.**

**Tier 4a — existence and abstract verified by the project owner; full text not read:**

- Mir MA & Panganiban AT 2006, *"Hantavirus N Protein Exhibits Genus-Specific Recognition
  of the Viral RNA Panhandle,"* J Virol, doi:10.1128/JVI.00820-06. Not in the project
  library. Supports §6.2's N-recognition footprint term (genus-specific panhandle
  recognition; SNV N requires the first 23 nt vs the terminal 9 for most hantavirus N).
  Read the full text before publication.
- *J Virol* 1995;69(12):8132–8136, PMC189768 — terminal 14 nt of S, M and L conserved
  across the genus. **Content verified; author list NOT verified.** Resolve the author list
  before it appears in any bibliography.
- "Viral reassortment as an information exchange between viral segments," *PNAS* 2012,
  doi:10.1073/pnas.1113300109. **Title and DOI verified; author list NOT verified.** See
  §13.3.

**Tier 4b — supplied but NOT verifiable during design; confirm before use:**

- Lorenz et al. 2011, *Algorithms Mol Biol* 6:26, doi:10.1186/1748-7188-6-26 (ViennaRNA 2.0)
- Bernhart et al. 2006, *Algorithms Mol Biol* 1:3 (RNAcofold)

  Neither is in the project library and web access was unavailable during design. Both are
  recorded from the project owner's memory and are explicitly **unverified**. They gate the
  RT constant in §6.1, so confirm them before that constant is trusted.

**Removed during review:** the "position-9 polymerase binding" claim previously in §6.2.
No primary source; not supported by the Mir & Panganiban abstract, which concerns N protein
rather than the L polymerase. See §6.2.

**Citation hygiene note.** `docs/Orthohantavirus Genomic Workflow and Reassortment Mapping -
NotebookML.csv` must not be used as a bibliography. An audit against the project's own
Zotero library found **16 references carrying a correct title with the wrong authors or
year**, including a fabricated hybrid that attaches McDonald et al.'s authorship to Lowen's
title. It also miscategorises DEmARC as a reassortment tool. Its *workflow structure* and
*marker list* are sound and may be reused; its reference list may not.

## 4. Architecture

```
rpi_generation = E × M_pack × M_pol × K
rpi_retention  = rpi_generation × F
```

**Multiplicative, not additive** — a deliberate departure from the additive house style of
`surveillance_priority` (`all_hanta.py:2182`). Three independent justifications:

1. Lowen (2018) Fig. 1 makes each stage a *necessary* gate: same host → same tissue →
   same cell → co-packaging. Necessary conditions compose by multiplication.
2. McDonald's cognate-competition result is inherently a ratio — non-cognate interactions
   are outcompeted by cognate ones, not summed against them.
3. An additive score would let molecular elegance compensate for ecological impossibility.
   A candidate whose parent lineages have never co-occurred in any host, place or decade
   must not be rescued by a good panhandle match.

Following house style, **every factor is emitted as its own column** so the product can be
audited, disputed, or re-weighted downstream — exactly the rationale documented in
`surveillance_priority`'s docstring.

## 5. Factor E — ecological opportunity

`E = h_host × g_geo × t_time`, each in [0, 1], over **all** implicated species pairs (not
just the first two — see §10.2).

- `h_host` — reservoir-host taxonomic proximity: same host species → same genus → same
  family → same order → different order.
- `g_geo` — geographic co-occurrence, country level, finer where metadata allows.
- `t_time` — overlap of sampling-year intervals.

**Data source and epistemic status — read carefully.**

The *co-occurrence frequencies* (which host taxa, countries and years actually co-occur)
are computed from the project's own 16,028-record dataset via `ortho_isolate_table.csv`.
That part is empirical.

The *decay weights* on the species→genus→family→order ladder **are a structural
assumption.** They cannot be estimated from 8 positives, and this spec does not pretend
otherwise. The functional form of the decay is therefore a first-class subject of the
sensitivity analysis (§12), not a fixed constant presented as measured.

**Supporting evidence for the axis itself** (not for the weights): Klempa (2018) describes
eastern Slovakia as a contact zone for three bank-vole phylogenetic clades — a "melting
pot" producing sympatric virus lineages, co-infection and reassortment. Lowen (2018)
identifies dual-exposure likelihood as a governing factor. Rodriguez's cross-host-species
arm (*P. maniculatus* vs *Sigmodon hispidus*) is the experimental correlate.

**Deliberately excluded: local prevalence.** Rivero report it trending *negative*, arguing
that high prevalence of one dominant reservoir amplifies a single viral background without
increasing exposure to divergent lineages. Lowen lists prevalence as a factor, so the two
sources are in tension. Since this is Tier 2 evidence and the direction is disputed,
prevalence is omitted rather than guessed at. Documented as an open question (§13).

## 6. Factor M — molecular compatibility

Split into two sub-factors because McDonald (2016) Fig. 3 identifies two mechanistically
distinct restrictions: incompatibility of **RNA–RNA** interactions (blocks co-packaging)
and incompatibility of **protein–RNA** interactions (blocks replication).

Both are **competition-normalised**:

```
M_x = affinity(non-cognate pairing) / affinity(cognate pairing)
```

so `M_x = 1.0` means the reassortant arrangement is thermodynamically indistinguishable
from the parental one, and `M_x → 0` means the cognate pairing dominates.

**Formulation note — do not implement this as a literal quotient of ΔG values.** Free
energies are negative and a raw ratio is both sign-unstable and meaningless when the
non-cognate duplex is unfavourable (ΔG > 0). Define it as a Boltzmann-style term on the
*difference*:

```
M_x = min(1, exp(-ΔΔG / RT)),  ΔΔG = ΔG_non-cognate - ΔG_cognate
```

which is bounded in (0, 1], is well defined for positive ΔG, and expresses the competition
as a relative occupancy rather than an arbitrary quotient.

**Sign convention — do not flip.** A weaker non-cognate duplex has a *less negative* ΔG, so
ΔΔG > 0 and `M_x < 1`. Reversing the subtraction inverts the index and would score
incompatible constellations as maximally plausible. Assert this direction in a unit test.

**Units pin.** RNAcofold reports ΔG in **kcal/mol** at a default temperature of **37 °C =
310.15 K**. `RT` must match those units:

```
RT = 0.0019872 kcal·mol⁻¹·K⁻¹ × 310.15 K ≈ 0.6163 kcal/mol
```

Both the temperature and the energy-parameter set are pinned alongside the ViennaRNA
version, because ΔG values are parameter-set dependent.

**Methods sentence (use verbatim).** *"The `exp(−ΔΔG/RT)` form is the relative Boltzmann
weight of the non-cognate versus cognate duplex — the equilibrium expression of the
competition argument of McDonald et al. (2016), in which low-affinity non-cognate RNA–RNA
interactions are outcompeted by higher-affinity cognate ones."* (McDonald et al. 2016, *Nat
Rev Microbiol* 14:448–460, doi:10.1038/nrmicro.2016.46 — Tier 1, read in full.)

**Citation status — VERIFY BEFORE PUBLICATION.** The ViennaRNA references supplied for this
section could not be checked during design: neither is in the project library, and web
access was unavailable in the design session. Recorded as given, both **unverified**:
Lorenz et al. 2011, *Algorithms Mol Biol* 6:26, doi:10.1186/1748-7188-6-26 (ViennaRNA
Package 2.0); Bernhart et al. 2006, *Algorithms Mol Biol* 1:3 (RNAcofold). Confirm both,
and confirm the installed version's documented default temperature, before either is cited
or the RT constant is trusted. This directly
implements McDonald's finding that reverse genetics can build reassortants that co-infection
almost never yields, because "low-affinity interactions between the non-cognate RNAs … are
readily outcompeted by the optimal, higher-affinity interactions between the cognate RNAs."

### 6.1 M_pack — RNA–RNA co-packaging

Scored by **RNAcofold ΔG ratios** (ViennaRNA, named and version-pinned) between the
exchanged segment's terminal region and the recipient background's, over the cognate
baseline. Physics, not fitting — this is what keeps the "no fitted coefficients" claim
honest.

**New dependency:** ViennaRNA is not currently installed and is absent from
`environment.yml`. Adding it is a real implementation cost, and the version must be pinned
because ΔG values are parameter-set dependent.

### 6.2 M_pol — protein–RNA polymerase recognition

**This sub-factor must not score the genus-conserved terminal tag.** The terminal 14 nt of
the S, M and L segments are highly conserved across the genus *and shared between the three
segments* (*J Virol* 1995;69(12):8132–8136, PMC189768 — **content verified by the project
owner; author list not verified**). Any motif-match metric computed there saturates near
1.0 for every pair and carries zero discriminating signal. Scoring it would manufacture a
factor that looks principled and measures nothing. That citation is the evidence for the
saturation argument and must accompany it.

`M_pol` is therefore scored on features **beyond** the invariant tag:

- **N-recognition footprint length.** Mir MA & Panganiban AT, *"Hantavirus N Protein
  Exhibits Genus-Specific Recognition of the Viral RNA Panhandle,"* J Virol 2006,
  doi:10.1128/JVI.00820-06 — **existence and abstract verified by the project owner; full
  text not read.** The abstract establishes two things this factor needs: hantavirus N
  recognises the panhandle **genus-specifically**, and **SNV N requires the first 23 nt
  where most hantavirus N proteins need only the terminal 9**. That 23-vs-9 asymmetry is
  the useful signal precisely because it *varies between lineages*, unlike the invariant
  14-nt tag.
- panhandle duplex length,
- NCR-internal secondary structure beyond the conserved terminus.

**Dropped: the position-9 polymerase claim.** An earlier draft scored "position-9 panhandle
mismatch" as a *polymerase*-binding feature. That claim could not be traced to a primary
source and is **not** what the Mir & Panganiban abstract supports — that abstract concerns
**N protein**, not the L polymerase. Rather than carry it as pending-primary-source, it is
removed. The verified N-recognition footprint above replaces it and is better suited to the
job. If a primary source for polymerase position-9 recognition is later found, it may be
reinstated as a separate term.

**Mandatory variance check.** Before `M_pol` is admitted to the product, compute it across
a representative sample of lineage pairs. **If it shows negligible variance, fix it at 1.0
and state so explicitly in the output schema and Methods.** A constant factor that is
declared constant is honest; a constant factor dressed as a measurement is not.

**Framing requirement.** The rotavirus precedent — polymerase recognising a seven-nucleotide
3′ motif, with groups A/C/D/F (UGUGACC/UGUGGCU) unable to replicate B/G/H templates
(AAAACCC/AAGACCC/UAUACCC) — is a **structurally motivated hypothesis** for hantavirus
L-polymerase terminal recognition, not an experimentally demonstrated hantavirus mechanism.
The Methods text must say so. Supporting but indirect: Lowen (2018) notes heterologous
combinations of bunyavirus polymerase and nucleoprotein genes are "poorly compatible."

## 7. Factor K — constellation prior

The empirical prior over which segment is the odd one out, conditioned on divergence
stratum.

**Raw observations:**

| Source | Pair | Stratum | M-odd | S-odd | L-odd | n |
|---|---|---|---|---|---|---|
| Rodriguez 1998 | SNV NMR11 × CC107 | close | 12 | 12 | 1 | 25 |
| Rizvanov 2004 | SNV × ANDV | distant | 10 | 0 | 0 | 10 |
| Klempa 2018 Tab. 2 | all published in vitro | mixed | all | 0 | 0 | — |

**Jeffreys-prior smoothing is mandatory.** Rizvanov's `S-odd = 0` and `L-odd = 0` are
*unobserved events in n = 10*, not measured impossibilities. A literal zero inside a
multiplicative model is a hidden hard gate — the identical cliff logic rejected for
`gate3`'s veto. Apply Dirichlet(½, ½, ½), i.e. add 0.5 to every category:

```
K_close   = (12.5, 12.5,  1.5) / 26.5 = (0.472, 0.472, 0.057)
K_distant = (10.5,  0.5,  0.5) / 11.5 = (0.913, 0.043, 0.043)
```

The smoothing constant is stated in the output schema and **varied in the sensitivity
analysis** (§12). It also makes the retrospective test (§11.1) non-circular: unsmoothed
zeros would trivially reproduce their own input.

**Why pooling across studies is legitimate here, given §0.4.** §0.4 forbids comparing
reassortment *rates* across studies because MOI differs (Rizvanov 3, Rodriguez 1.0–1.5). K
does not use rates. It uses the *conditional composition* of the reassortant pool — given
that a reassortant was recovered, which segment was the odd one out. MOI governs how often
co-infection occurs and therefore how many reassortants appear; there is **no known
packaging-step mechanism** by which it would favour M-odd over S-odd once co-infection has
happened. **A residual recovery-stage bias via trans-complementation is possible and
unquantifiable from the published methods.**

That caveat is *differentially applicable*, and the methods sections say why:

| Study | MOI | Genotyping | Amplification before assay |
|---|---|---|---|
| Rodriguez 1998 | 1.0–1.5 | Plaque-level; RNA extracted directly from agarose plugs | **None** — "without the need for further virus amplification in tissue culture… first-generation progeny" |
| Rizvanov 2004 | 3 | Plaque-level | **Yes** — "isolated plaques were passed in Vero E6 cells and the progeny assayed"; plus "two rounds of selection" |

Both are plaque-level, which bounds complementation's reach in each. But only Rizvanov's
design contains a passage/selection step at which recovery-stage bias could operate.
Rodriguez's first-generation design excludes it by construction. The trans-complementation
concern therefore applies asymmetrically — to the distant-stratum anchor, not the close
one. This is **analytical reasoning, not a measured result**; no direct citation for
trans-complementation at high MOI exists in the hantavirus literature.

This assumption is stated rather than proven, and it is what the leave-one-pair-out test in
§11.1 probes.

**Assumption register — pair-identity confound. UNTESTABLE WITH CURRENT DATA.**

The claim "distant pairs collapse to M-only" rests on a very thin evidence table:

- **Close stratum:** one pair — SNV NMR11 × CC107 (Rodriguez 1998), n = 25.
- **Distant stratum:** two pairs — SNV × ANDV (Rizvanov 2004), n = 10; and SNV × BCCV
  (Rodriguez 1998), n = 1. Total n = 11, all M-odd.
- **Aggregated pattern:** Klempa 2018 Table 2.

Divergence is therefore **confounded with pair identity**: every distant observation
involves a different virus pair from the close observation, so "distant" and "not-SNV×SNV"
cannot be separated. Worse, the strata are also confounded with MOI (close = MOI 1,
distant = MOI 3 for the n=10 arm). **No distant-pair MOI-1 experiment exists in the
literature**, so this cannot be resolved with published data. It is untestable with current
evidence and must be reported as such, not as a validated stratification.

**A correction to the aggregated pattern.** Klempa's abstract states that "only reassortants
with exchanged M segments could be generated." His own Table 2 does not support that as a
universal rule: for the close pair it lists S_CC107 M_CC107 L_NMR11 — an **L-odd**
reassortant — and in fact reproduces exactly Rodriguez's three observed genotypes (M-odd,
S-odd, L-odd). The M-only rule is **stratum-dependent, not genus-universal**. This spec
relies on the stratum-dependent reading, which is what both primary datasets actually show.

**Interpretation — this is the central molecular claim.** Divergence does not scale the
overall reassortment *rate* (§0.4 forbids that inference); it narrows which constellations
are *accessible*. Rizvanov states it directly: homologous L–M and L–S segments appeared in
equal proportion for closely related SNV strains, whereas distantly related viruses produced
"a single type of reassortant virus that contained S and L homologous segments." Klempa's
independent synthesis reaches the same place from the whole in vitro literature: "only
reassortants with exchanged M segments could be generated, suggesting that a high degree of
genetic compatibility is required for the S and L segments while the exchange of M segment
is better tolerated."

The pipeline already encodes this pattern in the `matches_canonical_m_exchange` column.

**Stratum boundary.** Close/distant strata are derived by fitting a mixture model to the
pairwise evolutionary distance (PED) distribution of the project's own dataset. This is a
**"DEmARC-style analysis"** — an adaptation of the distance-partitioning idea, **not** the
DEmARC tool. DEmARC proper is the ICTV species-demarcation criterion for *Hantaviridae*,
as correctly documented at `all_hanta.py:1866` and `README.md:487`. The distinction must be
maintained in code comments and Methods.

*Anchor caution:* the two Tier-1 sources disagree on SNV/BCCV divergence — Rodriguez's
Table 7 gives ~34% (N coding sequence), Rizvanov gives 70.7% S-segment nt similarity
(~29.3% divergence). Different regions and methods. Use the project's own PED values for
stratification and cite the published figures only as context.

## 8. Factor F — retention

`F` scales `rpi_generation` into `rpi_retention` using the measured generation→retention
shift:

- Generation (Rodriguez, in vitro): M : S : L = 12 : 12 : 1
- Retention (Henderson, natural, **Tier 3 secondhand**): M : S : L = 14 : 2 : 1

S-exchange is generated as often as M-exchange yet survives markedly less often in nature.
That ratio is the selection filter, measured rather than assumed. Gate 2's existing
structural checks (WAASA, TM/IS, pH switch) contribute as integrity terms.

**Scope, per §0.7:** `F` represents **immediate** post-reassortment fitness. McDonald's
Fig. 4 modes are *not* decomposed into three separate factors — mode (a) beneficial-allele
acquisition and mode (b) uncoupling of co-adapted allele sets are what `F` summarises;
mode (c), compensatory adaptation restoring fitness over time, is explicitly out of scope
and named as a limitation.

The small denominators (2 S-odd, 1 L-odd in the natural data) mean `F` is weakly determined.
Sensitivity analysis must show whether conclusions survive its variation.

## 9. Output schema

Appended to the `reassortment_detector.py` output row (all factors exposed, per house style):

| Column | Meaning |
|---|---|
| `rpi_generation` | E × M_pack × M_pol × K, rounded to 4 dp |
| `rpi_retention` | rpi_generation × F |
| `rpi_e_host`, `rpi_e_geo`, `rpi_e_time`, `rpi_e_total` | Factor E and components |
| `rpi_m_pack`, `rpi_m_pol` | Factor M sub-terms |
| `rpi_m_pol_status` | `scored` or `fixed_constant` (per §6.2 variance check) |
| `rpi_k_stratum` | `close` / `distant`, and which prior was applied |
| `rpi_k_value` | Smoothed constellation prior used |
| `rpi_f` | Retention factor |
| `rpi_smoothing` | Jeffreys constant in force (default 0.5) |
| `rpi_weights_version` | Version tag of the weight constant block |

Weights live in a single named constant block (mirroring `SURVEILLANCE_WEIGHTS`) so the
sensitivity analysis can sweep them without touching logic.

## 10. Integration

### 10.1 Parallel column; verdict unchanged

`combine_verdict` keeps its closed 10-verdict list **exactly as specified in the
2026-08-12 design**. No existing behaviour changes; that spec is not amended; the existing
`ortho_reassortment_table.csv` does not need regenerating for correctness.

Rationale, per the project owner's decision:

1. The verdict is an **observational report**; the RPI is a **mechanistic model**. Merging
   a belief into a detection report corrupts the audit trail.
2. Validation requires independence. Rank-consistency against the 8 supported calls is only
   a genuine test if those verdicts were computed without reference to the index.

### 10.2 RPI is computed for every candidate, including vetoed ones

Including candidates that `gate3` vetoed as `no_ecological_route`. This is not redundant —
it enables a **cross-check that validates both layers**: does the hard categorical veto
coincide with near-zero continuous `E`? Divergence between the two flags a defect in one of
them. This check is only possible because the index is computed in parallel rather than
gated by the verdict.

Note that `gate3_check` (`reassortment_detector.py:730`) currently examines only
`implicated_species[:2]` and has **no temporal dimension**, while factor E covers all
implicated species and includes `t_time`. Disagreement is therefore expected and
informative; it is the point of the check, not a bug.

## 11. Validation plan, in order of evidential weight

### 11.1 Primary — retrospective reproduction of the in vitro experiments

Feed the index the five published parental pairs and test whether it reproduces the
observed constellation distributions:

| Pair | Source | Observed |
|---|---|---|
| SNV NMR11 × CC107 | Rodriguez 1998 | M 12, S 12, L 1 |
| SNV × ANDV | Rizvanov 2004 | M 10, S 0, L 0 (monoploid) |
| SNV × BCCV | Rodriguez 1998 | 1 reassortant (M-odd) in 163 |
| PUUV × PHV | Handke 2010 (read the library PDF first) | M-exchanged |
| DOBV-Aa × DOBV-Af | Kirsanovs 2010 (Tier 3) | 65/207, M-exchanged |

This is the **primary evidence**. Ecological factor E is held at 1.0 for these, since the
experiments forced co-infection — which is exactly why the in vitro data calibrates the
molecular factors and cannot calibrate E.

Circularity guard: pairs used to *derive* K cannot serve as independent tests of K.
Report leave-one-pair-out reproduction.

### 11.2 Secondary — rank-consistency sanity check

Do the 8 `supported_reassortment_candidate` calls score above the 2,961 `no_conflict`
rows? **Present this as a sanity check with a wide confidence interval, not as
validation.** With 8 positives the interval on any rank statistic is very wide, and the
check must be reported with that interval attached rather than as a bare pass/fail.

### 11.3 Cross-layer check

Compare the `no_ecological_route` veto against near-zero `E`, per §10.2.

### 11.4 Optional external cross-check — GiRaF

Graph-incompatibility reassortment detection as an independent opinion on Gate 4, **scoped
to the 8 supported candidates only**. GiRaF consumes Bayesian posterior tree samples; the
pipeline runs IQ-TREE ML with UFBoot, so this requires a separate MrBayes/BEAST run that
does not compose with the pipeline at 2,995-row scale. Spot-check only.

The RDP5-style recombination screen (≥2 methods agreeing) remains as specified in the
2026-08-12 design §274 — a **control** ruling out recombination as an alternative
explanation for tree incongruence, not a reassortment detector. Reassortment and
recombination stay separately labelled.

## 12. Sensitivity analysis

Required before any result is reported. Vary each of the following ±50% and report rank
stability:

- E's species→genus→family→order decay weights **and the functional form of the decay**
  (§5 — these are assumptions, not measurements);
- the Jeffreys smoothing constant (§7);
- F's generation→retention ratio (§8), whose small denominators make it weakly determined;
- the close/distant stratum boundary (§7).

Report as rank-correlation against the baseline ordering. The headline claim the analysis
must support is: *do the top-ranked candidates stay top-ranked?*

## 13. Open questions

1. Does local prevalence belong in E? Rivero (Tier 2) says it trends negative; Lowen
   (Tier 1) lists it as governing. Currently omitted (§5).
2. Can directional asymmetry (§0.5) be captured? RNAcofold ΔG is computable in both
   directions — whether that recovers the φ6/φ13 phenomenon is untested.
3. Should cross-segment covariation (mutual information) be added as a fourth molecular
   term? It is the peer-reviewed ancestor of Rivero's linkage-disequilibrium finding.
   **If pursued, cite the primary paper directly, not through McDonald's reference chain:**
   "Viral reassortment as an information exchange between viral segments," *PNAS* 2012,
   doi:10.1073/pnas.1113300109 (Tier 4a — title and DOI verified, author list not).
   Zhang et al. 2014 (Tier 3, McDonald ref. 112) is a different, rotavirus-specific
   co-adaptation paper and is not a substitute for it.
4. Confirm Mir & Panganiban 2006 (Tier 4) before publication.
5. Read Handke et al. 2010 from the library PDF and promote it from Tier 3.

## 14. Definition of done

- [ ] `rpi_*` columns present for every candidate, vetoed ones included
- [ ] `combine_verdict` behaviour byte-identical to pre-change
- [ ] Retrospective reproduction (§11.1) run and reported, with leave-one-pair-out
- [ ] Rank-consistency (§11.2) reported **with confidence interval**
- [ ] `M_pol` variance check executed; status recorded in `rpi_m_pol_status`
- [ ] Sensitivity sweep (§12) complete
- [ ] ViennaRNA pinned to an exact version in `environment.yml` (bioconda `viennarna=<x.y.z>`),
      with the energy-parameter set and default temperature recorded alongside it, per §6.1
- [ ] ΔΔG sign convention asserted by a unit test (§6.1)
- [ ] Tier 4a/4b citations resolved: Mir & Panganiban full text read; author lists obtained
      for PMC189768 and the PNAS 2012 paper; Lorenz and Bernhart confirmed
- [ ] Limitations (§0) reproduced verbatim in README and any Methods text
