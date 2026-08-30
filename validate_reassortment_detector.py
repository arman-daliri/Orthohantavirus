#!/usr/bin/env python3
"""Validation harness for reassortment_detector.py.

WHAT THIS ANSWERS
-----------------
Does the detector actually fire on reassortants, and stay quiet on non-reassortants?

The unit tests in reassortment_detector.py check that each gate computes what it says
it computes. They do NOT answer the question above, because a detector can be
internally consistent and still be blind. This harness scores a panel of candidates
whose ground truth is known independently of the detector, and reports a confusion
matrix.

THE PANEL
---------
Four classes, each with a DIFFERENT expected outcome. A harness where everything is
expected to come back positive proves nothing -- see class C and class E.

  A. PUBLISHED INTER-SPECIES REASSORTANTS. Reconstructed in silico from the real
     genome sequences of the parental species named in the source paper. Expected:
     gate1 flags the exchanged segment.

  B. PUBLISHED INTRA-SPECIES REASSORTANTS. Same, but both parents are the same ICTV
     species. Expected: gate1 CANNOT see these -- see "known blind spot" below. This
     class exists to measure the blind spot, not to pass.

  C. TRUE NEGATIVES. Real, unmodified isolates: all three segments from one animal.
     Expected: no discordance. This is the class that can produce false positives, and
     without it a "100% detection rate" would be meaningless.

  D. SYNTHETIC INTER-SPECIES SWAPS. Constructed cross-species exchanges spanning every
     segment position (S, M and L donated in turn) and a range of parental distances.
     Not from any paper -- these exist to measure detection across a wider sweep than
     the published set covers, and to check the detector names the RIGHT segment.

  E. ADVERSARIAL. Junk, scrambled and truncated input. Expected: rejected at gate0 or
     otherwise handled, never silently scored as a confident reassortant.

PROVENANCE AND HONESTY
----------------------
Class A and B candidates are RECONSTRUCTIONS OF PUBLISHED CONSTELLATIONS, not the
published reassortant clones themselves. The papers report which parental strains were
crossed and which segment was exchanged; this harness rebuilds that genotype from the
parental sequences in all_hanta/. Where the paper's exact parental strain is present in
the dataset it is used and named in the output (`parents_are_exact_published_strains`);
where only the species is available, a representative isolate of that species is used
and the row says so. No sequence here is invented.

KNOWN BLIND SPOT (measured by class B, not hidden by it)
--------------------------------------------------------
gate1 detects reassortment by assigning each segment to its nearest reference SPECIES
by protein k-mer profile and looking for disagreement. When both parents belong to the
same species, all three segments assign to that same species and gate1 reports
"concordant" -- correctly, on its own terms, and uselessly for that event class.
Intra-species reassortment is what gate4's per-segment phylogenetics exists to catch.

GATE 4 IS NOT EXERCISED HERE
----------------------------
gate4 needs mafft and iqtree. Neither is installed in the environment this harness was
written in, and substituting fake binaries would produce trees from random numbers --
which would not validate anything. gate4 is therefore skipped (via a comparator
threshold that makes it return "insufficient_comparators" without running), and every
verdict below is a verdict-with-gate4-withheld. Run with --with-phylogeny once real
binaries are on PATH to lift this.

Usage:
    python validate_reassortment_detector.py
    python validate_reassortment_detector.py --json validation_report.json
"""

import argparse
import json
import random
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import reassortment_detector as rd

# Sentinel comparator threshold that makes gate4_check return early with
# "insufficient_comparators" instead of shelling out to mafft/iqtree. Chosen far above
# any real species' isolate count so it cannot accidentally admit a tree.
_PHYLO_DISABLED_MIN_TAXA = 10 ** 9


# ---------------------------------------------------------------------------
# Published crosses. `paper` and `reported` are quoted from the source; `expect` is
# this harness's prediction for gate1, derived from the cross being inter- or
# intra-species, NOT from the detector's output.
# ---------------------------------------------------------------------------

PUBLISHED_CROSSES: List[Dict[str, Any]] = [
    {
        "case_id": "RIZVANOV2004_SNVsl_x_ANDVm",
        "klass": "A_published_inter_species",
        "paper": "Rizvanov, Khaiboullina & St Jeor (2004) Virology 327:225-232, "
                 "doi:10.1016/j.virol.2004.07.012",
        "reported": "Reassortants generated in vitro between Sin Nombre virus and Andes "
                    "virus. Viable progeny carried the SNV S and L segments with the ANDV "
                    "M segment; no S- or L-exchanged reassortant was recovered.",
        "backbone_species": "Orthohantavirus sinnombreense",
        "donor_species": "Orthohantavirus andesense",
        "donated_segment": "M",
        "backbone_strain_hint": "nm r11",
        "donor_strain_hint": "chile",
        "expect_pattern": "m_only",
        "expect_discordant": "M",
    },
    {
        "case_id": "RIZVANOV2004_reciprocal_ANDVsl_x_SNVm",
        "klass": "A_published_inter_species",
        "paper": "Rizvanov, Khaiboullina & St Jeor (2004) Virology 327:225-232",
        "reported": "Reciprocal orientation of the same SNV x ANDV cross. Included as a "
                    "directionality control: the detector must not be sensitive to which "
                    "parent is called the backbone.",
        "backbone_species": "Orthohantavirus andesense",
        "donor_species": "Orthohantavirus sinnombreense",
        "donated_segment": "M",
        "backbone_strain_hint": "chile",
        "donor_strain_hint": "nm r11",
        "expect_pattern": "m_only",
        "expect_discordant": "M",
    },
    {
        "case_id": "HANDKE2010_PUUVsl_x_PHVm",
        "klass": "A_published_inter_species",
        "paper": "Handke, Oelschlegel, Franke, Wiedemann, Kruger & Rang (2010) "
                 "J Gen Virol 91:2351-2359, doi:10.1099/vir.0.021139-0",
        "reported": "Genetic reassortants generated between Puumala virus and Prospect "
                    "Hill virus. Only M-segment exchange yielded viable reassortants.",
        "backbone_species": "Orthohantavirus puumalaense",
        "donor_species": "Orthohantavirus prospectense",
        "donated_segment": "M",
        "expect_pattern": "m_only",
        "expect_discordant": "M",
    },
    {
        "case_id": "RODRIGUEZ1998_SNV_NMR11sl_x_CC107m",
        "klass": "B_published_intra_species",
        "paper": "Rodriguez, Owens, Peters & Nichol (1998) Virology 242:99-106, "
                 "doi:10.1006/viro.1997.8990",
        "reported": "Genetic reassortment among Sin Nombre viruses causing hantavirus "
                    "pulmonary syndrome. Coinfection of SNV strains NM R11 and Convict "
                    "Creek 107 (CC107) at MOI 1 yielded 25 reassortants; Table 7 records "
                    "12 M-exchanged and 12 S-exchanged progeny.",
        "backbone_species": "Orthohantavirus sinnombreense",
        "donor_species": "Orthohantavirus sinnombreense",
        "donated_segment": "M",
        "backbone_strain_hint": "nm r11",
        "donor_strain_hint": "cc107",
        "expect_pattern": "concordant",
        "expect_discordant": "",
        "blind_spot": True,
    },
    {
        "case_id": "RODRIGUEZ1998_SNV_NMR11ml_x_CC107s",
        "klass": "B_published_intra_species",
        "paper": "Rodriguez et al. (1998) Virology 242:99-106, Table 7",
        "reported": "The S-exchanged arm of the same experiment (12 of 25 reassortants). "
                    "Rodriguez's central result is that S exchange is as frequent as M "
                    "exchange between close parents.",
        "backbone_species": "Orthohantavirus sinnombreense",
        "donor_species": "Orthohantavirus sinnombreense",
        "donated_segment": "S",
        "backbone_strain_hint": "nm r11",
        "donor_strain_hint": "cc107",
        "expect_pattern": "concordant",
        "expect_discordant": "",
        "blind_spot": True,
    },
    {
        # D5 / Task 12. donor_species is "Orthohantavirus nigrorivense", NOT
        # "Orthohantavirus bayoui" -- an earlier lookup for this task assumed ICTV had
        # folded Black Creek Canal virus (BCCV) into the Bayou virus species, but
        # this project's own curated species table (all_hanta.py's SPECIES_INFO,
        # confirmed against output/ortho_species_table.csv) already resolves BCCV to
        # its own species, "Orthohantavirus nigrorivense" (abbrev "BCCV"), separate
        # from "Orthohantavirus bayoui". Verified directly: every isolate under
        # "Orthohantavirus nigrorivense" in all_hanta/ is a Black Creek Canal virus
        # deposit (NC_043073-075.1, L39949-51.1, GU997097.1, AB689163.1) -- no Bayou
        # virus sequence is lumped in under this name, so no donor_strain_hint is
        # needed to steer away from one (unlike the plan's original assumption).
        #
        # No BCCV isolate has a coding-complete S+M+L triplet in this dataset (every
        # deposit above is single-segment). That does NOT block this cross: the donor
        # only contributes the M segment, and find_donor does not require a triplet
        # (see its own docstring) -- NC_043073.1 (3668 nt) is a full-length,
        # ICTV-band M segment on its own.
        "case_id": "RODRIGUEZ1998_SNV_NMR11sl_x_BCCVm",
        "klass": "A_published_inter_species",
        "paper": "Rodriguez, Owens, Peters & Nichol (1998) Virology 242:99-106, "
                 "doi:10.1006/viro.1997.8990",
        "reported": "Coinfection of Sin Nombre virus strain NM R11 with Black Creek "
                    "Canal virus (a distinct ICTV species, not a same-species strain) "
                    "at MOI 1 yielded a single M-exchanged reassortant -- markedly "
                    "rarer than the same paper's intra-species NM R11 x CC107 cross "
                    "(12 of 25). RPI_INVITRO_BENCHMARK['SNV_x_BCCV'] "
                    "(reassortment_detector.py) already records these counts for the "
                    "K-stratum prior; this entry reconstructs the candidate itself.",
        "backbone_species": "Orthohantavirus sinnombreense",
        "donor_species": "Orthohantavirus nigrorivense",
        "donated_segment": "M",
        "backbone_strain_hint": "nm r11",
        # expect_pattern is "m_only", NOT the measured "undetectable_by_design" gate1
        # actually returns on this real pair (verified 2026-08-29) -- this module's own
        # docstring requires expect_* to be an independent prediction derived from the
        # cross being inter/intra-species, "NOT from the detector's output", and setting
        # it to whatever gate1 happens to emit would violate exactly that. The measured
        # mismatch is real and worth recording here rather than silently: BCCV's M
        # protein sits closer to the Orthohantavirus bayoui (Bayou virus) reference
        # profile than to any distinguishing nigrorivense-specific signal, and the gap
        # between top candidates falls inside undetectable_floor -- gate1 STILL names
        # the correct segment (discordant_segment == "M"), so this case passes the
        # PRIMARY criterion (`localized`) and lands in summarise()'s existing
        # "localized_but_parents_unresolvable" bucket, not "also_exact_pattern_match".
        # This is a genuinely new blind-spot flavour: a species pair that IS
        # taxonomically distinct in this project's curated table but is NOT reliably
        # distinguishable at the protein level from a third, related species -- gate1's
        # blind spot is not limited to literal same-species pairs (class B).
        "expect_pattern": "m_only",
        "expect_discordant": "M",
    },
    {
        "case_id": "KLEMPA2018_DOBV_intragenotype_M",
        "klass": "B_published_intra_species",
        "paper": "Klempa (2018) Virus Genes 54:638-646, doi:10.1007/s11262-018-1590-z",
        "reported": "Review: intra- and inter-lineage reassortment repeatedly reported "
                    "for Dobrava-Belgrade virus among its Apodemus-associated genotypes.",
        "backbone_species": "Orthohantavirus dobravaense",
        "donor_species": "Orthohantavirus dobravaense",
        "donated_segment": "M",
        "expect_pattern": "concordant",
        "expect_discordant": "",
        "blind_spot": True,
    },
    {
        "case_id": "PARK2025_HTNV_intraspecies_L",
        "klass": "B_published_intra_species",
        "paper": "Park, Kim, Kim, Kim & Song (2025) Sci Rep, "
                 "doi:10.1038/s41598-025-11003-z",
        "reported": "Longitudinal genomic surveillance of Orthohantavirus hantanense in "
                    "the Republic of Korea reports reassortant candidates with "
                    "heterogeneous L segments, i.e. S/L and M/L incongruence within the "
                    "species.",
        "backbone_species": "Orthohantavirus hantanense",
        "donor_species": "Orthohantavirus hantanense",
        "donated_segment": "L",
        "expect_pattern": "concordant",
        "expect_discordant": "",
        "blind_spot": True,
    },
    {
        "case_id": "KLEMPA2018_PUUV_intralineage_M",
        "klass": "B_published_intra_species",
        "paper": "Klempa (2018) Virus Genes 54:638-646",
        "reported": "Review: intra-lineage reassortment repeatedly reported for Puumala "
                    "virus among European lineages.",
        "backbone_species": "Orthohantavirus puumalaense",
        "donor_species": "Orthohantavirus puumalaense",
        "donated_segment": "M",
        "expect_pattern": "concordant",
        "expect_discordant": "",
        "blind_spot": True,
    },
]


class PanelBuilder:
    """Builds candidates from real sequences in all_hanta/."""

    def __init__(self, reference: rd.ReferenceData):
        self.reference = reference
        self.by_species: Dict[str, List[Any]] = {}
        for isolate in reference.isolates:
            self.by_species.setdefault(isolate.species_name, []).append(isolate)

    def genomic(self, accession: str) -> str:
        if not accession:
            return ""
        record = (self.reference.genomic_index.get(accession)
                  or self.reference.genomic_index.get(accession.split(".")[0]))
        return record.sequence if record else ""

    def _accession(self, isolate, segment: str) -> str:
        return {"S": isolate.s_accession, "M": isolate.m_accession,
                "L": isolate.l_accession}[segment]

    def _usable(self, isolate, segment: str) -> bool:
        """Is this isolate's `segment` a plausibly FULL-LENGTH sequence?

        Presence is not enough. The dataset is full of partial deposits, and an earlier
        version of this harness picked them: a PUUV backbone whose M was 1014 nt and L
        665 nt was assembled into a "reassortant" that gate0 then correctly rejected as
        junk. That measured this harness's isolate picking, not the detector. Requiring
        the ICTV segment band keeps the panel made of real, complete genomes so a
        junk_sequence verdict means something about the detector.
        """
        sequence = self.genomic(self._accession(isolate, segment))
        low, high = rd.ah.SEGMENT_LENGTH_RANGE[segment]
        return bool(sequence) and low <= len(sequence) <= high

    def find_donor(self, species: str, segment: str,
                   hint: Optional[str] = None, exclude=None) -> Optional[Any]:
        """An isolate of `species` carrying a full-length `segment`.

        A donor contributes ONE segment, so it does not need a complete triplet -- which
        matters for species like Prospect Hill that are in the dataset but have no
        complete triplet at all.
        """
        members = [m for m in self.by_species.get(species, [])
                   if m is not exclude and self._usable(m, segment)]
        if hint:
            hinted = [m for m in members
                      if hint.lower() in f"{m.isolate} {m.strain}".lower()]
            if hinted:
                return hinted[0]
        return members[0] if members else None

    def find_backbone(self, species: str, needed: Tuple[str, ...],
                      hint: Optional[str] = None) -> Optional[Any]:
        """An isolate of `species` carrying every segment in `needed`, full-length."""
        members = self.by_species.get(species, [])
        ok = [m for m in members if all(self._usable(m, seg) for seg in needed)]
        if hint:
            hinted = [m for m in ok if hint.lower() in f"{m.isolate} {m.strain}".lower()]
            if hinted:
                return hinted[0]
        return ok[0] if ok else None

    @staticmethod
    def label(isolate) -> str:
        name = (isolate.isolate or isolate.strain or "").strip()
        return name or isolate.isolate_key


def build_swap(builder: PanelBuilder, spec: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Assemble one reassortant: backbone contributes two segments, donor contributes one."""
    donated = spec["donated_segment"]
    kept = tuple(s for s in ("S", "M", "L") if s != donated)
    backbone = builder.find_backbone(
        spec["backbone_species"], kept, spec.get("backbone_strain_hint"))
    if backbone is None:
        return None
    donor = builder.find_donor(
        spec["donor_species"], donated, spec.get("donor_strain_hint"), exclude=backbone)
    if donor is None or donor is backbone:
        return None

    segments = {seg: builder.genomic(builder._accession(backbone, seg)) for seg in kept}
    segments[donated] = builder.genomic(builder._accession(donor, donated))

    hint_b, hint_d = spec.get("backbone_strain_hint"), spec.get("donor_strain_hint")
    exact = bool(hint_b and hint_d
                 and hint_b.lower() in f"{backbone.isolate} {backbone.strain}".lower()
                 and hint_d.lower() in f"{donor.isolate} {donor.strain}".lower())

    segment_sources = {seg: backbone for seg in kept}
    segment_sources[donated] = donor

    return {
        **spec,
        "candidate": rd.CandidateInput(
            candidate_id=spec["case_id"],
            s_genomic=segments["S"], m_genomic=segments["M"], l_genomic=segments["L"],
            known_segment_isolates=segment_sources),
        "backbone_strain": PanelBuilder.label(backbone),
        "donor_strain": PanelBuilder.label(donor),
        "backbone_accessions": {s: builder._accession(backbone, s) for s in kept},
        "donor_accession": builder._accession(donor, donated),
        "parents_are_exact_published_strains": exact,
    }


def build_panel(builder: PanelBuilder,
                seed: int = 20260824) -> Tuple[List[Dict[str, Any]], List[str]]:
    panel: List[Dict[str, Any]] = []
    skipped: List[str] = []

    # --- A and B: published crosses -------------------------------------------------
    for spec in PUBLISHED_CROSSES:
        case = build_swap(builder, spec)
        if case is None:
            skipped.append(f"{spec['case_id']} (parental sequences not in this dataset)")
            continue
        panel.append(case)

    # --- C: true negatives ----------------------------------------------------------
    # Real unmodified isolates, one per species, across as many species as the dataset
    # offers complete triplets for. These are the rows a false positive would show up in.
    negatives = []
    for species, members in sorted(builder.by_species.items()):
        triplets = [m for m in members if m.has_complete_triplet
                    and all(builder.genomic(builder._accession(m, s)) for s in "SML")]
        if triplets:
            negatives.append((species, triplets[0]))
    for species, isolate in negatives:
        panel.append({
            "case_id": f"NEGATIVE_{species.split()[-1]}_{PanelBuilder.label(isolate)}"[:70],
            "klass": "C_true_negative",
            "paper": "n/a -- real unmodified isolate from all_hanta/",
            "reported": "All three segments sequenced from one animal. Not a reassortant.",
            "expect_pattern": "concordant",
            "expect_discordant": "",
            "candidate": rd.CandidateInput(
                candidate_id=f"NEG_{PanelBuilder.label(isolate)}",
                s_genomic=builder.genomic(isolate.s_accession),
                m_genomic=builder.genomic(isolate.m_accession),
                l_genomic=builder.genomic(isolate.l_accession)),
            "backbone_species": species,
            "donor_species": species,
            "donated_segment": "",
            "backbone_strain": PanelBuilder.label(isolate),
            "donor_strain": "",
            "backbone_accessions": {s: builder._accession(isolate, s) for s in "SML"},
            "donor_accession": "",
            "parents_are_exact_published_strains": False,
        })

    # --- D: synthetic inter-species swaps -------------------------------------------
    # Every segment position donated in turn, over species pairs of varying distance.
    # The published set is all M-exchange (that is what the literature contains); this
    # is what shows the detector is not simply hardcoded to notice M.
    donor_pairs = [
        ("Orthohantavirus seoulense", "Orthohantavirus andesense"),
        ("Orthohantavirus hantanense", "Orthohantavirus puumalaense"),
        ("Orthohantavirus puumalaense", "Orthohantavirus seoulense"),
        ("Orthohantavirus andesense", "Orthohantavirus chocloense"),
        ("Orthohantavirus dobravaense", "Orthohantavirus hantanense"),
    ]
    for backbone_species, donor_species in donor_pairs:
        for segment in ("S", "M", "L"):
            spec = {
                "case_id": f"SYNTH_{backbone_species.split()[-1]}_x_"
                           f"{donor_species.split()[-1]}_{segment}",
                "klass": "D_synthetic_inter_species",
                "paper": "n/a -- constructed control, not a published event",
                "reported": f"{segment} segment of {donor_species} placed on a "
                            f"{backbone_species} backbone.",
                "backbone_species": backbone_species,
                "donor_species": donor_species,
                "donated_segment": segment,
                "expect_pattern": "m_only" if segment == "M" else "s_or_l_involving",
                "expect_discordant": segment,
            }
            case = build_swap(builder, spec)
            if case is None:
                skipped.append(f"{spec['case_id']} (sequences unavailable)")
                continue
            panel.append(case)

    # --- E: adversarial -------------------------------------------------------------
    rng = random.Random(seed)
    real = builder.find_backbone("Orthohantavirus seoulense", ("S", "M", "L"))
    if real is None:
        # Every other builder call in this function degrades to `skipped`; this one used
        # to dereference None and take the whole run down with an AttributeError AFTER
        # the reference load, losing classes A-D too. _usable is strict (exact ICTV
        # band), so this is reachable on any reduced --reference-input.
        skipped.append("E_adversarial (no full-length SEOV triplet to derive it from)")
        return panel, skipped
    real_s = builder.genomic(builder._accession(real, "S"))
    real_m = builder.genomic(builder._accession(real, "M"))
    real_l = builder.genomic(builder._accession(real, "L"))

    def adversarial(case_id, reported, s, m, l, expect_verdict):
        panel.append({
            "case_id": case_id, "klass": "E_adversarial",
            "paper": "n/a -- adversarial control",
            "reported": reported,
            "expect_pattern": None, "expect_discordant": None,
            "expect_verdict": expect_verdict,
            "candidate": rd.CandidateInput(case_id, s, m, l),
            "backbone_species": "", "donor_species": "", "donated_segment": "",
            "backbone_strain": "", "donor_strain": "",
            "backbone_accessions": {}, "donor_accession": "",
            "parents_are_exact_published_strains": False,
        })

    random_seq = "".join(rng.choice("ACGT") for _ in range(1800))
    adversarial("ADV_random_nucleotides",
                "Uniform random ACGT in all three segments. Contains no hantavirus ORF.",
                random_seq, random_seq, random_seq, "junk_sequence")

    def shuffled(seq):
        chars = list(seq)
        rng.shuffle(chars)
        return "".join(chars)

    # Either verdict is acceptable and they mean different things. A shuffle preserves
    # base composition and segment length, so it can land a long-enough ORF by chance
    # and clear gate0's structural check -- then fail the identity floor, because it
    # resembles nothing. "divergent_novel_candidate" is the honest label for that. What
    # matters is that it is never scored as a confident reassortant; it is flagged for
    # triage either way.
    adversarial("ADV_shuffled_real_sequence",
                "Real SEOV sequences with their bases shuffled -- identical base "
                "composition to a real genome, no preserved coding structure.",
                shuffled(real_s), shuffled(real_m), shuffled(real_l),
                ("junk_sequence", "divergent_novel_candidate"))
    adversarial("ADV_truncated_segments",
                "Real SEOV sequences truncated to 120 nt -- far below any hantavirus "
                "ORF length.",
                real_s[:120], real_m[:120], real_l[:120], "junk_sequence")
    # An explicit expected verdict, not None. With expect_verdict=None this case scored
    # `passed = verdict is not None`, which is true of EVERY verdict score_candidate can
    # return -- including one that had scored the missing-M row as a confident
    # reassortant. A control that cannot fail is not a control, and this branch has
    # spent its whole review history removing exactly that shape of assertion.
    adversarial("ADV_empty_M_segment",
                "Real SEOV S and L with an empty M -- a missing-data row, not a "
                "reassortant. Must not be scored as a confident call.",
                real_s, "", real_l, ("junk_sequence",))

    return panel, skipped


def score_case(case: Dict[str, Any], reference, precomputed, with_phylogeny: bool,
               output_dir: Path) -> Dict[str, Any]:
    args = argparse.Namespace(
        output=str(output_dir), mafft_bin="mafft", iqtree_bin="iqtree",
        bootstrap=1000, min_bootstrap_support=95.0,
        phylo_min_taxa=4 if with_phylogeny else _PHYLO_DISABLED_MIN_TAXA,
        phylo_seed=42, phylo_threads=1,
    )
    row = rd.score_candidate(case["candidate"], reference, precomputed, args)

    observed_pattern = row.get("gate1_pattern", "")
    observed_segment = row.get("gate1_discordant_segment", "")

    # Two levels, deliberately not conflated:
    #
    #   localized  -- gate1 named the segment that was actually exchanged. This is the
    #                 operationally useful result and the PRIMARY metric.
    #   exact      -- gate1 ALSO classified the pattern the way this harness predicted.
    #
    # They differ for a real reason. When the two candidate parent species are close
    # enough that gate1 cannot tell them apart it returns "undetectable_by_design"
    # while still naming the right segment -- an honest "I see the exchange, I cannot
    # resolve the parents", not a miss. Scoring only `exact` would count that as a
    # failure and understate the detector; scoring only `localized` would hide that the
    # verdict a user actually reads says "undetectable".
    localized = bool(case.get("expect_discordant")) and (
        observed_segment == case["expect_discordant"])
    exact = (observed_pattern == case["expect_pattern"]
             and observed_segment == case["expect_discordant"])

    if case["klass"] == "E_adversarial":
        expected = case.get("expect_verdict")
        # No permissive fallback. Every adversarial case must name the verdict(s) it
        # accepts; a missing expectation is a bug in the panel, not a pass.
        if not expected:
            raise ValueError(
                f"{case['case_id']}: adversarial cases must declare expect_verdict")
        accepted = expected if isinstance(expected, (list, tuple)) else (expected,)
        passed = row.get("verdict") in accepted
        localized = exact = passed
    elif case["klass"] in ("C_true_negative", "B_published_intra_species"):
        # Nothing to localize -- these must simply stay quiet.
        passed = (observed_pattern == case["expect_pattern"]
                  and observed_segment == case["expect_discordant"])
        localized = exact = passed
    else:
        passed = localized

    return {
        "case_id": case["case_id"],
        "class": case["klass"],
        "paper": case["paper"],
        "reported": case["reported"],
        "backbone_species": case.get("backbone_species", ""),
        "donor_species": case.get("donor_species", ""),
        "donated_segment": case.get("donated_segment", ""),
        "backbone_strain": case.get("backbone_strain", ""),
        "donor_strain": case.get("donor_strain", ""),
        "backbone_accessions": case.get("backbone_accessions", {}),
        "donor_accession": case.get("donor_accession", ""),
        "parents_are_exact_published_strains": case.get(
            "parents_are_exact_published_strains", False),
        "expected_pattern": case.get("expect_pattern"),
        "expected_discordant_segment": case.get("expect_discordant"),
        "observed_pattern": observed_pattern,
        "observed_discordant_segment": observed_segment,
        "observed_verdict": row.get("verdict"),
        # DERIVED, not read: score_candidate emits no gate0_passed_structure key, so
        # reading one gave a uniformly null column that looked like a recorded
        # diagnostic. The early junk return is taken on exactly `not
        # gate0.passed_structure`, so the verdict is an exact stand-in for it.
        "gate0_passed_structure": row.get("verdict") != "junk_sequence",
        "gate0_failures": row.get("gate0_failures", ""),
        "gate0_best_kmer_identity": row.get("gate0_best_kmer_identity"),
        "gate1_implicated_species": row.get("gate1_implicated_species"),
        "gate3_status": row.get("gate3_status"),
        "rpi_generation": row.get("rpi_generation"),
        "rpi_retention": row.get("rpi_retention"),
        "rpi_e_total": row.get("rpi_e_total"),
        "rpi_k_stratum": row.get("rpi_k_stratum"),
        "as_expected": bool(passed),
        "localized_correct_segment": bool(localized),
        "exact_pattern_match": bool(exact),
        "is_known_blind_spot": bool(case.get("blind_spot")),
    }


def summarise(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    def subset(klass):
        return [r for r in results if r["class"] == klass]

    detectable = subset("A_published_inter_species") + subset("D_synthetic_inter_species")
    negatives = subset("C_true_negative")
    blind = subset("B_published_intra_species")
    adversarial = subset("E_adversarial")

    detected = [r for r in detectable if r["localized_correct_segment"]]
    exact = [r for r in detectable if r["exact_pattern_match"]]
    unresolved = [r for r in detectable
                  if r["localized_correct_segment"] and not r["exact_pattern_match"]]
    false_positives = [r for r in negatives if not r["as_expected"]]

    return {
        "sensitivity_on_inter_species_reassortants": {
            "n": len(detectable), "correct": len(detected),
            "rate": (len(detected) / len(detectable)) if detectable else float("nan"),
            "definition": "PRIMARY. correct = gate1 named the exact segment that was "
                          "exchanged.",
            "also_exact_pattern_match": len(exact),
            "localized_but_parents_unresolvable": [r["case_id"] for r in unresolved],
            "unresolvable_note": "gate1 named the right segment but returned "
                                 "'undetectable_by_design' because the two candidate "
                                 "parent species are too close to tell apart. The "
                                 "exchange was found; the parents were not. Note that "
                                 "this is also the VERDICT the user reads.",
        },
        "specificity_on_true_negatives": {
            "n": len(negatives), "correct": len(negatives) - len(false_positives),
            "false_positives": [r["case_id"] for r in false_positives],
            "rate": ((len(negatives) - len(false_positives)) / len(negatives))
                    if negatives else float("nan"),
        },
        "known_blind_spot_intra_species": {
            "n": len(blind),
            "silent_as_predicted": sum(1 for r in blind if r["as_expected"]),
            "note": "These are REAL published reassortants that gate1 cannot see, "
                    "because both parents assign to the same reference species. This is "
                    "a measured limitation, not a passing test. gate4 phylogenetics is "
                    "the layer that would catch them.",
        },
        "adversarial": {
            "n": len(adversarial),
            "handled_as_expected": sum(1 for r in adversarial if r["as_expected"]),
            "verdicts": {r["case_id"]: r["observed_verdict"] for r in adversarial},
        },
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--reference-input", default="all_hanta")
    parser.add_argument("--json", help="write the full report to this path")
    parser.add_argument("--with-phylogeny", action="store_true",
                        help="run gate4 (needs mafft and iqtree on PATH)")
    parser.add_argument("--output", default="output/validation")
    args = parser.parse_args(argv)

    root = Path(args.reference_input)
    if not (root / "ncbi_dataset" / "data" / "data_report.jsonl").exists():
        print(f"ERROR: {root}/ NCBI package not found.", file=sys.stderr)
        return 2

    print(f"Loading reference from {root}/ ...", flush=True)
    reference = rd.load_reference_data(root)
    precomputed = rd._build_precomputed(reference)
    builder = PanelBuilder(reference)

    panel, skipped = build_panel(builder)
    print(f"Panel: {len(panel)} candidates"
          + (f" ({len(skipped)} skipped)" if skipped else ""), flush=True)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for case in panel:
        results.append(score_case(case, reference, precomputed,
                                  args.with_phylogeny, output_dir))
        print(".", end="", flush=True)
    print(flush=True)

    summary = summarise(results)

    width = max(len(r["case_id"]) for r in results) + 2
    current = None
    for row in sorted(results, key=lambda r: (r["class"], r["case_id"])):
        if row["class"] != current:
            current = row["class"]
            print(f"\n--- {current} ---")
        mark = "PASS" if row["as_expected"] else "FAIL"
        if row["is_known_blind_spot"]:
            mark = "BLIND" if row["as_expected"] else "FAIL"
        print(f"  [{mark:5s}] {row['case_id']:<{width}} "
              f"pattern={row['observed_pattern'] or '-':<22} "
              f"seg={row['observed_discordant_segment'] or '-':<3} "
              f"verdict={row['observed_verdict']}")

    print("\n=== SUMMARY ===")
    print(json.dumps(rd._json_safe(summary), indent=2, default=str, allow_nan=False))

    if skipped:
        print("\nSkipped (sequences not in this dataset):")
        for item in skipped:
            print("  -", item)

    if args.json:
        payload = {
            "gate4_exercised": args.with_phylogeny,
            "rpi_weights_version": rd.RPI_WEIGHTS_VERSION,
            "summary": summary,
            "results": results,
            "skipped": skipped,
        }
        # Same discipline as write_rpi_report, and it was missing here: Python emits
        # bare `NaN` tokens by default, which jq and JavaScript's JSON.parse reject
        # while Python's own json.loads accepts them. Not hypothetical -- summarise()
        # returns float("nan") for any class rate whose class is empty, which happens
        # on any reduced --reference-input. Non-finite values serialise as null; null
        # is "undefined", 0.0 would be a claim.
        Path(args.json).write_text(
            json.dumps(rd._json_safe(payload), indent=2, default=str, allow_nan=False),
            encoding="utf-8")
        print(f"\nFull report: {args.json}")

    sens = summary["sensitivity_on_inter_species_reassortants"]
    spec = summary["specificity_on_true_negatives"]
    adv = summary["adversarial"]
    blind = summary["known_blind_spot_intra_species"]

    # A class that failed to build must FAIL, not pass quietly. The criterion used to
    # read only sensitivity and specificity, so a run whose adversarial class was
    # skipped -- which happens whenever the reference lacks a full-length SEOV triplet,
    # since the None-guard above degrades to `skipped` -- reported PASS having tested
    # nothing adversarial at all. An absent measurement is not a passing measurement.
    present = {klass: sum(1 for r in results if r["class"] == klass)
               for klass in ("A_published_inter_species", "B_published_intra_species",
                             "C_true_negative", "D_synthetic_inter_species",
                             "E_adversarial")}
    empty = sorted(klass for klass, count in present.items() if count == 0)

    failures = []
    if empty:
        failures.append(f"class(es) built no candidates: {', '.join(empty)}")
    if sens["correct"] != sens["n"]:
        failures.append(f"sensitivity {sens['correct']}/{sens['n']}")
    if spec["false_positives"]:
        failures.append(f"false positives: {', '.join(map(str, spec['false_positives']))}")
    if adv["handled_as_expected"] != adv["n"]:
        failures.append(f"adversarial {adv['handled_as_expected']}/{adv['n']}")

    # Deliberately NOT a failure. The blind-spot class asserts that gate1 stays silent
    # on intra-species reassortment; if it ever stops being silent that is an
    # IMPROVEMENT, and failing the run for it would train whoever sees it to suppress
    # the improvement. It still has to be impossible to miss.
    if blind["silent_as_predicted"] != blind["n"]:
        print(f"\nNOTICE: {blind['n'] - blind['silent_as_predicted']} of {blind['n']} "
              f"intra-species cases were NOT silent. That is a change in gate1's reach "
              f"-- investigate before assuming it is a regression; it may be a gain.")

    if failures:
        print("\nRESULT: FAIL")
        for item in failures:
            print("  -", item)
        return 1
    print(f"\nRESULT: PASS  (classes: "
          + ", ".join(f"{k.split('_')[0]}={v}" for k, v in sorted(present.items())) + ")")
    return 0


if __name__ == "__main__":
    sys.exit(main())
