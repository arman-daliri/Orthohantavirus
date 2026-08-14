#!/usr/bin/env python3
"""Regression tests for all_hanta.py (Stage 1) functionality that remains in this file
after the phylogenetics code moved to phylogeny.py and its tests to test_phylogeny.py.

stdlib unittest, not pytest -- this repo has no test dependency installed. Run with:

    python3 -m unittest test_all_hanta -v

Three test classes, all fake-free (no faked Bio.Phylo/mafft/iqtree -- that lives in
test_phylogeny.py now):

- SegmentRankMetricsBackwardCompatTest: phylogeny.py's _segment_rank_metrics must reduce
  to all_hanta.py's own _sl_m_rank_metrics's exact values for the (M offset, S/L
  separation) pair it replaces -- not just a similar pattern.
- SpeciesEffortStatsTest: _species_effort_stats's effort_* covariates are computed per
  species from ALL accessions/isolates (not just the ML-usable subset), and count
  distinct years/countries correctly.
- SurveillancePriorityReassortmentTermTest: only surveillance_priority's
  'supported_reassortment_candidate' phylo_call may contribute priority_reassortment_term
  -- 'weak_signal_low_bootstrap_support' and 'no_conflict' must not.
"""

import unittest

import numpy as np

import all_hanta as ah
import phylogeny as ph


class SegmentRankMetricsBackwardCompatTest(unittest.TestCase):
    """_segment_rank_metrics must reduce to _sl_m_rank_metrics's exact values for the
    (M offset, S/L separation) pair it replaces -- not just a similar pattern. Task 1's
    acceptance criteria required species-scope output to stay unchanged, and this is the
    computation that guarantees it.
    """

    def test_matches_old_formula_on_random_orderings(self):
        rng = np.random.default_rng(12345)
        n = 8
        span = n - 1
        order = {}
        for segment in ("S", "M", "L"):
            dist = rng.random((n, n))
            dist = (dist + dist.T) / 2
            np.fill_diagonal(dist, np.inf)
            order[segment] = np.argsort(dist, axis=1)
        position = ah._rank_positions(order)

        for i in range(n):
            s_nn, m_nn, l_nn = order["S"][i, 0], order["M"][i, 0], order["L"][i, 0]
            old_sl, old_m = ah._sl_m_rank_metrics(position, span, i, s_nn, m_nn, l_nn)
            new = ph._segment_rank_metrics(position, span, i, {"S": s_nn, "M": m_nn, "L": l_nn})
            self.assertAlmostEqual(old_sl, new["tree_sl_rank_separation"], places=12)
            self.assertAlmostEqual(old_m, new["tree_m_rank_offset"], places=12)


def _make_isolate(key: str, species: str, s_acc: str, m_acc: str, l_acc: str,
                  species_status: str = "classified", host_order: str = "Rodentia") -> "ah.IsolateRecord":
    """Minimal but fully-populated IsolateRecord for phylogenetics-stage tests -- only the
    fields _analyze_group/collect_phylo_candidates actually read are given meaningful
    values; everything else is a valid placeholder of the right type.
    """
    return ah.IsolateRecord(
        isolate_key=key, species_name=species, species_taxid="999", species_status=species_status,
        biosample="", isolate=key, strain="", host_name="Test host", host_taxid="1",
        host_order=host_order, host_family="", host_role="reservoir_rodent",
        collection_date="", year="2020", location="", country="Testland", region="",
        grouping_tier="tier_2_species_isolate", grouping_confidence=1.0,
        s_accession=s_acc, m_accession=m_acc, l_accession=l_acc,
        s_completeness="complete", m_completeness="complete", l_completeness="complete",
        s_nt_length=200, m_nt_length=200, l_nt_length=200,
        s_protein_aa="M" * 60, m_protein_aa="M" * 60, l_protein_aa="M" * 60,
        concatenated_protein_aa="M" * 180,
        segments_present=["S", "M", "L"], coding_complete_segments=["S", "M", "L"],
        triplet_status="coding_complete_triplet", has_complete_triplet=True,
        is_lab_derived=False, is_vaccine_strain=False, has_nss_orf=False,
        metadata_conflict=False, taxonomic_incongruence=False, cross_species_partners=[],
        kmer_cluster_id="", s_protein_accessions=[], m_protein_accessions=[],
        l_protein_accessions=[], all_accessions=[s_acc, m_acc, l_acc], warnings=[],
    )


def _make_accession(accession: str, species: str, country: str, year: str) -> "ah.AccessionRecord":
    """Minimal but fully-populated AccessionRecord for _species_effort_stats tests."""
    return ah.AccessionRecord(
        accession=accession, sequence_hash="", species_name=species, species_taxid="999",
        species_status="classified", observed_organism=species, observed_taxid="999", ranks={},
        segment="S", segment_source="explicit", segment_confidence=1.0,
        completeness="complete", completeness_source="reported", is_coding_complete=True,
        nt_length=1500, biosample="", bioproject="", sra_accessions=[],
        isolate="", strain="", collection_date=year, year=year, location="", country=country,
        region="", host_name="Test host", host_taxid="1", host_order="Rodentia", host_family="",
        host_source="reported", host_role="reservoir_rodent", isolation_source="", tissue_source="",
        is_lab_host=False, lab_host="", is_vaccine_strain=False, protein_accessions=[],
        annotation_products=[], protein_aa="", protein_product="", has_nss_orf=False,
        group_key="", grouping_tier="", group_confidence=1.0, is_refseq=False,
        source_database="GenBank", warning_flags=[],
    )


class SpeciesEffortStatsTest(unittest.TestCase):
    """Task 4: effort_* covariates are computed per species from ALL accessions/isolates
    (not just the ML-usable subset), and count distinct years/countries correctly, so
    they measure surveillance intensity rather than accidentally re-deriving a biology
    signal from a filtered slice.
    """

    def setUp(self):
        # PUUV: 3 accessions, 2 countries, years 2001/2001/2010 -> span 9, 2 distinct years.
        self.accessions = [
            _make_accession("A1", "Orthohantavirus puumalaense", "Finland", "2001"),
            _make_accession("A2", "Orthohantavirus puumalaense", "Russia", "2001"),
            _make_accession("A3", "Orthohantavirus puumalaense", "Finland", "2010"),
            # PHV: single accession, blank year -- effort_first_year must stay "" (not a
            # crash, not a guessed 0).
            _make_accession("A4", "Orthohantavirus prospectense", "USA", ""),
        ]
        self.isolates = [
            _make_isolate("I1", "Orthohantavirus puumalaense", "A1", "A1", "A1"),
            _make_isolate("I2", "Orthohantavirus puumalaense", "A2", "A2", "A2"),
        ]
        # One isolate deliberately not a complete triplet, to check the triplet count.
        self.isolates[1].has_complete_triplet = False

    def test_counts_and_span_from_all_records(self):
        stats = ah._species_effort_stats(self.accessions, self.isolates)
        puuv = stats["Orthohantavirus puumalaense"]
        self.assertEqual(puuv["effort_n_accessions"], 3)
        self.assertEqual(puuv["effort_n_isolates"], 2)
        self.assertEqual(puuv["effort_n_countries"], 2)
        self.assertEqual(puuv["effort_n_years"], 2)
        self.assertEqual(puuv["effort_first_year"], 2001)
        self.assertEqual(puuv["effort_years_since_first"], 9)
        self.assertEqual(puuv["effort_n_complete_triplets"], 1)

    def test_blank_year_does_not_crash_and_stays_unset(self):
        stats = ah._species_effort_stats(self.accessions, self.isolates)
        phv = stats["Orthohantavirus prospectense"]
        self.assertEqual(phv["effort_first_year"], "")
        self.assertEqual(phv["effort_years_since_first"], 0)
        self.assertEqual(phv["effort_n_isolates"], 0)  # no isolate records for this species


class SurveillancePriorityReassortmentTermTest(unittest.TestCase):
    """Only 'supported_reassortment_candidate' may contribute priority_reassortment_term.
    'weak_signal_low_bootstrap_support' must not -- the whole point of the bootstrap gate
    is that an unsupported pattern is never reported as a finding, including in the risk
    score, which is exactly the mistake this task is not allowed to repeat.
    """

    def setUp(self):
        self.isolate = _make_isolate("iso1", "Orthohantavirus hantanense", "S1", "M1", "L1")
        self.risk = ah.species_risk_profile(self.isolate.species_name)

    def test_supported_call_contributes_the_weight(self):
        terms = ah.surveillance_priority(self.isolate, self.risk, None, "supported_reassortment_candidate")
        self.assertEqual(terms["priority_reassortment_term"], ah.SURVEILLANCE_WEIGHTS["supported_reassortment"])

    def test_weak_signal_contributes_nothing(self):
        terms = ah.surveillance_priority(self.isolate, self.risk, None, "weak_signal_low_bootstrap_support")
        self.assertEqual(terms["priority_reassortment_term"], 0.0)

    def test_no_conflict_contributes_nothing(self):
        terms = ah.surveillance_priority(self.isolate, self.risk, None, "no_conflict")
        self.assertEqual(terms["priority_reassortment_term"], 0.0)

    def test_no_phylo_call_contributes_nothing(self):
        terms = ah.surveillance_priority(self.isolate, self.risk, None)
        self.assertEqual(terms["priority_reassortment_term"], 0.0)

    def test_score_delta_equals_exactly_the_weight(self):
        without = ah.surveillance_priority(self.isolate, self.risk, None, None)
        with_support = ah.surveillance_priority(self.isolate, self.risk, None, "supported_reassortment_candidate")
        delta = with_support["surveillance_priority_score"] - without["surveillance_priority_score"]
        self.assertAlmostEqual(delta, ah.SURVEILLANCE_WEIGHTS["supported_reassortment"], places=6)


if __name__ == "__main__":
    unittest.main()
