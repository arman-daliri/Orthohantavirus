"""reassortment_detector.py -- thin shim.

Re-exports this module's original public surface from rpi_core.py, the
single consolidated module the RPI toolchain now lives in (branch
05-module-consolidation). See notes_consolidation.md for exactly what
moved and why; see rpi_core.py's own banner comments for this module's
section within it.
"""
import sys

from rpi_core import (
    PUUV_GPC_SEQUENCE, ANDV_GPC_SEQUENCE, HTNV_GPC_SEQUENCE, STRUCTURAL_ANCHORS,
    LINEAGE_MARKER_TABLE, NEW_WORLD_SPECIES, _ALIGNER, _cached_aligned_pair,
    find_aligned_residue, RPI_WEIGHTS, _weights, smoothed_prior,
    constellation_prior, ReferenceData, load_reference_data,
    protein_lengths_by_segment, genomic_sequences_by_segment,
    protein_profiles_by_segment, _STOP_CODONS, _CODON_TABLE, CandidateInput,
    _reverse_complement, _translate, _longest_orf_span_with_terminus,
    _longest_orf_span, longest_orf_protein, empirical_protein_length_range,
    empirical_panhandle_consensus, Gate0Result, _cached_reference_kmer_profile,
    _kmer_identity_to_reference, gate0_check, SpeciesHit,
    nearest_species_by_kmer, Gate1Result, gate1_check, Gate2TrackAResult,
    gate2_track_a, Gate2TrackBResult, gate2_track_b, Gate2TrackCResult,
    gate2_track_c, Gate3Result, gate3_check, candidate_cds_from_orf,
    build_candidate_isolate_record, build_comparator_members, Gate4Result,
    gate4_check, combine_verdict, divergence_stratum, host_genus_of,
    host_rank_proximity, _year_span, ecological_opportunity,
    _isolates_with_usable_protein, nearest_reference_isolate,
    resolve_e_parent_isolates, RPI_RT_KCAL_PER_MOL, boltzmann_compatibility,
    RPI_PANHANDLE_WINDOW_NT, panhandle_window, cofold_dg, viennarna_provenance,
    HANTAVIRUS_GENOMIC_5P_TAG, HANTAVIRUS_GENOMIC_3P_TAG,
    has_genomic_panhandle_termini, HANTAVIRUS_TERMINAL_TAG_NT,
    RPI_N_FOOTPRINT_NT, RPI_N_FOOTPRINT_DEFAULT_NT, RPI_M_POL_MIN_VARIANCE,
    RPI_MPOL_CDS_WINDOW_NT, n_recognition_footprint, _beyond_tag,
    polymerase_compatibility, m_pol_variance_check, _wild_isolates_of_species,
    _wild_species_genomic_sequences, _wild_species_cds_windows,
    species_consensus_panhandle_window, species_consensus_cds_window,
    rpi_species_assignment, score_m_pack, score_m_pol, retention_factor,
    RPI_WEIGHTS_VERSION, RPI_K_STRATUM_PROVISIONAL, _rpi_decimals, compute_rpi,
    RPI_INVITRO_BENCHMARK, RPI_BENCHMARK_CAVEAT, _composition_distance,
    _prior_counts_excluding, retrospective_reproduction, _make_ecology_isolate,
    _eco_isolate, _rpi_sensitivity_panel, _score_panel, sensitivity_sweep,
    _tied_ranks, _spearman, _worst_rank_correlation,
    RPI_HOST_DECAY_NOT_PERTURBED, _perturbed_host_decay, RPI_RANK_MIN_POSITIVES,
    rank_consistency, RPI_E_NEAR_ZERO_THRESHOLD, veto_vs_continuous_e,
    read_candidate, read_manifest, _build_precomputed, score_candidate,
    parse_args, RPI_REPORT_FILENAME, _json_safe, write_rpi_report,
    LongestOrfProteinTest, EmpiricalDerivationTest, Gate0CheckTest,
    NearestSpeciesByKmerTest, Gate1CheckTest, Gate2TrackATest, Gate2TrackBTest,
    Gate2TrackCTest, Gate3CheckTest, _make_isolate,
    _make_reference_isolate_with_cds, CandidateCdsFromOrfTest,
    BuildComparatorMembersTest, _write_fake_binary, Gate4CheckEndToEndTest,
    CombineVerdictTest, ProteinLengthsBySegmentTest,
    ProteinProfilesBySegmentTest, LoadReferenceDataIntegrationTest,
    _write_fasta, ReadCandidateTest, ReadManifestTest, ParseArgsTest,
    EndToEndRealDataTest, StructuralAnchorDataTest, FindAlignedResidueTest,
    ConstellationPriorTest, DivergenceStratumTest, _parent_map,
    EcologicalOpportunityTest, NearestReferenceIsolateTest,
    ResolveEParentIsolatesTest,
    ResolveEParentIsolatesSkipsExpensiveResolutionTest,
    BoltzmannCompatibilityTest, MPackTauRescalingTest, CofoldDgTest,
    PanhandleWindowTest, PolymeraseCompatibilityTest,
    HasGenomicPanhandleTerminiTest, _genomic_isolate, _genomic_reference,
    _TEST_PANHANDLE_WINDOW_NT, _tagged_seq, SpeciesConsensusPanhandleWindowTest,
    RpiSpeciesAssignmentTest, ScoreMPackTest, _orf_seq,
    SpeciesConsensusCdsWindowTest, ScoreMPolTest, RetentionFactorTest,
    RpiAssemblyTest, VerdictUnchangedRegressionTest,
    RetrospectiveReproductionTest, VetoVsContinuousETest, _parse_version,
    _read_requirement_pins, _read_conda_pins,
    ScorePanelDoesNotMutateGlobalWeightsTest, EnvironmentReconciliationTest,
    JsonSafeNumpyScalarTest, _reject_json_constant, WriteRpiReportTest,
    SensitivityAndRankTest
)
from rpi_core import _detector_main as main


if __name__ == "__main__":
    sys.exit(main())
