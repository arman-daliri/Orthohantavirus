"""rpi_discrimination.py -- thin shim.

Re-exports this module's original public surface from rpi_core.py, the
single consolidated module the RPI toolchain now lives in (branch
05-module-consolidation). See notes_consolidation.md for exactly what
moved and why; see rpi_core.py's own banner comments for this module's
section within it.
"""
import sys

from rpi_core import (
    RPI_DISCRIMINATION_SEED, RANDOM_FAKES_TARGET,
    WITHIN_SPECIES_SHUFFLES_PER_SPECIES, WITHIN_SPECIES_MIN_ISOLATES,
    COUNTRY_CONTINENT, ORTHO_REASSORTMENT_TABLE_PATH,
    reserved_calibration_isolate_keys, usable_pool, _candidate_from_isolates,
    _label, build_group_a_random, build_group_b1_within_species, _year_int,
    find_sympatric_species_pairs, build_group_b2_sympatric,
    load_supported_reassortment_candidates,
    load_published_inter_species_crosses, build_group_c_real_positives,
    _score_args, score_cases, build_discrimination_panel, summarise_groups,
    _fixture_isolate, _fixture_reference, CountryContinentTest,
    ReservedCalibrationIsolateKeysTest, FindSympatricSpeciesPairsTest,
    CandidateFromIsolatesTest, BuildGroupARandomTest,
    BuildGroupB1WithinSpeciesTest, BuildGroupB2SympatricTest,
    LoadSupportedReassortmentCandidatesTest, SummariseGroupsTest
)
from rpi_core import _discrimination_main as main


if __name__ == "__main__":
    sys.exit(main())
