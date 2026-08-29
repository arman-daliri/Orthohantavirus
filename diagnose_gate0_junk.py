"""Read-only diagnostic (roadmap Step 1): why REAL_ORTHO_sinnombreense_005
(isolate_key t2|Orthohantavirus sinnombreense|202003197) and
REAL_ORTHO_Orthohantavirus_006 (isolate_key t2|unclassified Orthohantavirus|msb89863)
hit gate0's junk_sequence verdict in the c_real_positive group.

Does not modify reassortment_detector.py or rpi_discrimination.py -- calls their public
gate0_check / _build_precomputed / load_reference_data only.
"""
import json
from pathlib import Path

import reassortment_detector as rd
import rpi_calibration as cal

ISOLATE_KEYS = [
    "t2|Orthohantavirus sinnombreense|202003197",
    "t2|unclassified Orthohantavirus|msb89863",
]


def main() -> None:
    reference = rd.load_reference_data(Path("all_hanta"))
    precomputed = rd._build_precomputed(reference)
    by_key = {i.isolate_key: i for i in reference.isolates}

    for key in ISOLATE_KEYS:
        isolate = by_key.get(key)
        print("=" * 100)
        print(f"isolate_key: {key}")
        if isolate is None:
            print("  NOT FOUND in reference.isolates")
            continue
        print(f"  species_name: {isolate.species_name}")
        print(f"  s_accession/m/l: {isolate.s_accession} / {isolate.m_accession} / {isolate.l_accession}")
        print(f"  s/m/l_completeness: {isolate.s_completeness} / {isolate.m_completeness} / {isolate.l_completeness}")
        print(f"  s/m/l_nt_length: {isolate.s_nt_length} / {isolate.m_nt_length} / {isolate.l_nt_length}")
        print(f"  triplet_status: {isolate.triplet_status}, has_complete_triplet: {isolate.has_complete_triplet}")
        print(f"  warnings: {isolate.warnings}")

        genomic = {
            seg: cal.genomic_sequence(reference, cal.segment_accession(isolate, seg))
            for seg in ("S", "M", "L")
        }
        for seg, seq in genomic.items():
            print(f"  {seg} genomic length actually loaded: {len(seq) if seq else 'EMPTY/None'}")

        candidate = rd.CandidateInput(
            candidate_id=key, s_genomic=genomic["S"], m_genomic=genomic["M"], l_genomic=genomic["L"])
        gate0 = rd.gate0_check(
            candidate, precomputed["protein_len_ranges"], precomputed["panhandle_consensus"],
            precomputed["protein_profiles"])
        print(f"  gate0.passed_structure: {gate0.passed_structure}")
        print(f"  gate0.passed_identity_floor: {gate0.passed_identity_floor}")
        print(f"  gate0.best_kmer_identity: {gate0.best_kmer_identity}")
        print("  gate0.failures:")
        for f in gate0.failures:
            print(f"    - {f}")
        print("  gate0.warnings:")
        for w in gate0.warnings:
            print(f"    - {w}")
        for seg, protein in gate0.orf_proteins.items():
            print(f"  {seg} longest_orf_protein length: {len(protein)}")


if __name__ == "__main__":
    main()
