# reassortment_detector.py
"""CLI entry point for the reassortment candidate detector.

Scores a brand-new candidate hantavirus genome (S+M+L genomic FASTA) against the
five-gate framework in docs/superpowers/specs/2026-08-12-reassortment-detector-design.md.
Imports all_hanta.py's ingest and phylogenetics machinery; makes no changes to it.
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Any, Dict, List

import all_hanta as ah
import detector_gates as dg
import detector_reference as dr


def read_candidate(candidate_id: str, s_path: Path, m_path: Path, l_path: Path) -> dg.CandidateInput:
    def read_one(path: Path) -> str:
        record = next(iter(ah.parse_fasta(Path(path)).values()), None)
        return record.sequence if record else ""

    return dg.CandidateInput(
        candidate_id=candidate_id,
        s_genomic=read_one(s_path), m_genomic=read_one(m_path), l_genomic=read_one(l_path),
    )


def read_manifest(path: Path) -> List[dg.CandidateInput]:
    candidates = []
    base_dir = Path(path).parent
    with Path(path).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            candidates.append(read_candidate(
                row["candidate_id"],
                base_dir / row["s_fasta"], base_dir / row["m_fasta"], base_dir / row["l_fasta"],
            ))
    return candidates


def _build_precomputed(reference: dr.ReferenceData) -> Dict[str, Any]:
    protein_lengths = dr.protein_lengths_by_segment(reference.isolates)
    genomic_sequences = dr.genomic_sequences_by_segment(reference.isolates, reference.genomic_index)
    return {
        "protein_len_ranges": {
            segment: dg.empirical_protein_length_range(lengths)
            for segment, lengths in protein_lengths.items()
        },
        "panhandle_consensus": {
            segment: dg.empirical_panhandle_consensus(sequences)
            for segment, sequences in genomic_sequences.items()
        },
        "protein_profiles": dr.protein_profiles_by_segment(reference.isolates),
    }


def score_candidate(
    candidate: dg.CandidateInput, reference: dr.ReferenceData, precomputed: Dict[str, Any],
    args: argparse.Namespace,
) -> Dict[str, Any]:
    gate0 = dg.gate0_check(
        candidate, precomputed["protein_len_ranges"], precomputed["panhandle_consensus"],
        precomputed["protein_profiles"])

    if not gate0.passed_structure:
        return {
            "candidate_id": candidate.candidate_id, "verdict": "junk_sequence",
            "gate0_failures": ";".join(gate0.failures), "gate0_warnings": ";".join(gate0.warnings),
        }

    gate1 = dg.gate1_check(gate0.orf_proteins, precomputed["protein_profiles"])
    gate2a = dg.gate2_track_a(gate0.orf_proteins.get("M", ""), gate1.implicated_species)
    gate2b = dg.gate2_track_b(gate0.orf_proteins, gate1.implicated_species)
    gate2c = dg.gate2_track_c(gate1)
    gate3 = dg.gate3_check(gate1.implicated_species, reference.isolates)

    candidate_cds = {
        "S": dg.candidate_cds_from_orf(candidate.s_genomic),
        "M": dg.candidate_cds_from_orf(candidate.m_genomic),
        "L": dg.candidate_cds_from_orf(candidate.l_genomic),
    }
    gate4 = dg.gate4_check(
        candidate, candidate_cds, gate1.implicated_species, reference, Path(args.output),
        args.mafft_bin, args.iqtree_bin, phylo_min_taxa=args.phylo_min_taxa,
        bootstrap=args.bootstrap, min_bootstrap_support=args.min_bootstrap_support,
        phylo_seed=args.phylo_seed, phylo_threads=args.phylo_threads,
    )

    verdict = dg.combine_verdict(gate0, gate1, gate2a, gate2b, gate2c, gate3, gate4)

    return {
        "candidate_id": candidate.candidate_id,
        "verdict": verdict,
        "gate0_warnings": ";".join(gate0.warnings),
        "gate0_best_kmer_identity": round(gate0.best_kmer_identity, 4),
        "gate1_pattern": gate1.pattern,
        "gate1_discordant_segment": gate1.discordant_segment,
        "gate1_implicated_species": ";".join(gate1.implicated_species),
        "gate2_waasa_intact": gate2a.waasa_intact,
        "gate2_tmis_conserved": gate2a.tmis_conserved,
        "gate2_ph_switch_conserved": gate2a.ph_switch_conserved,
        "gate2_structural_check_not_available": gate2a.structural_check_not_available,
        "gate2_lineage_marker_match": gate2b.lineage_marker_match,
        "gate2_matched_species": gate2b.matched_species,
        "gate2_sl_mismatch_high_risk": gate2c.sl_mismatch_high_risk,
        "gate2_receptor_class_note": gate2c.receptor_class_note,
        "gate3_status": gate3.status,
        "gate3_detail": gate3.detail,
        "gate4_phylo_status": gate4.phylo_status,
        "gate4_reassortment_call": gate4.reassortment_call,
        "gate4_discordant_segment": gate4.discordant_segment,
        "gate4_bootstrap_support": gate4.bootstrap_support,
    }


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-input", required=True,
                        help="Package root containing all_hanta/ncbi_dataset/data")
    parser.add_argument("--output", required=True, help="Output directory")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--candidate", help="candidate_id:s_fasta,m_fasta,l_fasta")
    group.add_argument("--candidates-manifest", help="CSV: candidate_id,s_fasta,m_fasta,l_fasta")
    parser.add_argument("--mafft-bin", default=None)
    parser.add_argument("--iqtree-bin", default=None)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--min-bootstrap-support", type=float, default=95.0)
    parser.add_argument("--phylo-min-taxa", type=int, default=4)
    parser.add_argument("--phylo-seed", type=int, default=42)
    parser.add_argument("--phylo-threads", type=int, default=1)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    ah.configure_logging(Path(args.output) / "reassortment_detector.log", verbose=args.verbose)

    if args.candidates_manifest:
        candidates = read_manifest(Path(args.candidates_manifest))
    else:
        candidate_id, rest = args.candidate.split(":", 1)
        s_fasta, m_fasta, l_fasta = rest.split(",")
        candidates = [read_candidate(candidate_id, Path(s_fasta), Path(m_fasta), Path(l_fasta))]

    reference = dr.load_reference_data(Path(args.reference_input))
    precomputed = _build_precomputed(reference)

    mafft_bin = ah.resolve_binary(args.mafft_bin, "mafft")
    iqtree_bin = ah.resolve_binary(args.iqtree_bin, "iqtree")
    args.mafft_bin, args.iqtree_bin = mafft_bin, iqtree_bin

    rows = [score_candidate(c, reference, precomputed, args) for c in candidates]

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / "candidate_reassortment_table.csv"
    fieldnames = sorted({key for row in rows for key in row})
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(f"{row['candidate_id']}: {row['verdict']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
