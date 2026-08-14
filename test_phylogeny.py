#!/usr/bin/env python3
"""Regression tests for phylogeny.py -- moved from test_all_hanta.py on 2026-08-13 when the
--phylogenetics code was split out of all_hanta.py into its own module. See
docs/superpowers/specs/2026-08-13-three-stage-pipeline-split-design.md for why.

stdlib unittest, matching the repo's house style. Run with:

    python3 -m unittest test_phylogeny -v

Fakes Bio.Phylo tree objects directly for the bootstrap/classify tests; fakes mafft and iqtree
as tiny throwaway scripts for the end-to-end test, since neither binary is installed in this
environment -- skipped automatically if Bio.Phylo/dendropy aren't importable.
"""

import random
import tempfile
import unittest
from io import StringIO
from pathlib import Path

import all_hanta as ah
import phylogeny as ph

try:
    import Bio.Phylo
    import dendropy  # noqa: F401
    import dendropy.calculate.treecompare  # noqa: F401
    _HAVE_PHYLO_DEPS = True
except ImportError:
    _HAVE_PHYLO_DEPS = False


class FakeConfidence:
    def __init__(self, confidence):
        self.confidence = confidence


class FakeTree:
    """Minimal stand-in for a Bio.Phylo tree: only what bootstrap_support_for_pair uses."""

    def __init__(self, support):
        self._support = support

    def find_clades(self, name):
        yield name

    def common_ancestor(self, clades):
        return FakeConfidence(self._support)


class ClassifyReassortmentCallTest(unittest.TestCase):
    """The generalized classifier must find an L-discordant (or S-discordant) pattern
    that the old M-only version structurally could not -- this is Task 1's headline
    acceptance criterion.
    """

    def _l_discordant_row(self):
        # S and M trees agree (low sm separation); L disagrees with both.
        return {
            "label": "I0000",
            "tree_s_nearest": "I0001", "tree_m_nearest": "I0001", "tree_l_nearest": "I0002",
            "tree_s_rank_offset": 0.1, "tree_m_rank_offset": 0.1, "tree_l_rank_offset": 0.9,
            "tree_sm_rank_separation": 0.05, "tree_sl_rank_separation": 0.9,
            "tree_ml_rank_separation": 0.9,
        }

    def test_l_discordant_supported(self):
        row = self._l_discordant_row()
        trees = {"S": None, "M": None, "L": FakeTree(99.0)}
        call, support, discordant, source = ph.classify_reassortment_call(row, trees, min_bootstrap_support=95.0)
        self.assertEqual(call, "supported_reassortment_candidate")
        self.assertEqual(discordant, "L")
        self.assertEqual(support, 99.0)
        self.assertEqual(source, "mrca")
        # The M-permissive prior from PMC6153690 must not fire for an L-discordant call.
        self.assertEqual(int(discordant == "M"), 0)

    def test_l_discordant_weak_bootstrap(self):
        row = self._l_discordant_row()
        trees = {"S": None, "M": None, "L": FakeTree(50.0)}
        call, support, discordant, source = ph.classify_reassortment_call(row, trees, min_bootstrap_support=95.0)
        self.assertEqual(call, "weak_signal_low_bootstrap_support")
        self.assertEqual(discordant, "L")
        self.assertEqual(support, 50.0)
        self.assertEqual(source, "mrca")

    def test_no_conflict_when_nothing_stands_out(self):
        row = {
            "label": "I0000",
            "tree_s_nearest": "I0001", "tree_m_nearest": "I0001", "tree_l_nearest": "I0002",
            "tree_s_rank_offset": 0.1, "tree_m_rank_offset": 0.1, "tree_l_rank_offset": 0.2,
            "tree_sm_rank_separation": 0.05, "tree_sl_rank_separation": 0.1,
            "tree_ml_rank_separation": 0.1,
        }
        trees = {"S": None, "M": None, "L": FakeTree(99.0)}
        call, support, discordant, source = ph.classify_reassortment_call(row, trees, min_bootstrap_support=95.0)
        self.assertEqual(call, "no_conflict")
        self.assertEqual(discordant, "")
        self.assertEqual(support, 0.0)
        self.assertEqual(source, "")

    def test_m_clade_bootstrap_support_is_none_when_m_not_discordant(self):
        """m_clade_bootstrap_support must stay None (not 0.0) when M isn't the discordant
        segment, so a genuine 0% support elsewhere is never confused with "not applicable".
        """
        row = self._l_discordant_row()
        trees = {"S": None, "M": None, "L": FakeTree(99.0)}
        _call, support, discordant, _source = ph.classify_reassortment_call(row, trees, min_bootstrap_support=95.0)
        m_clade_bootstrap_support = support if discordant == "M" else None
        discordant_clade_bootstrap_support = support if discordant else None
        self.assertIsNone(m_clade_bootstrap_support)
        self.assertEqual(discordant_clade_bootstrap_support, 99.0)


@unittest.skipUnless(_HAVE_PHYLO_DEPS, "biopython not installed")
class BootstrapSupportRootFallbackTest(unittest.TestCase):
    """Task 6 fix 1: when a pair's MRCA is the tree's root (unsupported by construction --
    IQ-TREE never writes a value there), bootstrap_support_for_pair must walk toward the
    tips for a real support value instead of silently reporting 0.0, which is
    indistinguishable from a genuine 0% UFBoot result.
    """

    def test_mrca_with_confidence_is_used_directly(self):
        tree = Bio.Phylo.read(StringIO("(((A:1,B:1)95:1,C:1)80:1,D:1);"), "newick")
        support, source = ph.bootstrap_support_for_pair(tree, "A", "B")
        self.assertEqual((support, source), (95.0, "mrca"))

    def test_unsupported_root_walks_toward_tips_for_nearest_support(self):
        # A vs D: MRCA is the unrooted tree's root (no confidence). Walking toward A
        # passes through the confidence=80 clade before reaching A itself; walking toward
        # D hits the tip directly with nothing in between.
        tree = Bio.Phylo.read(StringIO("(((A:1,B:1)95:1,C:1)80:1,D:1);"), "newick")
        support, source = ph.bootstrap_support_for_pair(tree, "A", "D")
        self.assertEqual((support, source), (80.0, "nearest_supported_ancestor"))

    def test_takes_the_more_conservative_side_when_both_are_found(self):
        # A vs E: walking toward A finds 90 one step down; walking toward E finds 60 one
        # step down. Must report the LOWER of the two, not whichever is found first.
        tree = Bio.Phylo.read(StringIO("((A:1,B:1)90:1,(D:1,E:1)60:1);"), "newick")
        support, source = ph.bootstrap_support_for_pair(tree, "A", "E")
        self.assertEqual((support, source), (60.0, "nearest_supported_ancestor"))

    def test_totally_unsupported_star_tree_reports_honestly(self):
        tree = Bio.Phylo.read(StringIO("(A:1,B:1,C:1);"), "newick")
        support, source = ph.bootstrap_support_for_pair(tree, "A", "B")
        self.assertEqual((support, source), (0.0, "unsupported_root"))


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


def _make_candidates(species_list, n_per_species=5, seed=1):
    rng = random.Random(seed)
    bases = "ACGT"

    def rand_seq(n):
        return "".join(rng.choice(bases) for _ in range(n))

    candidates = []
    for sp_i, species in enumerate(species_list):
        for i in range(n_per_species):
            key = f"{species}|{i}"
            s_acc, m_acc, l_acc = f"ACC{sp_i}{i}S", f"ACC{sp_i}{i}M", f"ACC{sp_i}{i}L"
            isolate = _make_isolate(key, species, s_acc, m_acc, l_acc)
            sequences = {"S": rand_seq(150), "M": rand_seq(150), "L": rand_seq(150)}
            candidates.append(ph.PhyloCandidate(isolate=isolate, tags=["complete_triplet"], sequences=sequences))
    return candidates


def _cds_index_for(s_acc, m_acc, l_acc, seq_len=150):
    """Minimal cds_index entries that select_primary_cds resolves to S/M/L via
    classify_protein_role: 'nucleocapsid' -> N (S), 'glyco' -> GPC (M), 'polymerase' ->
    RdRp (L).
    """
    rng = random.Random(hash((s_acc, m_acc, l_acc)) % (2**32))
    bases = "ACGT"

    def rand_seq(n):
        return "".join(rng.choice(bases) for _ in range(n))

    return {
        ah.accession_base(s_acc): [("nucleocapsid protein", rand_seq(seq_len))],
        ah.accession_base(m_acc): [("glycoprotein precursor", rand_seq(seq_len))],
        ah.accession_base(l_acc): [("RNA polymerase", rand_seq(seq_len))],
    }


class CollectPhyloCandidatesScopeTest(unittest.TestCase):
    """Task 3: unclassified isolates are excluded at scope='species' (unchanged) but
    admitted at scope='phylogroup' when observed host order resolves a phylogroup whose
    WHOLE order maps to one group (Eulipotyphla, Chiroptera) -- never for Rodentia, which
    spans three phylogroups, or for an isolate with no usable host order at all.
    """

    def _candidate(self, key, species_status, host_order):
        s_acc, m_acc, l_acc = f"{key}S", f"{key}M", f"{key}L"
        isolate = _make_isolate(key, "unclassified Orthohantavirus", s_acc, m_acc, l_acc,
                                species_status=species_status, host_order=host_order)
        cds_index = _cds_index_for(s_acc, m_acc, l_acc)
        return isolate, cds_index

    def _collect(self, isolates, cds_index, scope):
        return ph.collect_phylo_candidates(isolates, {}, cds_index, scope)

    def test_shrew_hosted_unclassified_admitted_only_at_phylogroup_scope(self):
        isolate, cds_index = self._candidate("shrew1", "unclassified", "Eulipotyphla")
        species_cands, species_excl = self._collect([isolate], cds_index, "species")
        self.assertEqual(species_cands, [])
        self.assertEqual(species_excl["species_not_classified"], 1)

        phylo_cands, phylo_excl = self._collect([isolate], cds_index, "phylogroup")
        self.assertEqual(len(phylo_cands), 1)
        self.assertEqual(dict(phylo_excl), {})
        basis = ph._assignment_basis(phylo_cands[0], "phylogroup")
        self.assertEqual(basis, "observed_host_order")
        self.assertEqual(ph._candidate_group_key(phylo_cands[0], "phylogroup"), "soricomorpha_borne")

    def test_bat_hosted_unclassified_admitted_as_chiroptera_borne(self):
        isolate, cds_index = self._candidate("bat1", "unclassified", "Chiroptera")
        phylo_cands, phylo_excl = self._collect([isolate], cds_index, "phylogroup")
        self.assertEqual(len(phylo_cands), 1)
        self.assertEqual(ph._candidate_group_key(phylo_cands[0], "phylogroup"), "chiroptera_borne")

    def test_rodent_hosted_unclassified_excluded_at_both_scopes(self):
        """Rodentia spans three phylogroups (Murinae/Arvicolinae/Sigmodontinae) -- order
        alone can't pick one, so this must NOT be silently guessed into any of them.
        """
        isolate, cds_index = self._candidate("rodent1", "unclassified", "Rodentia")
        _, species_excl = self._collect([isolate], cds_index, "species")
        _, phylo_excl = self._collect([isolate], cds_index, "phylogroup")
        self.assertEqual(species_excl["species_not_classified"], 1)
        self.assertEqual(phylo_excl["species_not_classified_and_no_resolvable_phylogroup"], 1)

    def test_no_host_order_unclassified_excluded_with_reason(self):
        isolate, cds_index = self._candidate("nohost1", "unclassified", "")
        _, phylo_excl = self._collect([isolate], cds_index, "phylogroup")
        self.assertEqual(phylo_excl["species_not_classified_and_no_resolvable_phylogroup"], 1)

    def test_classified_isolate_unaffected_by_scope(self):
        s_acc, m_acc, l_acc = "clsS", "clsM", "clsL"
        isolate = _make_isolate("cls1", "Orthohantavirus hantanense", s_acc, m_acc, l_acc,
                                species_status="classified", host_order="Rodentia")
        cds_index = _cds_index_for(s_acc, m_acc, l_acc)
        species_cands, _ = self._collect([isolate], cds_index, "species")
        phylo_cands, _ = self._collect([isolate], cds_index, "phylogroup")
        self.assertEqual(len(species_cands), 1)
        self.assertEqual(len(phylo_cands), 1)
        self.assertEqual(ph._assignment_basis(phylo_cands[0], "phylogroup"), "curated_reservoir")


class GroupCandidatesTest(unittest.TestCase):
    """group_candidates(scope='species') must be a drop-in replacement for the old
    group_candidates_by_species, and scope='phylogroup' must pool candidates that share a
    reservoir phylogroup even when their species differ -- the whole point of Task 2.
    """

    def setUp(self):
        # hantanense (HTNV) and dobravaense (DOBV) are both Apodemus-borne -> same
        # phylogroup (murinae_borne) via PHYLOGROUP, but different species.
        self.species_list = ["Orthohantavirus hantanense", "Orthohantavirus dobravaense"]
        self.candidates = _make_candidates(self.species_list)

    def test_species_scope_matches_old_wrapper(self):
        new, new_excl = ph.group_candidates(self.candidates, 4, "species")
        old, old_excl = ph.group_candidates_by_species(self.candidates, 4)
        self.assertEqual(new.keys(), old.keys())
        for key in new:
            self.assertEqual([c.isolate.isolate_key for c in new[key]],
                             [c.isolate.isolate_key for c in old[key]])
        self.assertEqual(dict(new_excl), dict(old_excl))
        self.assertEqual(set(new.keys()), set(self.species_list))

    def test_phylogroup_scope_pools_same_phylogroup_species(self):
        grouped, exclusions = ph.group_candidates(self.candidates, 4, "phylogroup")
        self.assertEqual(set(grouped.keys()), {"murinae_borne"})
        pooled_species = {c.isolate.species_name for c in grouped["murinae_borne"]}
        self.assertEqual(pooled_species, set(self.species_list))
        self.assertEqual(len(grouped["murinae_borne"]), len(self.candidates))
        self.assertEqual(dict(exclusions), {})


@unittest.skipUnless(_HAVE_PHYLO_DEPS, "biopython/dendropy not installed")
class PhylogeneticsScopeEndToEndTest(unittest.TestCase):
    """Runs _run_phylogenetics_scope for real, against fake mafft/iqtree scripts, to
    exercise directory layout, --phylo-scope orchestration and --phylo-outgroup routing
    end to end. The fake binaries produce structurally valid but biologically meaningless
    alignments/trees (padded sequences, a random caterpillar Newick with support values on
    internal nodes) -- enough for real Bio.Phylo/dendropy parsing to run, not a
    phylogenetics correctness check.
    """

    _FAKE_MAFFT = '''#!/usr/bin/env python3
import sys
path = sys.argv[-1]
seqs, label, buf = [], None, []
with open(path) as f:
    for line in f:
        line = line.rstrip("\\n")
        if line.startswith(">"):
            if label is not None:
                seqs.append((label, "".join(buf)))
            label, buf = line[1:], []
        else:
            buf.append(line)
    if label is not None:
        seqs.append((label, "".join(buf)))
maxlen = max(len(s) for _, s in seqs)
for label, s in seqs:
    print(f">{label}")
    print(s.ljust(maxlen, "-"))
'''

    _FAKE_IQTREE = '''#!/usr/bin/env python3
import random, sys
argv = sys.argv[1:]
def get_opt(flag):
    return argv[argv.index(flag) + 1] if flag in argv else None
alignment, prefix, outgroup, seed = get_opt("-s"), get_opt("--prefix"), get_opt("-o"), get_opt("--seed")
random.seed(int(seed) if seed else 0)
labels = []
with open(alignment) as f:
    for line in f:
        if line.startswith(">"):
            labels.append(line[1:].strip())
nodes = list(labels)
random.shuffle(nodes)
if outgroup and outgroup in nodes:
    nodes.remove(outgroup)
    nodes.insert(0, outgroup)
cur = f"{nodes[0]}:{random.random():.4f}"
for label in nodes[1:-1]:
    cur = f"({cur},{label}:{random.random():.4f}){random.randint(50, 100)}:{random.random():.4f}"
cur = f"({cur},{nodes[-1]}:{random.random():.4f})"
with open(f"{prefix}.treefile", "w") as f:
    f.write(cur + ";\\n")
'''

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        tmp_path = Path(self.tmpdir.name)
        self.mafft_bin = tmp_path / "fake_mafft.py"
        self.iqtree_bin = tmp_path / "fake_iqtree.py"
        self.mafft_bin.write_text(self._FAKE_MAFFT)
        self.iqtree_bin.write_text(self._FAKE_IQTREE)
        self.mafft_bin.chmod(0o755)
        self.iqtree_bin.chmod(0o755)
        self.output_dir = tmp_path / "output"
        self.output_dir.mkdir()

        self.species_list = ["Orthohantavirus hantanense", "Orthohantavirus dobravaense"]
        self.candidates = _make_candidates(self.species_list)

    def tearDown(self):
        self.tmpdir.cleanup()

    def _args(self, scope, outgroup=None):
        import argparse
        return argparse.Namespace(
            phylo_min_taxa=4, bootstrap=100, min_bootstrap_support=95.0, phylo_seed=42,
            phylo_threads=1, phylo_scope=scope, phylo_outgroup=outgroup,
        )

    def test_both_scopes_each_report_every_candidate_once(self):
        rows = []
        for scope in ("species", "phylogroup"):
            rows.extend(ph._run_phylogenetics_scope(
                scope, self.candidates, self.output_dir, self._args(scope),
                str(self.mafft_bin), str(self.iqtree_bin)))

        self.assertTrue(all(r["phylo_status"] == "analyzed" for r in rows))
        species_rows = [r for r in rows if r["analysis_scope"] == "species"]
        phylogroup_rows = [r for r in rows if r["analysis_scope"] == "phylogroup"]
        candidate_keys = sorted(c.isolate.isolate_key for c in self.candidates)
        self.assertEqual(sorted(r["isolate_key"] for r in species_rows), candidate_keys)
        self.assertEqual(sorted(r["isolate_key"] for r in phylogroup_rows), candidate_keys)

        # species_name is always the isolate's OWN species, even under phylogroup scope,
        # where the tree itself spans both species -- this is the field that lets a
        # consumer recover which phylogroup a row's species belongs to.
        phylogroup_species = {r["species_name"] for r in phylogroup_rows}
        self.assertEqual(phylogroup_species, set(self.species_list))

        # Species-scope output must live under phylogenetics/, unchanged from before
        # scope existed; phylogroup output must live somewhere else entirely.
        self.assertTrue((self.output_dir / "phylogenetics" / "orthohantavirus_hantanense").is_dir())
        self.assertTrue((self.output_dir / "phylogenetics" / "orthohantavirus_dobravaense").is_dir())
        self.assertTrue((self.output_dir / "phylogenetics_phylogroup" / "murinae_borne").is_dir())

    def test_outgroup_applied_only_where_present_and_logged_when_absent(self):
        outgroup_key = self.candidates[0].isolate.isolate_key  # a hantanense isolate
        with self.assertLogs(level="INFO") as captured:
            ph._run_phylogenetics_scope(
                "species", self.candidates, self.output_dir, self._args("species", outgroup_key),
                str(self.mafft_bin), str(self.iqtree_bin))
        messages = "\n".join(captured.output)
        # dobravaense's tree doesn't contain the outgroup isolate -- must be logged, not
        # silently ignored.
        self.assertIn("not a member of group Orthohantavirus dobravaense", messages)
        # hantanense's tree DOES contain it -- no such message for that group.
        self.assertNotIn("not a member of group Orthohantavirus hantanense", messages)


class PhyloCallByIsolateTest(unittest.TestCase):
    """phylo_call_by_isolate reduces an isolate's rows (up to one per --phylo-scope) to
    the single most informative call -- a species-scope 'no_conflict' must never mask a
    phylogroup-scope 'supported_reassortment_candidate' for the same isolate, since the
    two scopes answer different questions and can only ever add evidence, not contradict.
    """

    def test_picks_most_informative_call_across_scopes(self):
        rows = [
            {"isolate_key": "iso1", "reassortment_call": "no_conflict"},
            {"isolate_key": "iso1", "reassortment_call": "supported_reassortment_candidate"},
        ]
        self.assertEqual(ph.phylo_call_by_isolate(rows)["iso1"], "supported_reassortment_candidate")

    def test_not_analyzed_never_overrides_an_analyzed_call(self):
        rows = [
            {"isolate_key": "iso1", "reassortment_call": "weak_signal_low_bootstrap_support"},
            {"isolate_key": "iso1", "reassortment_call": "not_analyzed"},
        ]
        self.assertEqual(ph.phylo_call_by_isolate(rows)["iso1"], "weak_signal_low_bootstrap_support")

    def test_order_of_rows_does_not_matter(self):
        rows = [
            {"isolate_key": "iso1", "reassortment_call": "supported_reassortment_candidate"},
            {"isolate_key": "iso1", "reassortment_call": "no_conflict"},
        ]
        self.assertEqual(ph.phylo_call_by_isolate(rows), ph.phylo_call_by_isolate(list(reversed(rows))))


if __name__ == "__main__":
    unittest.main()
