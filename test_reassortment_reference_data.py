# test_reassortment_reference_data.py
import unittest
import reassortment_reference_data as rd


class StructuralAnchorDataTest(unittest.TestCase):
    def test_three_species_present_with_expected_basis(self):
        self.assertEqual(rd.STRUCTURAL_ANCHORS["Orthohantavirus andesense"]["basis"], "solved")
        self.assertEqual(rd.STRUCTURAL_ANCHORS["Orthohantavirus puumalaense"]["basis"], "solved")
        self.assertEqual(
            rd.STRUCTURAL_ANCHORS["Orthohantavirus hantanense"]["basis"],
            "homology_inferred_verified",
        )

    def test_andv_tmis_and_ph_switch_at_cited_literal_positions(self):
        seq = rd.STRUCTURAL_ANCHORS["Orthohantavirus andesense"]["sequence"]
        # 1-based positions from Guardado-Calvo et al. (PMC7572791), matching the raw
        # sequence directly -- ANDV is the numerical anchor, see spec section 15.
        self.assertEqual(seq[744], "Y")  # Y745
        self.assertEqual(seq[756], "E")  # E757
        self.assertEqual(seq[758], "D")  # D759
        self.assertEqual(seq[765], "W")  # W766
        self.assertEqual(seq[899], "F")  # F900

    def test_puuv_tmis_positions_are_offset_not_literal(self):
        anchor = rd.STRUCTURAL_ANCHORS["Orthohantavirus puumalaense"]
        seq = anchor["sequence"]
        # PUUV's own QOC68549.1 does NOT match Tscherne et al. 2025's literal integers
        # (confirmed empirically during design, see spec section 15) -- the table stores
        # the corrected, offset positions instead.
        self.assertEqual(seq[anchor["tmis"]["Y745"] - 1], "Y")
        self.assertEqual(seq[anchor["tmis"]["W766"] - 1], "W")
        self.assertEqual(seq[anchor["tmis"]["F900"] - 1], "F")
        self.assertEqual(seq[anchor["ph_switch"]["E757"] - 1], "E")
        self.assertEqual(seq[anchor["ph_switch"]["D759"] - 1], "D")

    def test_htnv_f900_equivalent_is_none_known_exception(self):
        anchor = rd.STRUCTURAL_ANCHORS["Orthohantavirus hantanense"]
        self.assertIsNone(anchor["tmis"]["F900"])
        self.assertIsNotNone(anchor["tmis"]["Y745"])
        self.assertIsNotNone(anchor["tmis"]["W766"])

    def test_waasa_motif_present_at_recorded_position_in_all_three(self):
        for species, anchor in rd.STRUCTURAL_ANCHORS.items():
            pos = anchor["waasa_position"]
            self.assertEqual(anchor["sequence"][pos - 1:pos - 1 + 5], "WAASA", species)

    def test_lineage_marker_table_ships_empty_by_design(self):
        self.assertEqual(rd.LINEAGE_MARKER_TABLE, {})

    def test_new_world_species_membership(self):
        self.assertIn("Orthohantavirus andesense", rd.NEW_WORLD_SPECIES)
        self.assertNotIn("Orthohantavirus hantanense", rd.NEW_WORLD_SPECIES)
        self.assertNotIn("Orthohantavirus puumalaense", rd.NEW_WORLD_SPECIES)


class FindAlignedResidueTest(unittest.TestCase):
    def test_identical_sequence_maps_to_same_position(self):
        seq = rd.STRUCTURAL_ANCHORS["Orthohantavirus andesense"]["sequence"]
        residue, aligned = rd.find_aligned_residue(seq, seq, 745)
        self.assertEqual(residue, "Y")
        self.assertTrue(aligned)

    def test_cross_species_alignment_recovers_the_conserved_residue(self):
        andv = rd.STRUCTURAL_ANCHORS["Orthohantavirus andesense"]["sequence"]
        puuv_anchor = rd.STRUCTURAL_ANCHORS["Orthohantavirus puumalaense"]
        residue, aligned = rd.find_aligned_residue(
            puuv_anchor["sequence"], andv, 745)
        self.assertEqual(residue, "Y")
        self.assertTrue(aligned)

    def test_short_unrelated_sequence_reports_not_aligned_or_wrong_residue(self):
        residue, aligned = rd.find_aligned_residue("MAAAAAAAAAA", "M" * 50, 745)
        self.assertFalse(aligned)


if __name__ == "__main__":
    unittest.main()
