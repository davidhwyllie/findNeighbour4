from unittest import TestCase
from unittest.mock import MagicMock

from findn.pca_covidtest_inc import VariantMatrix

class DummyClass:
    pass

class TestGetPositionCounts(TestCase):
    sample1 = {"A": [1,2], "M": [3], "invalid": 0}
    sample2 = {"A": [2], "M": [3, 4], "invalid": 0}
    def test_one_guid_correct_counts(self):
        PERSIST_object = DummyClass()
        PERSIST_object.refcompressedsequence_read = MagicMock(side_effect=[self.sample1])
        analysed_samples, vmodel, mmodel = VariantMatrix.get_position_counts({"s1"}, 100, PERSIST_object)
        self.assertEqual(analysed_samples, {"s1"})
        expected_vmodel = {1: 1, 2: 1}
        expected_mmodel = {3: 1}
        self.assertEqual(vmodel, expected_vmodel)
        self.assertEqual(mmodel, expected_mmodel)

    def test_two_guid_correct_counts(self):
        PERSIST_object = DummyClass()
        PERSIST_object.refcompressedsequence_read = MagicMock(
                side_effect=[self.sample1, self.sample2])
        analysed_samples, vmodel, mmodel = VariantMatrix.get_position_counts(["s1", "s2"], 100, PERSIST_object)
        self.assertEqual(analysed_samples, {"s1", "s2"})
        expected_vmodel = {1: 1, 2: 2}
        expected_mmodel = {3: 2, 4: 1}
        self.assertEqual(vmodel, expected_vmodel)
        self.assertEqual(mmodel, expected_mmodel)

class TestGetMissingness(TestCase):
    def test_get_missingness_cutoff(self):
        mmodel = {1: 2, 3: 12}
        selected_positions = [2,3,4,5]
        cutoff = VariantMatrix.get_missingness_cutoff(selected_positions, mmodel)
        self.assertEqual(cutoff, 6)

    def test_remove_positions_with_high_missingness(self):
        mmodel = {1: 2, 3: 12}
        selected_positions = {1,2,3,4,5}
        cutoff = 6
        num_removed = VariantMatrix.remove_high_missingness_positions(selected_positions, mmodel, cutoff)
        self.assertEqual(num_removed, 1)
        self.assertEqual(selected_positions, {1,2,4,5})
