from unittest import TestCase
from unittest.mock import MagicMock

from findn.pca_covidtest_inc import SNPMatrix

class DummyClass:
    pass

class TestGetBaseCounts(TestCase):
    sample1 = {"A": [1,2], "M": [3], "invalid": 0}
    sample2 = {"A": [2], "M": [3, 4], "invalid": 0}
    def test_one_guid_correct_counts(self):
        PERSIST_object = DummyClass()
        PERSIST_object.refcompressedsequence_read = MagicMock(side_effect=[self.sample1])
        vmodel, mmodel = SNPMatrix.get_base_counts(["sample"], 100, PERSIST_object)
        expected_vmodel = {1: 1, 2: 1}
        expected_mmodel = {3: 1}
        self.assertEqual(vmodel, expected_vmodel)
        self.assertEqual(mmodel, expected_mmodel)

    def test_two_guid_correct_counts(self):
        PERSIST_object = DummyClass()
        PERSIST_object.refcompressedsequence_read = MagicMock(
                side_effect=[self.sample1, self.sample2])
        vmodel, mmodel = SNPMatrix.get_base_counts(["sample", "sample"], 100, PERSIST_object)
        expected_vmodel = {1: 1, 2: 2}
        expected_mmodel = {3: 2, 4: 1}
        self.assertEqual(vmodel, expected_vmodel)
        self.assertEqual(mmodel, expected_mmodel)
