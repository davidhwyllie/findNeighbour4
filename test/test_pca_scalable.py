""" tests pca_scalable.py - software to do PCA
"""
import os
import shutil
import uuid
import unittest
import pandas as pd
from pca.pca_scalable import VariantMatrix, PCARunner
from findn.common_utils import ConfigManager
from localstore.localstoreutils import LocalStore


class Test_PCASequenceStore(unittest.TestCase):
    def setUp(self):
        # define a temporary directory
        self.tmpdir = os.path.join("unitTest_tmp", "pcastore", uuid.uuid4().hex)

        # start testing
        shutil.rmtree(
            self.tmpdir, ignore_errors=True
        )  # should not exist, as uuid is unique

    def tearDown(self):
        # remove any temporary directory
        shutil.rmtree(self.tmpdir, ignore_errors=True)


class Test_PCAScalable(unittest.TestCase):
    def setUp(self):
        # define a temporary directory
        self.tmpdir = os.path.join("unitTest_tmp", "pcastore", uuid.uuid4().hex)

        # start testing
        shutil.rmtree(
            self.tmpdir, ignore_errors=True
        )  # should not exist, as uuid is unique

    def tearDown(self):
        # remove any temporary directory
        shutil.rmtree(self.tmpdir, ignore_errors=True)


class Test_VariantMatrix_1a(Test_PCAScalable):
    """tests the VariantMatrix and PCARunner classes using a LocalStore object as a sequence provider"""

    def runTest(self):

        # local json test data n=6000
        LPERSIST = LocalStore("testdata/pca_scalable/test.tar")
        self.assertEqual(6000, len(LPERSIST.sequence_ids()))

        expected_sequences = []
        for sequence_id, rcs in LPERSIST.read_all():

            if rcs is not None:
                if rcs["invalid"] == 0:
                    expected_sequences.append(sequence_id)

        self.assertEqual(5710, len(set(expected_sequences)))
        self.assertEqual(
            len(expected_sequences), len(set(expected_sequences))
        )  # no duplicates

        cfm = ConfigManager("testdata/pca/config.json")
        CONFIG = cfm.read_config()

        # start testing
        v = VariantMatrix(CONFIG, LPERSIST, self.tmpdir, show_bar=False)

        # print("Preparing to analyse 5,710 valid sequences")
        v.prepare_to_analyse()

        # valid samples should be in v.PCASEQSTORE
        self.assertEqual(set(v.PCASEQSTORE.sequence_ids), set(expected_sequences))

        # now compute a variation model
        # print("Analysing variation")
        res = v.get_position_counts()
        self.assertIsInstance(res, dict)

        print("Building matrix")
        i = 0
        analysed_sequence_ids = set()
        guids = v.prepare(
            exclude_positions_with_missingness_fold_over_median=12,
            target_matrix_size=100,
            min_matrix_size=10,
        )

        self.assertEqual(len(guids), len(expected_sequences))

        for sequence_ids, mat in v.matrix_in_blocks(guids):
            i += 1
            for sequence_id in sequence_ids:
                analysed_sequence_ids.add(sequence_id)

        df = v.matrix_sequence_properties
        df = df[df["used_in_pca"]]
        analysed_sequence_ids = set(df.index)
        self.assertEqual(len(df.index), len(analysed_sequence_ids))
        # missing = set(expected_sequences) - set(analysed_sequence_ids)

        # analysed sequence ids will be smaller than expected sequences because no
        # i is the number of blocks
        self.assertEqual(i, 53)


class Test_VariantMatrix_1b(Test_PCAScalable):
    """tests the VariantMatrix and PCARunner classes using a LocalStore object as a sequence provider"""

    def runTest(self):

        # local json test data n=6000
        LPERSIST = LocalStore("testdata/pca_scalable/test.tar")

        self.assertEqual(6000, len(LPERSIST.sequence_ids()))

        cfm = ConfigManager("testdata/pca/config.json")
        CONFIG = cfm.read_config()

        # start testing
        v = VariantMatrix(CONFIG, LPERSIST, self.tmpdir, show_bar=False)

        # print("Preparing to analyse 6,000 sequences")
        v.prepare_to_analyse()

        # valid samples should be in v.PCASEQSTORE
        to_analyse = v.PCASEQSTORE.sequence_ids
        self.assertEqual(len(to_analyse), 5710)

        # now compute a variation model
        print("Analysing variation")
        res = v.get_position_counts()
        self.assertIsInstance(res, dict)

        print("Building matrix")

        # test run
        pcr = PCARunner(v)
        pcr.run(
            n_components=10,
            select_from=to_analyse,
            pca_parameters={},
            target_matrix_size=1000,
            min_matrix_size=1000,
        )

        # test cluster
        v = pcr.cluster()

        cl = pcr.vm.get_variationmodel_attribute("transformed_coordinate_categories")

        # 5300 samples * 10 coord
        self.assertIsInstance(cl, pd.DataFrame)  # it's a dataframe
        self.assertEqual(len(cl.index), 53000)


class Test_VariantMatrix_1c(Test_PCAScalable):
    """tests the VariantMatrix and PCARunner classes using a LocalStore object as a sequence provider
    with automatic batch size selection"""

    def runTest(self):

        # local json test data n=6000
        LPERSIST = LocalStore("testdata/pca_scalable/test.tar")
        self.assertEqual(6000, len(LPERSIST.sequence_ids()))

        expected_sequences = []
        for sequence_id, rcs in LPERSIST.read_all():
            if rcs is not None:
                if rcs["invalid"] == 0:
                    expected_sequences.append(sequence_id)

        self.assertEqual(5710, len(set(expected_sequences)))
        self.assertEqual(
            len(expected_sequences), len(set(expected_sequences))
        )  # no duplicates

        cfm = ConfigManager("testdata/pca/config.json")
        CONFIG = cfm.read_config()

        # start testing
        v = VariantMatrix(CONFIG, LPERSIST, self.tmpdir, show_bar=False)

        # print("Preparing to analyse 5,710 valid sequences")
        v.prepare_to_analyse()

        # valid samples should be in v.PCASEQSTORE
        self.assertEqual(set(v.PCASEQSTORE.sequence_ids), set(expected_sequences))

        # now compute a variation model
        # print("Analysing variation")
        res = v.get_position_counts()
        self.assertIsInstance(res, dict)

        print("Building matrix")
        i = 0
        analysed_sequence_ids = set()
        guids = v.prepare(exclude_positions_with_missingness_fold_over_median=12)

        self.assertEqual(len(guids), len(expected_sequences))

        for sequence_ids, mat in v.matrix_in_blocks(guids):
            i += 1
            for sequence_id in sequence_ids:
                analysed_sequence_ids.add(sequence_id)

        df = v.matrix_sequence_properties
        df = df[df["used_in_pca"]]
        analysed_sequence_ids = set(df.index)
        # missing = set(expected_sequences) - set(analysed_sequence_ids)

        # analysed sequence ids will be smaller than expected sequences because no
        self.assertEqual(len(analysed_sequence_ids), 5300)
        self.assertEqual(i, 1)
