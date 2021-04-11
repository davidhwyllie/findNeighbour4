""" integration test for findNeighbour4
assumes a findNeighbour4 server is running, with the connection string stated in demos/AC587/config/config_nocl.json.

An example command doing this would be 
pipenv run python3 findNeighbour4_server.py demos/AC587/config/config_nocl.json

The test loads the server with data from the AC587 test data
and compares the SNP distances and links with those from a seqComparer
instance.

"""

import unittest
import pandas as pd
from findn import DEFAULT_CONFIG_FILE
from pca.pca import PersistenceTest, MNStats, VariantMatrix, PCARunner 
from findn.common_utils import ConfigManager

class Test_PersistenceTest_1(unittest.TestCase):
    """ tests the PersistenceTest class"""
    def runTest(self):
        tp = PersistenceTest(connstring='thing',number_samples=251)
        tp.load_data(
            sample_ids_file = "/data/software/fn4dev/testdata/pca/seqs_5000test_ids.pickle",
            sequences_file = "/data/software/fn4dev/testdata/pca/seqs_5000test.pickle"
        )
        self.assertEqual(len(tp.seqs.keys()), 5000)
        self.assertEqual(len(tp.sample_ids), 5000)
        self.assertEqual(len(tp.refcompressedsequence_guids()), 5000)
        res = tp.refcompressedsequence_read('NOT THERE')
        self.assertEqual(res, None)

        guids = list(tp.sample_ids)
        res = tp.refcompressedsequence_read(guids[0])
        self.assertIsInstance(res, dict)


class Test_MNStats_1(unittest.TestCase):
    """ tests the MNStats class"""
    def runTest(self):
        
        rcs = {
            'A': {28281, 5387, 23603, 23270}, 
            'C': {6953, 15095, 16175, 24913, 28279}, 
            'T': {5985, 3266, 27971, 12069, 14407, 4299, 15278, 28047, 240, 912, 28976, 14675, 23708, 23062, 28280, 3036, 28094}, 
            'G': {24505, 23402, 28110}, 
            'N': {21764, 21765, 21766, 21767, 21768, 21769, 11287, 11288, 11289, 11290, 11291, 11292, 11293, 11294, 11295, 21990, 21991, 21992, 28270}, 
            'M': {}, 
            'invalid': 0}

        mns = MNStats(select_positions = list(range(29000)), analysed_reference_length = 29905)
        res = mns.examine(rcs)
        
        self.assertIsInstance(res, dict)
        self.assertEqual(res['N_total'],19)
        self.assertEqual(res['M_total'],0)
        
       

class Test_VariantMatrix_1(unittest.TestCase):
    """ tests the VariantMatrix and PCARunner classes"""
    def runTest(self):
        
        TPERSIST = PersistenceTest(connstring='thing',number_samples=251)
        TPERSIST.load_data(
            sample_ids_file = "/data/software/fn4dev/testdata/pca/seqs_5000test_ids.pickle",
            sequences_file = "/data/software/fn4dev/testdata/pca/seqs_5000test.pickle"
        )
        cfm = ConfigManager(DEFAULT_CONFIG_FILE)
        CONFIG = cfm.read_config()
        
        v = VariantMatrix(CONFIG, TPERSIST, show_bar=False)

        # test guids() method
        self.assertEqual(set(v.guids()), set(TPERSIST.refcompressedsequence_guids()))

        # test get_position_counts
        guids,vmodel,mmodel = v.get_position_counts()
        self.assertEqual(guids, set(TPERSIST.refcompressedsequence_guids()))
        self.assertIsInstance(vmodel, dict)
        self.assertIsInstance(mmodel, dict)
       
      
        # test get_missingness_cutoff
        m = v.get_missingness_cutoff(positions = vmodel.keys(), mmodel = mmodel)        # the missingness model
        self.assertEqual(m, 27)     

         # test build
        v.build()
        self.assertIsInstance(v.vm['variant_matrix'],pd.DataFrame)

        # test run
        pcr = PCARunner(v)
        pcr.run(n_components= 10, pca_parameters = {})

        # test cluster
        v= pcr.cluster()

        v.to_sqlite("unittest_tmp")