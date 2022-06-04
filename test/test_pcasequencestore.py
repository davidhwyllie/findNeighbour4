""" tests PCASequenceStore - a module to store sparse (matrix-like) representations of 
large numbers of sequences in a tar file
"""

import os
import numpy as np
import shutil
import uuid
import unittest

from localstore.pcasequencestore import PCASequenceStore


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


class Test_PCASequenceStore_1(Test_PCASequenceStore):
    """tests the PCASequenceStore class startup"""

    def runTest(self):

        reference = "ACTGACTG"
        pss1 = PCASequenceStore(reference, self.tmpdir)

        self.assertTrue(os.path.exists(self.tmpdir))

        pss2 = PCASequenceStore(reference, self.tmpdir)

        self.assertIsNone(pss1.mapping)
        self.assertIsNone(pss2.mapping)

        pss1._generate_mapping_file()
        pss2._generate_mapping_file()

        self.assertTrue(np.array_equal(pss1.mapping, pss2.mapping))

        # test encoding and decoding of variants
        for nt in ["A", "C", "G", "T"]:
            for pos in range(len(reference)):
                variant_int_id = pss1._variant_int_id(nt, pos)
                res = pss1.identify_variant(variant_int_id)
                self.assertEqual(nt, res["nt"])
                self.assertEqual(pos, res["pos"])


class Test_PCASequenceStore_1a(Test_PCASequenceStore):
    """tests the PCASequenceStore class startup"""

    def runTest(self):

        with self.assertRaises(TypeError):
            PCASequenceStore(24, self.tmpdir)


class Test_PCASequenceStore_2(Test_PCASequenceStore):
    """tests the PCASequenceStore class fasta encoding"""

    def runTest(self):

        # start testing
        reference = "ACTGACTG"
        PCASequenceStore(reference, self.tmpdir)


class Test_PCASequenceStore_3(Test_PCASequenceStore):
    """tests the PCASequenceStore class length verification"""

    def runTest(self):

        # start testing
        reference = "ACTGACTG"
        pss1 = PCASequenceStore(reference, self.tmpdir)

        with self.assertRaises(ValueError):
            pss1._encode_sequence("ACTGACTGN")


class Test_PCASequenceStore_4(Test_PCASequenceStore):
    """tests the PCASequenceStore class fasta encoding"""

    def runTest(self):

        # local json test data n=5000
        # LPERSIST = LocalStore("testdata/pca_scalable/test.tar")

        # start testing
        reference = "ACTG"
        pss1 = PCASequenceStore(reference, self.tmpdir)

        encoded = pss1._encode_sequence("ACTG")
        self.assertTrue(
            np.array_equal(
                encoded["variants"].toarray(),
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertTrue(np.array_equal(encoded["ns"].toarray(), [[0, 0, 0, 0]]))

        encoded = pss1._encode_sequence("CCTG")
        self.assertTrue(
            np.array_equal(
                encoded["variants"].toarray(),
                [[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertTrue(np.array_equal(encoded["ns"].toarray(), [[0, 0, 0, 0]]))

        encoded = pss1._encode_sequence("TCTG")
        self.assertTrue(
            np.array_equal(
                encoded["variants"].toarray(),
                [[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertTrue(np.array_equal(encoded["ns"].toarray(), [[0, 0, 0, 0]]))

        encoded = pss1._encode_sequence("GCTG")
        self.assertTrue(
            np.array_equal(
                encoded["variants"].toarray(),
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]],
            )
        )
        self.assertTrue(np.array_equal(encoded["ns"].toarray(), [[0, 0, 0, 0]]))

        encoded = pss1._encode_sequence("ACTA")
        self.assertTrue(
            np.array_equal(
                encoded["variants"].toarray(),
                [[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertTrue(np.array_equal(encoded["ns"].toarray(), [[0, 0, 0, 0]]))

        encoded = pss1._encode_sequence("ACTC")
        self.assertTrue(
            np.array_equal(
                encoded["variants"].toarray(),
                [[0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertTrue(np.array_equal(encoded["ns"].toarray(), [[0, 0, 0, 0]]))

        encoded = pss1._encode_sequence("ACTT")
        self.assertTrue(
            np.array_equal(
                encoded["variants"].toarray(),
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]],
            )
        )
        self.assertTrue(np.array_equal(encoded["ns"].toarray(), [[0, 0, 0, 0]]))

        encoded = pss1._encode_sequence("NCTG")
        self.assertTrue(
            np.array_equal(
                encoded["variants"].toarray(),
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertTrue(np.array_equal(encoded["ns"].toarray(), [[1, 0, 0, 0]]))

        encoded = pss1._encode_sequence("NCTN")
        self.assertTrue(
            np.array_equal(
                encoded["variants"].toarray(),
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertTrue(np.array_equal(encoded["ns"].toarray(), [[1, 0, 0, 1]]))


class Test_PCASequenceStore_5(Test_PCASequenceStore):
    """tests the PCASequenceStore class fasta encoding"""

    def runTest(self):

        # start testing
        reference = "ACTG"
        pss1 = PCASequenceStore(reference, self.tmpdir)

        # ACTG
        obj = {
            "A": {},
            "C": {},
            "T": {},
            "G": {},
            "N": {},
            "M": {},
            "invalid": 0,
            "U": [],
        }
        encoded = pss1._encode_rcs(obj)
        self.assertTrue(
            np.array_equal(
                encoded["variants"].toarray(),
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertTrue(np.array_equal(encoded["ns"].toarray(), [[0, 0, 0, 0]]))

        # "CCTG"
        obj = {
            "A": {},
            "C": {0},
            "T": {},
            "G": {},
            "N": {},
            "M": {},
            "invalid": 0,
            "U": [],
        }
        encoded = pss1._encode_rcs(obj)
        self.assertTrue(
            np.array_equal(
                encoded["variants"].toarray(),
                [[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertTrue(np.array_equal(encoded["ns"].toarray(), [[0, 0, 0, 0]]))

        # "TCTG"
        obj = {
            "A": {},
            "C": {},
            "T": {0},
            "G": {},
            "N": {},
            "M": {},
            "invalid": 0,
            "U": [],
        }
        encoded = pss1._encode_rcs(obj)
        self.assertTrue(
            np.array_equal(
                encoded["variants"].toarray(),
                [[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertTrue(np.array_equal(encoded["ns"].toarray(), [[0, 0, 0, 0]]))

        # "GCTG"
        obj = {
            "A": {},
            "C": {},
            "T": {},
            "G": {0},
            "N": {},
            "M": {},
            "invalid": 0,
            "U": [],
        }
        encoded = pss1._encode_rcs(obj)
        self.assertTrue(
            np.array_equal(
                encoded["variants"].toarray(),
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]],
            )
        )
        self.assertTrue(np.array_equal(encoded["ns"].toarray(), [[0, 0, 0, 0]]))

        # "ACTA"
        obj = {
            "A": {3},
            "C": {},
            "T": {},
            "G": {},
            "N": {},
            "M": {},
            "invalid": 0,
            "U": [],
        }
        encoded = pss1._encode_rcs(obj)
        self.assertTrue(
            np.array_equal(
                encoded["variants"].toarray(),
                [[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertTrue(np.array_equal(encoded["ns"].toarray(), [[0, 0, 0, 0]]))

        # "ACTC"
        obj = {
            "A": {},
            "C": {3},
            "T": {},
            "G": {},
            "N": {},
            "M": {},
            "invalid": 0,
            "U": [],
        }
        encoded = pss1._encode_rcs(obj)
        self.assertTrue(
            np.array_equal(
                encoded["variants"].toarray(),
                [[0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertTrue(np.array_equal(encoded["ns"].toarray(), [[0, 0, 0, 0]]))

        # "ACTT"
        obj = {
            "A": {},
            "C": {},
            "T": {3},
            "G": {},
            "N": {},
            "M": {},
            "invalid": 0,
            "U": [],
        }
        encoded = pss1._encode_rcs(obj)
        self.assertTrue(
            np.array_equal(
                encoded["variants"].toarray(),
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]],
            )
        )
        self.assertTrue(np.array_equal(encoded["ns"].toarray(), [[0, 0, 0, 0]]))

        # "NCTG"
        obj = {
            "A": {},
            "C": {},
            "T": {},
            "G": {},
            "N": {0},
            "M": {},
            "invalid": 0,
            "U": [0],
        }
        encoded = pss1._encode_rcs(obj)
        self.assertTrue(
            np.array_equal(
                encoded["variants"].toarray(),
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertTrue(np.array_equal(encoded["ns"].toarray(), [[1, 0, 0, 0]]))

        # "NCTN"
        obj = {
            "A": {},
            "C": {},
            "T": {},
            "G": {},
            "N": {0, 3},
            "M": {},
            "invalid": 0,
            "U": [0, 3],
        }
        encoded = pss1._encode_rcs(obj)
        self.assertTrue(
            np.array_equal(
                encoded["variants"].toarray(),
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertTrue(np.array_equal(encoded["ns"].toarray(), [[1, 0, 0, 1]]))


class Test_PCASequenceStore_6(Test_PCASequenceStore):
    """tests the PCASequenceStore class storage"""

    def runTest(self):

        # start testing
        reference = "ACTG"
        pss = PCASequenceStore(reference, self.tmpdir)

        pss.store_sequence("one", "ACTG")
        self.assertEqual(pss.sequence_ids, set(["one"]))
        pss.store_sequence("two", "ACTN")
        self.assertEqual(pss.sequence_ids, set(["one", "two"]))
        pss.store_sequence("three", "CCTN")
        self.assertEqual(pss.sequence_ids, set(["one", "two", "three"]))

        pss.flush()
        self.assertEqual(
            set(pss.localstore.sequence_ids()), set(["one", "two", "three"])
        )

        obj = {
            "A": {},
            "C": {},
            "T": {},
            "G": {},
            "N": {0, 3},
            "M": {},
            "invalid": 0,
            "U": [0, 3],
        }
        pss.store_rcs("four", obj)
        self.assertEqual(pss.sequence_ids, set(["one", "two", "three", "four"]))

        pss.flush()
        self.assertEqual(
            set(pss.localstore.sequence_ids()), set(["one", "two", "three", "four"])
        )


class Test_PCASequenceStore_7(Test_PCASequenceStore):
    """tests the PCASequenceStore class read_all and read_many methods"""

    def runTest(self):

        reference = "ACTG"
        pss = PCASequenceStore(reference, self.tmpdir)

        pss.store_sequence("one", "ACTG")
        self.assertEqual(pss.sequence_ids, set(["one"]))
        pss.store_sequence("two", "ACTN")
        self.assertEqual(pss.sequence_ids, set(["one", "two"]))
        pss.store_sequence("three", "CCTN")
        self.assertEqual(pss.sequence_ids, set(["one", "two", "three"]))

        pss.flush()

        result = set()
        for sequence_id, obj in pss.read_all():
            result.add(sequence_id)

        self.assertEqual(result, set(["one", "two", "three"]))

        result = set()
        for sequence_id, obj in pss.read_many(["one", "two", "three"]):
            result.add(sequence_id)

        self.assertEqual(result, set(["one", "two", "three"]))


class Test_PCASequenceStore_8(Test_PCASequenceStore):
    """tests the PCASequenceStore class summarise_all and summarise_many methods"""

    def runTest(self):

        reference = "ACTG"
        pss = PCASequenceStore(reference, self.tmpdir)

        pss.store_sequence("one", "ACTG")
        self.assertEqual(pss.sequence_ids, set(["one"]))
        pss.store_sequence("two", "ACTN")
        self.assertEqual(pss.sequence_ids, set(["one", "two"]))
        pss.store_sequence("three", "CCTN")
        self.assertEqual(pss.sequence_ids, set(["one", "two", "three"]))

        pss.flush()

        summary = pss.summarise_all()
        self.assertTrue(np.array_equal(summary["mmodel"], [[0, 0, 0, 2]]))
        self.assertTrue(np.array_equal(summary["vmodel"], [[1, 0, 0, 0]]))
        self.assertTrue(
            np.array_equal(
                summary["allelemodel"],
                [[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertEqual(summary["sequence_ids"], set(["one", "two", "three"]))

        summary = pss.summarise_many()
        self.assertTrue(np.array_equal(summary["mmodel"], [[0, 0, 0, 2]]))
        self.assertTrue(np.array_equal(summary["vmodel"], [[1, 0, 0, 0]]))
        self.assertTrue(
            np.array_equal(
                summary["allelemodel"],
                [[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertEqual(summary["sequence_ids"], set(["one", "two", "three"]))

        summary = pss.summarise_many(["one"])
        self.assertTrue(np.array_equal(summary["mmodel"], [[0, 0, 0, 0]]))
        self.assertTrue(np.array_equal(summary["vmodel"], [[0, 0, 0, 0]]))
        self.assertTrue(
            np.array_equal(
                summary["allelemodel"],
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertEqual(summary["sequence_ids"], set(["one"]))

        summary = pss.summarise_many(["one", "two", "three"])
        self.assertTrue(np.array_equal(summary["mmodel"], [[0, 0, 0, 2]]))
        self.assertTrue(np.array_equal(summary["vmodel"], [[1, 0, 0, 0]]))
        self.assertTrue(
            np.array_equal(
                summary["allelemodel"],
                [[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertEqual(summary["sequence_ids"], set(["one", "two", "three"]))

        summary = pss.summarise_many(["one"])
        self.assertTrue(np.array_equal(summary["mmodel"], [[0, 0, 0, 0]]))
        self.assertTrue(
            np.array_equal(
                summary["allelemodel"],
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertEqual(summary["sequence_ids"], set(["one"]))

        summary = pss.summarise_many(["two", "three"], starting_result=summary)
        self.assertTrue(np.array_equal(summary["mmodel"], [[0, 0, 0, 2]]))
        self.assertTrue(np.array_equal(summary["vmodel"], [[1, 0, 0, 0]]))
        self.assertTrue(
            np.array_equal(
                summary["allelemodel"],
                [[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertEqual(summary["sequence_ids"], set(["one", "two", "three"]))


class Test_PCASequenceStore_9(Test_PCASequenceStore):
    """tests the PCASequenceStore class update_summary methods"""

    def runTest(self):

        reference = "ACTG"
        pss = PCASequenceStore(reference, self.tmpdir)

        pss.store_sequence("one", "ACTG")
        self.assertEqual(pss.sequence_ids, set(["one"]))
        pss.store_sequence("two", "ACTN")
        self.assertEqual(pss.sequence_ids, set(["one", "two"]))
        pss.store_sequence("three", "CCTN")
        self.assertEqual(pss.sequence_ids, set(["one", "two", "three"]))

        pss.flush()

        summary = pss.update_summary()
        self.assertTrue(np.array_equal(summary["mmodel"], [[0, 0, 0, 2]]))
        self.assertTrue(
            np.array_equal(
                summary["allelemodel"],
                [[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertEqual(summary["sequence_ids"], set(["one", "two", "three"]))

        summary = pss.summarise_many(["one"])
        self.assertTrue(np.array_equal(summary["mmodel"], [[0, 0, 0, 0]]))
        self.assertTrue(
            np.array_equal(
                summary["allelemodel"],
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertEqual(summary["sequence_ids"], set(["one"]))

        summary = pss.update_summary()
        self.assertTrue(np.array_equal(summary["mmodel"], [[0, 0, 0, 2]]))
        self.assertTrue(
            np.array_equal(
                summary["allelemodel"],
                [[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
            )
        )
        self.assertEqual(summary["sequence_ids"], set(["one", "two", "three"]))


class Test_PCASequenceStore_10(Test_PCASequenceStore):
    """tests the PCASequenceStore class update_summary methods"""

    def runTest(self):

        reference = "ACTG"
        pss = PCASequenceStore(reference, self.tmpdir)

        allelemodel = [[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
        vmodel = pss.allelemodel2vmodel(allelemodel)
        self.assertTrue(np.array_equal(vmodel, [1, 0, 0, 0]))

        allelemodel = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
        vmodel = pss.allelemodel2vmodel(allelemodel)
        self.assertTrue(np.array_equal(vmodel, [0, 0, 0, 0]))

        allelemodel = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
        vmodel = pss.allelemodel2vmodel(allelemodel)
        self.assertTrue(np.array_equal(vmodel, [4, 4, 4, 4]))
