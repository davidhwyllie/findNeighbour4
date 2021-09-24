""" tests py_seqComparer.py

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

 

"""

import unittest
import json
from Bio import SeqIO
from findn.py_seqComparer import py_seqComparer


class test_py_seqComparer_51(unittest.TestCase):
    """tests mcompare"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = py_seqComparer(maxNs=1e8, reference=refSeq, snpCeiling=10)
        n = 0
        originals = [
            "AAACGN",
            "CCCCGN",
            "TTTCGN",
            "GGGGGN",
            "NNNCGN",
            "ACTCGN",
            "TCTNGN",
        ]
        guids = []
        for original in originals:
            n += 1
            c = sc.compress(original)
            guid = "{0}-{1}".format(original, n)
            guids.append(guid)
            sc.persist(c, guid=guid)

        res = sc.mcompare(guids[0])  # defaults to sample size 30
        self.assertEqual(len(res), len(originals) - 1)


class test_py_seqComparer_ec(unittest.TestCase):
    """tests exact comparison"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "G" * 30000
        sc = py_seqComparer(maxNs=1e8, reference=refSeq, snpCeiling=1000)

        obj1 = json.loads(
            """{
        "A": [], 
        "C": [], 
        "G": [], 
        "invalid": 0, 
        "M": [], 
        "N": [], 
        "T": [], 
        "U": []
        }"""
        )

        obj2 = json.loads(
            """{
        "A": [], 
        "C": [], 
        "G": [
        23402
        ], 
        "invalid": 0, 
        "M": {}, 
        "N": [
        385, 
        386, 
        387, 
        388, 
        389, 
        390, 
        391, 
        392, 
        393, 
        394
        ], 
        "T": [
        28931, 
        203, 
        29644, 
        6285, 
        21613, 
        240, 
        19184, 
        10448, 
        22226, 
        27768
        ], 
        "U": [
        385, 
        386, 
        387, 
        388, 
        389, 
        390, 
        391, 
        392, 
        393, 
        394
        ]
        }"""
        )

        sc.persist(obj1, "guid1")
        sc.persist(obj2, "guid2")
        dist = sc.compare("guid1", "guid2")

        self.assertEqual(dist, 11)


class test_py_seqComparer_49(unittest.TestCase):
    """tests reporting on stored contents"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = py_seqComparer(maxNs=1e8, reference=refSeq, snpCeiling=10)
        # need > 30 sequences
        originals = [
            "AAACGN",
            "CCCCGN",
            "TTTCGN",
            "GGGGGN",
            "NNNCGN",
            "ACTCGN",
            "TCTNGN",
            "AAACGN",
            "CCCCGN",
            "TTTCGN",
            "GGGGGN",
            "NNNCGN",
            "ACTCGN",
            "TCTNGN",
            "AAACGN",
            "CCCCGN",
            "TTTCGN",
            "GGGGGN",
            "NNNCGN",
            "ACTCGN",
            "TCTNGN",
            "AAACGN",
            "CCCCGN",
            "TTTCGN",
            "GGGGGN",
            "NNNCGN",
            "ACTCGN",
            "TCTNGN",
            "AAACGN",
            "CCCCGN",
            "TTTCGN",
            "GGGGGN",
            "NNNCGN",
            "ACTCGN",
            "TCTNGN",
            "AAACGN",
            "CCCCGN",
            "TTTCGN",
            "GGGGGN",
            "NNNCGN",
            "ACTCGN",
            "TCTNGN",
        ]
        guid_names = []
        n = 0
        for original in originals:
            n += 1
            c = sc.compress(original)
            this_guid = "{0}-{1}".format(original, n)
            sc.persist(c, guid=this_guid)
            guid_names.append(this_guid)

        res = sc.summarise_stored_items()
        self.assertTrue(isinstance(res, dict))
        self.assertEqual(set(res.keys()), set(["server|scstat|nSeqs"]))


class test_py_seqComparer_48(unittest.TestCase):
    """tests computations of p values from exact bionomial test"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = py_seqComparer(maxNs=1e8, reference=refSeq, snpCeiling=10)
        # need > 30 sequences
        originals = [
            "AAACGN",
            "CCCCGN",
            "TTTCGN",
            "GGGGGN",
            "NNNCGN",
            "ACTCGN",
            "TCTNGN",
            "AAACGN",
            "CCCCGN",
            "TTTCGN",
            "GGGGGN",
            "NNNCGN",
            "ACTCGN",
            "TCTNGN",
            "AAACGN",
            "CCCCGN",
            "TTTCGN",
            "GGGGGN",
            "NNNCGN",
            "ACTCGN",
            "TCTNGN",
            "AAACGN",
            "CCCCGN",
            "TTTCGN",
            "GGGGGN",
            "NNNCGN",
            "ACTCGN",
            "TCTNGN",
            "AAACGN",
            "CCCCGN",
            "TTTCGN",
            "GGGGGN",
            "NNNCGN",
            "ACTCGN",
            "TCTNGN",
            "AAACGN",
            "CCCCGN",
            "TTTCGN",
            "GGGGGN",
            "NNNCGN",
            "ACTCGN",
            "TCTNGN",
        ]
        guid_names = []
        n = 0
        for original in originals:
            n += 1
            c = sc.compress(original)
            this_guid = "{0}-{1}".format(original, n)
            sc.persist(c, guid=this_guid)
            guid_names.append(this_guid)


class test_py_seqComparer_46a(unittest.TestCase):
    """tests estimate_expected_unk, a function estimating the number of Ns in sequences
    by sampling"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = py_seqComparer(maxNs=1e8, reference=refSeq, snpCeiling=10)
        n = 0
        originals = [
            "AAACGN",
            "CCCCGN",
            "TTTCGN",
            "GGGGGN",
            "NNNCGN",
            "ACTCGN",
            "TCTNGN",
        ]
        guids = []
        for original in originals:
            n += 1
            c = sc.compress(original)
            guid = "{0}-{1}".format(original, n)
            guids.append(guid)
            sc.persist(c, guid=guid)

        res = sc.estimate_expected_unk()  # defaults to sample size 30
        self.assertEqual(res, None)

        # analyse the last two
        res = sc.estimate_expected_unk(sample_size=2, exclude_guids=guids[0:5])
        self.assertEqual(res, 1.5)

        # analyse the first two
        res = sc.estimate_expected_unk(sample_size=2, exclude_guids=guids[2:7])
        self.assertEqual(res, 1)


class test_py_seqComparer_46b(unittest.TestCase):
    """tests estimate_expected_unk, a function estimating the number of Ns in sequences
    by sampling"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = py_seqComparer(maxNs=3, reference=refSeq, snpCeiling=10)
        n = 0
        originals = [
            "AAACGN",
            "CCCCGN",
            "TTTCGN",
            "GGGGGN",
            "NNNCGN",
            "ACTCGN",
            "TCTGGN",
        ]
        guids = []
        for original in originals:
            n += 1
            c = sc.compress(original)
            guid = "{0}-{1}".format(original, n)
            guids.append(guid)
            sc.persist(c, guid=guid)

        res = sc.estimate_expected_unk()  # defaults to sample size 30
        self.assertEqual(res, None)

        # analyse them all
        res = sc.estimate_expected_unk(sample_size=7, exclude_guids=[])
        self.assertEqual(res, 1)

        # analyse them all
        res = sc.estimate_expected_unk(sample_size=6, exclude_guids=[])
        self.assertEqual(res, 1)


class test_py_seqComparer_46c(unittest.TestCase):
    """tests estimate_expected_unk_sites, a function estimating the number of Ns in sequences
    by sampling"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = py_seqComparer(maxNs=1e8, reference=refSeq, snpCeiling=10)
        n = 0
        originals = [
            "AAACGN",
            "CCCCGN",
            "TTTCGN",
            "GGGGGN",
            "NNNCGN",
            "ACTCGN",
            "TCTNGN",
        ]
        guids = []
        for original in originals:
            n += 1
            c = sc.compress(original)
            guid = "{0}-{1}".format(original, n)
            guids.append(guid)
            sc.persist(c, guid=guid)

        # analyse nothing
        res = sc.estimate_expected_unk_sites(
            sample_size=2, sites=set([]), exclude_guids=guids[0:5]
        )
        self.assertEqual(res, 0)

        # analyse the last two
        res = sc.estimate_expected_unk_sites(
            sample_size=2, sites=set([0, 1, 2, 3, 4, 5]), exclude_guids=guids[0:5]
        )
        self.assertEqual(res, 1.5)

        # analyse the first two
        res = sc.estimate_expected_unk_sites(
            sample_size=2, sites=set([0, 1, 2, 3, 4, 5]), exclude_guids=guids[2:7]
        )
        self.assertEqual(res, 1)


class test_py_seqComparer_45a(unittest.TestCase):
    """tests the generation of multiple alignments of variant sites."""

    def runTest(self):

        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = py_seqComparer(maxNs=1e8, reference=refSeq, snpCeiling=10)

        originals = [
            "AAACGN",
            "CCCCGN",
            "TTTCGN",
            "GGGGGN",
            "NNNCGN",
            "ACTCGN",
            "TCTNGN",
        ]
        guid_names = []
        n = 0
        for original in originals:
            n += 1
            c = sc.compress(original)
            this_guid = "{0}-{1}".format(original, n)
            sc.persist(c, guid=this_guid)
            guid_names.append(this_guid)

        res = sc.multi_sequence_alignment(guid_names)

        self.assertEqual(len(res.valid_guids), 7)
        self.assertEqual(res.variant_positions, [0, 1, 2, 3])


class test_py_seqComparer_45b(unittest.TestCase):
    """tests the generation of multiple alignments of variant sites."""

    def runTest(self):

        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = py_seqComparer(maxNs=6, reference=refSeq, snpCeiling=10)

        originals = [
            "AAACGN",
            "CCCCGN",
            "TTTCGN",
            "GGGGGN",
            "NNNCGN",
            "ACTCGN",
            "TCTGGN",
        ]
        guid_names = []
        n = 0
        for original in originals:
            n += 1
            c = sc.compress(original)
            this_guid = "{0}-{1}".format(original, n)
            sc.persist(c, guid=this_guid)
            guid_names.append(this_guid)

        res = sc.multi_sequence_alignment(guid_names)
        self.assertEqual(len(res.valid_guids), 7)
        self.assertEqual(res.variant_positions, [0, 1, 2, 3])


class test_py_seqComparer_1(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"
        sc = py_seqComparer(maxNs=1e8, snpCeiling=20, reference=refSeq)
        self.assertEqual(sc.reference, refSeq)


class test_py_seqComparer_2(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"
        sc = py_seqComparer(maxNs=1e8, snpCeiling=20, reference=refSeq)
        with self.assertRaises(TypeError):
            retVal = sc.compress(sequence="AC")
            self.assertTrue(retVal is not None)


class test_py_seqComparer_3(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"
        sc = py_seqComparer(maxNs=1e8, snpCeiling=20, reference=refSeq)
        retVal = sc.compress(sequence="ACTG")
        self.assertEqual(
            retVal,
            {
                "G": set([]),
                "A": set([]),
                "C": set([]),
                "T": set([]),
                "N": set([]),
                "U": set([]),
                "M": {},
                "invalid": 0,
            },
        )


class test_py_seqComparer_3b(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"
        sc = py_seqComparer(maxNs=1e8, snpCeiling=20, reference=refSeq)
        retVal = sc.compress(sequence="ACTQ")
        self.assertEqual(
            retVal,
            {
                "G": set([]),
                "A": set([]),
                "C": set([]),
                "T": set([]),
                "N": set([]),
                "U": set([3]),
                "M": {3: "Q"},
                "invalid": 0,
            },
        )


class test_py_seqComparer_3c(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"
        sc = py_seqComparer(maxNs=1e8, snpCeiling=20, reference=refSeq)
        retVal = sc.compress(sequence="NYTQ")
        self.assertEqual(
            retVal,
            {
                "G": set([]),
                "A": set([]),
                "C": set([]),
                "T": set([]),
                "N": set([0]),
                "U": set([0, 1, 3]),
                "M": {1: "Y", 3: "Q"},
                "invalid": 0,
            },
        )


class test_py_seqComparer_4(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"
        sc = py_seqComparer(maxNs=1e8, snpCeiling=20, reference=refSeq)

        retVal = sc.compress(sequence="ACTN")
        self.assertEqual(
            retVal,
            {
                "G": set([]),
                "A": set([]),
                "C": set([]),
                "T": set([]),
                "N": set([3]),
                "M": {},
                "U": set([3]),
                "invalid": 0,
            },
        )


class test_py_seqComparer_5(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"
        sc = py_seqComparer(maxNs=1e8, snpCeiling=20, reference=refSeq)
        retVal = sc.compress(sequence="ACT-")
        self.assertEqual(
            retVal,
            {
                "G": set([]),
                "A": set([]),
                "C": set([]),
                "T": set([]),
                "N": set([3]),
                "M": {},
                "U": set([3]),
                "invalid": 0,
            },
        )


class test_py_seqComparer_6(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"

        sc = py_seqComparer(maxNs=1e8, snpCeiling=20, reference=refSeq)

        retVal = sc.compress(sequence="TCT-")
        self.assertEqual(
            retVal,
            {
                "G": set([]),
                "A": set([]),
                "C": set([]),
                "T": set([0]),
                "N": set([3]),
                "M": {},
                "U": set([3]),
                "invalid": 0,
            },
        )


class test_py_seqComparer_7(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"

        sc = py_seqComparer(maxNs=1e8, snpCeiling=20, reference=refSeq)
        retVal = sc.compress(sequence="ATT-")
        self.assertEqual(
            retVal,
            {
                "G": set([]),
                "A": set([]),
                "C": set([]),
                "T": set([1]),
                "N": set([3]),
                "M": {},
                "U": set([3]),
                "invalid": 0,
            },
        )


class test_py_seqComparer_6b(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"

        sc = py_seqComparer(maxNs=1e8, snpCeiling=20, reference=refSeq)
        originals = [
            "AAAA",
            "CCCC",
            "TTTT",
            "GGGG",
            "NNNN",
            "ACTG",
            "ACTC",
            "TCTN",
            "NYTQ",
            "QRST",
        ]
        for original in originals:

            compressed_sequence = sc.compress(sequence=original)

            roundtrip = sc.uncompress(compressed_sequence)
            self.assertEqual(original, roundtrip)


class test_py_seqComparer_6c(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"

        sc = py_seqComparer(maxNs=1e8, snpCeiling=20, reference=refSeq)
        originals = ["NNNN"]
        for original in originals:

            compressed_sequence = sc.compress(sequence=original)
            roundtrip = sc.uncompress(compressed_sequence)
            self.assertEqual(original, roundtrip)


class test_py_seqComparer_6d(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"

        sc = py_seqComparer(maxNs=3, snpCeiling=20, reference=refSeq)
        originals = ["NNNN"]
        for original in originals:

            compressed_sequence = sc.compress(sequence=original)
            with self.assertRaises(ValueError):
                sc.uncompress(compressed_sequence)


class test_py_seqComparer_16(unittest.TestCase):
    """tests the comparison of two sequences where both differ from the reference."""

    def runTest(self):
        # generate compressed sequences
        refSeq = "ACTG"
        sc = py_seqComparer(maxNs=1e8, reference=refSeq, snpCeiling=10)

        seq1 = sc.compress("AAAA")
        seq2 = sc.compress("CCCC")
        self.assertEqual(sc.countDifferences(seq1, seq2), 4)


class test_py_seqComparer_16b(unittest.TestCase):
    """tests the comparison of two sequences where both differ from the reference."""

    def runTest(self):
        # generate compressed sequences
        refSeq = "ACTG"
        sc = py_seqComparer(maxNs=1e8, reference=refSeq, snpCeiling=10)

        seq1 = sc.compress("AAAA")
        seq2 = sc.compress("RRCC")
        self.assertEqual(sc.countDifferences(seq1, seq2), 2)


class test_py_seqComparer_16c(unittest.TestCase):
    """tests the comparison of two sequences where both differ from the reference."""

    def runTest(self):
        # generate compressed sequences
        refSeq = "ACTG"
        sc = py_seqComparer(maxNs=1e8, reference=refSeq, snpCeiling=10)

        seq1 = sc.compress("AAAA")
        seq2 = sc.compress("RRNN")
        self.assertEqual(sc.countDifferences(seq1, seq2), 0)


class test_py_seqComparer_17(unittest.TestCase):
    """tests the comparison of two sequences where one is invalid"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "ACTG"
        sc = py_seqComparer(maxNs=3, reference=refSeq, snpCeiling=10)

        seq1 = sc.compress("AAAA")
        seq2 = sc.compress("NNNN")
        self.assertEqual(sc.countDifferences(seq1, seq2), None)


class test_py_seqComparer_cmp(unittest.TestCase):
    """tests the comparison of two sequences where both differ from the reference."""

    def runTest(self):
        # generate compressed sequences
        refSeq = "ACTG"
        sc = py_seqComparer(maxNs=1e8, reference=refSeq, snpCeiling=10)

        seq1 = sc.compress("AAAA")
        seq2 = sc.compress("CCCC")
        sc.persist(seq1, "s1")
        sc.persist(seq2, "s2")

        self.assertEqual(sc.compare("s1", "s2"), 4)

        with self.assertRaises(KeyError):
            sc.compare("s1", "not_there")


class test_py_seqComparer_saveload3(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"
        sc = py_seqComparer(maxNs=1e8, snpCeiling=20, reference=refSeq)
        compressedObj = sc.compress(sequence="ACTT")
        sc.persist(compressedObj, "one")
        retVal = sc.load(guid="one")
        self.assertEqual(compressedObj, retVal)


class test_py_seqComparer_save_remove(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"
        sc = py_seqComparer(maxNs=1e8, snpCeiling=20, reference=refSeq)
        compressedObj = sc.compress(sequence="ACTT")
        sc.persist(compressedObj, "one")
        retVal = sc.iscachedinram(guid="one")
        self.assertEqual(True, retVal)
        sc.remove("one")
        retVal = sc.iscachedinram(guid="one")
        self.assertEqual(False, retVal)


class test_py_seqComparer_24(unittest.TestCase):
    """tests N compression"""

    def runTest(self):

        refSeq = "ACTGTTAATTTTTTTTTGGGGGGGGGGGGAA"
        sc = py_seqComparer(maxNs=1e8, snpCeiling=20, reference=refSeq)

        retVal = sc.compress(sequence="ACTGTTAANNNNNNNNTGGGGGGGGGGGGAA")
        self.assertEqual(
            retVal,
            {
                "G": set([]),
                "A": set([]),
                "C": set([]),
                "T": set([]),
                "M": {},
                "N": set([8, 9, 10, 11, 12, 13, 14, 15]),
                "U": set([8, 9, 10, 11, 12, 13, 14, 15]),
                "invalid": 0,
            },
        )
        retVal = sc.compress(sequence="NNTGTTAANNNNNNNNTGGGGGGGGGGGGAA")
        self.assertEqual(
            retVal,
            {
                "G": set([]),
                "A": set([]),
                "C": set([]),
                "T": set([]),
                "M": {},
                "N": set([0, 1, 8, 9, 10, 11, 12, 13, 14, 15]),
                "U": set([0, 1, 8, 9, 10, 11, 12, 13, 14, 15]),
                "invalid": 0,
            },
        )


class test_py_seqComparer_29(unittest.TestCase):
    """tests _setStats"""

    def runTest(self):

        refSeq = "ACTGTTAATTTTTTTTTGGGGGGGGGGGGAA"
        sc = py_seqComparer(maxNs=1e8, snpCeiling=20, reference=refSeq)
        compressedObj1 = sc.compress(sequence="GGGGTTAANNNNNNNNNGGGGGAAAAGGGAA")
        compressedObj2 = sc.compress(sequence="ACTGTTAATTTTTTTTTNNNNNNNNNNNNNN")
        (n1, n2, nall, rv1, rv2, retVal) = sc._setStats(
            compressedObj1["N"], compressedObj2["N"]
        )
        self.assertEqual(
            retVal,
            set(
                [
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    16,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                    24,
                    25,
                    26,
                    27,
                    28,
                    29,
                    30,
                ]
            ),
        )

        compressedObj1 = sc.compress(sequence="GGGGTTAANNNNNNNNTGGGGGAAAAGGGAA")
        compressedObj2 = sc.compress(sequence="ACTGTTAATTTTTTTTTNNNNNNNNNNNNNN")
        (n1, n2, nall, rv1, rv2, retVal) = sc._setStats(
            compressedObj1["N"], compressedObj2["N"]
        )
        self.assertEqual(
            retVal,
            set(
                [
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                    24,
                    25,
                    26,
                    27,
                    28,
                    29,
                    30,
                ]
            ),
        )

        compressedObj1 = sc.compress(sequence="NNNGTTAANNNNNNNNTGGGGGAAAAGGGAA")
        compressedObj2 = sc.compress(sequence="ACTGTTAATTTTTTTTTNNNNNNNNNNNNNN")
        (n1, n2, nall, rv1, rv2, retVal) = sc._setStats(
            compressedObj1["N"], compressedObj2["N"]
        )
        self.assertEqual(
            retVal,
            set(
                [
                    0,
                    1,
                    2,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                    24,
                    25,
                    26,
                    27,
                    28,
                    29,
                    30,
                ]
            ),
        )

        compressedObj1 = sc.compress(sequence="NNNGTTAANNNNNNNNTGGGGGAAAAGGGAA")
        compressedObj2 = sc.compress(sequence="ACTNNNNNTTTTTTTTTNNNNNNNNNNNNNN")
        (n1, n2, nall, rv1, rv2, retVal) = sc._setStats(
            compressedObj1["N"], compressedObj2["N"]
        )
        self.assertEqual(
            retVal,
            set(
                [
                    0,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                    24,
                    25,
                    26,
                    27,
                    28,
                    29,
                    30,
                ]
            ),
        )

        compressedObj1 = sc.compress(sequence="NNNGTTAANNNNNNNNTGGGGGAAAAGGGAA")
        compressedObj2 = sc.compress(sequence="ACTNNNNNTTTTTTTTTQQQQQQQQQQQQQQ")
        (n1, n2, nall, rv1, rv2, retVal) = sc._setStats(
            compressedObj1["N"], compressedObj2["N"]
        )
        self.assertEqual(
            retVal, set([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])
        )
        (n1, n2, nall, rv1, rv2, retVal) = sc._setStats(
            compressedObj1["M"], compressedObj2["M"]
        )
        self.assertEqual(
            retVal, set([17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30])
        )

        compressedObj1 = sc.compress(sequence="qqqGTTAAqqqqqqqqTGGGGGAAAAGGGAA")
        compressedObj2 = sc.compress(sequence="ACTqqqqqTTTTTTTTTqqqqqqqqqqqqqq")
        (n1, n2, nall, rv1, rv2, retVal) = sc._setStats(
            compressedObj1["M"], compressedObj2["M"]
        )
        self.assertEqual(
            retVal,
            set(
                [
                    0,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                    24,
                    25,
                    26,
                    27,
                    28,
                    29,
                    30,
                ]
            ),
        )


class test_py_seqComparer_37(unittest.TestCase):
    """tests the loading of an exclusion file"""

    def runTest(self):

        # default exclusion file
        refSeq = "ACTG"
        sc = py_seqComparer(maxNs=1e8, reference=refSeq, snpCeiling=1)
        self.assertEqual(
            sc.excluded_hash(), "Excl 0 nt [d751713988987e9331980363e24189ce]"
        )


class test_py_seqComparer_38(unittest.TestCase):
    """tests the loading of an exclusion file"""

    def runTest(self):

        # no exclusion file
        refSeq = "ACTG"
        sc = py_seqComparer(maxNs=1e8, reference=refSeq, snpCeiling=1)
        self.assertEqual(
            sc.excluded_hash(), "Excl 0 nt [d751713988987e9331980363e24189ce]"
        )


class test_py_seqComparer_40(unittest.TestCase):
    """tests the computation of a hash of a compressed object"""

    def runTest(self):

        # generate compressed sequences
        refSeq = "ACTG"
        sc = py_seqComparer(maxNs=1e8, reference=refSeq, snpCeiling=10)
        compressed_sequence = sc.compress(sequence="TTAA")

        res = sc.compressed_sequence_hash(compressed_sequence)
        self.assertEqual(res, "da8785691df5858b0b847db59bdefd11")


class test_py_seqComparer_45(unittest.TestCase):
    """tests insertion of large sequences"""

    def runTest(self):
        inputfile = "reference/NC_000962.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                goodseq = str(record.seq)
                badseq = "".join("N" * len(goodseq))
                originalseq = list(str(record.seq))
        sc = py_seqComparer(maxNs=1e8, reference=record.seq, snpCeiling=100)
        n_pre = 0
        guids_inserted = list()
        for i in range(1, 4):  # 40

            seq = originalseq
            if i % 5 == 0:
                is_mixed = True
                guid_to_insert = "mixed_{0}".format(n_pre + i)
            else:
                is_mixed = False
                guid_to_insert = "nomix_{0}".format(n_pre + i)
            # make i mutations at position 500,000

            offset = 500000
            nVariants = 0
            for j in range(i):
                mutbase = offset + j
                ref = seq[mutbase]
                if is_mixed is False:
                    nVariants += 1
                    if not ref == "T":
                        seq[mutbase] = "T"
                    if not ref == "A":
                        seq[mutbase] = "A"
                if is_mixed is True:
                    seq[mutbase] = "N"
            seq = "".join(seq)

            if i % 11 == 0:
                seq = badseq  # invalid

            guids_inserted.append(guid_to_insert)
            if not is_mixed:
                # print("Adding TB sequence {2} of {0} bytes with {1} Ns and {3} variants relative to ref.".format(len(seq), seq.count('N'), guid_to_insert, nVariants))
                pass
            else:
                # print("Adding mixed TB sequence {2} of {0} bytes with {1} Ns relative to ref.".format(len(seq), seq.count('N'), guid_to_insert))
                pass
            self.assertEqual(len(seq), 4411532)  # check it's the right sequence

            c = sc.compress(seq)
            sc.persist(c, guid=guid_to_insert)


class test_py_seqComparer_47(unittest.TestCase):
    """tests raise_error"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGGGGGGGG"
        sc = py_seqComparer(maxNs=1e8, reference=refSeq, snpCeiling=10)
        with self.assertRaises(ZeroDivisionError):
            sc.raise_error("token")


class test_py_seqComparer_47dist(unittest.TestCase):
    """tests distmat, a function yielding a distance matrix."""

    def runTest(self):

        # generate compressed sequences
        refSeq = "GGGGGGGGGGGG"
        sc = py_seqComparer(maxNs=1e8, reference=refSeq, snpCeiling=10)

        originals = ["AAACACTGACTG", "CCCCACTGACTG", "TTTCACTGACTG"]
        for original in originals:
            c = sc.compress(original)
            sc.persist(c, guid=original)

        n = 0
        for item in sc.distmat(half=False, diagonal=True):
            n += 1
        l_originals = len(originals)

        self.assertEqual(n, l_originals * l_originals)

        n = 0
        for item in sc.distmat(half=False, diagonal=False):
            n += 1

        l_originals = len(originals)

        self.assertEqual(n, (l_originals * l_originals) - l_originals)

        n = 0
        for item in sc.distmat(half=True, diagonal=False):
            n += 1

        l_originals = len(originals)

        self.assertEqual(n, (l_originals * (l_originals - 1) / 2))


class test_py_seqComparer_50(unittest.TestCase):
    """tests estimate_expected_proportion, a function computing the proportion of Ns expected based on the median
    Ns in a list of sequences"""

    def runTest(self):
        refSeq = "GGGGGGGGGGGG"
        sc = py_seqComparer(maxNs=1e8, reference=refSeq, snpCeiling=10)

        res = sc.estimate_expected_proportion([])
        self.assertTrue(res is None)

        res = sc.estimate_expected_proportion(["AA", "AA"])
        self.assertTrue(res is None)

        res = sc.estimate_expected_proportion(["AA", "AA", "AA"])
        self.assertTrue(res is not None)
        self.assertTrue(res == 0)

        res = sc.estimate_expected_proportion(["AAN", "AAN", "AAN"])
        self.assertTrue(res is not None)
        self.assertAlmostEqual(res, 1 / 3)
