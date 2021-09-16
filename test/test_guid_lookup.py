""" runs unittest for guidLookup

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

 

"""

import unittest
import datetime
import uuid
import time
from findn.guidLookup import guidSearcher, guidDbSearcher
from findn.persistence import Persistence
from findn.NucleicAcid import NucleicAcid

## persistence unit tests
UNITTEST_MONGOCONN = "mongodb://localhost"
UNITTEST_RDBMSCONN = "sqlite://"


class test_gm_1(unittest.TestCase):
    def runTest(self):

        gs = guidSearcher()
        gs.add("b1")
        self.assertEqual(gs.guids, ["b1"])
        gs.add("b3")
        self.assertEqual(gs.guids, ["b1", "b3"])
        gs.add("b3")
        self.assertEqual(gs.guids, ["b1", "b3"])
        gs.add("b2")
        self.assertEqual(gs.guids, ["b1", "b2", "b3"])
        gs.add("a1")
        self.assertEqual(gs.guids, ["a1", "b1", "b2", "b3"])
        gs.add("c1")
        self.assertEqual(gs.guids, ["a1", "b1", "b2", "b3", "c1"])


class test_gm_2(unittest.TestCase):
    def runTest(self):

        gs = guidSearcher()
        gs.add("b1")
        gs.add("b3")
        gs.add("b3")
        gs.add("b2")
        gs.add("a1")
        gs.add("c1")

        retVal = gs.search("b")
        self.assertEqual(retVal, ["b1", "b2", "b3"])
        retVal = gs.search("b", max_returned=30)
        self.assertEqual(retVal, ["b1", "b2", "b3"])
        retVal = gs.search("b", max_returned=1)

        self.assertEqual(retVal, ["b1"])
        retVal = gs.search("b", max_returned=1, return_subset=True)
        self.assertEqual(len(retVal), 1)
        retVal = gs.search("b", max_returned=1, return_subset=False)
        self.assertEqual(len(retVal), 0)

        retVal = gs.search("b", max_returned=2, return_subset=True)
        self.assertEqual(len(retVal), 2)
        retVal = gs.search("b", max_returned=2, return_subset=False)
        self.assertEqual(len(retVal), 0)
        retVal = gs.search("z")
        self.assertEqual(retVal, [])


class test_gm_3(unittest.TestCase):
    def runTest(self):

        gs = guidSearcher()
        gs.add("b1")

        retVal = gs.search("b1")
        self.assertEqual(len(retVal), 1)


@unittest.skip("benchmark disabled")
class test_gm_benchmark(unittest.TestCase):
    def runTest(self):
        """tests time to search a list of 1M guids

        median addition time : 0.00026 sec (0.26 ms)
        median search time   : 0.00002 sec (0.02 ms)"""

        gs = guidSearcher()

        print("generating list of guids")

        added = []
        for i in range(1000000):
            new_guid = str(uuid.uuid4())
            added.append(new_guid)

        t2 = datetime.datetime.now()
        print("adding")
        for item in added:
            gs.add(item)

        t3 = datetime.datetime.now()
        print("searching")
        for item in added:
            gs.search(item)

        t4 = datetime.datetime.now()
        print("SEARCH PER SAMPLE", (t4 - t3) / len(added))
        print("ADD PER SAMPLE", (t3 - t2) / len(added))


# class for testing guidDbsearcher
class Test_Base1t(unittest.TestCase):
    """initialise FN persistence and adds data, 0.1 secs apart.  Used for testing queries examining order of recovery of samples."""

    def setUp(self):

        pm = Persistence()
        
        self.t = pm.get_storage_object(connString=UNITTEST_MONGOCONN, debug=2)
        self.t._delete_existing_data()
        dna = NucleicAcid()

        # add some sequences
        seqs = {"b1": "ACGT", "b2": "NACT", "b3": "TTTT", "a1": "CCCC", "c1": "TTTT"}
        for guid in seqs.keys():
            time.sleep(0.1)
            seq = seqs[guid]
            dna.examine(seq)
            self.t.guid_annotate(
                guid=guid, nameSpace="DNAQuality", annotDict=dna.composition
            )
        self.seqs = seqs


class test_gdm_1(Test_Base1t):
    def runTest(self):

        gs = guidDbSearcher(self.t)
        self.assertEqual(gs.guids, ["a1", "b1", "b2", "b3", "c1"])


class test_gdm_2(Test_Base1t):
    def runTest(self):

        gs = guidDbSearcher(self.t)

        retVal = gs.search("b")
        self.assertEqual(retVal, ["b1", "b2", "b3"])
        retVal = gs.search("b", max_returned=30)
        self.assertEqual(retVal, ["b1", "b2", "b3"])
        retVal = gs.search("b", max_returned=1)

        self.assertEqual(retVal, ["b1"])
        retVal = gs.search("b", max_returned=1, return_subset=True)
        self.assertEqual(len(retVal), 1)
        retVal = gs.search("b", max_returned=1, return_subset=False)
        self.assertEqual(len(retVal), 0)

        retVal = gs.search("b", max_returned=2, return_subset=True)
        self.assertEqual(len(retVal), 2)
        retVal = gs.search("b", max_returned=2, return_subset=False)
        self.assertEqual(len(retVal), 0)
        retVal = gs.search("z")
        self.assertEqual(retVal, [])


class test_gdm_3(Test_Base1t):
    def runTest(self):

        gs = guidDbSearcher(self.t, recheck_interval_seconds=1)

        retVal = gs.search("b")
        self.assertEqual(retVal, ["b1", "b2", "b3"])

        dna = NucleicAcid()

        # add more
        # add some sequences
        seqs = {"b4": "ACGT", "b5": "NACT", "b6": "TTTT", "c1": "TTTT"}
        for guid in seqs.keys():
            time.sleep(0.1)
            seq = seqs[guid]
            dna.examine(seq)
            self.t.guid_annotate(
                guid=guid, nameSpace="DNAQuality", annotDict=dna.composition
            )

        # has not searched recently - won't find the new ones
        retVal = gs.search("b")
        self.assertEqual(retVal, ["b1", "b2", "b3"])

        time.sleep(1)

        # should now find the new ones
        retVal = gs.search("b")
        self.assertEqual(retVal, ["b1", "b2", "b3", "b4", "b5", "b6"])


class test_gdm_4(Test_Base1t):
    def runTest(self):
        """Tests add method in guidDbSearcher, which should return a NotImplementedError"""
        gs = guidDbSearcher(self.t)
        with self.assertRaises(NotImplementedError):
            gs.add("b1")


# class for testing guidDbsearcher using rdbms
class Test_Base1tr(unittest.TestCase):
    """initialise FN persistence and adds data, 0.1 secs apart.  Used for testing queries examining order of recovery of samples."""

    def setUp(self):

        pm = Persistence()

        self.t = pm.get_storage_object(connString=UNITTEST_RDBMSCONN, debug=2)

        dna = NucleicAcid()

        # add some sequences
        seqs = {"b1": "ACGT", "b2": "NACT", "b3": "TTTT", "a1": "CCCC", "c1": "TTTT"}
        for guid in seqs.keys():
            time.sleep(0.1)
            seq = seqs[guid]
            dna.examine(seq)
            self.t.refcompressedseq_store(guid, {"seq": seq, "invalid": 0})
            self.t.guid_annotate(
                guid=guid, nameSpace="DNAQuality", annotDict=dna.composition
            )
        self.seqs = seqs


class test_gdm_1r(Test_Base1t):
    def runTest(self):

        gs = guidDbSearcher(self.t)
        self.assertEqual(gs.guids, ["a1", "b1", "b2", "b3", "c1"])


class test_gdm_2r(Test_Base1tr):
    def runTest(self):

        gs = guidDbSearcher(self.t)

        retVal = gs.search("b")
        self.assertEqual(retVal, ["b1", "b2", "b3"])
        retVal = gs.search("b", max_returned=30)
        self.assertEqual(retVal, ["b1", "b2", "b3"])
        retVal = gs.search("b", max_returned=1)

        self.assertEqual(retVal, ["b1"])
        retVal = gs.search("b", max_returned=1, return_subset=True)
        self.assertEqual(len(retVal), 1)
        retVal = gs.search("b", max_returned=1, return_subset=False)
        self.assertEqual(len(retVal), 0)

        retVal = gs.search("b", max_returned=2, return_subset=True)
        self.assertEqual(len(retVal), 2)
        retVal = gs.search("b", max_returned=2, return_subset=False)
        self.assertEqual(len(retVal), 0)
        retVal = gs.search("z")
        self.assertEqual(retVal, [])


class test_gdm_3r(Test_Base1tr):
    def runTest(self):

        gs = guidDbSearcher(self.t, recheck_interval_seconds=1)

        retVal = gs.search("b")
        self.assertEqual(retVal, ["b1", "b2", "b3"])

        dna = NucleicAcid()

        # add more
        # add some sequences
        seqs = {"b4": "ACGT", "b5": "NACT", "b6": "TTTT"}
        for guid in seqs.keys():
            time.sleep(0.1)
            seq = seqs[guid]
            dna.examine(seq)
            self.t.refcompressedseq_store(guid, {"seq": seq, "invalid": 0})
            self.t.guid_annotate(
                guid=guid, nameSpace="DNAQuality", annotDict=dna.composition
            )

        # has not searched recently - won't find the new ones
        retVal = gs.search("b")
        self.assertEqual(retVal, ["b1", "b2", "b3"])

        time.sleep(1)

        # should now find the new ones
        retVal = gs.search("b")
        self.assertEqual(retVal, ["b1", "b2", "b3", "b4", "b5", "b6"])


class test_gdm_4r(Test_Base1tr):
    def runTest(self):
        """Tests add method in guidDbSearcher, which should return a NotImplementedError"""
        gs = guidDbSearcher(self.t)
        with self.assertRaises(NotImplementedError):
            gs.add("b1")
