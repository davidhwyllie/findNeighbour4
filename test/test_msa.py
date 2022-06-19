""" tests msa.py

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@ukhsa.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

 

"""
import os
import unittest
import time
import pandas as pd
from findn.msa import MSAResult, MSAStore
from findn.persistence import Persistence

# unittests
UNITTEST_MONGOCONN: str = "mongodb://localhost"
UNITTEST_RDBMSCONN: str = "sqlite://"


class Test_MSA(unittest.TestCase):
    """tests the MSAResult class"""

    def runTest(self):
        inputdict = {
            "fconst": {},
            "variant_positions": [0, 1, 2, 3],
            "invalid_guids": [],
            "valid_guids": [
                "AAACGY-1",
                "CCCCGY-2",
                "TTTCGY-3",
                "GGGGGY-4",
                "NNNCGY-5",
                "ACTCGY-6",
                "TCTQGY-7",
                "AAACGY-8",
            ],
            "expected_p1": 0.16666666666666666,
            "sample_size": 30,
            "df_dict": {
                "AAACGY-1": {
                    "aligned_seq": "AAAC",
                    "aligned_mseq": "AAAC",
                    "aligned_seq_len": 4,
                    "allN": 0,
                    "alignN": 0,
                    "allM": 1,
                    "alignM": 0,
                    "allN_or_M": 1,
                    "alignN_or_M": 0,
                    "p_value1": 1.0,
                    "p_value2": 1.0,
                    "p_value3": 1.0,
                    "p_value4": 1.0,
                    "observed_proportion": 0.0,
                    "expected_proportion1": 0.16666666666666666,
                    "expected_proportion2": 0.5,
                    "expected_proportion3": 0.5,
                    "expected_proportion4": 0.0,
                    "what_tested": "M",
                },
                "CCCCGY-2": {
                    "aligned_seq": "CCCC",
                    "aligned_mseq": "CCCC",
                    "aligned_seq_len": 4,
                    "allN": 0,
                    "alignN": 0,
                    "allM": 1,
                    "alignM": 0,
                    "allN_or_M": 1,
                    "alignN_or_M": 0,
                    "p_value1": 1.0,
                    "p_value2": 1.0,
                    "p_value3": 1.0,
                    "p_value4": 1.0,
                    "observed_proportion": 0.0,
                    "expected_proportion1": 0.16666666666666666,
                    "expected_proportion2": 0.5,
                    "expected_proportion3": 0.5,
                    "expected_proportion4": 0.0,
                    "what_tested": "M",
                },
            },
            "what_tested": "M",
            "outgroup": None,
            "creation_time": "2019-11-17T23:46:00.098151",
        }

        m = MSAResult(**inputdict)
        self.assertEqual(type(m.msa_fasta()), str)
        self.assertEqual(type(m.msa_html()), str)
        self.assertEqual(type(m.msa_dict()), dict)
        self.assertEqual(type(m.serialise()), dict)
        self.assertEqual(type(m.msa_interactive_depiction()), str)
        self.assertEqual(type(m.df), pd.DataFrame)
        self.assertEqual(
            m.token, "msa|M|no_og|ddc4781ec984b66b0b5bf006a71b29cf1f523740"
        )


# skip these tests if the NO_MONGO_TESTS variable exists
mongo_test = unittest.skipIf(
    os.environ.get("NO_MONGO_TESTS", False), "no mongo tests performed"
)


@mongo_test
class Test_MSAStore_mongodb(unittest.TestCase):
    """tests the MSAStore class"""

    def runTest(self):

        inputdict1 = {
            "fconst": {},
            "variant_positions": [0, 1, 2, 3],
            "invalid_guids": [],
            "valid_guids": ["AAACGY-1", "CCCCGY-2"],
            "expected_p1": 0.16666666666666666,
            "sample_size": 30,
            "df_dict": {
                "AAACGY-1": {
                    "aligned_seq": "AAMM",
                    "aligned_mseq": "AAYR",
                    "aligned_seq_len": 4,
                    "allN": 0,
                    "alignN": 0,
                    "allM": 3,
                    "alignM": 2,
                    "allN_or_M": 1,
                    "alignN_or_M": 0,
                    "p_value1": 1.0,
                    "p_value2": 1.0,
                    "p_value3": 1.0,
                    "p_value4": 1.0,
                    "observed_proportion": 0.0,
                    "expected_proportion1": 0.16666666666666666,
                    "expected_proportion2": 0.5,
                    "expected_proportion3": 0.5,
                    "expected_proportion4": 0.0,
                    "what_tested": "M",
                },
                "CCCCGY-2": {
                    "aligned_seq": "CCCC",
                    "aligned_mseq": "CCCC",
                    "aligned_seq_len": 4,
                    "allN": 0,
                    "alignN": 0,
                    "allM": 1,
                    "alignM": 0,
                    "allN_or_M": 1,
                    "alignN_or_M": 0,
                    "p_value1": 1.0,
                    "p_value2": 1.0,
                    "p_value3": 1.0,
                    "p_value4": 1.0,
                    "observed_proportion": 0.0,
                    "expected_proportion1": 0.16666666666666666,
                    "expected_proportion2": 0.5,
                    "expected_proportion3": 0.5,
                    "expected_proportion4": 0.0,
                    "what_tested": "M",
                },
            },
            "what_tested": "M",
            "outgroup": None,
            "creation_time": "2019-11-17T23:46:00.098151",
        }

        inputdict2 = {
            "fconst": {},
            "variant_positions": [0, 1, 2, 3],
            "invalid_guids": [],
            "valid_guids": ["AAACGY-1", "CCCCGY-2"],
            "expected_p1": 0.16666666666666666,
            "sample_size": 30,
            "df_dict": {
                "AAACGY-1": {
                    "aligned_seq": "AAAM",
                    "aligned_mseq": "AAAR",
                    "aligned_seq_len": 4,
                    "allN": 0,
                    "alignN": 0,
                    "allM": 1,
                    "alignM": 1,
                    "allN_or_M": 1,
                    "alignN_or_M": 0,
                    "p_value1": 1.0,
                    "p_value2": 1.0,
                    "p_value3": 1.0,
                    "p_value4": 1.0,
                    "observed_proportion": 0.0,
                    "expected_proportion1": 0.16666666666666666,
                    "expected_proportion2": 0.5,
                    "expected_proportion3": 0.5,
                    "expected_proportion4": 0.0,
                    "what_tested": "M",
                },
                "CCCCGY-2": {
                    "aligned_seq": "CCCC",
                    "aligned_mseq": "CCCC",
                    "aligned_seq_len": 4,
                    "allN": 0,
                    "alignN": 0,
                    "allM": 1,
                    "alignM": 0,
                    "allN_or_M": 1,
                    "alignN_or_M": 0,
                    "p_value1": 1.0,
                    "p_value2": 1.0,
                    "p_value3": 1.0,
                    "p_value4": 1.0,
                    "observed_proportion": 0.0,
                    "expected_proportion1": 0.16666666666666666,
                    "expected_proportion2": 0.5,
                    "expected_proportion3": 0.5,
                    "expected_proportion4": 0.0,
                    "what_tested": "M",
                },
            },
            "what_tested": "M",
            "outgroup": None,
            "creation_time": "2019-11-17T23:46:00.098151",
        }

        m1 = MSAResult(**inputdict1)
        m2 = MSAResult(**inputdict2)

        guids1 = m1.valid_guids
        guids2 = m2.valid_guids

        pm = Persistence()
        p = pm.get_storage_object(connString=UNITTEST_MONGOCONN, debug=2, verbose=True)
        ms = MSAStore(p, in_ram_persistence_time=1)

        t1 = ms.get_token("M", False, guids1)
        t2 = ms.get_token("M", False, guids2)

        self.assertFalse(ms.is_in_ram(t1))
        self.assertFalse(ms.is_in_ram(t2))

        ms.cache_in_ram(token=t1, msa_result=m1)
        ms.cache_in_ram(token=t2, msa_result=m2)
        self.assertTrue(ms.is_in_ram(t1))
        self.assertTrue(ms.is_in_ram(t2))

        # test purging of expired samples
        time.sleep(2)
        ms.purge_ram()

        self.assertFalse(ms.is_in_ram(t1))
        self.assertFalse(ms.is_in_ram(t2))

        # store on disc
        ms.persist(t1, m1)
        ms.persist(t2, m2)

        m1r = ms.load(t1)
        m2r = ms.load(t2)

        self.assertEqual(m1r.valid_guids, m1.valid_guids)
        self.assertEqual(m2r.valid_guids, m2.valid_guids)

        ms.unpersist([])

        m1r = ms.load(t1)
        m2r = ms.load(t2)

        self.assertIsNone(m1r)
        self.assertIsNone(m2r)


rdbms_test = unittest.skipIf(os.environ.get("NO_RDBMS_TESTS", False), "No rdbms tests")


@rdbms_test
class Test_MSAStore_rdbms(unittest.TestCase):
    """tests the MSAStore class"""

    def runTest(self):

        inputdict1 = {
            "fconst": {},
            "variant_positions": [0, 1, 2, 3],
            "invalid_guids": [],
            "valid_guids": ["AAACGY-1", "CCCCGY-2"],
            "expected_p1": 0.16666666666666666,
            "sample_size": 30,
            "df_dict": {
                "AAACGY-1": {
                    "aligned_seq": "AAMM",
                    "aligned_mseq": "AAYR",
                    "aligned_seq_len": 4,
                    "allN": 0,
                    "alignN": 0,
                    "allM": 3,
                    "alignM": 2,
                    "allN_or_M": 1,
                    "alignN_or_M": 0,
                    "p_value1": 1.0,
                    "p_value2": 1.0,
                    "p_value3": 1.0,
                    "p_value4": 1.0,
                    "observed_proportion": 0.0,
                    "expected_proportion1": 0.16666666666666666,
                    "expected_proportion2": 0.5,
                    "expected_proportion3": 0.5,
                    "expected_proportion4": 0.0,
                    "what_tested": "M",
                },
                "CCCCGY-2": {
                    "aligned_seq": "CCCC",
                    "aligned_mseq": "CCCC",
                    "aligned_seq_len": 4,
                    "allN": 0,
                    "alignN": 0,
                    "allM": 1,
                    "alignM": 0,
                    "allN_or_M": 1,
                    "alignN_or_M": 0,
                    "p_value1": 1.0,
                    "p_value2": 1.0,
                    "p_value3": 1.0,
                    "p_value4": 1.0,
                    "observed_proportion": 0.0,
                    "expected_proportion1": 0.16666666666666666,
                    "expected_proportion2": 0.5,
                    "expected_proportion3": 0.5,
                    "expected_proportion4": 0.0,
                    "what_tested": "M",
                },
            },
            "what_tested": "M",
            "outgroup": None,
            "creation_time": "2019-11-17T23:46:00.098151",
        }

        inputdict2 = {
            "fconst": {},
            "variant_positions": [0, 1, 2, 3],
            "invalid_guids": [],
            "valid_guids": ["AAACGY-1", "CCCCGY-2"],
            "expected_p1": 0.16666666666666666,
            "sample_size": 30,
            "df_dict": {
                "AAACGY-1": {
                    "aligned_seq": "AAAM",
                    "aligned_mseq": "AAAR",
                    "aligned_seq_len": 4,
                    "allN": 0,
                    "alignN": 0,
                    "allM": 1,
                    "alignM": 1,
                    "allN_or_M": 1,
                    "alignN_or_M": 0,
                    "p_value1": 1.0,
                    "p_value2": 1.0,
                    "p_value3": 1.0,
                    "p_value4": 1.0,
                    "observed_proportion": 0.0,
                    "expected_proportion1": 0.16666666666666666,
                    "expected_proportion2": 0.5,
                    "expected_proportion3": 0.5,
                    "expected_proportion4": 0.0,
                    "what_tested": "M",
                },
                "CCCCGY-2": {
                    "aligned_seq": "CCCC",
                    "aligned_mseq": "CCCC",
                    "aligned_seq_len": 4,
                    "allN": 0,
                    "alignN": 0,
                    "allM": 1,
                    "alignM": 0,
                    "allN_or_M": 1,
                    "alignN_or_M": 0,
                    "p_value1": 1.0,
                    "p_value2": 1.0,
                    "p_value3": 1.0,
                    "p_value4": 1.0,
                    "observed_proportion": 0.0,
                    "expected_proportion1": 0.16666666666666666,
                    "expected_proportion2": 0.5,
                    "expected_proportion3": 0.5,
                    "expected_proportion4": 0.0,
                    "what_tested": "M",
                },
            },
            "what_tested": "M",
            "outgroup": None,
            "creation_time": "2019-11-17T23:46:00.098151",
        }

        m1 = MSAResult(**inputdict1)
        m2 = MSAResult(**inputdict2)

        guids1 = m1.valid_guids
        guids2 = m2.valid_guids

        pm = Persistence()
        p = pm.get_storage_object(connString=UNITTEST_RDBMSCONN, debug=2, verbose=True)
        ms = MSAStore(p, in_ram_persistence_time=1)

        t1 = ms.get_token("M", False, guids1)
        t2 = ms.get_token("M", False, guids2)

        self.assertFalse(ms.is_in_ram(t1))
        self.assertFalse(ms.is_in_ram(t2))

        ms.cache_in_ram(token=t1, msa_result=m1)
        ms.cache_in_ram(token=t2, msa_result=m2)
        self.assertTrue(ms.is_in_ram(t1))
        self.assertTrue(ms.is_in_ram(t2))

        # test purging of expired samples
        time.sleep(2)
        ms.purge_ram()

        self.assertFalse(ms.is_in_ram(t1))
        self.assertFalse(ms.is_in_ram(t2))

        # store on disc
        ms.persist(t1, m1)
        ms.persist(t2, m2)

        m1r = ms.load(t1)
        m2r = ms.load(t2)

        self.assertEqual(m1r.valid_guids, m1.valid_guids)
        self.assertEqual(m2r.valid_guids, m2.valid_guids)

        ms.unpersist([])

        m1r = ms.load(t1)
        m2r = ms.load(t2)

        self.assertIsNone(m1r)
        self.assertIsNone(m2r)
