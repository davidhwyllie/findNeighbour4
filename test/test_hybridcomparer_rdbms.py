""" tests hybridcomparer.py using rdbms

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.  See see <https://www.gnu.org/licenses/>.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

"""

import unittest
from findn.hybridComparer import hybridComparer
from catwalk.pycw_client import CatWalk
from findn.msa import MSAResult
import os
from Bio import SeqIO

## persistence unit tests
UNITTEST_RDBMSCONN = "sqlite://"  # "unittest_oracle"

# skip these tests if the NORDBMSTESTS variable exists
rdbms_test = unittest.skipIf(
    os.environ.get("NO_RDBMS_TESTS", False), "Non graphical tests only"
)


@rdbms_test
class test_hybridComparer_update_preComparer_settings(unittest.TestCase):
    """tests update preComparer settings"""

    def runTest(self):
        # generate compressed sequences

        refSeq = "GGGGGG"

        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
            unittesting=True,
        )

        n = 0
        originals = ["AAACGN", "CCCCGN"]
        guids = []
        for original in originals:
            n += 1
            c = sc.compress(original)

            guid = "{0}-{1}".format(original, n)
            guids.append(guid)
            sc.persist(c, guid=guid)

        preComparer_settings = {
            "selection_cutoff": 20,
            "over_selection_cutoff_ignore_factor": 10,
            "uncertain_base": "M",
        }

        sc.PERSIST.config_store("preComparer", preComparer_settings)


@rdbms_test
class test_hybridComparer_mcompare(unittest.TestCase):
    """tests mcompare"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
            unittesting=True,
        )

        sc.PERSIST._delete_existing_data()
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
        res = res["neighbours"]
        self.assertEqual(len(res), len(originals) - 1)


@rdbms_test
class test_hybridComparer_summarise_stored_items(unittest.TestCase):
    """tests reporting on stored contents"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
            unittesting=True,
        )
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
        self.assertTrue("server|pcstat|nSeqs" in set(res.keys()))


@rdbms_test
class test_hybridComparer_48(unittest.TestCase):
    """tests computations of p values from exact bionomial test"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
            unittesting=True,
        )
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


@rdbms_test
class test_hybridComparer_47b3(unittest.TestCase):
    """tests generation of a multisequence alignment with
    testing for the correct selection of non-variant bases"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=6,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
            unittesting=True,
        )
        originals = [
            "AAACGR",
            "CCCCGY",
            "TTTCGN",
            "RGGGGN",
            "YNNCGN",
            "ACTCGN",
            "TCTNGN",
            "AAACNN",
            "CCCCNN",
        ]

        # there's variation at positions 0,1,2,3
        #'AAACGR'
        #'CCCCGY',
        #'TTTCGN'
        #'GGGGGN'
        #'RNNCGN'
        #'YCTCGN'
        #'TCTNGN'
        #'AAACNN'
        #'CCCCNN',
        guid_names = []
        n = 0
        for original in originals:
            n += 1
            c = sc.compress(original)
            this_guid = "{0}".format(original)
            sc.persist(c, guid=this_guid)
            guid_names.append(this_guid)

        msa = sc.multi_sequence_alignment(guid_names[0:8])

        self.assertEqual(msa.variant_positions, [0, 1, 2, 3])
        for i, ix in enumerate(msa.df.index):
            self.assertEqual(ix[:6], originals[i])  # the sequence id is the sequence
            self.assertEqual(
                msa.df.loc[ix, "aligned_mseq"], originals[i][:4]
            )  # first 4


@rdbms_test
class test_hybridComparer_47c(unittest.TestCase):
    """tests generation of a multisequence alignment with
    testing for the proportion of Ns.
    Tests situation with externally supplied _p1"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=8,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
            unittesting=True,
        )
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

        # but with expected_N supplied;
        msa = sc.multi_sequence_alignment(guid_names[0:8], expected_p1=0.995)
        # there's variation at positions 0,1,2,3
        self.assertTrue(isinstance(msa, MSAResult))
        expected_cols = set(
            [
                "what_tested",
                "aligned_seq",
                "aligned_mseq",
                "aligned_seq_len",
                "aligned_seq_len",
                "allN",
                "alignN",
                "allM",
                "alignM",
                "allN_or_M",
                "alignN_or_M",
                "p_value1",
                "p_value2",
                "p_value3",
                "p_value4",
                "observed_proportion",
                "expected_proportion1",
                "expected_proportion2",
                "expected_proportion3",
                "expected_proportion4",
            ]
        )
        self.assertEqual(set(msa.df.columns.values), expected_cols)

        self.assertEqual(len(msa.df.index), 8)

        msa = sc.multi_sequence_alignment(guid_names[0:8], expected_p1=0.995)

        self.assertEqual(set(msa.df.columns.values), expected_cols)

        self.assertEqual(
            set(msa.df.index.tolist()),
            set(
                [
                    "AAACGN-1",
                    "NNNCGN-5",
                    "CCCCGN-2",
                    "TTTCGN-3",
                    "GGGGGN-4",
                    "ACTCGN-6",
                    "TCTNGN-7",
                    "AAACGN-8",
                ]
            ),
        )
        self.assertTrue(
            msa.df.loc["AAACGN-1", "expected_proportion1"] is not None
        )  # check it computed a value
        self.assertEqual(
            msa.df.loc["AAACGN-1", "expected_proportion1"], 0.995
        )  # check is used the value passed


@rdbms_test
class test_hybridComparer_47b2(unittest.TestCase):
    """tests generation of a multisequence alignment with
    testing for the proportion of Ms.
    Tests all three outputs."""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=6,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
            unittesting=True,
        )
        # need > 30 sequences
        originals = [
            "AAACGY",
            "CCCCGY",
            "TTTCGY",
            "GGGGGY",
            "NNNCGY",
            "ACTCGY",
            "TCTQGY",
            "AAACGY",
            "CCCCGY",
            "TTTCGY",
            "GGGGGY",
            "NNNCGY",
            "ACTCGY",
            "TCTNGY",
            "AAACGY",
            "CCCCGY",
            "TTTCGY",
            "GGGGGY",
            "NNNCGY",
            "ACTCGY",
            "TCTNGY",
            "AAACGY",
            "CCCCGY",
            "TTTCGY",
            "GGGGGY",
            "NNNCGY",
            "ACTCGY",
            "TCTNGY",
            "AAACGY",
            "CCCCGY",
            "TTTCGY",
            "GGGGGY",
            "NNNCGY",
            "ACTCGY",
            "TCTNGY",
            "AAACGY",
            "CCCCGY",
            "TTTCGY",
            "GGGGGY",
            "NNNCGY",
            "ACTCGY",
            "TCTNGY",
        ]

        guid_names = []
        n = 0
        for original in originals:
            n += 1
            c = sc.compress(original)
            this_guid = "{0}-{1}".format(original, n)
            sc.persist(c, guid=this_guid)
            guid_names.append(this_guid)

        msa = sc.multi_sequence_alignment(guid_names[0:8], uncertain_base_type="M")

        # there's variation at positions 0,1,2,3
        #'AAACGY'
        #'CCCCGY',
        #'TTTCGY'
        #'GGGGGY'
        #'NNNCGY'
        #'ACTCGY'
        #'TCTNGY'
        #'AAACGY'
        #'CCCCGY',
        self.assertEqual(msa.variant_positions, [0, 1, 2, 3])

        # there's variation at positions 0,1,2,3
        expected_cols = set(
            [
                "what_tested",
                "aligned_seq",
                "aligned_mseq",
                "aligned_seq_len",
                "aligned_seq_len",
                "allN",
                "alignN",
                "allM",
                "alignM",
                "allN_or_M",
                "alignN_or_M",
                "p_value1",
                "p_value2",
                "p_value3",
                "p_value4",
                "observed_proportion",
                "expected_proportion1",
                "expected_proportion2",
                "expected_proportion3",
                "expected_proportion4",
            ]
        )

        self.assertEqual(set(msa.df.columns.values), expected_cols)

        self.assertEqual(len(msa.df.index), 8)
        self.assertTrue(
            msa.df.loc["AAACGY-1", "expected_proportion1"] is not None
        )  # check it computed a value
        self.assertEqual(
            set(msa.df.index.tolist()),
            set(
                [
                    "AAACGY-1",
                    "NNNCGY-5",
                    "CCCCGY-2",
                    "TTTCGY-3",
                    "GGGGGY-4",
                    "ACTCGY-6",
                    "TCTQGY-7",
                    "AAACGY-8",
                ]
            ),
        )


@rdbms_test
class test_hybridComparer_47b(unittest.TestCase):
    """tests generation of a multisequence alignment with
    testing for the proportion of Ns.
    Tests all three outputs."""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=6,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
            unittesting=True,
        )
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

        msa = sc.multi_sequence_alignment(guid_names[0:8])

        # there's variation at positions 0,1,2,3
        #'AAACGN'
        #'CCCCGN',
        #'TTTCGN'
        #'GGGGGN'
        #'NNNCGN'
        #'ACTCGN'
        #'TCTNGN'
        #'AAACGN'
        #'CCCCGN',
        self.assertEqual(msa.variant_positions, [0, 1, 2, 3])

        # there's variation at positions 0,1,2,3
        expected_cols = set(
            [
                "what_tested",
                "aligned_seq",
                "aligned_mseq",
                "aligned_seq_len",
                "aligned_seq_len",
                "allN",
                "alignN",
                "allM",
                "alignM",
                "allN_or_M",
                "alignN_or_M",
                "p_value1",
                "p_value2",
                "p_value3",
                "p_value4",
                "observed_proportion",
                "expected_proportion1",
                "expected_proportion2",
                "expected_proportion3",
                "expected_proportion4",
            ]
        )

        self.assertEqual(set(msa.df.columns.values), expected_cols)
        self.assertEqual(len(msa.df.index), 8)
        self.assertTrue(
            msa.df.loc["AAACGN-1", "expected_proportion1"] is not None
        )  # check it computed a value
        self.assertEqual(
            set(msa.df.index.tolist()),
            set(
                [
                    "AAACGN-1",
                    "CCCCGN-2",
                    "TTTCGN-3",
                    "GGGGGN-4",
                    "NNNCGN-5",
                    "ACTCGN-6",
                    "TCTNGN-7",
                    "AAACGN-8",
                ]
            ),
        )


@rdbms_test
class test_hybridComparer_estimate_expected_unk(unittest.TestCase):
    """tests estimate_expected_unk, a function estimating the number of Ns in sequences
    by sampling"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
            unittesting=True,
        )
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
        res = sc.estimate_expected_unk(
            sample_size=2, unk_type="N", exclude_guids=guids[0:5]
        )
        self.assertEqual(res, 1.5)

        # analyse the first two
        res = sc.estimate_expected_unk(
            sample_size=2, unk_type="N", exclude_guids=guids[2:7]
        )
        self.assertEqual(res, 1)


@rdbms_test
class test_hybridComparer_estimate_expected_unk_sites(unittest.TestCase):
    """tests estimate_expected_unk_sites, a function estimating the number of Ns or Ms in sequences
    by sampling"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
            unittesting=True,
        )
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

        # evaluate Ms
        res = sc.estimate_expected_unk_sites(
            unk_type="N", sites=set([]), sample_size=30
        )
        self.assertEqual(res, None)
        res = sc.estimate_expected_unk_sites(unk_type="N", sites=set([]), sample_size=7)
        self.assertEqual(res, 0)
        res = sc.estimate_expected_unk_sites(
            unk_type="N", sites=set([0, 1, 2, 3, 4, 5]), sample_size=7
        )
        self.assertEqual(res, 1)
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
            unittesting=True,
        )
        n = 0
        originals = [
            "AAACGM",
            "CCCCGM",
            "TTTCGM",
            "GGGGGM",
            "MMMCGM",
            "ACTCGM",
            "TCTMGM",
        ]
        guids = []
        for original in originals:
            n += 1
            c = sc.compress(original)
            guid = "{0}-{1}".format(original, n)
            guids.append(guid)
            sc.persist(c, guid=guid)

        # analyse
        res = sc.estimate_expected_unk_sites(unk_type="M", sites=set([]), sample_size=7)
        self.assertEqual(res, 0)
        res = sc.estimate_expected_unk_sites(
            unk_type="M", sites=set([0, 1, 2, 3, 4, 5]), sample_size=7
        )
        self.assertEqual(res, 1)


@rdbms_test
class test_hybridComparer_45a(unittest.TestCase):
    """tests the generation of multiple alignments of variant sites."""

    def runTest(self):

        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
            unittesting=True,
        )
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

        self.assertEqual(len(sc.pc.seqProfile.keys()), 7)

        msa = sc.multi_sequence_alignment(guid_names)

        # there's variation at positions 0,1,2,3
        self.assertEqual(msa.alignment_length, 4)
        self.assertEqual(msa.variant_positions, [0, 1, 2, 3])


@rdbms_test
class test_hybridComparer_45b(unittest.TestCase):
    """tests the generation of multiple alignments of variant sites."""

    def runTest(self):

        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=6,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
            unittesting=True,
        )

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

        msa = sc.multi_sequence_alignment(guid_names)

        # there's variation at positions 0,1,2,3
        self.assertEqual(msa.alignment_length, 4)
        self.assertEqual(msa.variant_positions, [0, 1, 2, 3])


@rdbms_test
class test_hybridComparer_1(unittest.TestCase):
    """test init"""

    def runTest(self):
        refSeq = "ACTG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )
        self.assertEqual(sc.reference, refSeq)


@rdbms_test
class test_hybridComparer_2(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )
        with self.assertRaises(TypeError):
            sc.compress(sequence="AC")


@rdbms_test
class test_hybridComparer_3(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )
        retVal = sc.compress(sequence="ACTG")
        self.assertEqual(
            retVal,
            {
                "G": set([]),
                "A": set([]),
                "C": set([]),
                "T": set([]),
                "N": set([]),
                "M": {},
                "invalid": 0,
            },
        )


@rdbms_test
class test_hybridComparer_3b(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )
        retVal = sc.compress(sequence="ACTQ")
        self.assertEqual(
            retVal,
            {
                "G": set([]),
                "A": set([]),
                "C": set([]),
                "T": set([]),
                "N": set([]),
                "M": {3: "Q"},
                "invalid": 0,
            },
        )


@rdbms_test
class test_hybridComparer_3c(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )
        retVal = sc.compress(sequence="NYTQ")
        self.assertEqual(
            retVal,
            {
                "G": set([]),
                "A": set([]),
                "C": set([]),
                "T": set([]),
                "N": set([0]),
                "M": {1: "Y", 3: "Q"},
                "invalid": 0,
            },
        )


@rdbms_test
class test_hybridComparer_4(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )

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
                "invalid": 0,
            },
        )


@rdbms_test
class test_hybridComparer_5(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )
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
                "invalid": 0,
            },
        )


@rdbms_test
class test_hybridComparer_6(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"

        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )

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
                "invalid": 0,
            },
        )


@rdbms_test
class test_hybridComparer_7(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"

        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )
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
                "invalid": 0,
            },
        )


@rdbms_test
class test_hybridComparer_6b(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"

        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )
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


@rdbms_test
class test_hybridComparer_6c(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"

        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
            unittesting=True,
        )
        originals = ["NNNN"]
        for original in originals:

            compressed_sequence = sc.compress(sequence=original)
            roundtrip = sc.uncompress(compressed_sequence)
            self.assertEqual(original, roundtrip)


@rdbms_test
class test_hybridComparer_6d(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"

        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=3,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
            unittesting=True,
        )
        originals = ["NNNN"]
        for original in originals:

            compressed_sequence = sc.compress(sequence=original)
            with self.assertRaises(ValueError):
                sc.uncompress(compressed_sequence)


@rdbms_test
class test_hybridComparer_16(unittest.TestCase):
    """tests the comparison of two sequences where both differ from the reference."""

    def runTest(self):
        # generate compressed sequences
        refSeq = "ACTG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
            unittesting=True,
        )

        seq1 = sc.compress("AAAA")
        seq2 = sc.compress("CCCC")
        self.assertEqual(sc.countDifferences("k1", "k2", seq1, seq2), ("k1", "k2", 4))


@rdbms_test
class test_hybridComparer_16b(unittest.TestCase):
    """tests the comparison of two sequences where both differ from the reference."""

    def runTest(self):
        # generate compressed sequences
        refSeq = "ACTG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
            unittesting=True,
        )

        seq1 = sc.compress("AAAA")
        seq2 = sc.compress("RRCC")
        self.assertEqual(sc.countDifferences("k1", "k2", seq1, seq2), ("k1", "k2", 2))


@rdbms_test
class test_hybridComparer_16c(unittest.TestCase):
    """tests the comparison of two sequences where both differ from the reference."""

    def runTest(self):
        # generate compressed sequences
        refSeq = "ACTG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            unittesting=True,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )

        seq1 = sc.compress("AAAA")
        seq2 = sc.compress("RRNN")
        self.assertEqual(sc.countDifferences("k1", "k2", seq1, seq2), ("k1", "k2", 0))


@rdbms_test
class test_hybridComparer_17(unittest.TestCase):
    """tests the comparison of two sequences where one is invalid"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "ACTG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=3,
            reference=refSeq,
            snpCeiling=10,
            unittesting=True,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )

        seq1 = sc.compress("AAAA")
        seq2 = sc.compress("NNNN")
        self.assertEqual(sc.countDifferences("k1", "k2", seq1, seq2), None)


@rdbms_test
class test_hybridComparer_18(unittest.TestCase):
    """tests the comparison of two sequences where one is invalid"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "ACTG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=3,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
            unittesting=True,
        )

        seq1 = sc.compress("AAAA")
        seq2 = sc.compress("NNNN")

        sc.persist(seq2, "k2")
        sc.persist(seq1, "k1")

        res = sc.mcompare("k1")
        res = res["neighbours"]
        self.assertEqual(len(res), 0)

        res = sc.mcompare("k2")
        res = res["neighbours"]
        self.assertEqual(len(res), 0)


@rdbms_test
class test_hybridComparer_saveload3(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTGGG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
            },
            unittesting=True,
        )
        compressedObj = sc.compress(sequence="ACTTMN")
        sc.persist(compressedObj, "one")
        retVal = sc.load(guid="one")

        self.assertEqual(compressedObj, retVal)
        self.assertEqual(
            len(sc.pc.seqProfile.keys()), 1
        )  # one entry in the preComparer

        compressedObj = sc.compress(sequence="ACTTTT")
        sc.persist(compressedObj, "one")  # should succeed, but add nothing
        self.assertEqual(
            len(sc.pc.seqProfile.keys()), 1
        )  # one entry in the preComparer

        compressedObj = sc.compress(sequence="ACTTTA")
        sc.persist(compressedObj, "two")
        retVal = sc.load(guid="two")
        self.assertEqual(compressedObj, retVal)
        self.assertEqual(
            len(sc.pc.seqProfile.keys()), 2
        )  # one entry in the preComparer


@rdbms_test
class test_hybridComparer_remove_all(unittest.TestCase):
    def runTest(self):
        refSeq = "ACTG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )
        compressedObj = sc.compress(sequence="ACTT")
        sc.persist(compressedObj, "one")
        sc_obj = sc.load(guid="one")
        sc.seqProfile["one"] = sc_obj
        self.assertEqual(
            len(sc.seqProfile.keys()), 1
        )  # 1 entry in the temporary  dictionary  either
        self.assertEqual(
            len(sc.pc.seqProfile.keys()), 1
        )  # one entry in the preComparer

        compressedObj = sc.compress(sequence="ACTT")
        sc.persist(compressedObj, "two")
        sc.load(guid="two")
        self.assertEqual(
            len(sc.pc.seqProfile.keys()), 2
        )  # two entry in the preComparer

        sc.remove_all_temporary_seqs()
        self.assertEqual(len(sc.pc.seqProfile.keys()), 2)  # 2 in preComparer
        self.assertEqual(
            len(sc.seqProfile.keys()), 0
        )  # no entry in the temporary dictioanry


@rdbms_test
class test_hybridComparer_24(unittest.TestCase):
    """tests N compression"""

    def runTest(self):

        refSeq = "ACTGTTAATTTTTTTTTGGGGGGGGGGGGAA"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )

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
                "invalid": 0,
            },
        )


@rdbms_test
class test_hybridComparer_29(unittest.TestCase):
    """tests _setStats"""

    def runTest(self):

        refSeq = "ACTGTTAATTTTTTTTTGGGGGGGGGGGGAA"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )
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


@rdbms_test
class test_hybridComparer_37(unittest.TestCase):
    """tests the loading of an exclusion file"""

    def runTest(self):

        # default exclusion file
        refSeq = "ACTG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=1,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )
        self.assertEqual(
            sc.excluded_hash(), "Excl 0 nt [d751713988987e9331980363e24189ce]"
        )


@rdbms_test
class test_hybridComparer_38(unittest.TestCase):
    """tests the loading of an exclusion file"""

    def runTest(self):

        # no exclusion file
        refSeq = "ACTG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=1,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )
        self.assertEqual(
            sc.excluded_hash(), "Excl 0 nt [d751713988987e9331980363e24189ce]"
        )


@rdbms_test
class test_hybridComparer_40(unittest.TestCase):
    """tests the computation of a hash of a compressed object"""

    def runTest(self):

        # generate compressed sequences
        refSeq = "ACTG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )
        compressed_sequence = sc.compress(sequence="TTAA")

        res = sc.compressed_sequence_hash(compressed_sequence)
        self.assertEqual(res, "6ce0e55c4ab092f560e03c5d2de53098")


@rdbms_test
class test_hybridComparer_45(unittest.TestCase):
    """tests insertion of large sequences"""

    def runTest(self):
        inputfile = "reference/NC_000962.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                goodseq = str(record.seq)
                badseq = "".join("N" * len(goodseq))
                originalseq = list(str(record.seq))
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            reference=record.seq,
            snpCeiling=100,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )
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


@rdbms_test
class test_hybridComparer_47(unittest.TestCase):
    """tests raise_error"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGGGGGGGG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )
        with self.assertRaises(ZeroDivisionError):
            sc.raise_error("token")


@rdbms_test
class test_hybridComparer_50(unittest.TestCase):
    """tests estimate_expected_proportion, a function computing the proportion of Ns expected based on the median
    Ns in a list of sequences"""

    def runTest(self):
        refSeq = "GGGGGGGGGGGG"
        sc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            unittesting=True,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )

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


@rdbms_test
class test_hybridComparer_51(unittest.TestCase):
    """tests repopulate_sample."""

    def runTest(self):

        # remove any temporary database

        sqlite_file = "unitTest_tmp/test_51.sqlite"
        if os.path.exists(sqlite_file):
            os.unlink(sqlite_file)

        # remove any temporary database
        refSeq = "GGGGGG"
        sc = hybridComparer(
            PERSIST="sqlite:///{0}".format(sqlite_file),
            maxNs=1e8,
            reference=refSeq,
            unittesting=False,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )

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

        self.assertEqual(len(sc.pc.seqProfile.keys()), 7)
        self.assertEqual(len(sc.PERSIST.guids()), 7)

        ## try to repopulate it
        refSeq = "GGGGGG"
        sc = hybridComparer(
            PERSIST="sqlite:///{0}".format(sqlite_file),
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            unittesting=False,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {},
            },
        )

        self.assertEqual(len(sc.pc.seqProfile.keys()), 0)
        self.assertEqual(len(sc.PERSIST.guids()), 7)

        sc.repopulate_sample(n=7)  # defaults to 100
        self.assertEqual(len(sc.pc.seqProfile.keys()), 7)


@rdbms_test
class test_hybridComparer_52(unittest.TestCase):
    """tests catwalk with uncertain_base = M"""

    def runTest(self):
        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                refSeq = str(record.seq)
                originalseq = list(str(record.seq))
        hc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
            },
            unittesting=True,
        )

        hc.PERSIST._delete_existing_data()

        inserted_guids = ["guid_ref"]
        obj = hc.compress(refSeq)
        hc.persist(obj, "guid_ref")

        for k in [0, 1]:
            for i in [0, 1, 2]:
                guid_to_insert = "msa2_{1}_guid_{0}".format(k * 100 + i, k)
                inserted_guids.append(guid_to_insert)
                muts = 0
                seq = originalseq.copy()
                # make  mutations
                if k == 1:
                    for j in range(1000000, 1000010):  # make 10 mutants at position 1m
                        mutbase = j
                        ref = seq[mutbase]
                        if not ref == "T":
                            seq[mutbase] = "T"
                        if not ref == "A":
                            seq[mutbase] = "A"
                        muts += 1

                seq = "".join(seq)
                obj = hc.compress(seq)
                hc.persist(obj, guid_to_insert)

        self.assertEqual(hc.pc.guids(), set(inserted_guids))

        for guid in inserted_guids:
            res = hc.mcompare(guid)
            self.assertEqual(res["timings"]["catWalk_enabled"], True)
            self.assertEqual(res["timings"]["preComparer_distances_are_exact"], False)
            self.assertTrue(res["timings"]["candidates"] > 0)


@rdbms_test
class test_hybridComparer_53(unittest.TestCase):
    """tests catwalk with uncertain_base = N_or_M.
    cf. hybridComparer test 11, which tests python vs. catwalk methods.

    functions correctly if catwalk is not already running and prepopulated"""

    def runTest(self):

        cw = CatWalk(
            cw_binary_filepath=None,
            reference_name="H37RV",
            reference_filepath="reference/TB-ref.fasta",
            mask_filepath="reference/TB-exclude-adaptive.txt",
            max_distance=20,
            bind_host="localhost",
            bind_port=5999,
        )

        # stop the server if it is running
        cw.stop()
        self.assertFalse(cw.server_is_running())

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                refSeq = str(record.seq)
                originalseq = list(str(record.seq))
        hc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=20,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "N_or_M",
                "over_selection_cutoff_ignore_factor": 1,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
            },
            unittesting=True,
        )

        hc.PERSIST._delete_existing_data()

        inserted_guids = ["guid_ref"]
        obj = hc.compress(refSeq)
        hc.persist(obj, "guid_ref")

        for k in [0, 1]:
            for i in [0, 1, 2]:
                guid_to_insert = "msa2_{1}_guid_{0}".format(k * 100 + i, k)
                inserted_guids.append(guid_to_insert)
                muts = 0
                seq = originalseq.copy()
                # make  mutations
                if k == 1:
                    for j in range(
                        500000, 500005
                    ):  # make 5 mutants at position 500000k
                        mutbase = j
                        ref = seq[mutbase]
                        if not ref == "T":
                            seq[mutbase] = "T"
                        else:
                            seq[mutbase] = "A"
                        muts += 1

                seq = "".join(seq)
                obj = hc.compress(seq)
                hc.persist(obj, guid_to_insert)

        self.assertEqual(hc.pc.guids(), set(inserted_guids))

        for guid in inserted_guids:
            res = hc.mcompare(guid)
            # print("UNITTEST 53: MCOMPARE TIMINGS", guid, res['timings'])
            # print("UNITTEST 53: MCOMPARE NEIGHBOURS", res['neighbours'])
            self.assertEqual(
                len(res["neighbours"]), len(inserted_guids) - 1
            )  # doesn't match self, but does match guid_ref
            self.assertEqual(res["timings"]["catWalk_enabled"], True)
            self.assertEqual(res["timings"]["preComparer_distances_are_exact"], True)
            self.assertTrue(res["timings"]["candidates"] == 0)  # doesn't report these
            self.assertEqual(res["timings"]["preCompared"], len(inserted_guids))
            self.assertTrue(
                res["timings"]["matches"] == len(inserted_guids) - 1
            )  # doesn't report self self matches

            # test everything was stored OK
            stored_guids = hc.PERSIST.guid2neighbours(
                guid, returned_format=3, cutoff=20
            )["neighbours"]
            # print("UNITTEST 53: STORED GUIDS LINKED TO ", guid, 'are', stored_guids)
            expected = set(inserted_guids) - set([guid])  # everything except itself
            # print("UNITTEST 53: EXPECTED_LINKS to ", guid, 'are', expected)
            # print("UNITTEST 53: STORED GUIDS to ", guid, 'are', stored_guids)
            self.assertEqual(expected, set(stored_guids))


@rdbms_test
class test_hybridComparer_54(unittest.TestCase):
    """tests catwalk with insertion disabled (catwalk should not start)"""

    def runTest(self):
        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                refSeq = str(record.seq)

        hc = hybridComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "N_or_M",
                "over_selection_cutoff_ignore_factor": 5,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
            },
            unittesting=True,
            disable_insertion=True,
        )
        self.assertFalse(hc.pc.catWalk_enabled)

        obj = hc.compress(refSeq)
        with self.assertRaises(NotImplementedError):
            hc.persist(obj, "guid_ref")
