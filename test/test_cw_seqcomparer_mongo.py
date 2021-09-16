""" tests cw_seqComparer.py

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

"""

import unittest
from findn.cw_seqComparer import NoCWParametersProvidedError, cw_seqComparer
from Bio import SeqIO

## persistence unit tests
UNITTEST_MONGOCONN = "mongodb://localhost"
UNITTEST_RDBMSCONN = "sqlite://"


class setup_ref(unittest.TestCase):
    def setUp(self):

        self.inputfile = "reference/NC_000962.fasta"
        with open(self.inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                self.goodseq = str(record.seq)
                self.badseq = "".join("N" * len(self.goodseq))
                self.originalseq = list(str(record.seq))


class test_cw_seqComparer_estimate_expected_unk_1(setup_ref):
    """tests estimate_expected_unk, a function estimating the number of Ns in sequences
    by sampling"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = cw_seqComparer(
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
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
            },
            unittesting=True,
        )

        res = sc.estimate_expected_unk(
            sample_size=2, unk_type="N"
        )

        self.assertIsNone(res)
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

        res = sc.estimate_expected_unk(sample_size=30)
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

class test_cw_seqComparer_estimate_expected_unk_sites(setup_ref):
    """tests estimate_expected_unk_sites, a function estimating the number of Ns or Ms in sequences
    by sampling"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = cw_seqComparer(
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
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
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
        self.assertEqual(res, 1)
        res = sc.estimate_expected_unk_sites(
            unk_type="N", sites=set([0, 1, 2, 3, 4, 5]), sample_size=7
        )
        self.assertEqual(res, 1)
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = cw_seqComparer(
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
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
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
        self.assertEqual(res, 1)
        res = sc.estimate_expected_unk_sites(
            unk_type="M", sites=set([0, 1, 2, 3]), sample_size=7
        )
        self.assertEqual(res, 0)


class test_cw_seqComparer_msa_1(setup_ref):
    """tests the generation of multiple alignments of variant sites."""

    def runTest(self):

        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = cw_seqComparer(
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
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
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

        self.assertEqual(len(guid_names), 7)

        msa = sc.multi_sequence_alignment(guid_names)

        # there's variation at positions 0,1,2,3
        self.assertEqual(msa.alignment_length, 4)
        self.assertEqual(msa.variant_positions, [0, 1, 2, 3])


class test_cw_seqComparer_msa_b(setup_ref):
    """tests the generation of multiple alignments of variant sites."""

    def runTest(self):

        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = cw_seqComparer(
            maxNs=6,
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
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
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


class test_cw_seqComparer_repopulate(setup_ref):
    """tests persist, guids, mcompare and methods"""

    def runTest(self):
        # load catwalk and database
        refSeq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGA"
        sc = cw_seqComparer(
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 1,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
            },
            unittesting=True,
        )

        sc.PERSIST._delete_existing_data()

        s1 = sc.summarise_stored_items()
        self.assertEqual(s1["server|catwalk|n_samples"], 0)

        n = 0
        originals = [
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGN",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCGN",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTCGN",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGN",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNCGN",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTCGN",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATCTNGN",
        ]

        guids = []
        for original in originals:
            n += 1
            c = sc.compress(original)
            guid = "{0}-{1}".format(original, n)
            guids.append(guid)

            sc.persist(c, guid=guid)

        # check that all have inserted
        self.assertEqual(set(guids), set(sc.guids()))
        s2 = sc.summarise_stored_items()
        self.assertEqual(s2["server|catwalk|n_samples"], len(originals))
        self.assertEqual(len(sc.PERSIST.guids()), len(guids))
        self.assertEqual(len(sc.PERSIST.guids_valid()), len(guids))
        sc.catWalk.stop()  # terminate the catwalk instance
        sc.catWalk.start()
        sc.repopulate_all()
        s3 = sc.summarise_stored_items()
        self.assertEqual(s3["server|catwalk|n_samples"], len(originals))


class test_cw_seqComparer_mcompare(setup_ref):
    """tests persist, guids, mcompare and methods"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "GGGGGG"
        sc = cw_seqComparer(
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 1,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
            },
            unittesting=True,
        )

        sc.PERSIST._delete_existing_data()

        s1 = sc.summarise_stored_items()
        self.assertEqual(s1["server|catwalk|n_samples"], 0)

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

        # check that all have inserted
        self.assertEqual(set(guids), set(sc.guids()))
        s2 = sc.summarise_stored_items()
        self.assertEqual(s2["server|catwalk|n_samples"], len(originals))

        for test_guid in guids:
            res = sc.mcompare(test_guid)  # defaults to sample size 30
            res = res["neighbours"]
            self.assertEqual(len(res), len(originals) - 1)


class test_no_catwalk_settings(setup_ref):
    """tests initialisation with no catwalk"""

    def runTest(self):
        # generate compressed sequences

        refSeq = "GGGGGG"

        with self.assertRaises(NoCWParametersProvidedError):
            cw_seqComparer(
                maxNs=1e8,
                reference=refSeq,
                snpCeiling=10,
                preComparer_parameters={},
                unittesting=True,
            )


class test_invalid_inputs(setup_ref):
    """tests initialisation with invalid inputs"""

    def runTest(self):
        # generate compressed sequences

        refSeq = "GGGGGG"

        with self.assertRaises(TypeError):
            cw_seqComparer(
                maxNs=1e8,
                reference=refSeq,
                snpCeiling=10,
                preComparer_parameters={"catWalk_parameters": {}},
                unittesting=True,
                PERSIST=-1,
            )

        with self.assertRaises(TypeError):
            cw_seqComparer(
                maxNs=1e8,
                reference=refSeq,
                snpCeiling=10,
                preComparer_parameters={"catWalk_parameters": {}},
                unittesting=True,
                PERSIST=-1,
            )

        with self.assertRaises(TypeError):
            cw_seqComparer(
                maxNs=1e8,
                reference=refSeq,
                snpCeiling=10,
                preComparer_parameters={"catWalk_parameters": {}},
                unittesting=True,
                PERSIST=None,
            )


class test_cw_seqComparer_eh(setup_ref):
    """tests the loading of an exclusion file"""

    def runTest(self):

        # default exclusion file
        refSeq = "ACTG"
        sc = cw_seqComparer(
            unittesting=True,
            maxNs=1e8,
            reference=refSeq,
            snpCeiling=1,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 1,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
            },
        )
        self.assertEqual(
            sc.excluded_hash(), "Excl 0 nt [d751713988987e9331980363e24189ce]"
        )


class test_cw_seqComparer_compress_uncompress(setup_ref):
    def runTest(self):
        refSeq = "ACTG"

        sc = cw_seqComparer(
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 1,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
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


class test_cw_seqComparer_load_refseq(setup_ref):
    """test init"""

    def runTest(self):
        refSeq = "ACTG"
        sc = cw_seqComparer(
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 1,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
            },
        )
        self.assertEqual(sc.reference, refSeq)


class test_cw_seqComparer_wrong_length(setup_ref):
    def runTest(self):
        refSeq = "ACTG"
        sc = cw_seqComparer(
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 1,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
            },
        )
        with self.assertRaises(TypeError):
            sc.compress(sequence="AC")


class test_cw_seqComparer_is_ref(setup_ref):
    def runTest(self):
        refSeq = "ACTG"
        sc = cw_seqComparer(
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 1,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
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


class test_cw_seqComparer_store_m(setup_ref):
    def runTest(self):
        refSeq = "ACTG"
        sc = cw_seqComparer(
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 1,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
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


class test_cw_seqComparer_store_ms(setup_ref):
    def runTest(self):
        refSeq = "ACTG"
        sc = cw_seqComparer(
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 1,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
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


class test_cw_seqComparer_store_n(setup_ref):
    def runTest(self):
        refSeq = "ACTG"
        sc = cw_seqComparer(
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 1,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
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


class test_cw_seqComparer_store_dash_as_n(setup_ref):
    def runTest(self):
        refSeq = "ACTG"
        sc = cw_seqComparer(
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 1,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
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


class test_cw_seqComparer_store_n_and_var(setup_ref):
    def runTest(self):
        refSeq = "ACTG"

        sc = cw_seqComparer(
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 1,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
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


class test_cw_seqComparer_store_n_and_othervar(setup_ref):
    def runTest(self):
        refSeq = "ACTG"

        sc = cw_seqComparer(
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 1,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
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


class test_cw_seqComparer_compress_uncompress_roundtrip(setup_ref):
    def runTest(self):
        refSeq = "ACTG"

        sc = cw_seqComparer(
            unittesting=True,
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 1,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
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


class test_cw_seqComparer_compress_uncompress_roundtrip_2(setup_ref):
    def runTest(self):
        refSeq = "ACTG"

        sc = cw_seqComparer(
            maxNs=1e8,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 1,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
            },
            unittesting=True,
        )
        originals = ["NNNN"]
        for original in originals:

            compressed_sequence = sc.compress(sequence=original)
            roundtrip = sc.uncompress(compressed_sequence)
            self.assertEqual(original, roundtrip)


class test_cw_seqComparer_compress_uncompress_roundtrip_3(setup_ref):
    def runTest(self):
        refSeq = "ACTG"

        sc = cw_seqComparer(
            maxNs=3,
            snpCeiling=20,
            reference=refSeq,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 1,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
            },
            unittesting=True,
        )
        originals = ["NNNN"]
        for original in originals:

            compressed_sequence = sc.compress(sequence=original)
            with self.assertRaises(ValueError):
                sc.uncompress(compressed_sequence)


class test_cw_seqComparer_one_invalid_sequence(setup_ref):
    """tests the comparison of two sequences where one is invalid"""

    def runTest(self):
        # generate compressed sequences
        refSeq = "ACTG"
        sc = cw_seqComparer(
            maxNs=3,
            reference=refSeq,
            snpCeiling=10,
            preComparer_parameters={
                "selection_cutoff": 20,
                "uncertain_base": "M",
                "over_selection_cutoff_ignore_factor": 1,
                "catWalk_parameters": {
                    "bind_port": 5999,
                    "bind_host": "localhost",
                    "cw_binary_filepath": None,
                    "reference_name": "H37RV",
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
            },
            unittesting=True,
        )

        seq1 = sc.compress("AAAA")
        seq2 = sc.compress("NNNN")

        sc.persist(seq2, "k2")
        sc.persist(seq1, "k1")

        # should be in persistence
        self.assertEqual(sc.PERSIST.guid_valid("k1"), 0)
        self.assertEqual(sc.PERSIST.guid_valid("k2"), 1)

        res = sc.mcompare("k1")
        self.assertEqual(res["invalid"], 0)
        res = res["neighbours"]
        self.assertEqual(len(res), 0)

        res = sc.mcompare("k2")
        self.assertEqual(res["invalid"], 1)
        res = res["neighbours"]
        self.assertEqual(len(res), 0)


class test_cw_seqComparer_csh(setup_ref):
    """tests the computation of a hash of a compressed object"""

    def runTest(self):

        # generate compressed sequences
        refSeq = "ACTG"
        sc = cw_seqComparer(
            unittesting=True,
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
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
            },
        )
        compressed_sequence = sc.compress(sequence="TTAA")

        res = sc.compressed_sequence_hash(compressed_sequence)
        self.assertEqual(res, "6ce0e55c4ab092f560e03c5d2de53098")


class test_cw_seqComparer_epp(setup_ref):
    """tests estimate_expected_proportion, a function computing the proportion of Ns expected based on the median
    Ns in a list of sequences"""

    def runTest(self):
        refSeq = "GGGGGGGGGGGG"
        sc = cw_seqComparer(
            unittesting=True,
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
                    "reference_filepath": self.inputfile,
                    "mask_filepath": "reference/TB-exclude-adaptive.txt",
                },
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
