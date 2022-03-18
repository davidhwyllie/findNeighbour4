""" tests cw_seqComparer.py

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

"""
import os
import unittest
from localstore.localstoreutils import LocalStore
from findn.cw_seqComparer import cw_seqComparer
from Bio import SeqIO

## persistence unit tests
UNITTEST_MONGOCONN = "mongodb://localhost"
UNITTEST_RDBMSCONN = "sqlite://"


class setup_ref(unittest.TestCase):
    def setUp(self):

        self.refSeq = None
        self.inputfile = "reference/NC_000962.fasta"
        with open(self.inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                self.refSeq = str(record.seq)

        if self.refSeq is None:
            raise ValueError("did not read reference sequence")

class test_cw_seqComparer_fr_0(setup_ref):
    """test type checking by cw_seqComparer """

    def runTest(self):

        # create a list of 5 sequences 
        seq_list = []
        expected_sequence_ids = []
        for i in range(5):
            sequence_id = "SEQ_{0}".format(i)
            expected_sequence_ids.append(sequence_id)
            seq_list.append(
                (sequence_id, self.refSeq)
            )
        # create a tar archive
        tarfile_name = "unittest_tmp/test.tar"
        try:
            os.unlink(tarfile_name)
        except FileNotFoundError:
            pass

        self.assertEqual(False, os.path.exists(tarfile_name))

        # create the local data store
        ls = LocalStore(
            "unittest_tmp/test.tar",
            write_batch_size=3
        )

        # test different localstore information
        cw_seqComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            reference=self.refSeq,
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

        cw_seqComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            reference=self.refSeq,
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
            localstore =None,
            unittesting=True,
        )

        cw_seqComparer(
            PERSIST=UNITTEST_RDBMSCONN,
            maxNs=1e8,
            reference=self.refSeq,
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
            localstore =ls,
            unittesting=True,
        )

        with self.assertRaises(TypeError):
            cw_seqComparer(
                PERSIST=UNITTEST_RDBMSCONN,
                maxNs=1e8,
                reference=self.refSeq,
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
                localstore ="a string",
                unittesting=True,
            )
        

class test_cw_seqComparer_fr_1(setup_ref):
    """loads 5 sequences into cw_seqComparer """

    def runTest(self):

        # create a list of 5 sequences 
        seq_list = []
        expected_sequence_ids = []

        test_size = 10
        for i in range(test_size):
            sequence_id = "SEQ_{0}".format(i)
            expected_sequence_ids.append(sequence_id)
            seq_list.append(
                (sequence_id, self.refSeq)
            )
        # create a tar archive
        tarfile_name = "unittest_tmp/test.tar"
        try:
            os.unlink(tarfile_name)
        except FileNotFoundError:
            pass

        self.assertEqual(False, os.path.exists(tarfile_name))

        # create the local data store
        ls = LocalStore(
            "unittest_tmp/test.tar",
            write_batch_size=3
        )

        # empty database
        sc = cw_seqComparer(
            PERSIST=UNITTEST_MONGOCONN,
            maxNs=1e8,
            reference=self.refSeq,
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
            localstore =ls,
            unittesting=True,
        )

        # generate compressed sequences
        sc = cw_seqComparer(
            PERSIST=UNITTEST_MONGOCONN,
            maxNs=1e8,
            reference=self.refSeq,
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
            localstore =ls,
            unittesting=False,
        )

        # check nothing is in the database
        self.assertEqual(sc.PERSIST.guids(), set([]))

        # check there is nothing in the local store
        self.assertEqual(len(ls.sequence_ids()), 0)

        # add samples to sc
        print("Adding samples")
        i = 0
        for (sequence_id, sequence) in seq_list:
            i += 1 
            if i % 50 == 0:
                print("Loaded ", i)
            
            cs = sc.compress(self.refSeq)
            sc.persist(cs, sequence_id)

        # check the five samples have been added
        self.assertEqual(sc.PERSIST.guids(), set(expected_sequence_ids))

        # check there is nothing in the local store
        self.assertEqual(len(ls.sequence_ids()), 0)

        # synchronise localstore
        sc.update_localstore()

        # check there are test_size samples in the localstore
        self.assertEqual(len(ls.sequence_ids()), test_size)

        # forcibly terminate catwalk
        sc.catWalk.stop()

        # open another localstore
        ls = LocalStore(
            "unittest_tmp/test.tar"
        )

        # check there are samples in the localstore
        self.assertEqual(len(ls.sequence_ids()), test_size)

        # open another seqcomparer.  Do not wipe existing data (unittesting = False)
        sc = cw_seqComparer(
            PERSIST=UNITTEST_MONGOCONN,
            maxNs=1e8,
            reference=self.refSeq,
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
            localstore =ls,
            unittesting=False,
        )

        # check there are samples in the database
        self.assertEqual(sc.PERSIST.guids(), set(expected_sequence_ids))

        # check there are samples in the catwalk
        self.assertEqual(set(sc.catWalk.sample_names()), set(expected_sequence_ids))
