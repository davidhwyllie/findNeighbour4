""" tests msaviewer.py

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@ukhsa.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

 

"""

import unittest
import pandas as pd
import collections
from findn.msaviewer import SimulateSequenceData, DepictMSA


class Test_SimulateSequenceData_1(unittest.TestCase):
    """tests make_seq"""

    def runTest(self):
        ssd = SimulateSequenceData()

        seq = ssd.make_seq(30)
        self.assertIsInstance(seq, str)
        self.assertEqual(30, len(seq))


class Test_SimulateSequenceData_2(unittest.TestCase):
    """tests mutate_seq"""

    def runTest(self):
        ssd = SimulateSequenceData()

        seq = ssd.make_seq(100)
        mseq = ssd.mutate_seq(seq)

        self.assertIsInstance(mseq, str)
        self.assertEqual(100, len(mseq))

        seq = ssd.make_seq(100)
        mseq = ssd.mutate_seq(seq, iupac=True)

        self.assertIsInstance(mseq, str)
        self.assertEqual(100, len(mseq))


class Test_SimulateSequenceData_3(unittest.TestCase):
    """tests msa generation"""

    def runTest(self):
        ssd = SimulateSequenceData()
        msa = ssd.make_msa(10, 40)
        self.assertIsInstance(msa, pd.DataFrame)
        self.assertEqual(len(msa.index), 10)


class Test_DepictMSA_1(unittest.TestCase):
    """tests msa object"""

    def runTest(self):
        ssd = SimulateSequenceData()
        nSeqs = 20
        alignLen = 200

        x_labels = ["seq{0}".format(x) for x in range(alignLen)]
        x_labels_int = [x for x in range(alignLen)]

        msa = ssd.make_msa(nSeqs, alignLen)
        dep_msa = DepictMSA(msa)
        self.assertEqual(dep_msa.align_width, alignLen)
        self.assertEqual(dep_msa.nSeqs, nSeqs)
        self.assertIsInstance(dep_msa.composition, collections.Counter)
        self.assertIsInstance(dep_msa.rectangles, pd.DataFrame)
        retVal = dep_msa.render_msa()

        self.assertIsInstance(retVal, str)

        retVal = dep_msa.render_msa()
        self.assertIsInstance(retVal, str)

        msa = ssd.make_msa(nSeqs, alignLen)
        dep_msa = DepictMSA(msa)

        self.assertEqual(dep_msa.align_width, alignLen)
        self.assertEqual(dep_msa.nSeqs, nSeqs)
        self.assertIsInstance(dep_msa.composition, collections.Counter)
        self.assertIsInstance(dep_msa.rectangles, pd.DataFrame)
        retVal = dep_msa.render_msa()
        self.assertIsInstance(retVal, str)

        dep_msa = DepictMSA(msa, positions_analysed=x_labels)
        retVal = dep_msa.render_msa()
        self.assertIsInstance(retVal, str)

        dep_msa = DepictMSA(msa, positions_analysed=x_labels_int)
        retVal = dep_msa.render_msa()
        self.assertIsInstance(retVal, str)

        with self.assertRaises(ValueError):
            dep_msa = DepictMSA(msa, identify_sequence_by=["no field"])
            retVal = dep_msa.render_msa()
            self.assertIsInstance(retVal, str)

        dep_msa = DepictMSA(msa, identify_sequence_by=["Surname", "Forename"])
        retVal = dep_msa.render_msa()
        self.assertIsInstance(retVal, str)
        self.assertTrue("</html>" in retVal)

        dep_msa = DepictMSA(
            msa, identify_sequence_by=["Surname", "Forename"], max_elements_in_plot=20
        )
        retVal = dep_msa.render_msa()
        self.assertIsInstance(retVal, str)
        self.assertTrue("</html>" in retVal)
