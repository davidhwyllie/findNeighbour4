""" tests cw_seqComparer.py

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@ukhsa.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

 

"""

import unittest

from findn.identify_sequence_set import IdentifySequenceSet


class Test_iss_1(unittest.TestCase):
    """tests the identifier() method"""

    def runTest(self):
        m = IdentifySequenceSet()
        self.assertEqual(
            m.make_identifier("msa", "N", True, [1, 2, 3]),
            "msa|N|og|cfd68e36493b6db0bc44c859e1e49290e07efb00",
        )
        self.assertEqual(
            m.make_identifier("msa", "N", False, [1, 2, 3]),
            "msa|N|no_og|cfd68e36493b6db0bc44c859e1e49290e07efb00",
        )


class Test_iss(unittest.TestCase):
    """tests the _hashComponents() method"""

    def runTest(self):
        m = IdentifySequenceSet()
        self.assertEqual(
            m._hashComponents([1, 2, 3]), "cfd68e36493b6db0bc44c859e1e49290e07efb00"
        )
        self.assertEqual(
            m._hashComponents(["a", "b", "c"]),
            "a58bef80a0b3c42054baeea0edc2308989c7562c",
        )
