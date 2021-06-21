""" tests nucleicacid.py

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

from findn.NucleicAcid import NucleicAcid


class Test_NucleicAcid_Base1(unittest.TestCase):
    """base class used for tests of NucleicAcid class"""

    def setUp(self):
        pass

    def tearDown(self):
        pass


class Test_NucleicAcid_1(Test_NucleicAcid_Base1):
    """when initialised, are all outputs appropriately set."""

    def runTest(self):
        na = NucleicAcid()

        self.assertEqual(na.isValid, None)
        self.assertEqual(na.composition["MD5"], None)


class Test_NucleicAcid_2(Test_NucleicAcid_Base1):
    """checks number of nucleic acids correctly."""

    def runTest(self):
        na = NucleicAcid()
        na.examine("ACGT")
        self.assertEqual(na.isValid, True)
        self.assertEqual(na.composition["MD5"], "f1f8f4bf413b16ad135722aa4591043e")


class Test_NucleicAcid_3(Test_NucleicAcid_Base1):
    """checks number of nucleic acids"""

    def runTest(self):
        na = NucleicAcid()
        na.examine("ACGTN-")
        self.assertEqual(na.isValid, True)


class Test_NucleicAcid_4(Test_NucleicAcid_Base1):
    """checks composition of nucleic acids"""

    def runTest(self):
        na = NucleicAcid()
        self.assertRaises(ValueError, na.examine, nucleicAcidString="ACGTN-X")


class Test_NucleicAcid_5(Test_NucleicAcid_Base1):
    """checks that what is passed is a string"""

    def runTest(self):
        na = NucleicAcid()
        self.assertRaises(TypeError, na.examine, nucleicAcidString=35)


class Test_NucleicAcid_6(Test_NucleicAcid_Base1):
    """checks that what is passed is a string"""

    def runTest(self):
        na = NucleicAcid()
        self.assertRaises(ValueError, na.examine, nucleicAcidString="")


class Test_NucleicAcid_7(Test_NucleicAcid_Base1):
    """checks number of nucleic acids"""

    def runTest(self):
        na = NucleicAcid()
        na.examine("AACCCGT")
        self.assertEqual(na.composition["A"], 2)
        self.assertEqual(na.composition["C"], 3)


class Test_NucleicAcid_8(Test_NucleicAcid_Base1):
    """checks number of nucleic acids"""

    def runTest(self):
        na = NucleicAcid()
        na.examine("AAAN")
        self.assertEqual(na.composition["propACTG"], 0.75)


class Test_NucleicAcid_9(Test_NucleicAcid_Base1):
    """checks number of nucleic acids"""

    def runTest(self):
        na = NucleicAcid()
        na.examine("AAAR")
        self.assertEqual(na.composition["mixed"], 1)
        self.assertEqual(na.composition["R"], 1)


class Test_NucleicAcid_10(Test_NucleicAcid_Base1):
    """checks trinucleotide code is acceptable"""

    def runTest(self):
        na = NucleicAcid()
        na.examine("AAAB")
        self.assertEqual(na.isValid, True)
        self.assertEqual(na.composition["mixed"], 1)
        self.assertEqual(na.composition["B"], 1)
