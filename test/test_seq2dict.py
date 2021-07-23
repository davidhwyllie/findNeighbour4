import unittest
from findn.seq2json import SeqDictConverter


class Test_SDC_1(unittest.TestCase):
    """tests storage of json data for testing"""

    def runTest(self):

        # test whether a reference compressed data structure can be recovered from json
        input = {
            "A": set([1, 2, 3]),
            "C": set([6]),
            "T": set([4]),
            "G": set([5]),
            "M": {11: "Y", 12: "k"},
            "invalid": 1,
        }
        self.assertIsInstance(input, dict)
        m = SeqDictConverter()
        j1 = m.to_json(input)
        self.assertIsInstance(j1, str)
        j2 = m.from_json(j1)
        self.assertIsInstance(j2, dict)
        self.assertDictEqual(input, j2)
