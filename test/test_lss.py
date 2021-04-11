import unittest
from make_large_sequence_set import LargeSequenceSetGenerator

class test_lss_generator_1(unittest.TestCase):
    """ tests large sequence generator """
    def runTest(self):
       # generate compressed sequences
        refSeq='GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG'
        prefix = 'A'
        n_sequences = 10
        k1 = 10
        s = 2
        max_c = 10
        p = 0.1
        lsg = LargeSequenceSetGenerator(refSeq)
        nGen = 0
        for guid, seq in lsg.make_sequences(prefix, n_sequences, k1, s, max_c, p):
            nGen+=1
            self.assertEqual(len(seq) , len(refSeq))
        self.assertEqual(nGen, n_sequences * max_c)