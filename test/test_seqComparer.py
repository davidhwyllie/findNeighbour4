""" tests seqComparer.py

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
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


from findn.seqComparer import seqComparer



class test_seqComparer_51(unittest.TestCase):
    """ tests mcompare """
    def runTest(self):
        # generate compressed sequences
        refSeq='GGGGGG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        n=0
        originals = [ 'AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN' ]
        guids = []
        for original in originals:
            n+=1
            c = sc.compress(original)
            guid = "{0}-{1}".format(original,n )
            guids.append(guid)
            sc.persist(c, guid=guid)

        res = sc.mcompare(guids[0])      # defaults to sample size 30
        self.assertEqual(len(res), len(originals)-1)
        print("completed")
	
  
class test_seqComparer_49(unittest.TestCase):
    """ tests reporting on stored contents """
    def runTest(self):
        # generate compressed sequences
        refSeq='GGGGGG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        # need > 30 sequences
        originals = ['AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN','AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN',
                     'AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN','AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN',
                     'AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN','AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN']
        guid_names = []
        n=0
        for original in originals:
            n+=1
            c = sc.compress(original)
            this_guid = "{0}-{1}".format(original,n )
            sc.persist(c, guid=this_guid)
            guid_names.append(this_guid)

        res = sc.summarise_stored_items()
        self.assertTrue(isinstance(res, dict))
        self.assertEqual(set(res.keys()), set(['server|scstat|nSeqs']))
class test_seqComparer_48(unittest.TestCase):
    """ tests computations of p values from exact bionomial test """
    def runTest(self):
        # generate compressed sequences
        refSeq='GGGGGG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        # need > 30 sequences
        originals = ['AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN','AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN',
                     'AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN','AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN',
                     'AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN','AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN']
        guid_names = []
        n=0
        for original in originals:
            n+=1
            c = sc.compress(original)
            this_guid = "{0}-{1}".format(original,n )
            sc.persist(c, guid=this_guid)
            guid_names.append(this_guid)

class test_seqComparer_47c(unittest.TestCase):
    """ tests generation of a multisequence alignment with
        testing for the proportion of Ns.
        Tests situation with externally supplied _p1"""
    def runTest(self):
        # generate compressed sequences
        refSeq='GGGGGG'
        sc=seqComparer( maxNs = 8,
                       reference=refSeq,
                       snpCeiling =10)
        # need > 30 sequences
        originals = ['AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN','AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN',
                     'AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN','AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN',
                     'AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN','AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN']

        guid_names = []
        n=0
        for original in originals:
            n+=1
            c = sc.compress(original)
            this_guid = "{0}-{1}".format(original,n )
            sc.persist(c, guid=this_guid)
            guid_names.append(this_guid)

        # but with expected_N supplied;
        df= sc.multi_sequence_alignment(guid_names[0:8], output='df', expected_p1=0.995)      
        # there's variation at positions 0,1,2,3
        self.assertTrue(isinstance(df, pd.DataFrame))
        expected_cols = set(['what_tested','aligned_seq','aligned_seq_len','aligned_seq_len','allN','alignN','allM','alignM','allN_or_M','alignN_or_M','p_value1','p_value2','p_value3', 'p_value4', 'observed_proportion','expected_proportion1','expected_proportion2','expected_proportion3','expected_proportion4'])
        self.assertEqual(set(df.columns.values),expected_cols)

        self.assertEqual(len(df.index),8)
        res= sc.multi_sequence_alignment(guid_names[0:8], output='df_dict', expected_p1=0.995)
        #print(res)
        df = pd.DataFrame.from_dict(res,orient='index')

        self.assertEqual(set(df.columns.values),expected_cols)
    
        self.assertEqual(set(df.index.tolist()), set(['AAACGN-1','NNNCGN-5','CCCCGN-2','TTTCGN-3','GGGGGN-4','ACTCGN-6', 'TCTNGN-7','AAACGN-8']))
        self.assertTrue(df.loc['AAACGN-1','expected_proportion1'] is not None)        # check it computed a value
        self.assertEqual(df.loc['AAACGN-1','expected_proportion1'], 0.995)        # check is used the value passed

class test_seqComparer_47b2(unittest.TestCase):
    """ tests generation of a multisequence alignment with
        testing for the proportion of Ms.
        Tests all three outputs."""
    def runTest(self):
        # generate compressed sequences
        refSeq='GGGGGG'
        sc=seqComparer( maxNs = 6,
                       reference=refSeq,
                       snpCeiling =10)
        # need > 30 sequences
        originals = ['AAACGY','CCCCGY','TTTCGY','GGGGGY','NNNCGY','ACTCGY', 'TCTQGY','AAACGY','CCCCGY','TTTCGY','GGGGGY','NNNCGY','ACTCGY', 'TCTNGY',
                     'AAACGY','CCCCGY','TTTCGY','GGGGGY','NNNCGY','ACTCGY', 'TCTNGY','AAACGY','CCCCGY','TTTCGY','GGGGGY','NNNCGY','ACTCGY', 'TCTNGY',
                     'AAACGY','CCCCGY','TTTCGY','GGGGGY','NNNCGY','ACTCGY', 'TCTNGY','AAACGY','CCCCGY','TTTCGY','GGGGGY','NNNCGY','ACTCGY', 'TCTNGY']

        guid_names = []
        n=0
        for original in originals:
            n+=1
            c = sc.compress(original)
            this_guid = "{0}-{1}".format(original,n )
            sc.persist(c, guid=this_guid)
            guid_names.append(this_guid)

        res= sc.multi_sequence_alignment(guid_names[0:8], output='dict', uncertain_base_type='M')

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
        self.assertEqual(res['variant_positions'],[0,1,2,3])
        df= sc.multi_sequence_alignment(guid_names[0:8], output='df', uncertain_base_type='M')
        
        # there's variation at positions 0,1,2,3
        self.assertTrue(isinstance(df, pd.DataFrame))
        expected_cols = set(['what_tested','aligned_seq','aligned_seq_len','aligned_seq_len','allN','alignN','allM','alignM','allN_or_M','alignN_or_M','p_value1','p_value2','p_value3', 'p_value4', 'observed_proportion','expected_proportion1','expected_proportion2','expected_proportion3','expected_proportion4'])
 

        self.assertEqual(set(df.columns.values),expected_cols)

        self.assertEqual(len(df.index),8)
        res= sc.multi_sequence_alignment(guid_names[0:8], output='df_dict')
        df = pd.DataFrame.from_dict(res,orient='index')

        self.assertTrue(df.loc['AAACGY-1','expected_proportion1'] is not None)        # check it computed a value
        self.assertEqual(set(df.index.tolist()), set(['AAACGY-1','NNNCGY-5','CCCCGY-2','TTTCGY-3','GGGGGY-4','ACTCGY-6', 'TCTQGY-7','AAACGY-8']))

class test_seqComparer_47b(unittest.TestCase):
    """ tests generation of a multisequence alignment with
        testing for the proportion of Ns.
        Tests all three outputs."""
    def runTest(self):
        # generate compressed sequences
        refSeq='GGGGGG'
        sc=seqComparer( maxNs = 6,
                       reference=refSeq,
                       snpCeiling =10)
        # need > 30 sequences
        originals = ['AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN','AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN',
                     'AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN','AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN',
                     'AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN','AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN']

        guid_names = []
        n=0
        for original in originals:
            n+=1
            c = sc.compress(original)
            this_guid = "{0}-{1}".format(original,n )
            sc.persist(c, guid=this_guid)
            guid_names.append(this_guid)

        res= sc.multi_sequence_alignment(guid_names[0:8], output='dict')

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
        self.assertEqual(res['variant_positions'],[0,1,2,3])
        df= sc.multi_sequence_alignment(guid_names[0:8], output='df')
        
        # there's variation at positions 0,1,2,3
        self.assertTrue(isinstance(df, pd.DataFrame))
        expected_cols = set(['what_tested','aligned_seq','aligned_seq_len','aligned_seq_len','allN','alignN','allM','alignM','allN_or_M','alignN_or_M','p_value1','p_value2','p_value3', 'p_value4', 'observed_proportion','expected_proportion1','expected_proportion2','expected_proportion3','expected_proportion4'])
 
        self.assertEqual(set(df.columns.values),expected_cols)
        self.assertEqual(len(df.index),8)
        res= sc.multi_sequence_alignment(guid_names[0:8], output='df_dict')
        df = pd.DataFrame.from_dict(res,orient='index')
        self.assertTrue(df.loc['AAACGN-1','expected_proportion1'] is not None)        # check it computed a value
        self.assertEqual(set(df.index.tolist()), set(['AAACGN-1','CCCCGN-2','TTTCGN-3','GGGGGN-4','NNNCGN-5','ACTCGN-6', 'TCTNGN-7','AAACGN-8']))

class test_seqComparer_47a(unittest.TestCase):
    """ tests generation of a multisequence alignment with
        testing for the proportion of Ns.
        Tests all three outputs."""
    def runTest(self):
        # generate compressed sequences
        refSeq='GGGGGG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        # need > 30 sequences
        originals = ['AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN','AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN',
                     'AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN','AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN',
                     'AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN','AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN']
        guid_names = []
        n=0
        for original in originals:
            n+=1
            c = sc.compress(original)
            this_guid = "{0}-{1}".format(original,n )
            sc.persist(c, guid=this_guid)
            guid_names.append(this_guid)

        res= sc.multi_sequence_alignment(guid_names[0:8], output='dict')
        # there's variation at positions 0,1,2,3
        self.assertEqual(res['variant_positions'],[0,1,2,3])

        df= sc.multi_sequence_alignment(guid_names[0:8], output='df')
        # there's variation at positions 0,1,2,3
        self.assertTrue(isinstance(df, pd.DataFrame))
        expected_cols = set(['what_tested','aligned_seq','aligned_seq_len','aligned_seq_len','allN','alignN','allM','alignM','allN_or_M','alignN_or_M','p_value1','p_value2','p_value3', 'p_value4', 'observed_proportion','expected_proportion1','expected_proportion2','expected_proportion3','expected_proportion4'])
 
        self.assertEqual(set(df.columns.values),expected_cols)
        self.assertEqual(len(df.index),8)
        res= sc.multi_sequence_alignment(guid_names[0:8], output='df_dict')
        df = pd.DataFrame.from_dict(res,orient='index')
    
        self.assertEqual(set(df.index.tolist()), set(guid_names[0:8]))
 
class test_seqComparer_46a(unittest.TestCase):
    """ tests estimate_expected_unk, a function estimating the number of Ns in sequences
        by sampling """
    def runTest(self):
        # generate compressed sequences
        refSeq='GGGGGG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        n=0
        originals = [ 'AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN' ]
        guids = []
        for original in originals:
            n+=1
            c = sc.compress(original)
            guid = "{0}-{1}".format(original,n )
            guids.append(guid)
            sc.persist(c, guid=guid)
          
        res = sc.estimate_expected_unk()      # defaults to sample size 30
        self.assertEqual(res, None)
        
        # analyse the last two
        res = sc.estimate_expected_unk(sample_size=2, exclude_guids = guids[0:5])      
        self.assertEqual(res, 1.5)

        # analyse the first two
        res = sc.estimate_expected_unk(sample_size=2, exclude_guids = guids[2:7])      
        self.assertEqual(res, 1)
class test_seqComparer_46b(unittest.TestCase):
    """ tests estimate_expected_unk, a function estimating the number of Ns in sequences
        by sampling """
    def runTest(self):
        # generate compressed sequences
        refSeq='GGGGGG'
        sc=seqComparer( maxNs = 3,
                       reference=refSeq,
                       snpCeiling =10)
        n=0
        originals = [ 'AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTGGN' ]
        guids = []
        for original in originals:
            n+=1
            c = sc.compress(original)
            guid = "{0}-{1}".format(original,n )
            guids.append(guid)
            sc.persist(c, guid=guid)
          
        res = sc.estimate_expected_unk()      # defaults to sample size 30
        self.assertEqual(res, None)
        
        # analyse them all
        res = sc.estimate_expected_unk(sample_size=7, exclude_guids = [])      
        self.assertEqual(res, None)

        # analyse them all
        res = sc.estimate_expected_unk(sample_size=6, exclude_guids = [])      
        self.assertEqual(res, 1)
class test_seqComparer_46c(unittest.TestCase):
    """ tests estimate_expected_unk_sites, a function estimating the number of Ns in sequences
        by sampling """
    def runTest(self):
        # generate compressed sequences
        refSeq='GGGGGG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        n=0
        originals = [ 'AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN' ]
        guids = []
        for original in originals:
            n+=1
            c = sc.compress(original)
            guid = "{0}-{1}".format(original,n )
            guids.append(guid)
            sc.persist(c, guid=guid)
                 
        # analyse nothing
        res = sc.estimate_expected_unk_sites(sample_size=2, sites=set([]), exclude_guids = guids[0:5])      
        self.assertEqual(res, 0)

        # analyse the last two
        res = sc.estimate_expected_unk_sites(sample_size=2, sites=set([0,1,2,3,4,5]), exclude_guids = guids[0:5])      
        self.assertEqual(res, 1.5)

        # analyse the first two
        res = sc.estimate_expected_unk_sites(sample_size=2, sites=set([0,1,2,3,4,5]), exclude_guids = guids[2:7])      
        self.assertEqual(res, 1)       
class test_seqComparer_45a(unittest.TestCase):
    """ tests the generation of multiple alignments of variant sites."""
    def runTest(self):
        
        # generate compressed sequences
        refSeq='GGGGGG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        
        originals = [ 'AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTNGN' ]
        guid_names = []
        n=0
        for original in originals:
            n+=1
            c = sc.compress(original)
            this_guid = "{0}-{1}".format(original,n )
            sc.persist(c, guid=this_guid)
            guid_names.append(this_guid)

        res= sc.multi_sequence_alignment(guid_names)
        # there's variation at positions 0,1,2,3
        df = pd.DataFrame.from_dict(res['guid2sequence'], orient='index')
        df.columns=res['variant_positions']
        self.assertEqual(len(df.index), 7)
        self.assertEqual(res['variant_positions'],[0,1,2,3])

class test_seqComparer_45b(unittest.TestCase):
    """ tests the generation of multiple alignments of variant sites."""
    def runTest(self):
        
        # generate compressed sequences
        refSeq='GGGGGG'
        sc=seqComparer( maxNs = 6,
                       reference=refSeq,
                       snpCeiling =10)
        
        originals = [ 'AAACGN','CCCCGN','TTTCGN','GGGGGN','NNNCGN','ACTCGN', 'TCTGGN' ]
        guid_names = []
        n=0
        for original in originals:
            n+=1
            c = sc.compress(original)
            this_guid = "{0}-{1}".format(original,n )
            sc.persist(c, guid=this_guid)
            guid_names.append(this_guid)

        res= sc.multi_sequence_alignment(guid_names)
        # there's variation at positions 0,1,2,3
        df = pd.DataFrame.from_dict(res['guid2sequence'], orient='index')
        df.columns=res['variant_positions']
        self.assertEqual(len(df.index), 7)
        self.assertEqual(res['variant_positions'],[0,1,2,3])

class test_seqComparer_1(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        self.assertEqual(sc.reference,refSeq)     
class test_seqComparer_2(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        with self.assertRaises(TypeError):
            retVal=sc.compress(sequence='AC')
class test_seqComparer_3(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq )
        retVal=sc.compress(sequence='ACTG')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([]), 'U': set([]), 'M':{}, 'invalid':0})
class test_seqComparer_3b(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq )
        retVal=sc.compress(sequence='ACTQ')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([]), 'U':set([3]), 'M':{3:'Q'}, 'invalid':0})
class test_seqComparer_3c(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq )
        retVal=sc.compress(sequence='NYTQ')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([0]), 'U': set([0,1,3]), 'M':{1:'Y',3:'Q'}, 'invalid':0})
class test_seqComparer_4(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)

        retVal=sc.compress(sequence='ACTN')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([3]),  'M':{}, 'U':set([3]), 'invalid':0})
class test_seqComparer_5(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        retVal=sc.compress(sequence='ACT-')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([3]), 'M':{}, 'U':set([3]), 'invalid':0})         
class test_seqComparer_6(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'

        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)

        retVal=sc.compress(sequence='TCT-')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([0]), 'N': set([3]), 'M':{}, 'U':set([3]), 'invalid':0})
class test_seqComparer_7(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'

        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        retVal=sc.compress(sequence='ATT-')
        self.assertEqual(retVal,{ 'G': set([]), 'A': set([]), 'C': set([]), 'T': set([1]), 'N': set([3]), 'M':{}, 'U':set([3]), 'invalid':0})

class test_seqComparer_6b(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'

        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        originals = [ 'AAAA','CCCC','TTTT','GGGG','NNNN','ACTG','ACTC', 'TCTN','NYTQ','QRST']
        for original in originals:

            compressed_sequence=sc.compress(sequence=original)
          
            roundtrip = sc.uncompress(compressed_sequence)
            self.assertEqual(original, roundtrip)

class test_seqComparer_6c(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'

        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        originals = [ 'NNNN']
        for original in originals:

            compressed_sequence=sc.compress(sequence=original)
            roundtrip = sc.uncompress(compressed_sequence)
            self.assertEqual(original, roundtrip)

class test_seqComparer_6d(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'

        sc=seqComparer( maxNs = 3, snpCeiling = 20,reference=refSeq)
        originals = [ 'NNNN']
        for original in originals:

            compressed_sequence=sc.compress(sequence=original)
            with self.assertRaises(ValueError):
                roundtrip = sc.uncompress(compressed_sequence)
           
class test_seqComparer_16(unittest.TestCase):
    """ tests the comparison of two sequences where both differ from the reference. """
    def runTest(self):   
        # generate compressed sequences
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        
        seq1 = sc.compress('AAAA')
        seq2 = sc.compress('CCCC')
        self.assertEqual(sc.countDifferences(seq1, seq2),4)
class test_seqComparer_16b(unittest.TestCase):
    """ tests the comparison of two sequences where both differ from the reference. """
    def runTest(self):   
        # generate compressed sequences
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        
        seq1 = sc.compress('AAAA')
        seq2 = sc.compress('RRCC')
        self.assertEqual(sc.countDifferences(seq1,seq2),2)
class test_seqComparer_16c(unittest.TestCase):
    """ tests the comparison of two sequences where both differ from the reference. """
    def runTest(self):   
        # generate compressed sequences
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        
        seq1 = sc.compress('AAAA')
        seq2 = sc.compress('RRNN')
        self.assertEqual(sc.countDifferences(seq1,seq2),0)
class test_seqComparer_17(unittest.TestCase):
    """ tests the comparison of two sequences where one is invalid """
    def runTest(self):   
        # generate compressed sequences
        refSeq='ACTG'
        sc=seqComparer( maxNs = 3,
                       reference=refSeq,
                       snpCeiling =10)
        
        seq1 = sc.compress('AAAA')
        seq2 = sc.compress('NNNN')
        self.assertEqual(sc.countDifferences(seq1,seq2),None)
class test_seqComparer_saveload3(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        compressedObj =sc.compress(sequence='ACTT')
        sc.persist(compressedObj, 'one' )     
        retVal=sc.load(guid='one' )
        self.assertEqual(compressedObj,retVal)        
class test_seqComparer_save_remove(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        compressedObj =sc.compress(sequence='ACTT')
        sc.persist(compressedObj, 'one' )     
        retVal=sc.iscachedinram(guid='one' )
        self.assertEqual(True,retVal)        
        sc.remove('one')
        retVal=sc.iscachedinram(guid='one' )
        self.assertEqual(False,retVal)  
class test_seqComparer_24(unittest.TestCase):
    """ tests N compression """
    def runTest(self):
        
        refSeq=                     'ACTGTTAATTTTTTTTTGGGGGGGGGGGGAA'
        sc=seqComparer(maxNs = 1e8, snpCeiling = 20,reference=refSeq)

        retVal=sc.compress(sequence='ACTGTTAANNNNNNNNTGGGGGGGGGGGGAA')
        self.assertEqual(retVal,{ 'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'M':{}, 'N': set([8,9,10,11,12,13,14,15]),'U': set([8,9,10,11,12,13,14,15]), 'invalid':0})
        retVal=sc.compress(sequence='NNTGTTAANNNNNNNNTGGGGGGGGGGGGAA')
        self.assertEqual(retVal,{ 'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'M':{}, 'N': set([0,1,8,9,10,11,12,13,14,15]),  'U': set([0,1,8,9,10,11,12,13,14,15]),'invalid':0})
       
class test_seqComparer_29(unittest.TestCase):
    """ tests _setStats """
    def runTest(self):
        
        refSeq=                             'ACTGTTAATTTTTTTTTGGGGGGGGGGGGAA'
        sc=seqComparer(maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        compressedObj1=sc.compress(sequence='GGGGTTAANNNNNNNNNGGGGGAAAAGGGAA')
        compressedObj2=sc.compress(sequence='ACTGTTAATTTTTTTTTNNNNNNNNNNNNNN')
        (n1,n2,nall,rv1,rv2,retVal) =sc._setStats(compressedObj1['N'],compressedObj2['N'])
        self.assertEqual(retVal, set([8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]))
 
        compressedObj1=sc.compress(sequence='GGGGTTAANNNNNNNNTGGGGGAAAAGGGAA')
        compressedObj2=sc.compress(sequence='ACTGTTAATTTTTTTTTNNNNNNNNNNNNNN')
        (n1,n2,nall,rv1,rv2,retVal)=sc._setStats(compressedObj1['N'],compressedObj2['N'])
        self.assertEqual(retVal,set([8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30]))
        
        compressedObj1=sc.compress(sequence='NNNGTTAANNNNNNNNTGGGGGAAAAGGGAA')
        compressedObj2=sc.compress(sequence='ACTGTTAATTTTTTTTTNNNNNNNNNNNNNN')
        (n1,n2,nall,rv1,rv2,retVal)=sc._setStats(compressedObj1['N'],compressedObj2['N'])
        self.assertEqual(retVal,set([0,1,2,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30]))


        compressedObj1=sc.compress(sequence='NNNGTTAANNNNNNNNTGGGGGAAAAGGGAA')
        compressedObj2=sc.compress(sequence='ACTNNNNNTTTTTTTTTNNNNNNNNNNNNNN')
        (n1,n2,nall,rv1,rv2,retVal)=sc._setStats(compressedObj1['N'],compressedObj2['N'])
        self.assertEqual(retVal,set([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30]))
 
        compressedObj1=sc.compress(sequence='NNNGTTAANNNNNNNNTGGGGGAAAAGGGAA')
        compressedObj2=sc.compress(sequence='ACTNNNNNTTTTTTTTTQQQQQQQQQQQQQQ')
        (n1,n2,nall,rv1,rv2,retVal)=sc._setStats(compressedObj1['N'],compressedObj2['N'])
        self.assertEqual(retVal,set([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]))
        (n1,n2,nall,rv1,rv2,retVal)=sc._setStats(compressedObj1['M'],compressedObj2['M'])
        self.assertEqual(retVal,set([17,18,19,20,21,22,23,24,25,26,27,28,29,30]))
                
        compressedObj1=sc.compress(sequence='qqqGTTAAqqqqqqqqTGGGGGAAAAGGGAA')
        compressedObj2=sc.compress(sequence='ACTqqqqqTTTTTTTTTqqqqqqqqqqqqqq')
        (n1,n2,nall,rv1,rv2,retVal)=sc._setStats(compressedObj1['M'],compressedObj2['M'])
        self.assertEqual(retVal,set([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30]))
 
class test_seqComparer_37(unittest.TestCase):
    """ tests the loading of an exclusion file """
    def runTest(self):
        
        # default exclusion file
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, reference=refSeq, snpCeiling =1)
        self.assertEqual( sc.excluded_hash(), 'Excl 0 nt [d751713988987e9331980363e24189ce]')

class test_seqComparer_38(unittest.TestCase):
    """ tests the loading of an exclusion file """
    def runTest(self):
        
        # no exclusion file
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, reference=refSeq, snpCeiling =1)
        self.assertEqual( sc.excluded_hash(), 'Excl 0 nt [d751713988987e9331980363e24189ce]')
 
class test_seqComparer_40(unittest.TestCase):
    """ tests the computation of a hash of a compressed object """
    def runTest(self):

        # generate compressed sequences
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, reference=refSeq, snpCeiling =10)
        compressed_sequence = sc.compress(sequence='TTAA')

        res = sc.compressed_sequence_hash(compressed_sequence)
        self.assertEqual(res, "da8785691df5858b0b847db59bdefd11")
        
        

class test_seqComparer_45(unittest.TestCase):
    """ tests insertion of large sequences """
    def runTest(self):
        inputfile = "reference/NC_000962.fasta"
        with open(inputfile, 'rt') as f:
            for record in SeqIO.parse(f,'fasta'):
                    goodseq = str(record.seq)
                    badseq = ''.join('N'*len(goodseq))
                    originalseq = list(str(record.seq))
        sc=seqComparer( maxNs = 1e8,
                           reference=record.seq,
                           snpCeiling =100)
        n_pre =  0          
        guids_inserted = list()			
        for i in range(1,4):        #40
            
            seq = originalseq
            if i % 5 ==0:
                is_mixed = True
                guid_to_insert = "mixed_{0}".format(n_pre+i)
            else:
                is_mixed = False
                guid_to_insert = "nomix_{0}".format(n_pre+i)	
            # make i mutations at position 500,000
            
            offset = 500000
            nVariants = 0
            for j in range(i):
                mutbase = offset+j
                ref = seq[mutbase]
                if is_mixed == False:
                    nVariants +=1
                    if not ref == 'T':
                        seq[mutbase] = 'T'
                    if not ref == 'A':
                        seq[mutbase] = 'A'
                if is_mixed == True:
                        seq[mutbase] = 'N'					
            seq = ''.join(seq)
            
            if i % 11 == 0:
                seq = badseq        # invalid
                
            guids_inserted.append(guid_to_insert)			
            if not is_mixed:
                    #print("Adding TB sequence {2} of {0} bytes with {1} Ns and {3} variants relative to ref.".format(len(seq), seq.count('N'), guid_to_insert, nVariants))
                    pass
            else:
                    #print("Adding mixed TB sequence {2} of {0} bytes with {1} Ns relative to ref.".format(len(seq), seq.count('N'), guid_to_insert))
                    pass    
            self.assertEqual(len(seq), 4411532)		# check it's the right sequence
    
            c = sc.compress(seq)
            sc.persist(c, guid=guid_to_insert )
 
class test_seqComparer_47(unittest.TestCase):
    """ tests raise_error"""
    def runTest(self):
                # generate compressed sequences
        refSeq='GGGGGGGGGGGG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        with self.assertRaises(ZeroDivisionError):
            sc.raise_error("token")

class test_seqComparer_48(unittest.TestCase):
    """ tests distmat, a function yielding a distance matrix."""
    def runTest(self):
        
        # generate compressed sequences
        refSeq='GGGGGGGGGGGG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        
        originals = [ 'AAACACTGACTG','CCCCACTGACTG','TTTCACTGACTG' ]
        for original in originals:   
            c = sc.compress(original)
            sc.persist(c, guid=original )
        
        n = 0
        for item in sc.distmat(half=False, diagonal=True):
            n +=1
        l = len(originals)
        
        self.assertEqual(n, l*l)

        n = 0
        for item in sc.distmat(half=False, diagonal=False):
            n +=1

        l = len(originals)
        
        self.assertEqual(n, (l*l)-l)


        n = 0
        for item in sc.distmat(half=True, diagonal=False):
            n +=1

        l = len(originals)
        
        self.assertEqual(n, (l*(l-1)/2))
        
class test_seqComparer_50(unittest.TestCase):
    """ tests estimate_expected_proportion, a function computing the proportion of Ns expected based on the median
    Ns in a list of sequences"""
    def runTest(self):
        refSeq='GGGGGGGGGGGG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
  
        res = sc.estimate_expected_proportion([])
        self.assertTrue(res is None)

        res = sc.estimate_expected_proportion(['AA','AA'])
        self.assertTrue(res is None)
        
        res = sc.estimate_expected_proportion(['AA','AA', 'AA'])
        self.assertTrue(res is not None)
        self.assertTrue(res == 0)
        
        res = sc.estimate_expected_proportion(['AAN','AAN', 'AAN'])
        self.assertTrue(res is not None)
        self.assertAlmostEqual(res, 1/3)
        
