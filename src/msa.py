#!/usr/bin/env python3
""" manages multisequence alignments """

import unittest
import pandas as pd
import datetime
import time
from mongoStore import fn3persistence
from identify_sequence_set import IdentifySequenceSet

class MSAStore():
    """ stores multisequence alignments, both in ram and in a database """
    def __init__(self, PERSIST, in_ram_persistence_time=3600):
        """ sets up an MSA cache, which stores MSAs on disc and in ram 
            PERSIST is an fn3persistence object
            in_ram_persistence_time the time MSAs are cached in ram, in seconds"""
        self.cache = {}
        self.PERSIST = PERSIST
        self.iss = IdentifySequenceSet()
        self.in_ram_persistence_time = in_ram_persistence_time

    def get_token(self, what, has_outgroup, guids):
        """ computes a token representing a set of guids with or without an outgroup """
        return self.iss.make_identifier('msa', what, has_outgroup is not None, guids)

    def existing_tokens(self):
        """ returns a set of all tokens on disc """
        return set(self.PERSIST.msa_stored_ids())

    def is_in_ram(self,token):
        """ whether the msa identifed by token is cached in ram """
        return token in self.cache.keys()

        return token in self.cache.keys()
    def persist(self, token, msa_result):
        """ stores a MSAResult object in mongo identified by token """
        self.PERSIST.msa_store(token, msa_result.serialise())
        self.cache_in_ram(token, msa_result)

    def load(self, token):
        """ recovers an MSAResult identified by token.  Caches it in ram"""
        persisted_dict = self.PERSIST.msa_read(token)
        if persisted_dict is not None:
            msa_result=MSAResult(**persisted_dict)
            self.cache_in_ram(token, msa_result)
            return msa_result
        else:
            return None

    def cache_in_ram(self, token, msa_result):
        """ stores a MSAResult object in ram identified by token """
        self.purge_ram()
        now = datetime.datetime.now()
        expiry_time = now + datetime.timedelta(seconds = self.in_ram_persistence_time)
        if not token in self.cache.keys():
            self.cache[token]={'expiry_time':expiry_time, 'MSA':msa_result  }
        else:         
            self.cache[token]['expiry_time']  = expiry_time     # update expiry time

    def purge_ram(self):
        """ removes any in-ram objects which are time expired """
        now = datetime.datetime.now()
        to_delete = set()
        for token in self.cache.keys():
            if now > self.cache[token]['expiry_time']:
                to_delete.add(token)
        for token in to_delete:
                del(self.cache[token])

    def unpersist(self, whitelist):
        """ removes any on-disc msa results the ids of which are not in the whitelist """
        self.PERSIST.msa_delete_unless_whitelisted(whitelist)

class MSAResult():
    """ a representation of a multisequence alignment """
    def __init__(self, 
                    variant_positions,
                    invalid_guids,
                    valid_guids,
                    expected_p1,
                    sample_size,
                    df_dict,
                    what_tested,
                    outgroup,
                    creation_time
                ):

        """ a representation of a multisequence alignment  
        
        Parameters:
                    'variant_positions': positions of variation in msa
                    'invalid_guids': invalid guids associated with, but not included in, msa
                    'valid_guids': guids in msa
                    'expected_p1': expected_p1 (see ._msa in hybridComparer for discussion)
                    'sample_size': sample size used to estimate p1
                    'df_dict': a dictionary serialising the pandas dataframe containing the msa
                    'what_tested': what (M/N/M_or_N) was tested
                    'outgroup': either none, or a guid used as the outgroup.
                    'creation_time': timestamp the msa was created
        """

        self.variant_positions = variant_positions
        self.alignment_length=len(variant_positions)
        self.invalid_guids = invalid_guids
        self.valid_guids = valid_guids
        self.expected_p1 = expected_p1
        self.sample_size = sample_size
        self.df = pd.DataFrame.from_dict(df_dict, orient='index')

        self.what_tested = what_tested
        self.outgroup = outgroup
        self.creation_time= creation_time

        iss = IdentifySequenceSet()
        self.token = iss.make_identifier('msa', self.what_tested, outgroup is not None, valid_guids)


    def msa_df(self):
        """ return the msa as a dataframe """
        return self.df
    def msa_html(self):
        """ returns the msa as an html table """
        return self.df.to_html()
    def msa_fasta(self):    
        """ returns the msa as fasta """
        fasta= ""
        for guid in self.df.index:
            fasta=fasta + ">{0}\n{1}\n".format(guid, self.df.loc[guid,'aligned_seq'])
        return fasta	    
    def msa_dict(self):
        """ return the msa as a dictionary """
        return(self.df.to_dict(orient='index'))
    def serialise(self):
        """ return the entire object in a serialisable format (dictionary) """
        return {
                    'variant_positions':self.variant_positions,
                    'invalid_guids': self.invalid_guids,
                    'valid_guids':self.valid_guids,
                    'expected_p1':self.expected_p1,
                    'sample_size':self.sample_size,
                    'df_dict':self.df.to_dict(orient='index'),
                    'what_tested':self.what_tested,
                    'outgroup':self.outgroup,
                    'creation_time':self.creation_time  
                    }
    
class Test_MSA(unittest.TestCase):
    """ tests the MSAResult class"""
    def runTest(self):
        inputdict = {'variant_positions': [0, 1, 2, 3], 'invalid_guids': [], 'valid_guids': ['AAACGY-1', 'CCCCGY-2', 'TTTCGY-3', 'GGGGGY-4', 'NNNCGY-5', 'ACTCGY-6', 'TCTQGY-7', 'AAACGY-8'], 'expected_p1': 0.16666666666666666, 'sample_size': 30, 'df_dict': {'AAACGY-1': {'aligned_seq': 'AAAC', 'aligned_seq_len': 4, 'allN': 0, 'alignN': 0, 'allM': 1, 'alignM': 0, 'allN_or_M': 1, 'alignN_or_M': 0, 'p_value1': 1.0, 'p_value2': 1.0, 'p_value3': 1.0, 'p_value4': 1.0, 'observed_proportion': 0.0, 'expected_proportion1': 0.16666666666666666, 'expected_proportion2': 0.5, 'expected_proportion3': 0.5, 'expected_proportion4': 0.0, 'what_tested': 'M'}, 'CCCCGY-2': {'aligned_seq': 'CCCC', 'aligned_seq_len': 4, 'allN': 0, 'alignN': 0, 'allM': 1, 'alignM': 0, 'allN_or_M': 1, 'alignN_or_M': 0, 'p_value1': 1.0, 'p_value2': 1.0, 'p_value3': 1.0, 'p_value4': 1.0, 'observed_proportion': 0.0, 'expected_proportion1': 0.16666666666666666, 'expected_proportion2': 0.5, 'expected_proportion3': 0.5, 'expected_proportion4': 0.0, 'what_tested': 'M'}, 'TTTCGY-3': {'aligned_seq': 'TTTC', 'aligned_seq_len': 4, 'allN': 0, 'alignN': 0, 'allM': 1, 'alignM': 0, 'allN_or_M': 1, 'alignN_or_M': 0, 'p_value1': 1.0, 'p_value2': 1.0, 'p_value3': 1.0, 'p_value4': 1.0, 'observed_proportion': 0.0, 'expected_proportion1': 0.16666666666666666, 'expected_proportion2': 0.5, 'expected_proportion3': 0.5, 'expected_proportion4': 0.0, 'what_tested': 'M'}, 'GGGGGY-4': {'aligned_seq': 'GGGG', 'aligned_seq_len': 4, 'allN': 0, 'alignN': 0, 'allM': 1, 'alignM': 0, 'allN_or_M': 1, 'alignN_or_M': 0, 'p_value1': 1.0, 'p_value2': 1.0, 'p_value3': 1.0, 'p_value4': 1.0, 'observed_proportion': 0.0, 'expected_proportion1': 0.16666666666666666, 'expected_proportion2': 0.5, 'expected_proportion3': 0.5, 'expected_proportion4': 0.0, 'what_tested': 'M'}, 'NNNCGY-5': {'aligned_seq': 'NNNC', 'aligned_seq_len': 4, 'allN': 3, 'alignN': 3, 'allM': 1, 'alignM': 0, 'allN_or_M': 4, 'alignN_or_M': 3, 'p_value1': 1.0, 'p_value2': 1.0, 'p_value3': 1.0, 'p_value4': 1.0, 'observed_proportion': 0.0, 'expected_proportion1': 0.16666666666666666, 'expected_proportion2': 0.5, 'expected_proportion3': 0.5, 'expected_proportion4': 0.0, 'what_tested': 'M'}, 'ACTCGY-6': {'aligned_seq': 'ACTC', 'aligned_seq_len': 4, 'allN': 0, 'alignN': 0, 'allM': 1, 'alignM': 0, 'allN_or_M': 1, 'alignN_or_M': 0, 'p_value1': 1.0, 'p_value2': 1.0, 'p_value3': 1.0, 'p_value4': 1.0, 'observed_proportion': 0.0, 'expected_proportion1': 0.16666666666666666, 'expected_proportion2': 0.5, 'expected_proportion3': 0.5, 'expected_proportion4': 0.0, 'what_tested': 'M'}, 'TCTQGY-7': {'aligned_seq': 'TCTM', 'aligned_seq_len': 4, 'allN': 0, 'alignN': 0, 'allM': 2, 'alignM': 1, 'allN_or_M': 2, 'alignN_or_M': 1, 'p_value1': 0.5177469135802468, 'p_value2': 0.0, 'p_value3': 0.9375, 'p_value4': 0.0, 'observed_proportion': 0.25, 'expected_proportion1': 0.16666666666666666, 'expected_proportion2': 0.5, 'expected_proportion3': 0.5, 'expected_proportion4': 0.0, 'what_tested': 'M'}, 'AAACGY-8': {'aligned_seq': 'AAAC', 'aligned_seq_len': 4, 'allN': 0, 'alignN': 0, 'allM': 1, 'alignM': 0, 'allN_or_M': 1, 'alignN_or_M': 0, 'p_value1': 1.0, 'p_value2': 1.0, 'p_value3': 1.0, 'p_value4': 1.0, 'observed_proportion': 0.0, 'expected_proportion1': 0.16666666666666666, 'expected_proportion2': 0.5, 'expected_proportion3': 0.5, 'expected_proportion4': 0.0, 'what_tested': 'M'}}, 'what_tested': 'M', 'outgroup': None, 'creation_time': '2019-11-17T23:46:00.098151'}

        m = MSAResult(**inputdict)
        self.assertEqual(type(m.msa_fasta()), str)
        self.assertEqual(type(m.msa_html()), str)
        self.assertEqual(type(m.msa_dict()), dict)
        self.assertEqual(type(m.serialise()), dict)
        self.assertEqual(type(m.df), pd.DataFrame)
        self.assertEqual(m.token, "msa|M|no_og|ddc4781ec984b66b0b5bf006a71b29cf1f523740")

UNITTEST_MONGOCONN = "mongodb://localhost"

class Test_MSAStore(unittest.TestCase):
    """ tests the MSAStore class"""
    def runTest(self):

        inputdict1 = {'variant_positions': [0, 1, 2, 3], 'invalid_guids': [], 'valid_guids': ['AAACGY-1', 'CCCCGY-2', 'TTTCGY-3', 'GGGGGY-4', 'NNNCGY-5', 'ACTCGY-6', 'TCTQGY-7', 'AAACGY-8'], 'expected_p1': 0.16666666666666666, 'sample_size': 30, 'df_dict': {'AAACGY-1': {'aligned_seq': 'AAAC', 'aligned_seq_len': 4, 'allN': 0, 'alignN': 0, 'allM': 1, 'alignM': 0, 'allN_or_M': 1, 'alignN_or_M': 0, 'p_value1': 1.0, 'p_value2': 1.0, 'p_value3': 1.0, 'p_value4': 1.0, 'observed_proportion': 0.0, 'expected_proportion1': 0.16666666666666666, 'expected_proportion2': 0.5, 'expected_proportion3': 0.5, 'expected_proportion4': 0.0, 'what_tested': 'M'}, 'CCCCGY-2': {'aligned_seq': 'CCCC', 'aligned_seq_len': 4, 'allN': 0, 'alignN': 0, 'allM': 1, 'alignM': 0, 'allN_or_M': 1, 'alignN_or_M': 0, 'p_value1': 1.0, 'p_value2': 1.0, 'p_value3': 1.0, 'p_value4': 1.0, 'observed_proportion': 0.0, 'expected_proportion1': 0.16666666666666666, 'expected_proportion2': 0.5, 'expected_proportion3': 0.5, 'expected_proportion4': 0.0, 'what_tested': 'M'}, 'TTTCGY-3': {'aligned_seq': 'TTTC', 'aligned_seq_len': 4, 'allN': 0, 'alignN': 0, 'allM': 1, 'alignM': 0, 'allN_or_M': 1, 'alignN_or_M': 0, 'p_value1': 1.0, 'p_value2': 1.0, 'p_value3': 1.0, 'p_value4': 1.0, 'observed_proportion': 0.0, 'expected_proportion1': 0.16666666666666666, 'expected_proportion2': 0.5, 'expected_proportion3': 0.5, 'expected_proportion4': 0.0, 'what_tested': 'M'}, 'GGGGGY-4': {'aligned_seq': 'GGGG', 'aligned_seq_len': 4, 'allN': 0, 'alignN': 0, 'allM': 1, 'alignM': 0, 'allN_or_M': 1, 'alignN_or_M': 0, 'p_value1': 1.0, 'p_value2': 1.0, 'p_value3': 1.0, 'p_value4': 1.0, 'observed_proportion': 0.0, 'expected_proportion1': 0.16666666666666666, 'expected_proportion2': 0.5, 'expected_proportion3': 0.5, 'expected_proportion4': 0.0, 'what_tested': 'M'}, 'NNNCGY-5': {'aligned_seq': 'NNNC', 'aligned_seq_len': 4, 'allN': 3, 'alignN': 3, 'allM': 1, 'alignM': 0, 'allN_or_M': 4, 'alignN_or_M': 3, 'p_value1': 1.0, 'p_value2': 1.0, 'p_value3': 1.0, 'p_value4': 1.0, 'observed_proportion': 0.0, 'expected_proportion1': 0.16666666666666666, 'expected_proportion2': 0.5, 'expected_proportion3': 0.5, 'expected_proportion4': 0.0, 'what_tested': 'M'}, 'ACTCGY-6': {'aligned_seq': 'ACTC', 'aligned_seq_len': 4, 'allN': 0, 'alignN': 0, 'allM': 1, 'alignM': 0, 'allN_or_M': 1, 'alignN_or_M': 0, 'p_value1': 1.0, 'p_value2': 1.0, 'p_value3': 1.0, 'p_value4': 1.0, 'observed_proportion': 0.0, 'expected_proportion1': 0.16666666666666666, 'expected_proportion2': 0.5, 'expected_proportion3': 0.5, 'expected_proportion4': 0.0, 'what_tested': 'M'}, 'TCTQGY-7': {'aligned_seq': 'TCTM', 'aligned_seq_len': 4, 'allN': 0, 'alignN': 0, 'allM': 2, 'alignM': 1, 'allN_or_M': 2, 'alignN_or_M': 1, 'p_value1': 0.5177469135802468, 'p_value2': 0.0, 'p_value3': 0.9375, 'p_value4': 0.0, 'observed_proportion': 0.25, 'expected_proportion1': 0.16666666666666666, 'expected_proportion2': 0.5, 'expected_proportion3': 0.5, 'expected_proportion4': 0.0, 'what_tested': 'M'}, 'AAACGY-8': {'aligned_seq': 'AAAC', 'aligned_seq_len': 4, 'allN': 0, 'alignN': 0, 'allM': 1, 'alignM': 0, 'allN_or_M': 1, 'alignN_or_M': 0, 'p_value1': 1.0, 'p_value2': 1.0, 'p_value3': 1.0, 'p_value4': 1.0, 'observed_proportion': 0.0, 'expected_proportion1': 0.16666666666666666, 'expected_proportion2': 0.5, 'expected_proportion3': 0.5, 'expected_proportion4': 0.0, 'what_tested': 'M'}}, 'what_tested': 'M', 'outgroup': None, 'creation_time': '2019-11-17T23:46:00.098151'}
        inputdict2 = {'variant_positions': [0, 1, 2, 3], 'invalid_guids': [], 'valid_guids': ['AAACGY-1', 'CCCCGY-2'], 'expected_p1': 0.16666666666666666, 'sample_size': 30, 'df_dict': {'AAACGY-1': {'aligned_seq': 'AAAC', 'aligned_seq_len': 4, 'allN': 0, 'alignN': 0, 'allM': 1, 'alignM': 0, 'allN_or_M': 1, 'alignN_or_M': 0, 'p_value1': 1.0, 'p_value2': 1.0, 'p_value3': 1.0, 'p_value4': 1.0, 'observed_proportion': 0.0, 'expected_proportion1': 0.16666666666666666, 'expected_proportion2': 0.5, 'expected_proportion3': 0.5, 'expected_proportion4': 0.0, 'what_tested': 'M'}, 'CCCCGY-2': {'aligned_seq': 'CCCC', 'aligned_seq_len': 4, 'allN': 0, 'alignN': 0, 'allM': 1, 'alignM': 0, 'allN_or_M': 1, 'alignN_or_M': 0, 'p_value1': 1.0, 'p_value2': 1.0, 'p_value3': 1.0, 'p_value4': 1.0, 'observed_proportion': 0.0, 'expected_proportion1': 0.16666666666666666, 'expected_proportion2': 0.5, 'expected_proportion3': 0.5, 'expected_proportion4': 0.0, 'what_tested': 'M'} }, 'what_tested': 'M', 'outgroup': None, 'creation_time': '2019-11-17T23:46:00.098151'}

        m1 = MSAResult(**inputdict1)
        m2 = MSAResult(**inputdict2)

        guids1 = m1.valid_guids
        guids2 = m2.valid_guids

        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
        ms = MSAStore(p, in_ram_persistence_time=1)

        t1 = ms.get_token('M',None, guids1)
        t2 = ms.get_token('M',None, guids2)
 
        self.assertFalse(ms.is_in_ram(t1))
        self.assertFalse(ms.is_in_ram(t2))

        ms.cache_in_ram(token=t1,msa_result=m1)
        ms.cache_in_ram(token=t2,msa_result= m2)
        self.assertTrue(ms.is_in_ram(t1))
        self.assertTrue(ms.is_in_ram(t2))

        # test purging of expired samples
        time.sleep(2)
        ms.purge_ram()
  
        self.assertFalse(ms.is_in_ram(t1))
        self.assertFalse(ms.is_in_ram(t2))

        # store on disc
        ms.persist(t1,m1)
        ms.persist(t2,m2)

        m1r = ms.load(t1)
        m2r = ms.load(t2)

        self.assertEqual(m1r.valid_guids, m1.valid_guids)
        self.assertEqual(m2r.valid_guids, m2.valid_guids)

        ms.unpersist([])
        
        m1r = ms.load(t1)
        m2r = ms.load(t2)

        self.assertIsNone(m1r)
        self.assertIsNone(m2r)

   
