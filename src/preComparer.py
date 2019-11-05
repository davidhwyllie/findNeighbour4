#!/usr/bin/env python3

# python3 code to compare fasta sequences
import unittest
import os
import glob
import sys
import copy
import datetime
import pickle
import hashlib
import math
import multiprocessing
import uuid
import json
import psutil
from gzip import GzipFile
import random
import itertools
import numpy as np
from scipy.stats import binom_test, binom, median_absolute_deviation
import pandas as pd
from collections import Counter

# only used for unit testing
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide

class preComparer():
    """ compares reference compressed sequences, using an approximate method which does not store Ns.
        it does not compare sequences yielding SNP distances: rather, it provides and indication of whether
        an exact SNP distance should be computed.

        This approach is much faster and uses much less RAM than previous alternatives.
        However, it requires that parameters for the function are carefully set and validated to ensure the 
        approximate computation delivers fully sensitive results.

    """
    def __init__(self,
                    
                    selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor = 5,
                    mixed_reporting_cutoff=0,
                    N_mean=39420,
                    N_sd = 9210,
                    highN_z_reporting_cutoff = 2,
                    alpha = 1e-14,
                    probN_inflation_factor = 3,           
                    n_positions_examined = 3854241
                ):

        """ instantiates the sequence comparer, an object which manages in-memory reference compressed sequences.
        
        It does not manage persistence, nor does it automatically load sequences.
        
        Parameters:
        ===========
        selection_cutoff: (int)
            snp distances more than this are not of interest epidemiologically

        over_selection_cutoff_ignore_factor:
            estimated SNP distances, ignoring Ns, more than over_selection_cutoff_ignore_factor * selection_cutoff do not need to be considered.  For example, if a SNP cutoff was 20, and over_selection_cutoff_ignore_factor is 5, we can safely consider isolates 
        mixed_reporting_cutoff:
            regard sample pairs with more Ms than this in the alignment as mixed ## ignored 

        N_mean, n_sd:
            mean and standard deviation of the number of Ns per sample in the population.  Used to compute z score for each sample  

	##### the below are currently ignored -not used
        highN_z_reporting_cutoff:  
            samples with Z-scores for numbers of Ns higher than this are regarded as having high Ns
        alpha:
            the type I error of the estimation performed.  For example, if alpha = 1x10^12, we expect a sample with true pairwise difference less than selection_cutoff to be called as unrelated once every 10^12 comparisons.
        probN_inflation_factor:
            used to estimate the expected number of Ns in any comparison.
            probN is estimated as n_positions_examined/number of Ns per sample.
            the expected number of Ns in a comparison (or rather the mean of the distribution of Ns in a comparison)
            is computed as probN * probN_inflation_factor.  This is an important parameter and has to be estimated from data.  See preComparer_validator.py for details.
	#### end of section which may be removed

	Algorithmically, the compare() method 
	* examines whether either sequence is tagged as invalid.  if so, reports no need to retest; otherwise
    * compares numbers of Ns by Z-score.  If either of these exceeds the Z-score, 
        if penalised_destim > self.over_selection_cutoff_ignore_factor * self.selection_cutoff [essentially: the SNP distance estimate is very high]
            reports no need to retest
        otherwise
	    reports a need to retest

    For optimal performance, these cutoffs need to be tuned.  Code doing so is available in preComparer_calibrator.py.
        Returns:
        ========
        None

        - to run unit tests, do
        python3 -m unittest preComparer
        """
         
        self.set_operating_parameters( selection_cutoff,
                    over_selection_cutoff_ignore_factor,
                    mixed_reporting_cutoff,
                    N_mean,
                    N_sd,
                    highN_z_reporting_cutoff,
                    alpha,
                    probN_inflation_factor,           
                    n_positions_examined)
  
        
        # initialise pairwise sequences for comparison.
        self._refresh()

    def set_operating_parameters(self, selection_cutoff,
                    over_selection_cutoff_ignore_factor,
                    mixed_reporting_cutoff,
                    N_mean,
                    N_sd,
                    highN_z_reporting_cutoff,
                    alpha,
                    probN_inflation_factor,           
                    n_positions_examined):
        """ sets the parameters used by the algorithm """
        self.selection_cutoff = selection_cutoff 
        self.over_selection_cutoff_ignore_factor = over_selection_cutoff_ignore_factor
        self.mixed_reporting_cutoff= mixed_reporting_cutoff
        self.N_mean= N_mean

        self.N_sd = N_sd 
        self.highN_z_reporting_cutoff = highN_z_reporting_cutoff
        self.alpha = alpha
        self.oneminusalpha = 1-alpha
        self.probN_inflation_factor = probN_inflation_factor
        
        self.n_positions_examined = n_positions_examined

    def check_operating_parameters(self, selection_cutoff,
                    over_selection_cutoff_ignore_factor,
                    mixed_reporting_cutoff,
                    N_mean,
                    N_sd,
                    highN_z_reporting_cutoff,
                    alpha,
                    probN_inflation_factor,           
                    n_positions_examined):
        """ returns true if the existing operating parameters are those passed, otherwise returns false """
  
        if  (self.selection_cutoff == selection_cutoff  and 
            self.over_selection_cutoff_ignore_factor == over_selection_cutoff_ignore_factor and 
            self.mixed_reporting_cutoff == mixed_reporting_cutoff and 
            self.N_mean == N_mean and 
            self.N_sd == N_sd and
            self.highN_z_reporting_cutoff == highN_z_reporting_cutoff and
            self.alpha == alpha and
            self.probN_inflation_factor == probN_inflation_factor and     
            self.n_positions_examined == n_positions_examined):
            return True
        else:
            return False
      
    def _refresh(self):
        """ initialise in ram stores """
        self.seqProfile={}      # the positions which differ from the reference
        self.composition = {}   # composition statistics for each sequence
        self.N_positions = dict()	# number of Ns at positions where there are Ns
        self.M_positions = dict()	# number of Ns at positions where there are Ns
        self.binom_results = {} # results of binomial tests: they are expensive to compute, so we cache the results

    def raise_error(self,token):
        """ raises a ZeroDivisionError, with token as the message.
            useful for unit tests of error logging """
        raise ZeroDivisionError(token)
    
    def persist(self, obj, guid):
        """ keeps a reference compressed object into RAM.
            Note: the sequences are stored on disc/db relative to the reference.
            Compression relative to each other is carried out post-hoc in ram

            Parameters:
                obj: 
                a reference compressed object, a dictionary of the form
                {'A':set([1,2,3,4]} which would represent non-reference A at positions 1-4.

                guid:
                a unique identifier for a sequence

            Returns:
                invalid: 0 if valid, 1 if invalid
            """

        # check it does not exist.  Does not allow overwriting
        if guid in self.seqProfile.keys():
            raise KeyError("Duplicate guid supplied; already exists in preComparer: {0}".format(guid))
        # check it is not invalid
        isinvalid = False
        try:
            isinvalid = (obj['invalid']==1)
        except KeyError:
            pass

        # store composition
        # store the positions of Ns, and their number
        try:
            for pos in obj['N']:        # the Ns
                try:
                    n = self.N_positions[pos]
                except KeyError:        # no N at that position
                    n = 0
                self.N_positions[pos] = n+1

        except KeyError:                # no Ns
            pass

        # store the positions of Ms, and their number
        try:
            for pos in obj['M'].keys():        # the Ms
                try:
                    n = self.M_positions[pos]
                except KeyError:        # no N at that position
                    n = 0
                self.M_positions[pos] = n+1

        except KeyError:                # no Ns
            pass
        obj_composition= {'A':0,'C':0,'G':0,'T':0,'M':0,'N':0,'propN':0, 'Z':0, 'invalid':0}

        if isinvalid:
            obj_composition['invalid']=1
        else:
            obj['invalid']=0

        for key in set(obj.keys()).intersection(['A','C','G','T','M','N']):
                obj_composition[key] = len(obj[key])
        for key in set(['A','C','G','T']) - set(obj.keys()):        # what is missing
                obj[key]=set()      # add empty set
                    
        try:
            obj_composition['propN']=len(obj['N'])/self.n_positions_examined
            obj_composition['Z']=(len(obj['N'])-self.N_mean)/self.N_sd

        except KeyError:
            pass        # no N
        self.composition[guid] = obj_composition

        # create a smaller object to store.
        smaller_obj = {}
        for item in ['A','C','G','T','M','invalid']:
            try:
                smaller_obj[item] = copy.deepcopy(obj[item])
            except KeyError:
                pass        # if it doesn't exist, that's OK	

        self.seqProfile[guid]=smaller_obj
        return obj['invalid']

    def remove(self, guid):
        """ removes a reference compressed object into RAM.
            If compression relative to other sequences has been carried out post-hoc in ram,
            only the sequence is removed; any consensus linked to it (and potentially to other sequences)
            remain unaltered.
            """
        try:
               del self.seqProfile[guid]
        except KeyError:
               pass 	# we permit attempts to delete things which don't exist

    def load(self, guid):
        """ returns a reference compressed object into RAM.
            Note: this function loads stored on disc/db relative to the reference.
            Compression relative to each other is carried out post-hoc in ram
            """
        return self.seqProfile[guid]
      
    def M_in_variant_positions(self, guids=None, variant_positions=set()):
        """ reports on the number of M bases in, and outside, the positions of variation, where:
        M bases = mixed bases
        positions of variation = a set consisting of positions of interest.  This version of the preComparer does not compute these automatically, as doing so slowed up operation.

        Parameters:
        guids: a list on which to report whether Ms are in variant positions.  If a string is passed, it is assumed to be a single guid; if None is passed, or nothing provided, all guids will be analysed.
        variant_positions: a set of positions of interest, e.g. positions of ancetral variation.  These have to be supplied by the user.
        Returns:
        a dictionary, keyed by guid, with parameters about the number of Ms in each sequence."""

        if isinstance(guids, str):
            guids = [guids]
        
        if guids is None:
            guids = set(self.seqProfile.keys())
        else:
            guids = set(guids)
            
            # check all exist
            missing = guids - set(self.seqProfile.keys())
            if len(missing)>0:
                raise KeyError("Asked to report on missing data {0}".format(missing))

        retVal = {}
        
        for guid in guids:
            if self.seqProfile[guid]['invalid'] == 1:
                retVal[guid] = {'M_total':None,
                      'M_in_vpos':None,
                      'M_not_in_vpos':None,
                      'number_variant_positions':len(variant_positions),
                      'expected_p':None,
                      'observed_p':None,
                      'p_value':None}
    
            else:

                try:
                    Ms = set(self.seqProfile[guid]['M'].keys())
                except KeyError:
                    # there are no Ms
                    Ms = set()
                nMs = len(Ms)

                Ms_in_vpos = Ms.intersection(variant_positions)
                nMs_in_vpos = len(Ms_in_vpos)
                nMs_not_in_vpos= nMs - nMs_in_vpos

                if len(variant_positions)>0:
                    observed_p = nMs_in_vpos/len(variant_positions)
                else:
                    observed_p = None
                if self.n_positions_examined-len(variant_positions)>0:
                    expected_p = nMs_not_in_vpos/(self.n_positions_examined-len(variant_positions))
                else:
                    expected_p = None

                if expected_p is None:     # we don't have an expectation, so we can't assess the binomial test;
                    p_value = None
                elif len(variant_positions)==0:      # we don't have any information to work with
                    p_value = None                
                else:  
                    if nMs_in_vpos== 0:
                        p_value = 1     # no need to compute
                    else:
                        p_value = binom_test(nMs_in_vpos,len(variant_positions), expected_p, alternative='greater')                    
                                                     
                retVal[guid]={ 
                          'M_total':nMs, 
                          'M_in_vpos':nMs_in_vpos, 
                          'M_not_in_vpos':nMs_not_in_vpos, 
                          'number_variant_positions':len(variant_positions), 
                          'expected_p':expected_p, 
                          'observed_p':observed_p,  
                          'p_value':p_value 
                            }

        return retVal

    def mcompare(self, guid, guids=None):
        """ performs comparison of one guid with 
        all guids, which are all stored samples stored with .persist().
        """

        # if guids are not specified, we do all vs all
        if guids is None:
            guids = set(self.seqProfile.keys())
        
        if not guid in self.seqProfile.keys():
            raise KeyError("Asked to compare {0}  but guid requested has not been stored.  call .persist() on the sample to be added before using mcompare.")
        
        guids = list(set(guids))       
        sampleCount = len(guids)
        neighbours = []
        
        for key2 in guids:
            if not guid==key2:
             res=self.compare(guid,key2)            
             neighbours.append(res)     
        return(neighbours)
       
    def summarise_stored_items(self):
        """ counts how many sequences exist of various types """
        retVal = {}
        retVal['server|pcstat|nSeqs'] = len(self.seqProfile.keys())
        #retVal['server|pcstat|nBinomialResults'] = len(self.binom_results.keys())
   
        return(retVal)

    def iscachedinram(self,guid):
        """ returns true or false depending whether we have a local copy of the refCompressed representation of a sequence (name=guid) in this machine """
        if guid in self.seqProfile.keys():
            return(True)
        else:
            return(False)

    def guidscachedinram(self):
        """ returns all guids with sequence profiles currently in this machine """
        retVal=set()
        for item in self.seqProfile.keys():
            retVal.add(item)
        return(retVal)

   
    def compare(self,key1, key2):
        """ compares seqProfile[key1] and seqProfile[key2]. 
        """

        # check keys are present
        if not key1 in self.seqProfile.keys():
            raise KeyError("{0} not present in preComparer".format(key1))
        if not key2 in self.seqProfile.keys():
            raise KeyError("{0} not present in preComparer".format(key2))

	# compute Ms in non-reference positions
        try:	
            m1 = set(list(self.seqProfile[key1]['M'].keys()))
        except KeyError:
            m1 = set()
        try:
    	    m2 = set(list(self.seqProfile[key2]['M'].keys()))
        except KeyError:
            m2 = set()

        ms = m1 | m2
	    
	# compute positions which differ;
        nDiff=0

        differing_positions = set()
        for nucleotide in ['C','G','A','T']:

            # we do not consider differences relative to the reference if the other nucleotide is an M
            nonM_seq1=self.seqProfile[key1][nucleotide]-ms
            nonM_seq2=self.seqProfile[key2][nucleotide]-ms
            differing_positions = differing_positions | (nonM_seq1 ^ nonM_seq2)
        
        ms= differing_positions.intersection(ms)
        n_ms = len(ms)

        destim = len(differing_positions)       # number of positions varying

        # determine whether we can 'exit fast' - stop computation at an early stage
        # optimisation discussed with Denis Volk
        # assign default return values
        res = {'guid1':key1, 
			'guid2':key2, 
			'destim':destim, 
			'upperCI':None, 
			'guid1_invalid':self.seqProfile[key1]['invalid'],
			'guid2_invalid':self.seqProfile[key2]['invalid'],
			'penalised_destim':None, 
			'guid1_Ms':self.composition[key1]['M'],
			'guid2_Ms':self.composition[key2]['M'],
			'guid1_Ns':self.composition[key1]['N'],
			'guid2_Ns':self.composition[key2]['N'],
			'guid1_Z_Ns':self.composition[key1]['Z'],
			'guid2_Z_Ns':self.composition[key2]['Z'],
			'mixed_in_cmp':n_ms,
            		'no_retest':False,
			'reported_category':'not assigned'}

        # if it is invalid, we can exit fast
        if self.seqProfile[key1]['invalid']+self.seqProfile[key2]['invalid']>0:
            res['reported_category'] = 'No retest - invalid'
            res['no_retest'] = True
            return res

        # if destim is very large, we can stop.
        if destim > self.over_selection_cutoff_ignore_factor * self.selection_cutoff:		# very large
           
                res['reported_category'] = 'No retest - High estimated distance'
                res['no_retest'] = True
        else:
                res['reported_category'] = 'Retest - Estimated distance small'
                res['no_retest'] = False
        return res


class test_preComparer_1(unittest.TestCase):
    """ tests __init__ method"""
    def runTest(self):
        # initialise comparer
        sc=preComparer(  selection_cutoff = 20,
                    mixed_reporting_cutoff=0,
                    N_mean=40000,
                    N_sd = 3000,
                    highN_z_reporting_cutoff = 2,
                    alpha = 1e-14,
                    probN_inflation_factor = 3,           
                    n_positions_examined = 3812800)

class test_preComparer_1b(unittest.TestCase):
    """ tests check_operating_parameters method"""
    def runTest(self):
        # initialise comparer
        sc=preComparer( selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor = 5,
                    mixed_reporting_cutoff=0,
                    N_mean=40000,
                    N_sd = 3000,
                    highN_z_reporting_cutoff = 2,
                    alpha = 1e-14,
                    probN_inflation_factor = 3,           
                    n_positions_examined = 3812800)

        sc.set_operating_parameters(selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor = 5,
                    mixed_reporting_cutoff=0,
                    N_mean=40000,
                    N_sd = 3000,
                    highN_z_reporting_cutoff = 2,
                    alpha = 1e-14,
                    probN_inflation_factor = 3,           
                    n_positions_examined = 3812800)

        self.assertTrue(sc.check_operating_parameters(selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor = 5,
                    mixed_reporting_cutoff=0,
                    N_mean=40000,
                    N_sd = 3000,
                    highN_z_reporting_cutoff = 2,
                    alpha = 1e-14,
                    probN_inflation_factor = 3,           
                    n_positions_examined = 3812800))

        sc.set_operating_parameters(  selection_cutoff = 10,
                    mixed_reporting_cutoff=0,
                    over_selection_cutoff_ignore_factor = 5,
                    N_mean=40000,
                    N_sd = 3000,
                    highN_z_reporting_cutoff = 2,
                    alpha = 1e-14,
                    probN_inflation_factor = 3,           
                    n_positions_examined = 3812800)

        self.assertFalse(sc.check_operating_parameters(selection_cutoff = 20,
                    mixed_reporting_cutoff=0,
                    N_mean=40000,
                    over_selection_cutoff_ignore_factor = 5,
                    N_sd = 3000,
                    highN_z_reporting_cutoff = 2,
                    alpha = 1e-14,
                    probN_inflation_factor = 3,           
                    n_positions_examined = 3812800))

class test_preComparer_2(unittest.TestCase):
    """ tests storage """
    def runTest(self):
        # initialise comparer
        sc=preComparer(  selection_cutoff = 20,
                    mixed_reporting_cutoff=0,
                    N_mean=40000,
                    N_sd = 3000,
                    highN_z_reporting_cutoff = 2,
                    alpha = 1e-14,
                    probN_inflation_factor = 3,           
                    n_positions_examined = 3812800)
        obj = {'A':set([1,2,3,4])}
        sc.persist(obj,'guid1')
        self.assertEqual(sc.composition['guid1']['A'],4)
        self.assertEqual(sc.composition['guid1']['C'],0)
        self.assertEqual(sc.composition['guid1']['T'],0)
        self.assertEqual(sc.composition['guid1']['G'],0)
        self.assertEqual(sc.composition['guid1']['N'],0)
        self.assertEqual(sc.composition['guid1']['M'],0)
        self.assertEqual(sc.composition['guid1']['invalid'],0)

class test_preComparer_3(unittest.TestCase):
    """ tests storage of invalid samples """
    def runTest(self):
        
        sc=preComparer(  selection_cutoff = 20,
                    mixed_reporting_cutoff=0,
                    N_mean=40000,
                    N_sd = 3000,
                    highN_z_reporting_cutoff = 2,
                    alpha = 1e-14,
                    probN_inflation_factor = 3,           
                    n_positions_examined = 3812800)
        obj = {'A':set([1,2,3,4]), 'invalid':1}
        sc.persist(obj,'guid1')
        self.assertEqual(sc.composition['guid1']['A'],4)
        self.assertEqual(sc.composition['guid1']['C'],0)
        self.assertEqual(sc.composition['guid1']['T'],0)
        self.assertEqual(sc.composition['guid1']['G'],0)
        self.assertEqual(sc.composition['guid1']['N'],0)
        self.assertEqual(sc.composition['guid1']['M'],0)
        self.assertEqual(sc.composition['guid1']['invalid'],1)

        obj = {'A':set([1,2,3,4]), 'invalid':0}
        sc.persist(obj,'guid2')
        self.assertEqual(sc.composition['guid2']['A'],4)
        self.assertEqual(sc.composition['guid2']['C'],0)
        self.assertEqual(sc.composition['guid2']['T'],0)
        self.assertEqual(sc.composition['guid2']['G'],0)
        self.assertEqual(sc.composition['guid2']['N'],0)
        self.assertEqual(sc.composition['guid2']['M'],0)
        self.assertEqual(sc.composition['guid2']['invalid'],0)

        obj = {'A':set([1,2,3,4]), 'N':set([10,11]), 'invalid':0}
        sc.persist(obj,'guid3')
        self.assertEqual(sc.composition['guid3']['A'],4)
        self.assertEqual(sc.composition['guid3']['C'],0)
        self.assertEqual(sc.composition['guid3']['T'],0)
        self.assertEqual(sc.composition['guid3']['G'],0)
        self.assertEqual(sc.composition['guid3']['N'],2)
        self.assertEqual(sc.composition['guid3']['M'],0)
        self.assertEqual(sc.composition['guid3']['invalid'],0)

        obj = {'A':set([1,2,3,4]), 'N':set([10,11,12]), 'M':{13:'Y'}, 'invalid':0}
        sc.persist(obj,'guid4')
        self.assertEqual(sc.composition['guid4']['A'],4)
        self.assertEqual(sc.composition['guid4']['C'],0)
        self.assertEqual(sc.composition['guid4']['T'],0)
        self.assertEqual(sc.composition['guid4']['G'],0)
        self.assertEqual(sc.composition['guid4']['N'],3)
        self.assertEqual(sc.composition['guid4']['M'],1)
        self.assertEqual(sc.composition['guid4']['invalid'],0)

        self.assertEqual(sc.guidscachedinram(),set(['guid1','guid2','guid3','guid4']))

        self.assertEqual(sc.N_positions, {10:2,11:2,12:1})
        self.assertEqual(sc.M_positions, {13:1})
class test_preComparer_4(unittest.TestCase):
    """ tests reporting of server status """
    def runTest(self):
        # initialise comparer
        sc=preComparer(  selection_cutoff = 20,
                    mixed_reporting_cutoff=0,
                    N_mean=40000,
                    N_sd = 3000,
                    highN_z_reporting_cutoff = 2,
                    alpha = 1e-14,
                    probN_inflation_factor = 3,           
                    n_positions_examined = 3812800)
        res1 = sc.summarise_stored_items()
        self.assertEqual({'server|pcstat|nSeqs': 0, 'server|pcstat|nBinomialResults': 0}, res1)

        obj = {'A':set([1,2,3,4])}
        sc.persist(obj,'guid1')
        res2 = sc.summarise_stored_items()
        self.assertEqual({'server|pcstat|nSeqs': 1, 'server|pcstat|nBinomialResults': 0}, res2)

class test_preComparer_5(unittest.TestCase):
    """ tests comparison """
    def runTest(self):
        
        sc=preComparer(  selection_cutoff = 20,
                    mixed_reporting_cutoff=0,
                    N_mean=40000,
                    N_sd = 3000,
                    highN_z_reporting_cutoff = 2,
                    alpha = 1e-14,
                    probN_inflation_factor = 3,           
                    n_positions_examined = 3812800)
       
        obj = {'A':set([1,2,3,4]), 'invalid':0}
        sc.persist(obj,'guid2')
      
        obj = {'A':set([1,2,3,4]), 'N':set([10,11]), 'invalid':0}
        sc.persist(obj,'guid3')

        res = sc.compare('guid2','guid3')
        self.assertEqual(res['reported_category'],'Retest - Estimated distance small')


class test_preComparer_6(unittest.TestCase):
    """ tests comparison """
    def runTest(self):
        
        sc=preComparer(  selection_cutoff = 20,
                    mixed_reporting_cutoff=0,
                    N_mean=40000,
                    N_sd = 3000,
                    highN_z_reporting_cutoff = 2,
                    alpha = 1e-14,
                    probN_inflation_factor = 3,           
                    n_positions_examined = 3812800)
       
        obj = {'A':set([1,2,3,4]), 'invalid':0}
        sc.persist(obj,'guid2')
      
        obj = {'A':set([3,4,5,6]), 'N':set([10,11]), 'invalid':0}
        sc.persist(obj,'guid3')

        res = sc.compare('guid2','guid3')
        self.assertEqual(res['reported_category'],'Retest - Estimated distance small')

class test_preComparer_7(unittest.TestCase):
    """ tests comparison """
    def runTest(self):
        
        sc=preComparer(  selection_cutoff = 20,
                    mixed_reporting_cutoff=0,
                    N_mean=40000,
                    N_sd = 3000,
                    highN_z_reporting_cutoff = 2,
                    alpha = 1e-14,
                    probN_inflation_factor = 3,           
                    n_positions_examined = 3812800)
       
        obj = {'A':set([1,2,3,4]), 'invalid':0}
        sc.persist(obj,'guid2')
      
        obj = {'A':set([3,4,5,6]), 'N':set([10,11]), 'invalid':1}
        sc.persist(obj,'guid3')
 
        # check invalid set
        self.assertEqual(sc.seqProfile['guid3']['invalid'], 1)
        self.assertEqual(sc.seqProfile['guid2']['invalid'], 0)
              
        res = sc.compare('guid2','guid3')
        self.assertEqual(res['reported_category'],'No retest - invalid')
        res = sc.summarise_stored_items()
        self.assertEqual(res['server|pcstat|nBinomialResults'],0)

class test_preComparer_8(unittest.TestCase):
    """ tests comparison """
    def runTest(self):
        
        sc=preComparer(  selection_cutoff = 20,
                    mixed_reporting_cutoff=0,
                    N_mean=40000,
                    N_sd = 3000,
                    highN_z_reporting_cutoff = 2,
                    alpha = 1e-14,
                    probN_inflation_factor = 3,           
                    n_positions_examined = 3812800)
       
        obj = {'A':set([1,2,3,4]), 'invalid':0}
        sc.persist(obj,'guid2')
      
        # check all keys are present
        self.assertEqual(set(sc.seqProfile['guid2'].keys()), set(['A','C','G','T','invalid']))
       
        obj = {'A':set([1,2,3,4]), 'invalid':1}
        sc.persist(obj,'guid3')
      
        # check all keys are present
        self.assertEqual(set(sc.seqProfile['guid3'].keys()), set(['A','C','G','T','invalid']))

        # check invalid set
        self.assertEqual(sc.seqProfile['guid3']['invalid'], 1)
 
class test_preComparer_9(unittest.TestCase):
    """ tests mcompare """
    def runTest(self):
        
        sc=preComparer(  selection_cutoff = 20,
                    mixed_reporting_cutoff=0,
                    N_mean=40000,
                    N_sd = 9000,
                    highN_z_reporting_cutoff = 2,
                    alpha = 1e-14,
                    probN_inflation_factor = 3,           

                    n_positions_examined = 3812800)
       
        obj = {'A':set([1,2,3,4]), 'invalid':0}
        sc.persist(obj,'guid2')
      
        # check all keys are present
        self.assertEqual(set(sc.seqProfile['guid2'].keys()), set(['A','C','G','T','invalid']))
       
        obj = {'A':set([1,2,3,4,5]), 'invalid':0}
        sc.persist(obj,'guid3')
      
        res = sc.mcompare('guid2')
        self.assertEqual(len(res),1)

class test_preComparer_10(unittest.TestCase):
    """ tests M_in_variant_positions """
    def runTest(self):
        
        sc=preComparer(  selection_cutoff = 20,
                    mixed_reporting_cutoff=0,
                    N_mean=40000,
                    N_sd = 3000,
                    highN_z_reporting_cutoff = 2,
                    alpha = 1e-14,
                    probN_inflation_factor = 3,           

                    n_positions_examined = 3812800)
     
        res = sc.M_in_variant_positions()
        self.assertEqual(res, {})
 
        obj = {'A':set([1,2,3,4]), 'invalid':0}
        sc.persist(obj,'guid2')
 
        res = sc.M_in_variant_positions(variant_positions=set([1,2,3,4]))
        self.assertEqual(res, {'guid2':{'M_total':0, 'M_in_vpos':0, 'M_not_in_vpos':0, 'number_variant_positions':4, 'expected_p':0, 'observed_p':0, 'p_value':1}})

        obj = {'A':set([1,2,3,4]), 'M':{5:'Y',6:'Y'}, 'invalid':0}
        sc.persist(obj,'guid3')
        res = sc.M_in_variant_positions(variant_positions=set([1,2,3,4]))
        self.assertEqual(len(res.keys()), 2)
        self.assertEqual(res['guid3']['M_not_in_vpos'],2)



        obj = {'M':{1:'Y',2:'Y',3:'Y', 4:'Y'}, 'invalid':0}
        sc.persist(obj,'guid4')
        res = sc.M_in_variant_positions(variant_positions=set([1,2,3,4]))
        self.assertEqual(res['guid4']['M_not_in_vpos'],0)
        self.assertEqual(res['guid4']['M_in_vpos'],4)
        self.assertEqual(res['guid4']['p_value'],0)
