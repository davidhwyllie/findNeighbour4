#!/usr/bin/env python3

# python3 code to compare fasta sequences
import unittest
import os
import glob
import sys

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
from scipy.stats import binom_test
import pandas as pd
from collections import Counter

# only used for unit testing
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide

class seqComparer():
    def __init__(self,
                    reference,
                    maxNs,
                    snpCeiling,

                    excludePositions=set(),
                    return_none_if_high_snp_distance=True
                ):

        """ instantiates the sequence comparer, an object which manages reference compressed sequences.
        
        It does not manage persistence, nor does it automatically load sequences.
        
        reference is a string consisting of the reference sequence.
        This is required because, as a data compression technique,
        only differences from the reference are stored.
               
        excludePositions contains a zero indexed set of bases which should not be considered at all in the sequence comparisons.
        Any bases which are always N should be added to this set.
        Not doing so will substantially degrade the algorithm's performance.
     
        self.return_none_if_high_snp_distance: if True (default), does not report high SNP distances; returns None if the distance is higher than a SNP cutoff

        If the number of Ns are more than maxNs, no data from the sequence is stored.
         
       
        Results > snpCeiling are not returned or stored.
            
        unknown_base_type is either N or M, and is used for computation of mixture statistics
        David Wyllie, Nov 2018
        
        - to run unit tests, do
        python3 -m unittest seqComparer
        """
        
     
        # we support three kinds of sequences.
        # sequence in strings;
        # reference based compression relative to reference 'compressed_sequence';
        self.return_none_if_high_snp_distance=return_none_if_high_snp_distance 
        self.compressed_sequence_keys = set(['invalid','A','C','G','T', 'N', 'M', 'U'])  
        self.snpCeiling = snpCeiling
              
        # sequences with more than maxNs Ns will be considered invalid and their details (apart from their invalidity) will not be stored.
        self.maxNs=maxNs

        # check composition of the reference.
        self.reference=str(reference)           # if passed a Bio.Seq object, coerce to string.
        letters=Counter(self.reference)   
        if len(set(letters.keys())-set(['A','C','G','T']) )>0:
            raise TypeError("Reference sequence supplied contains characters other than ACTG: {0}".format(letters))
                       
        # load the excluded bases
        self.excluded=excludePositions
        
        # define what is included
        self.included=set(range(len(self.reference)))-self.excluded
        
        # initialise pairwise sequences for comparison.
        self._refresh()

 
    def persist(self, object, guid):
        """ keeps a reference compressed object into RAM.
            Note: the sequences are stored on disc/db relative to the reference.
            Compression relative to each other is carried out post-hoc in ram
            """

        # older databases don't store the M/N combination positions in the 'U' key
        # we create it on load into RAM
        if not 'U' in object.keys():
            object['U'] = set()
            if 'N' in object.keys():
                object['U']=object['N']     # if no Ms, then this is a reference taking ~ no memory
            if 'M' in object.keys():
                if len(object['M'])>0:
                    object['U']=object['N'].intersection(object['M'])
            
        self.seqProfile[guid]=object

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

    def remove_all(self):
        """ empties any sequence data from ram """
        self._refresh()
               
    def _refresh(self):
        """ empties any sequence data from ram """
        self.seqProfile={}
         
    def mcompare(self, guid, guids=None):
        """ performs comparison of one guid with 
        all guids, which are also stored samples.
        """

        # if guids are not specified, we do all vs all
        if guids is None:
            guids = set(self.seqProfile.keys())
        
        if not guid in self.seqProfile.keys():
            raise KeyError("Asked to compare {0}  but guid requested has not been stored.  call .persist() on the sample to be added before using mcompare.".format(guid))
        
        guids = list(set(guids))       
        sampleCount = len(guids)
        neighbours = []
        
        for key2 in guids:
            if not guid==key2:
                (guid1,guid2,dist,n1,n2,nboth, N1pos, N2pos, Nbothpos)=self.countDifferences_byKey(keyPair=(guid,key2),
                                                                                                      cutoff = self.snpCeiling)            
                neighbours.append([guid1,guid2,dist,n1,n2,nboth,N1pos, N2pos, Nbothpos])

       
        return(neighbours)
 
    def distmat(self, half=True, diagonal=False):
        """ returns a distance matrix.
        parameters
        If half=True, returns only half the matrix.
        If diagonal=True, return the diagonal (all of which are zeros)
        
        returns:
        A generator yielding neighbours, in format [guid1,guid2,dist,n1,n2,nboth,N1pos, N2pos, Nbothpos]
        """
        
        for guid1 in self.seqProfile.keys():
            for guid2 in self.seqProfile.keys():
                include = False
                if not half and not guid1==guid2:
                    include = True
                if guid1==guid2 and diagonal:
                    include = True
                if guid1<guid2 and half:
                    include = True

                if include:
                    yield self.countDifferences_byKey(keyPair=(guid1,guid2))            
        
    def summarise_stored_items(self):
        """ counts how many sequences exist of various types """
        retVal = {}
        retVal['server|scstat|nSeqs'] = len(self.seqProfile.keys())
              
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

    def _delta(self,x):
        """ returns the difference between two numbers in a tuple x """
        return(x[1]-x[0])

    def excluded_hash(self):
        """ returns a string containing the number of nt excluded, and a hash of their positions.
        This is useful for version tracking & storing patterns of masking. """
        l = sorted(list(self.excluded))
        len_l = len(l)
        h = hashlib.md5()
        h.update(json.dumps(l).encode('utf-8'))
        md5_l = h.hexdigest()
        return("Excl {0} nt [{1}]".format(len_l, md5_l))
    
    def uncompress(self, compressed_sequence):
        """ returns a sequence from a compressed_sequence """
        if 'invalid' in compressed_sequence.keys():
            if compressed_sequence['invalid']==1:
                raise ValueError("Cannot uncompress an invalid sequence, because the sequence it is not stored {0}".format(compressed_sequence.keys()))
           
        seq = list(self.reference)
        
        # mark all positions excluded as N
        for x in self.excluded:
            seq[x]='N'
        
        # mark any bases definitively called as whatever they are called as
        for item in ['A','C','T','G','N']:
            for x in compressed_sequence[item]:     # these are different from reference;
                seq[x]=item
                
        # add any mixed bases
        for x in compressed_sequence['M'].keys():   # the 'M' option includes iupac coding
            seq[x]= compressed_sequence['M'][x]
        return(''.join(seq))
    
    def compress(self, sequence):
        """ reads a string sequence and extracts position - genome information from it.
        returns a dictionary consisting of zero-indexed positions of non-reference bases.
        
        """
        if not len(sequence)==len(self.reference):
            raise TypeError("sequence must of the same length as reference; seq is {0} and ref is {1}".format(len(sequence),len(self.reference)))
        if len(self.reference)==0:
            raise TypeError("reference cannot be of zero length")
               
        # we consider - characters to be the same as N
        sequence=sequence.replace('-','N')
        
        # we only record differences relative to to refSeq.
        # anything the same as the refSeq is not recorded.
        # a dictionary, M, records the mixed base calls.
        # we also store a set containing either N or M.  This is wasteful of RAM, and the need for it could be removed but
        # it is included here as it massively increases computational speed; the process of combining sets of Ns and Ms is very expensive.

        diffDict={ 'A':set([]),'C':set([]),'T':set([]),'G':set([]),'N':set([]), 'M':{}, 'U':set([])}        

        for i in self.included:     # for the bases we need to compress
            if not sequence[i]==self.reference[i]:      # if it's not reference
                if sequence[i] in ['A','C','T','G','N']:
                    diffDict[sequence[i]].add(i)        # if it's a definitively called base
                else:
                    # we regard it as a code representing a mixed base.  we store the results in a dictionary
                    diffDict['M'][i] = sequence[i]
                    diffDict['U'].add(i)
            if sequence[i] == 'N':
                diffDict['U'].add(i)
			                 
        # check how many Ns or Ms  
        if len(diffDict['U'])>self.maxNs:
            # we store it, but not with sequence details if is invalid
            diffDict={'invalid':1}
        else:
            diffDict['invalid']=0
            
        return(diffDict)
            
    def _setStats(self, i1, i2):
        """ compares either:
        * two sets (if i1 or i2 is a set)
        OR
        * the keys of the dictionaries i1 or i2 is a dictionary

        returns 
        * the number of elements in i1
        * the number of elements in i2
        * the number of elements in the union of i1 and i2
        * i1
        * i2
        * the union of set1 and set2   
        
        """
        if isinstance(i1,set):
            retVal1= i1
        elif isinstance(i1,dict):
            retVal1= set(i1.keys())
        if isinstance(i2,set):
            retVal2= i2
        elif isinstance(i2,dict):
            retVal2= set(i2.keys())
               
        retVal=retVal2 | retVal1

        return(len(retVal1), len(retVal2), len(retVal), retVal1, retVal2, retVal)   

    def countDifferences_byKey(self, keyPair, cutoff=None):
        """ compares the in memory refCompressed sequences at
        self.seqProfile[key1] and self.seqProfile[key2]

        Returns the number of SNPs between self._seq1 and self._seq2, and,
        if the pairwise SNP distance is less than cutoff,
        the number of N and Ms in the two sequences and the union of their positions.
        """

        if not type(keyPair) is tuple:
            raise TypeError("Wanted tuple keyPair, but got keyPair={0} with type {1}".format(keyPair, type(keyPair)))
        if not len(keyPair)==2:
            raise TypeError("Wanted a keyPair with two elements, but got {0}".format(keyPair))
        
        ## test the keys exist
        (key1,key2)=keyPair
        if not key1 in self.seqProfile.keys():
            raise KeyError("Key1={0} does not exist in the in-memory store.".format(key1))
        if not key2 in self.seqProfile.keys():
            raise KeyError("Key1={0} does not exist in the in-memory store.".format(key1))
         
        # if cutoff is not specified, we use snpCeiling
        if cutoff is None:
            cutoff = self.snpCeiling
            
        ## do the computation  
        nDiff=self.countDifferences(self.seqProfile[key1],self.seqProfile[key2],cutoff=cutoff)

        if nDiff is None:
            return((key1, key2, nDiff, None, None, None, None, None, None))
        elif nDiff<=cutoff:
            seq1_uncertain = self.seqProfile[key1]['N'] 
            seq2_uncertain = self.seqProfile[key2]['N'] 
            (n1, n2, nboth, N1pos, N2pos, Nbothpos) = self._setStats(seq1_uncertain, seq2_uncertain)
            return((key1, key2, nDiff, n1,n2,nboth, N1pos, N2pos, Nbothpos))
        else:
            return((key1, key2, nDiff, None, None, None, None, None, None))

    
    def countDifferences(self, seq1, seq2, cutoff=None):
        """ compares seq1 with seq2.
        
        Ns and Ms (uncertain bases) are ignored in snp computations.

	"""
        #  if cutoff is not specified, we use snpCeiling
        if cutoff is None:
            cutoff = self.snpCeiling
     
        nDiff=0
        if seq1 is None or seq2 is None:
            return None 
                 
        if seq1['invalid']==1 or seq2['invalid']==1:
            return None 
         
        # compute positions which differ;
        differing_positions = set()
        for nucleotide in ['C','G','A','T']:
       
            # we do not consider differences relative to the reference if the other nucleotide is an N or M
            nonN_seq1=seq1[nucleotide]-seq2['U']
            nonN_seq2=seq2[nucleotide]-seq1['U']
            differing_positions = differing_positions | (nonN_seq1 ^ nonN_seq2)
        
            # if the number of differences is already larger than the cutoff,
            # then we do not need to perform additional computations; we return None,
            # which is indicated there is no link less than or equal to cutoff
            # unless we are told to return exact differences
            if self.return_none_if_high_snp_distance and len(differing_positions) > cutoff:
                return None

        nDiff = len(differing_positions)
        
        if self.return_none_if_high_snp_distance and nDiff>cutoff:
            return(None)
        else:
            return(nDiff)
       
    def compressed_sequence_hash(self, compressed_sequence):
        """ returns a string containing a hash of a compressed object.
        Used for identifying compressed objects, including consensus sequences.
        """
        keys = sorted(compressed_sequence.keys())
        serialised_compressed_sequence = ""
        for key in keys:
            if isinstance(compressed_sequence[key], set):
                l = sorted(list(compressed_sequence[key]))
            else:
                l = compressed_sequence[key]
            serialised_compressed_sequence = serialised_compressed_sequence + key + ":" + str(l) + ';'
        h = hashlib.md5()
        h.update(serialised_compressed_sequence.encode('utf-8'))
        md5 = h.hexdigest()
        return(md5)
    
    
    def estimate_expected_proportion(self, seqs):
        """ computes the median Ns for seqs, a list.
        Returns None if
        * the length of the sequences in seqs are 0, or
        * there are <= 3 seqs
        """
        if len(seqs) < 3:
            return None
        if len(seqs[0]) == 0:
            return None
        Ns = []
        for seq in seqs:
            Ns.append(seq.count('N'))
        return np.median(Ns)/len(seqs[0])
    
    def estimate_expected_unk(self, sample_size=30, exclude_guids=set(), unk_type='N'):
        """ computes the median unk_type for sample_size guids, randomly selected from all guids except for exclude_guids.
        Used to estimate the expected number of unk_type bases in an alignment.
        unk_type can be one of 'N' 'M' 'N_or_M'.
        """
        
        guids = list(set(self.seqProfile.keys())-set(exclude_guids))
        np.random.shuffle(list(guids))
        
        retVal = None       # cannot compute 
        unks = []
        for guid in guids:
            this_unk = 0
            try:

                if unk_type in ['N','N_or_M']:
                    this_unk = this_unk + len(self.seqProfile[guid]['N'])
                if unk_type in ['M','N_or_M']:
                    this_unk = this_unk + len(self.seqProfile[guid]['M'])
                unks.append(this_unk)
            except KeyError:
                # it is invalid
                pass
            if len(unks)>=sample_size:
                break
        if len(unks)>=sample_size:     
            return np.median(unks)
        else: 
            return None
  
    def estimate_expected_unk_sites(self, sample_size=30, sites = set(), exclude_guids=set(), unk_type='N'):
        """ computes the median unk_type for sample_size guids, randomly selected from all guids except for exclude_guids.
        Only reports unk_type bases at sites().
        Used to estimate the expected number of unk_type (one of 'N','M',or 'N_or_M' in an alignment """
        
        guids = list(set(self.seqProfile.keys())-set(exclude_guids))
        np.random.shuffle(list(guids))
  
        retVal = None       # cannot compute 
        unks = []
        for guid in guids:
            try:
                this_unk = 0
                if unk_type in ['N','N_or_M']:
                    this_unk = this_unk + len(self.seqProfile[guid]['N'].intersection(sites))
                if unk_type in ['M','N_or_M']:
                    this_unk = this_unk + len(set(self.seqProfile[guid]['M'].keys()).intersection(sites))
                unks.append(this_unk)
            except ValueError:
                # it is invalid
                pass
            if len(unks)>=sample_size:
                break
        if len(unks)>=sample_size:     
            return np.median(unks)
        else: 
            return None
    
    def multi_sequence_alignment(self, guids, output='dict', sample_size=30, expected_p1=None, uncertain_base_type='N'):
        """ computes a multiple sequence alignment containing only sites which vary between guids.
        
        sample_size is the number of samples to randomly sample to estimate the expected number of Ns in
        the population of sequences currently in the server.  From this, the routine computes expected_p1,
        which is expected_expected_N/ the length of sequence.
        if expected_p1 is supplied, then such sampling does not occur.
        
        output can be either
        'dict', in which case the output is presented as dictionaries mapping guid to results; or
        'df' in which case the results is a pandas data frame like the below, where the index consists of the
        guids identifying the sequences, or

            (index)      aligned_seq  allN  alignN   p_value
            AAACGN-1        AAAC     1       0  0.250000
            CCCCGN-2        CCCC     1       0  0.250000
            TTTCGN-3        TTTC     1       0  0.250000
            GGGGGN-4        GGGG     1       0  0.250000
            NNNCGN-5        NNNC     4       3  0.003906
            ACTCGN-6        ACTC     1       0  0.250000
            TCTNGN-7        TCTN     2       1  0.062500
            AAACGN-8        AAAC     1       0  0.250000
            
        'df_dict'.  This is a serialisation of the above, which correctly json serialised.  It can be turned back into a
        pandas DataFrame as follows:
        
        res= sc.multi_sequence_alignment(guid_names[0:8], output='df_dict')     # make the dictionary, see unit test _47
        df = pd.DataFrame.from_dict(res,orient='index')                         # turn it back.
        
        The p values reported are derived from exact, one-sided binomial tests as implemented in python's scipy.stats.binom_test().

        TEST 1:
        This tests the hypothesis that the number of Ns in the *alignment*
        is GREATER than those expected from the expected_N in the population of whole sequences.
 
        Does so by comparing the observed number of Ns in the alignment (alignN),
        given the alignment length (4 in the above case) and an expectation of the proportion of bases which will be N.
        The expected number of Ns is estimated by
        i) randomly sampling sample_size guids from those stored in the server and
        observing the number of Ns per base across the genome.  The estimate_expected_unk() function performs this.
        ii) randomly sampling sample_size guids from those stored in the server and
        observing the number of Ns per base across the relevant  genome.  The estimate_expected_unk() function performs this.
          
        This approach determines the median number of Ns in valid sequences, which (if bad samples with large Ns are rare)
        is a relatively unbiased estimate of the median number of Ns in the good quality samples.
        
        If there  are not enough samples in the server to obtain an estimate, p_value is not computed, being
        reported as None.
  
        TEST 2:
        This tests the hypothesis that the number of Ns in the *alignment*
        is GREATER than those expected from the expected_N in the population of whole sequences
        *at the bases examined in the alignment*.
        This might be relevant if these particular bases are generally hard to call.
 
        Does so by comparing the observed number of Ns in the alignment (alignN),
        given the alignment length (4 in the above case) and an expectation of the proportion of bases which will be N.
        The expected number of Ns is estimated by randomly sampling sample_size guids from those stored in the server and
        observing the number of Ns per base at the relevant sites.  The estimate_expected_unk_sites() function performs this.
   
        This approach determines the median number of Ns in valid sequences, which (if bad samples with large Ns are rare)
        is a relatively unbiased estimate of the median number of Ns in the good quality samples.
        
        If there  are not enough samples in the server to obtain an estimate, p_value is not computed, being
        reported as None.
             
        TEST 3: tests whether the proportion of Ns in the alignment is greater
        than in the bases not in the alignment, for this sequence.
        
        ## TEST 4:  
        tests whether the proportion of Ns in the alignment  for this sequence
        is greater than the median proportion of Ns in the alignment for all other sequences.
        This is a sensible option if the test is performed per-cluster, and the proportion of Ns in the
        aligned sequences differs markedly by cluster.
        
        This test is computed if there are four or more samples in the cluster, and if the alignment length is non-zero.

        """
        
        # -1 validate input
        if expected_p1 is not None:
            if expected_p1 < 0 or expected_p1 > 1:
                raise ValueError("Expected_p1 must lie between 0 and 1")
        if sample_size is None:
            sample_size = 30


        # step 0: find all valid guids

        valid_guids = []
        invalid_guids = []

        comparatorSeq = {}
        for guid in guids:
            try:
                comparatorSeq[guid] = self.seqProfile[guid]
                valid_guids.append(guid)
            except ValueError:
                invalid_guids.append(guid)
        # Estimate expected N as median(observed Ns),
        # which is a valid thing to do if the proportion of mixed samples is low.
        
        if expected_p1 is None:
            expected_N1 = self.estimate_expected_unk(sample_size=sample_size, exclude_guids= invalid_guids)
            if expected_N1 is None:
                expected_p1 = None
            else:
                expected_p1 = expected_N1 / len(self.reference)
        else:
            expected_N1 = np.floor(expected_p1 * len(self.reference))
            
        return self._msa(valid_guids, invalid_guids, expected_p1, output, sample_size, uncertain_base_type)
    
    def _msa(self, valid_guids, invalid_guids, expected_p1, output, sample_size, uncertain_base_type='N'):
        """ perform multisequence alignment and significance tests.
        
        Parameters:
        valid_guids:  the guids in valid_guids are those on which the msa is computed.
        
        Additionally, the software performs four significance tests
        on the number of uncertain_bases in an alignment containing only the variant bases between variant_guids.
        
        If uncertain_base_type == 'N' analyses Ns.
        If uncertain_base_type == 'M' analyses mixed bases.
        If uncertain_base_type == 'N_or_M' analyses both combined
          
        expected_p1:  The expected proportion of Ns or Ms (as specified by uncertain_base_type) is expected_p1.
        
        invalid_guids: these are not analysed when computed expected numbers of uncertain_base_type.
        
        output:   can be either
        'dict', in which case the output is presented as dictionaries mapping guid to results; or
        'df' in which case the results is a pandas data frame like the below, where the index consists of the
        guids identifying the sequences, or

            (index)      aligned_seq  allN  alignN   p_value
            AAACGN-1        AAAC     1       0  0.250000
            CCCCGN-2        CCCC     1       0  0.250000
            TTTCGN-3        TTTC     1       0  0.250000
            GGGGGN-4        GGGG     1       0  0.250000
            NNNCGN-5        NNNC     4       3  0.003906
            ACTCGN-6        ACTC     1       0  0.250000
            TCTNGN-7        TCTN     2       1  0.062500
            AAACGN-8        AAAC     1       0  0.250000
            
        'df_dict'.  This is a serialisation of the above, which correctly json serialised.  It can be turned back into a
        pandas DataFrame as follows:
        
        res= sc.multi_sequence_alignment(guid_names[0:8], output='df_dict')     # make the dictionary, see unit test _47
        df = pd.DataFrame.from_dict(res,orient='index')                         # turn it back.
        
        The p values reported are derived from exact, one-sided binomial tests as implemented in pythons scipy.stats.binom_test().
        
        TEST 1:
        This tests the hypothesis that the number of Ns in the *alignment*
        is GREATER than those expected from the expected_N in the population of whole sequences.
 
        Does so by comparing the observed number of Ns in the alignment (alignN),
        given the alignment length (4 in the above case) and an expectation of the proportion of bases which will be N.
        The expected number of Ns is estimated by
        i) randomly sampling sample_size guids from those stored in the server and
        observing the number of Ns per base across the genome.  The estimate_expected_unk() function performs this.
        ii) randomly sampling sample_size guids from those stored in the server and
        observing the number of Ns per base across the relevant  genome.  The estimate_expected_unk() function performs this.
          
        This approach determines the median number of Ns in valid sequences, which (if bad samples with large Ns are rare)
        is a relatively unbiased estimate of the median number of Ns in the good quality samples.
        
        If there  are not enough samples in the server to obtain an estimate, p_value is not computed, being
        reported as None.
  
        TEST 2:
        This tests the hypothesis that the number of Ns in the *alignment*
        is GREATER than those expected from the expected_N in the population of whole sequences
        *at the bases examined in the alignment*.
        This might be relevant if these particular bases are generally hard to call.
 
        Does so by comparing the observed number of Ns in the alignment (alignN),
        given the alignment length (4 in the above case) and an expectation of the proportion of bases which will be N.
        The expected number of Ns is estimated by randomly sampling sample_size guids from those stored in the server and
        observing the number of Ns per base at the relevant sites.  The estimate_expected_unk_sites() function performs this.
   
        This approach determines the median number of Ns in valid sequences, which (if bad samples with large Ns are rare)
        is a relatively unbiased estimate of the median number of Ns in the good quality samples.
        
        If there  are not enough samples in the server to obtain an estimate, p_value is not computed, being
        reported as None.
                  
        TEST 3: tests whether the proportion of Ns in the alignment is greater
        than in the bases not in the alignment, for this sequence.
        
        TEST 4: tests whether the proportion of Ns in the alignment  for this sequence
        is greater than the median proportion of Ns in the alignment for all other sequences.
        This is a sensible option if the test is performed per-cluster, and the proportion of Ns in the
        aligned sequences differs markedly by cluster.
        
        This test is computed if there are two or more samples in the cluster.

        """
        

        # step 1: define all positions in these sequences for which there is a non-reference position
        # we ignore Ms and Ns in this analysis
        # we also compute the total number of Ms or Ns in each guid
        nrps = {}
        guid2all = {'N':{},'M':{}, 'N_or_M':{}}
        guid2align = {'N':{},'M':{}, 'N_or_M':{}}

        for guid in valid_guids:
            if self.seqProfile[guid]['invalid'] ==1:
                raise TypeError("Invalid sequence {0} passed in valid_guids".format(guid))    

            for unk_type in ['N','M']:

                unks = len(self.seqProfile[guid][unk_type])
                guid2all[unk_type][guid] = unks

            guid2all['N_or_M'][guid] = guid2all['N'][guid]+guid2all['M'][guid]
            
            for base in ['A','C','T','G']:
                positions = self.seqProfile[guid][base]
                for position in positions:
                  if not position in nrps.keys():     # if it's non-reference, and we've got no record of this position
                     nrps[position]=set()             # then we generate a set of bases at this position                   
                  nrps[position].add(base)            # either way add the current non-reference base there
                  
        # step 2: for the non-reference called positions, check if there's a reference base there.
        for guid in valid_guids:
            for position in nrps.keys():
                psn_accounted_for  = 0
                for base in ['A','C','T','G']:
                    if position in self.seqProfile[guid][base]:
                        psn_accounted_for = 1
                if psn_accounted_for ==0 :
                    # it is reference; this guid has no record of a variant base at this position, so it must be reference.
                    nrps[position].add(self.reference[position])
                 
        # step 3: find those which have multiple bases at a position
        variant_positions = set()
        for position in nrps.keys():
            if len(nrps[position])>1:
                variant_positions.add(position)

        # step 4: determine the sequences of all bases.
        ordered_variant_positions = sorted(list(variant_positions))
        guid2seq = {}
        guid2msa_seq={}
        for guid in valid_guids:
            guid2seq[guid]=[]
            for position in ordered_variant_positions:
                this_base = self.reference[position]
                for base in ['A','C','T','G','N','M']:
                    if not base == 'M':
                        positions = self.seqProfile[guid][base]
                    else:
                        positions = self.seqProfile[guid][base].keys()

                    if position in positions:
                        this_base = base
                guid2seq[guid].append(this_base)
            guid2msa_seq[guid] = ''.join(guid2seq[guid])
        
        # step 5: determine the expected_p2 at the ordered_variant_positions:

        expected_N2 = self.estimate_expected_unk_sites(sample_size=sample_size, exclude_guids= invalid_guids, sites=set(ordered_variant_positions), unk_type = uncertain_base_type)
        if expected_N2 is None:
            expected_p2 = None
        elif len(ordered_variant_positions) is 0:
            expected_p2 = None
        else:
            expected_p2 = expected_N2 / len(ordered_variant_positions)
        
        # step 6: perform Binomial tests on all samples
        if len(valid_guids)>0:
            guid2pvalue1 = {}
            guid2pvalue2 = {}
            guid2pvalue3 = {}
            guid2pvalue4 = {}
            guid2observed_p = {}
            guid2expected_p1 = {}
            guid2expected_p2 = {}
            guid2expected_p3 = {}
            guid2expected_p4 = {}
            for guid in valid_guids:
                
                # compute p value 1.  This tests the hypothesis that the number of Ns in the *alignment*
                # is GREATER than those expected from the expected_N in the population of whole sequences.
                for unk_type in ['N','M']:
                    guid2align[unk_type][guid]= guid2msa_seq[guid].count(unk_type)
                guid2align['N_or_M'][guid] = guid2align['N'][guid] + guid2align['M'][guid]
                
                if expected_p1 is None:     # we don't have an expectation, so we can't assess the first binomial test;
                    p_value1 = None
                    observed_p = None
                elif len(guid2msa_seq[guid])==0:      # we don't have any information to work with
                    p_value1 = None
                    observed_p = None                    
                else:  
                    observed_p = guid2align[uncertain_base_type][guid]/len(guid2msa_seq[guid])
                    p_value1 = binom_test(guid2align[uncertain_base_type][guid],len(guid2msa_seq[guid]), expected_p1, alternative='greater')
                    #print(uncertain_base_type, guid2msa_seq[guid],guid2align[uncertain_base_type][guid],len(guid2msa_seq[guid]), observed_p, expected_p1, p_value1)
                    
                guid2pvalue1[guid]=p_value1
                guid2observed_p[guid]=observed_p
                guid2expected_p1[guid]=expected_p1
                
                # compute p value 2.  This tests the hypothesis that the number of Ns in the *alignment*
                # is GREATER than those expected from the expected_N in the population of whole sequences
                # at these sites.
                if expected_p2 is None:     # we don't have an expectation, so we can't assess the binomial test;
                    p_value2 = None
                elif len(guid2msa_seq[guid])==0:      # we don't have any information to work with
                    p_value2 = None                
                else:  
                    p_value2 = binom_test(guid2align[uncertain_base_type][guid],len(guid2msa_seq[guid]), expected_p2, alternative='greater')                    
                guid2pvalue2[guid]=p_value2
                guid2expected_p2[guid]=expected_p2
                
                                
                # compute p value 3.  This tests the hypothesis that the number of Ns in the alignment of THIS SEQUENCE
                # is GREATER than the number of Ns not in the alignment  IN THIS SEQUENCE
                # based on sequences not in the alignment

                expected_p3 = (guid2all[uncertain_base_type][guid]-guid2align[uncertain_base_type][guid])/(len(self.reference)-len(guid2msa_seq[guid]))
                p_value = binom_test(guid2align[uncertain_base_type][guid],len(guid2msa_seq[guid]), expected_p3, alternative='greater')
                guid2pvalue3[guid]=p_value
                guid2expected_p3[guid]=expected_p3
                
                # compute p value 4.
                # tests whether the proportion of Ns in the alignment  for this sequence
                # is greater than the median proportion of Ns in the alignment for all other sequences.
                # This is a sensible option if the test is performed per-cluster, and the proportion of Ns in the
                # aligned sequences differs markedly by cluster.
                # This test is computed if there are two or more samples in the cluster.
               
                # get the sequences which are not guid;
                other_seqs = []
                for this_guid in guid2msa_seq.keys():
                    if not guid == this_guid:
                        other_seqs.append(guid2msa_seq[this_guid])
                expected_p4 = self.estimate_expected_proportion(other_seqs)
                if expected_p4 is None:
                    p_value = None
                else:
                    p_value = binom_test(guid2align[uncertain_base_type][guid],len(guid2msa_seq[guid]), expected_p4, alternative='greater')
                guid2pvalue4[guid]=p_value
                guid2expected_p4[guid]=expected_p4
                 
                
            # assemble dataframe
            df1 = pd.DataFrame.from_dict(guid2msa_seq, orient='index')
            df1.columns=['aligned_seq']
            df1['aligned_seq_len'] = len(ordered_variant_positions)
            
            df2n = pd.DataFrame.from_dict(guid2all['N'], orient='index')
            df2n.columns=['allN']
            df3n = pd.DataFrame.from_dict(guid2align['N'], orient='index')
            df3n.columns=['alignN']

            df2m = pd.DataFrame.from_dict(guid2all['M'], orient='index')
            df2m.columns=['allM']
            df3m = pd.DataFrame.from_dict(guid2align['M'], orient='index')
            df3m.columns=['alignM']

            df2mn = pd.DataFrame.from_dict(guid2all['N_or_M'], orient='index')
            df2mn.columns=['allN_or_M']
            df3mn = pd.DataFrame.from_dict(guid2align['N_or_M'], orient='index')
            df3mn.columns=['alignN_or_M']
            
            df4 = pd.DataFrame.from_dict(guid2pvalue1, orient='index')
            df4.columns=['p_value1']
            df5 = pd.DataFrame.from_dict(guid2pvalue2, orient='index')
            df5.columns=['p_value2']
            df6 = pd.DataFrame.from_dict(guid2pvalue3, orient='index')
            df6.columns=['p_value3']
            df7 = pd.DataFrame.from_dict(guid2pvalue4, orient='index')
            df7.columns=['p_value4']
            df8 = pd.DataFrame.from_dict(guid2observed_p, orient='index')
            df8.columns=['observed_proportion']
            df9 = pd.DataFrame.from_dict(guid2expected_p1, orient='index')
            df9.columns=['expected_proportion1']
            df10 = pd.DataFrame.from_dict(guid2expected_p3, orient='index')
            df10.columns=['expected_proportion2']
            df11 = pd.DataFrame.from_dict(guid2expected_p3, orient='index')
            df11.columns=['expected_proportion3']
            df12 = pd.DataFrame.from_dict(guid2expected_p4, orient='index')
            df12.columns=['expected_proportion4']
            df12['what_tested'] = uncertain_base_type
            
            df = df1.merge(df2n, left_index=True, right_index=True)
            df = df.merge(df3n, left_index=True, right_index=True)
            df = df.merge(df2m, left_index=True, right_index=True)
            df = df.merge(df3m, left_index=True, right_index=True)
            df = df.merge(df2mn, left_index=True, right_index=True)
            df = df.merge(df3mn, left_index=True, right_index=True)
            df = df.merge(df4, left_index=True, right_index=True)
            df = df.merge(df5, left_index=True, right_index=True)
            df = df.merge(df6, left_index=True, right_index=True)
            df = df.merge(df7, left_index=True, right_index=True)
            df = df.merge(df8, left_index=True, right_index=True)
            df = df.merge(df9, left_index=True, right_index=True)         
            df = df.merge(df10, left_index=True, right_index=True)
            df = df.merge(df11, left_index=True, right_index=True)
            df = df.merge(df12, left_index=True, right_index=True)
                  
            retDict = {
                'variant_positions':ordered_variant_positions,
                    'invalid_guids': invalid_guids,
                    'what_tested':[uncertain_base_type]*len(df.index),
                    'guid2sequence':guid2seq,
                    'guid2allN':guid2all['N'],
                    'guid2allM':guid2all['M'],
                    'guid2allN_or_M':guid2all['N_or_M'],
                    'guid2msa_seq':guid2msa_seq,
                    'guid2observed_proportion':guid2observed_p,
                    'guid2expected_p1':guid2expected_p1,
                    'guid2expected_p2':guid2expected_p2,
                    'guid2expected_p3':guid2expected_p3,
                    'guid2expected_p4':guid2expected_p4,
                    'guid2pvalue1':guid2pvalue1,
                    'guid2pvalue2':guid2pvalue2,
                    'guid2pvalue3':guid2pvalue3,
                    'guid2pvalue4':guid2pvalue4,
                    'guid2alignN':guid2align['N'],
                    'guid2alignM':guid2align['M'],
                    'guid2alignN_or_M':guid2align['N_or_M']
                    }
        
        else:
            return None
                
        if output=='dict':    
            return(retDict)
        elif output=='df':
            return(df)
        elif output=='df_dict':
            return(df.to_dict(orient='index'))
        else:
            raise ValueError("Don't know how to format {0}.  Valid options are {'df','df_dict', 'dict'}".format(output))
    def raise_error(self,token):
        """ raises a ZeroDivisionError, with token as the message.
            useful for unit tests of error logging """
        raise ZeroDivisionError(token)
 


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
        inputfile = "../reference/NC_000962.fasta"
        with open(inputfile, 'rt') as f:
            for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):
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
        
