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
                    debugMode=False,
                    excludePositions=set(),
                    snpCompressionCeiling = 250,
                    cpuCount = 1
                ):

        """ instantiates the sequence comparer, an object which manages in-memory reference compressed sequences.
        
        It does not manage persistence, nor does it automatically load sequences.
        
        reference is a string consisting of the reference sequence.
        This is required because, as a data compression technique,
        only differences from the reference are stored.
               
        excludePositions contains a zero indexed set of bases which should not be considered at all in the sequence comparisons.
        Any bases which are always N should be added to this set.
        Not doing so will substantially degrade the algorithm's performance.
        
        If debugMode==True, the server will only load 500 samples.
        
        If the number of Ns are more than maxNs, no data from the sequence is stored.
        If the number of Ns exceeds nCompressionCutoff, Ns are stored not as single base positions but as ranges.  This markedly reduces memory
        usage if the Ns are in long blocks, but slows the comparison rate from about 5,000 per second to about 150/second.
        
        Results > snpCeiling are not returned or stored.
        
        If snpCompressionCeiling is not None, then will consider samples up to snpCompressionCeiling
        when performing deltas compression relative to close neighbours.
        
        cpuCount is ignored in this version.
        
        unknown_base_type is either N or M, and is used for computation of mixture statistics
        David Wyllie, Nov 2018
        
        - to run unit tests, do
        python3 -m unittest seqComparer
        """
        
     
        # we support three kinds of sequences.
        # sequence in strings;
        # reference based compression relative to reference 'compressed_sequence';
        # reference based compression relative to a consensus 'patch_and_consensus'.
        # we detected the latter two by their keys.
        
        self.compressed_sequence_keys = set(['invalid','A','C','G','T', 'N', 'M'])  
        self.patch_and_consensus_keys=set(['consensus_md5','patch'])
        self.patch_keys = set(['+','-', 'M'])
        self.consensus_keys = set(['A','C','G','T', 'N'])
        # store snpCeilings.
        self.snpCeiling = snpCeiling
        self.snpCompressionCeiling = snpCompressionCeiling
              
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

        # prepare to load signatures into memory if directed to do so.
        self.seqProfile={}
        self.consensi = {}      # where consensus sequences are stored in ram
        
        # this version only uses one thread
        self.cpuCount = 1
 
    def raise_error(self,token):
        """ raises a ZeroDivisionError, with token as the message.
            useful for unit tests of error logging """
        raise ZeroDivisionError(token)
    
    def persist(self, object, guid):
        """ keeps a reference compressed object into RAM.
            Note: the sequences are stored on disc/db relative to the reference.
            Compression relative to each other is carried out post-hoc in ram
            """
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
      
    def _refresh(self):
        self._seq1=None
        self._seq2=None
        self.seq1md5=None
        self.seq2md5=None
    
    def mcompare(self, guid, guids=None):
        """ performs comparison of one guid with 
        all guids, which are also stored samples.
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
                (guid1,guid2,dist,n1,n2,nboth, N1pos, N2pos, Nbothpos)=self.countDifferences_byKey(keyPair=(guid,key2),
                                                                                                      cutoff = self.snpCeiling)            
                neighbours.append([guid1,guid2,dist,n1,n2,nboth,N1pos, N2pos, Nbothpos])

       
        return(neighbours)
 
    def distmat(self, half=True, diagonal=False):
        """ returns a distance matrix.
        parameters
        If half=True, returns only half the matrix.
        If diagonal=True, return the diagonal (all of which are zeros)
        
        could refactor to support multithreading.
        
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
                    yield self.countDifferences_byKey(keyPair=(guid1,guid2), cutoff = self.snpCompressionCeiling)            
        
    def summarise_stored_items(self):
        """ counts how many sequences exist of various types """
        retVal = {}
        retVal['server|scstat|nSeqs'] = len(self.seqProfile.keys())
        retVal['server|scstat|nConsensi'] = len(self.consensi.keys())
        retVal['server|scstat|nInvalid'] = 0
        retVal['server|scstat|nCompressed'] =0
        retVal['server|scstat|nRecompressed'] =0
        
        if len(self.seqProfile.keys())==0:
            return(retVal)

        for guid in self.seqProfile.keys():
            if 'invalid' in self.seqProfile[guid]:
                if self.seqProfile[guid]['invalid'] == 1:
                    retVal['server|scstat|nInvalid'] +=1
            if set(self.seqProfile[guid].keys())==self.patch_and_consensus_keys:
                retVal['server|scstat|nRecompressed'] +=1
            else:
                retVal['server|scstat|nCompressed'] +=1
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
    def _guid(self):
        """ returns a new guid, generated de novo """
        return(str(uuid.uuid1()))

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
          
        compressed_sequence = self._computeComparator(compressed_sequence)    # decompress if it is a patch_consensus
        
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
        diffDict={ 'A':set([]),'C':set([]),'T':set([]),'G':set([]),'N':set([]), 'M':{}}        

        for i in self.included:     # for the bases we need to compress
            if not sequence[i]==self.reference[i]:      # if it's not reference
                if sequence[i] in ['A','C','T','G','N']:
                    diffDict[sequence[i]].add(i)        # if it's a definitively called base
                else:
                    # we regard it as a code representing a mixed base.  we store the results in a dictionary
                    diffDict['M'][i] = sequence[i]
                 
        # check how many Ns   
        if len(diffDict['N'])+len(diffDict['M'].keys())>self.maxNs:
            # we store it, but not with sequence details if is invalid
            diffDict={'invalid':1}
        else:
            diffDict['invalid']=0
            
        return(diffDict)
            
    def _computeComparator(self, sequence):
        """ generates a reference compressed version of sequence.
        Acceptable inputs are :
        i) a string containing sequence
        ii) a reference compressed version of the sequence
        iii) a reference compressed version relative to a consensus
        """
                   
        if isinstance(sequence, str):
            return(self.compress(sequence))
        elif isinstance(sequence, dict):
            try:
                if sequence['invalid']==1:
                    raise ValueError("Cannot uncompress an invalid sequence, as it is not stored. {0}".format(sequence.keys()))
            except KeyError:
                pass
            
            if set(sequence.keys())==self.compressed_sequence_keys:
                return(sequence)
            elif set(sequence.keys())==self.patch_and_consensus_keys:
                return(
                    self.apply_patch(sequence['patch'], self.consensi[sequence['consensus_md5']])
                    ) #decompress relative to a patch
            else:
                raise KeyError("Was passed a dictionary with keys {0} but cannot handle this".format(sequence.keys()))
        else:
            raise TypeError("Cannot use object of class {0} as a sequence".format(type(sequence)))

    
    def setComparator1(self,sequence):
        """ stores a reference compressed sequence (no patch) in self._seq1. If the sequence is invalid, stores None"""
        try:
            self._seq1=self._computeComparator(sequence)
        except ValueError:
            # it's invalid
            self._seq1 = None
            
    def setComparator2(self,sequence):
        """ stores a reference compressed sequence (no patch) in self._seq2. If the sequence is invalid, stores None. """
        try:
            self._seq2=self._computeComparator(sequence)
        except ValueError:
            # it's invalid
            self._seq2 = None    
        
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
        # if either sequence is considered invalid (e.g. high Ns) then we report no neighbours.
        self.setComparator1(self.seqProfile[key1])      # need to investigate - is this making large numbers of copies, not just one
        self.setComparator2(self.seqProfile[key2])      # which eat up RAM during (for example) mCompare operations
        nDiff=self.countDifferences(cutoff=cutoff)

        if nDiff is None:
            return((key1, key2, nDiff, None, None, None, None, None, None))
        elif nDiff<=cutoff:
            seq1_uncertain = self._seq1['N'] | set(self._seq1['M'].keys())
            seq2_uncertain = self._seq2['N'] | set(self._seq1['M'].keys())
            (n1, n2, nboth, N1pos, N2pos, Nbothpos) = self._setStats(seq1_uncertain, seq2_uncertain)
            return((key1, key2, nDiff, n1,n2,nboth, N1pos, N2pos, Nbothpos))
        else:
            return((key1, key2, nDiff, None, None, None, None, None, None))

    
    def countDifferences(self,cutoff=None):
        """ compares self._seq1 with self._seq2;
        these are set with self.setComparator1 and 2 respectively.
        Returns the number of SNPs between self._seq1 and self._seq2.
        
        Ns and Ms (uncertain bases) are ignored.
        
        Transparently decompresses any sequences stored as deltas relative to a consensus
        operation rate about 25000 per second."""
        #  if cutoff is not specified, we use snpCeiling
        if cutoff is None:
            cutoff = self.snpCeiling
     
        nDiff=0
        if self._seq1 is None or self._seq2 is None:
            return(None)
                 
        if self._seq1['invalid']==1 or self._seq2['invalid']==1:
            return(None)
         
        # compute positions which differ;
        differing_positions = set()
        for nucleotide in ['C','G','A','T']:
       
            # we do not consider differences relative to the reference if the other nucleotide is an N or M
            nonN_seq1=self._seq1[nucleotide]-(self._seq2['N']|set(self._seq2['M'].keys()))
            nonN_seq2=self._seq2[nucleotide]-(self._seq1['N']|set(self._seq1['M'].keys()))
            differing_positions = differing_positions | (nonN_seq1 ^ nonN_seq2)
        
        nDiff = len(differing_positions)
        
        if nDiff>cutoff:
            return(None)
        else:
            return(nDiff)
    
    def consensus(self, compressed_sequences, cutoff_proportion):
        """ from a list of compressed sequences (as generated by compress())
        generate a consensus consisting of the variation present in at least cutoff_proportion of sequences.
        
        returns the consensus object, which is in the same format as that generated by compress()
        """
        
        # for the compressed sequences in the iterable compressed_sequences, compute a frequency distribution of all variants.
        
        # exclude invalid compressed sequences
        valid_compressed_sequences= []
        for compressed_sequence in compressed_sequences:
            compressed_sequence = self._computeComparator(compressed_sequence)
            if compressed_sequence['invalid']==0:
                valid_compressed_sequences.append(compressed_sequence)
                
        # if there are no valid compressed sequences, we return no consensus.
        if len(valid_compressed_sequences)==0:
            return({'A':set(), 'C':set(),'T':set(), 'G':set(), 'N':set()})

        # otherwise we compute the consensus
        counter = dict()
        for item in ['A','C','T','G','N']:
            if not item in counter.keys():
                counter[item]=dict()
            for result in compressed_sequences:
                if item in result.keys():
                    for position in result[item]:
                        if not position in counter[item].keys():
                            counter[item][position]=0
                        counter[item][position]+=1
        
        # next create a diff object reflecting any variants present in at least cutoff_proportion of the time
        cutoff_number = len(compressed_sequences)*cutoff_proportion
        delta = dict()
        for item in ['A','C','T','G','N']:
            if not item in delta.keys():
                delta[item]=set()
                for position in counter[item]:
                    if counter[item][position] >= cutoff_number:
                        delta[item].add(position)
        return(delta)
    
    def generate_patch(self, compressed_sequence, consensus):
        """ generates a 'patch' or difference between a compressed sequence and a consensus.

        anything which is in consensus and compressed_sequence does not need to be in patch;
        anything which is in consensus and not in  compressed_sequence needs are the 'subtract positions';
        anything which is in compressed_sequence and not consensus in  are the 'add positions'
        
        """
        
        add_positions = {'A':set(),'C':set(),'T':set(),'G':set(),'N':set()}
        subtract_positions = {'A':set(),'C':set(),'T':set(),'G':set(),'N':set()}
        for item in ['A','C','T','G','N']:
            add_positions[item] = compressed_sequence[item]-consensus[item]
            subtract_positions[item]= consensus[item]-compressed_sequence[item]
        for item in ['A','C','T','G','N']:
            if len(add_positions[item])==0:
                del add_positions[item]  # don't store empty sets; they cost ~ 120 bytes each
            if len(subtract_positions[item])==0:
                del subtract_positions[item]  # don't store empty sets; they cost ~ 120 bytes each
          
        retVal = {'+':add_positions, '-':subtract_positions, 'M':compressed_sequence['M']}
        return(retVal)
    
    def apply_patch(self, patch, consensus):
        """ generates a compressed_sequence from a patch and a consensus.
        """
        # sanity check
        if not patch.keys() == self.patch_keys:
            raise TypeError("Patch passed has wrong keys {0}".format(patch.keys))
        if not consensus.keys() == self.consensus_keys:
            raise TypeError("Consensus passed has wrong keys {0}".format(consensus.keys))
        compressed_sequence = {'invalid':0, 'A':set(),'C':set(),'T':set(),'G':set(),'N':set()}
        for item in ['A','C','T','G','N']:
            # empty sets are not stored in a patch
            if item in patch['+']:
                add_these = patch['+'][item]
            else:
                add_these = set()
            if item in patch['-']:
                subtract_these = patch['-'][item]
            else:
                subtract_these = set()
                
            compressed_sequence[item]=  (consensus[item]|add_these)-subtract_these
        compressed_sequence['M']=patch['M']
        return(compressed_sequence)
    
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
    
    def remove_unused_consensi(self):
        """ identifies and removes any consensi which are not used """
        
        # determine all the consensi which are referred to
        used_consensi_md5 = set()
        for guid in self.seqProfile.keys():
            if 'consensus_md5' in self.seqProfile[guid].keys(): 
                used_consensi_md5.add(self.seqProfile[guid]['consensus_md5'])
        initial_consensi = set(self.consensi.keys())
        for consensus_md5 in initial_consensi:
            if not consensus_md5 in used_consensi_md5:
                del self.consensi[consensus_md5]               
    def compress_relative_to_consensus(self, guid, cutoff_proportion=0.8):
        """ identifies sequences similar to the sequence identified by guid.
        Returns any guids which have been compressed as part of the operation"""
        visited_guids = [guid]
        visited_sequences = [self.seqProfile[guid]]
        for compare_with in self.seqProfile.keys():
            result = self.countDifferences_byKey((guid,compare_with),
                                                 self.snpCompressionCeiling)
            # work outward, finding neighbours of seed_sequence up to self.snpCompressionCeiling
            if result[1] is not guid and result[2] is not None:      # we have a close neighbour
               if result[2]<self.snpCompressionCeiling:
                    # we have found something similar, with which we should compress;
                    visited_sequences.append(self.seqProfile[result[1]])
                    visited_guids.append(result[1])
        
        # compute the consensus for these  and store in consensi
        if len(visited_sequences)>1:    # we can compute a consensus
            consensus = self.consensus(visited_sequences, cutoff_proportion)
            consensus_md5 = self.compressed_sequence_hash(consensus)
            self.consensi[consensus_md5]= consensus
            
            # compress the in-memory instances of these samples
            for guid in visited_guids:
                # decompress the in-memory sequence if it is compressed, and re-compress
                this_seqProfile = self.seqProfile[guid]
                self.seqProfile[guid] = {
                    'patch':self.generate_patch(
                            self._computeComparator(this_seqProfile),
                            consensus),
                    'consensus_md5':consensus_md5
                }
            
            # cleanup; remove any consensi which are not needed
            self.remove_unused_consensi()
        
        # return visited_guids
        return(visited_guids)
    
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
            try:
                seq = self._computeComparator(self.seqProfile[guid])
                this_unk = 0
                if unk_type in ['N','N_or_M']:
                    this_unk = this_unk + len(seq['N'])
                if unk_type in ['M','N_or_M']:
                    this_unk = this_unk + len(seq['M'])
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
                seq = self._computeComparator(self.seqProfile[guid])
                this_unk = 0
                if unk_type in ['N','N_or_M']:
                    this_unk = this_unk + len(seq['N'].intersection(sites))
                if unk_type in ['M','N_or_M']:
                    this_unk = this_unk + len(set(seq['M'].keys()).intersection(sites))
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
                comparatorSeq[guid] = self._computeComparator(self.seqProfile[guid])
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
        comparatorSeq={}
        guid2all = {'N':{},'M':{}, 'N_or_M':{}}
        guid2align = {'N':{},'M':{}, 'N_or_M':{}}

        for guid in valid_guids:
            comparatorSeq[guid] = self._computeComparator(self.seqProfile[guid])    
            seq = comparatorSeq[guid]
            for unk_type in ['N','M']:
                guid2all[unk_type][guid] = len(comparatorSeq[guid][unk_type])
            guid2all['N_or_M'][guid] = guid2all['N'][guid]+guid2all['M'][guid]
            
            for base in ['A','C','T','G']:
                positions = seq[base]
                for position in positions:
                  if not position in nrps.keys():     # if it's non-reference, and we've got no record of this position
                     nrps[position]=set()             # then we generate a set of bases at this position                   
                  nrps[position].add(base)            # either way add the current non-reference base there
                  
        # step 2: for the non-reference called positions, check if there's a reference base there.
        for guid in valid_guids:
            seq = comparatorSeq[guid]
            for position in nrps.keys():
                psn_accounted_for  = 0
                for base in ['A','C','T','G']:
                    if position in seq[base]:
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
            seq = comparatorSeq[guid]
            for position in ordered_variant_positions:
                this_base = self.reference[position]
                for base in ['A','C','T','G','N','M']:
                    if not base == 'M':
                        positions = seq[base]
                    else:
                        positions = seq[base].keys()

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
        self.assertEqual(set(res.keys()), set(['server|scstat|nSeqs', 'server|scstat|nConsensi', 'server|scstat|nInvalid', 'server|scstat|nCompressed', 'server|scstat|nRecompressed']))
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
        sc=seqComparer( maxNs = 3,
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

        self.assertEqual(len(df.index),7)
        res= sc.multi_sequence_alignment(guid_names[0:8], output='df_dict', expected_p1=0.995)
        #print(res)
        df = pd.DataFrame.from_dict(res,orient='index')

        self.assertEqual(set(df.columns.values),expected_cols)
    
        self.assertEqual(set(df.index.tolist()), set(['AAACGN-1','CCCCGN-2','TTTCGN-3','GGGGGN-4','ACTCGN-6', 'TCTNGN-7','AAACGN-8']))
        self.assertTrue(df.loc['AAACGN-1','expected_proportion1'] is not None)        # check it computed a value
        self.assertEqual(df.loc['AAACGN-1','expected_proportion1'], 0.995)        # check is used the value passed
        #print(df['p_value1'])
        #print(df['allN'])
        #print(df['alignN'])
class test_seqComparer_47b2(unittest.TestCase):
    """ tests generation of a multisequence alignment with
        testing for the proportion of Ms.
        Tests all three outputs."""
    def runTest(self):
        # generate compressed sequences
        refSeq='GGGGGG'
        sc=seqComparer( maxNs = 3,
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

        self.assertEqual(len(df.index),7)
        res= sc.multi_sequence_alignment(guid_names[0:8], output='df_dict')
        df = pd.DataFrame.from_dict(res,orient='index')

        self.assertTrue(df.loc['AAACGY-1','expected_proportion1'] is not None)        # check it computed a value
        self.assertEqual(set(df.index.tolist()), set(['AAACGY-1','CCCCGY-2','TTTCGY-3','GGGGGY-4','ACTCGY-6', 'TCTQGY-7','AAACGY-8']))

class test_seqComparer_47b(unittest.TestCase):
    """ tests generation of a multisequence alignment with
        testing for the proportion of Ns.
        Tests all three outputs."""
    def runTest(self):
        # generate compressed sequences
        refSeq='GGGGGG'
        sc=seqComparer( maxNs = 3,
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
        self.assertEqual(len(df.index),7)
        res= sc.multi_sequence_alignment(guid_names[0:8], output='df_dict')
        df = pd.DataFrame.from_dict(res,orient='index')
        self.assertTrue(df.loc['AAACGN-1','expected_proportion1'] is not None)        # check it computed a value
        self.assertEqual(set(df.index.tolist()), set(['AAACGN-1','CCCCGN-2','TTTCGN-3','GGGGGN-4','ACTCGN-6', 'TCTNGN-7','AAACGN-8']))

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
        sc=seqComparer( maxNs = 3,
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
        self.assertEqual(len(df.index), 6)
        self.assertEqual(res['variant_positions'],[0,1,2,3])

class test_seqComparer_45c(unittest.TestCase):
    """ tests the generation of multiple alignments of variant sites."""
    def runTest(self):
        
        # generate compressed sequences
        refSeq='GGGGGG'
        sc=seqComparer( maxNs = 3,
                       reference=refSeq,
                       snpCeiling =10)
        
        originals = ['NNNCGN' ]
        guid_names = []
        n=0
        for original in originals:
            n+=1
            c = sc.compress(original)
            this_guid = "{0}-{1}".format(original,n )
            sc.persist(c, guid=this_guid)
            guid_names.append(this_guid)

        res= sc.multi_sequence_alignment(guid_names)
        self.assertTrue(res is None)
   
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
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([]), 'M':{}, 'invalid':0})
class test_seqComparer_3b(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq )
        retVal=sc.compress(sequence='ACTQ')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([]), 'M':{3:'Q'}, 'invalid':0})
class test_seqComparer_3c(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq )
        retVal=sc.compress(sequence='NYTQ')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([0]), 'M':{1:'Y',3:'Q'}, 'invalid':0})
class test_seqComparer_4(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)

        retVal=sc.compress(sequence='ACTN')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([3]),  'M':{}, 'invalid':0})
class test_seqComparer_5(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        retVal=sc.compress(sequence='ACT-')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([3]), 'M':{}, 'invalid':0})         
class test_seqComparer_6(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'

        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)

        retVal=sc.compress(sequence='TCT-')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([0]), 'N': set([3]), 'M':{}, 'invalid':0})
class test_seqComparer_7(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'

        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        retVal=sc.compress(sequence='ATT-')
        self.assertEqual(retVal,{ 'G': set([]), 'A': set([]), 'C': set([]), 'T': set([1]), 'N': set([3]), 'M':{}, 'invalid':0})

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
          
class test_seqComparer_8(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(maxNs = 1e8, snpCeiling = 20,reference=refSeq)

        sc.setComparator1(sequence='ACTG')
        sc.setComparator2(sequence='ACTG')

class test_seqComparer_9(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        sc.setComparator1(sequence='ACTG')
        sc.setComparator2(sequence='ACTG')
        self.assertEqual(sc.countDifferences(),0)
class test_seqComparer_10(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        sc.setComparator1(sequence='TTTG')
        sc.setComparator2(sequence='ACTG')
        self.assertEqual(sc.countDifferences(),2)
class test_seqComparer_11(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        sc.setComparator1(sequence='TTTG')
        sc.setComparator2(sequence='NNTG')
        self.assertEqual(sc.countDifferences(),0)
class test_seqComparer_11b(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        sc.setComparator1(sequence='TTTG')
        sc.setComparator2(sequence='MMTG')
        self.assertEqual(sc.countDifferences(),0)
class test_seqComparer_12(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)

        sc.setComparator2(sequence='TTTG')
        sc.setComparator1(sequence='NNTG')
        self.assertEqual(sc.countDifferences(),0)
class test_seqComparer_13(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        sc.setComparator2(sequence='TTTG')
        sc.setComparator1(sequence='--TG')
        self.assertEqual(sc.countDifferences(),0)
class test_seqComparer_13b(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        sc.setComparator2(sequence='TTTG')
        sc.setComparator1(sequence='RRRR')
        self.assertEqual(sc.countDifferences(),0)
class test_seqComparer_13c(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        sc.setComparator1(sequence='TTTG')
        sc.setComparator2(sequence='RRRR')
        self.assertEqual(sc.countDifferences(),0)
class test_seqComparer_14(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        sc.setComparator2(sequence='TTAA')
        sc.setComparator1(sequence='--AG')
        self.assertEqual(sc.countDifferences(),1)
class test_seqComparer_15(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        sc.setComparator1(sequence='TTAA')
        sc.setComparator2(sequence='--AG')
        self.assertEqual(sc.countDifferences(),1)
        
class test_seqComparer_16(unittest.TestCase):
    """ tests the comparison of two sequences where both differ from the reference. """
    def runTest(self):   
        # generate compressed sequences
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        
        sc._seq1 = sc.compress('AAAA')
        sc._seq2 = sc.compress('CCCC')
        self.assertEqual(sc.countDifferences(),4)
class test_seqComparer_16b(unittest.TestCase):
    """ tests the comparison of two sequences where both differ from the reference. """
    def runTest(self):   
        # generate compressed sequences
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        
        sc._seq1 = sc.compress('AAAA')
        sc._seq2 = sc.compress('RRCC')
        self.assertEqual(sc.countDifferences(),2)
class test_seqComparer_16c(unittest.TestCase):
    """ tests the comparison of two sequences where both differ from the reference. """
    def runTest(self):   
        # generate compressed sequences
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        
        sc._seq1 = sc.compress('AAAA')
        sc._seq2 = sc.compress('RRNN')
        self.assertEqual(sc.countDifferences(),0)
class test_seqComparer_17(unittest.TestCase):
    """ tests the comparison of two sequences where one is invalid """
    def runTest(self):   
        # generate compressed sequences
        refSeq='ACTG'
        sc=seqComparer( maxNs = 3,
                       reference=refSeq,
                       snpCeiling =10)
        
        sc._seq1 = sc.compress('AAAA')
        sc._seq2 = sc.compress('NNNN')
        self.assertEqual(sc.countDifferences(),None)
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
        self.assertEqual(retVal,{ 'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'M':{}, 'N': set([8,9,10,11,12,13,14,15]), 'invalid':0})
        retVal=sc.compress(sequence='NNTGTTAANNNNNNNNTGGGGGGGGGGGGAA')
        self.assertEqual(retVal,{ 'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'M':{}, 'N': set([0,1,8,9,10,11,12,13,14,15]), 'invalid':0})
       
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
 
class test_seqComparer_30(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, reference=refSeq, snpCeiling= 1)
        sc.setComparator1(sequence='ACTG')
        sc.setComparator2(sequence='ACTG')
class test_seqComparer_31(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, reference=refSeq, snpCeiling =1 )
        sc.setComparator1(sequence='ACTG')
        sc.setComparator2(sequence='ACTG')
        self.assertEqual(sc.countDifferences(),0)
class test_seqComparer_32(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, reference=refSeq, snpCeiling =1)
        sc.setComparator1(sequence='TTTG')
        sc.setComparator2(sequence='ACTG')
        self.assertEqual(sc.countDifferences(),None)
class test_seqComparer_33(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, reference=refSeq, snpCeiling =1)
        sc.setComparator1(sequence='TTTG')
        sc.setComparator2(sequence='NNTG')
        self.assertEqual(sc.countDifferences(),0)
class test_seqComparer_34(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, reference=refSeq, snpCeiling =1)
        sc.setComparator2(sequence='TTTG')
        sc.setComparator1(sequence='NNTG')
        self.assertEqual(sc.countDifferences(),0)
class test_seqComparer_35(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 2, reference=refSeq, snpCeiling =1)
        sc.setComparator2(sequence='TTTG')
        sc.setComparator1(sequence='NNNG')
        self.assertEqual(sc.countDifferences(),0)
class test_seqComparer_13(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, reference=refSeq, snpCeiling =1)
        sc.setComparator2(sequence='TTTG')
        sc.setComparator1(sequence='--TG')
        self.assertEqual(sc.countDifferences(),0)
class test_seqComparer_35(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, reference=refSeq, snpCeiling =1)
        sc.setComparator2(sequence='TTAA')
        sc.setComparator1(sequence='--AG')
        self.assertEqual(sc.countDifferences(),1)
class test_seqComparer_36(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, reference=refSeq, snpCeiling =1)
        sc.setComparator1(sequence='TTAA')
        sc.setComparator2(sequence='--AG')
        self.assertEqual(sc.countDifferences(),1)
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


class test_seqComparer_39a(unittest.TestCase):
    """ tests the computation of a consensus sequence """
    def runTest(self):
        
        # generate compressed sequences
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, reference=refSeq, snpCeiling =10)
        compressed_sequences = []
        compressed_sequences.append(sc.compress(sequence='TTAA'))
        compressed_sequences.append(sc.compress(sequence='TTTA'))
        compressed_sequences.append(sc.compress(sequence='TTGA'))
        compressed_sequences.append(sc.compress(sequence='TTAA'))

        cutoff_proportion = 0.5
        consensus = sc.consensus(compressed_sequences, cutoff_proportion)

        expected_consensus = { 'T': {0, 1}, 'N': set(), 'A': {2, 3}, 'C': set(), 'G': set()}       
        self.assertEqual(consensus, expected_consensus)


class test_seqComparer_39b(unittest.TestCase):
    """ tests the computation of a consensus sequence """
    def runTest(self):
        
        # generate compressed sequences
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, reference=refSeq, snpCeiling =10)
        compressed_sequences = []
        
        cutoff_proportion = 0.5
        delta = sc.consensus(compressed_sequences, cutoff_proportion)

        expected_delta = {'A':set(), 'C':set(),'T':set(), 'G':set(), 'N':set()} 
        
        self.assertEqual(delta, expected_delta)
          
class test_seqComparer_39c(unittest.TestCase):
    """ tests the computation of a consensus sequence """
    def runTest(self):
        
        # generate compressed sequences
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, reference=refSeq, snpCeiling =10)
        compressed_sequences = []
        compressed_sequences.append(sc.compress(sequence='TTAA'))
        compressed_sequences.append(sc.compress(sequence='TTTA'))
        compressed_sequences.append(sc.compress(sequence='TTGA'))
        compressed_sequences.append(sc.compress(sequence='TTAN'))

        cutoff_proportion = 0.5
        consensus = sc.consensus(compressed_sequences, cutoff_proportion)

        expected_consensus = { 'T': {0, 1}, 'N': set(), 'A': {2, 3}, 'C': set(), 'G': set()}       
        self.assertEqual(consensus, expected_consensus)


class test_seqComparer_39d(unittest.TestCase):
    """ tests the computation of a consensus sequence """
    def runTest(self):
        
        # generate compressed sequences
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, reference=refSeq, snpCeiling =10)
        compressed_sequences = []
        compressed_sequences.append(sc.compress(sequence='TTAA'))
        compressed_sequences.append(sc.compress(sequence='TTTA'))
        compressed_sequences.append(sc.compress(sequence='TTGM'))
        compressed_sequences.append(sc.compress(sequence='TTAN'))

        cutoff_proportion = 0.5
        consensus = sc.consensus(compressed_sequences, cutoff_proportion)

        expected_consensus = { 'T': {0, 1}, 'N': set(), 'A': {2, 3}, 'C': set(), 'G': set()}       
        self.assertEqual(consensus, expected_consensus)
   
class test_seqComparer_40(unittest.TestCase):
    """ tests the computation of a hash of a compressed object """
    def runTest(self):

        # generate compressed sequences
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, reference=refSeq, snpCeiling =10)
        compressed_sequence = sc.compress(sequence='TTAA')

        res = sc.compressed_sequence_hash(compressed_sequence)
        self.assertEqual(res, "6ce0e55c4ab092f560e03c5d2de53098")
        
        
class test_seqComparer_41(unittest.TestCase):
    """ tests the computation of a difference relative to a reference + delta """
    def runTest(self):
        
        # generate compressed sequences
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, reference=refSeq, snpCeiling =10)
        compressed_sequences = []
        compressed_sequences.append(sc.compress(sequence='TTAA'))
        compressed_sequences.append(sc.compress(sequence='TTTA'))
        compressed_sequences.append(sc.compress(sequence='TTGA'))
        compressed_sequences.append(sc.compress(sequence='TTAA'))

        cutoff_proportion = 0.5
        consensus = sc.consensus(compressed_sequences, cutoff_proportion)

        originals = [ 'AAAA','CCCC','TTTT','GGGG','NNNN','ACTG','ACTC', 'TCTN', 'MCTN', 'MMTT', 'QQQQ']
        for original in originals:

            compressed_sequence = sc.compress(sequence=original)
            patch = sc.generate_patch(compressed_sequence, consensus)
            roundtrip_compressed_sequence = sc.apply_patch(patch, consensus)
            self.assertEqual(compressed_sequence, roundtrip_compressed_sequence)
            roundtrip = sc.uncompress(roundtrip_compressed_sequence)
            self.assertEqual(roundtrip, original)

class test_seqComparer_42(unittest.TestCase):
    """ tests the compression relative to a consensus """
    def runTest(self):
        
        # generate compressed sequences
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        
        originals = [ 'AAAC','CCCC','TTTC','GGGC','NNNC','ACTC','ACTC', 'TCTN' ]
        for original in originals:   
            c = sc.compress(original)
            sc.persist(c, guid=original )

        sc.compress_relative_to_consensus(guid = 'AAAC', cutoff_proportion = 0.5)
        
        for original in originals:
            self.assertEqual(original, sc.uncompress(sc.seqProfile[original]))

class test_seqComparer_43(unittest.TestCase):
    """ tests the compression relative to a consensus with a consensus present"""
    def runTest(self):
        
        # generate compressed sequences
        refSeq='GGGGGGGGGGGG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        
        originals = [ 'AAACACTGACTG','CCCCACTGACTG','TTTCACTGACTG','GGGCACTGACTG','NNNCACTGACTG','ACTCACTGACTG','ACTCACTGACTG', 'TCTNACTGACTG' ]
        for original in originals:   
            c = sc.compress(original)
            sc.persist(c, guid=original )

        sc.compress_relative_to_consensus(guid = 'AAACACTGACTG', cutoff_proportion = 0.8)
        
        for original in originals:
            self.assertEqual(original, sc.uncompress(sc.seqProfile[original]))

class test_seqComparer_44(unittest.TestCase):
    """ tests the compression relative to a consensus with a consensus present"""
    def runTest(self):
        
        # generate compressed sequences
        refSeq='GGGG'
        sc=seqComparer( maxNs = 2,
                       reference=refSeq,
                       snpCeiling =10)

        with self.assertRaises(KeyError):
            res = sc.countDifferences_byKey(('AAAC','NNNC'))
  
        originals = [ 'AAAC', 'NNNC' ]
        for original in originals:   
            c = sc.compress(original)
            sc.persist(c, guid=original )


        # use a Tuple
        res = sc.countDifferences_byKey(('AAAC','NNNC'))
        self.assertEqual(res[2], None)
        
        refSeq='GGGG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        originals = [ 'AAAC', 'NNNC' ]
        for original in originals:   
            c = sc.compress(original)
            sc.persist(c, guid=original )

        res = sc.countDifferences_byKey(('AAAC','NNNC'))
        self.assertEqual(res[2], 0)

   
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
            if i % 5 == 0:
                sc.compress_relative_to_consensus(guid_to_insert)

class test_seqComparer_46(unittest.TestCase):
    """ tests the compression relative to a consensus with a consensus present.
    then adds more sequences, changing the consensus."""
    def runTest(self):
        
        # generate compressed sequences
        refSeq='GGGGGGGGGGGG'
        sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       snpCeiling =10)
        
        originals = [ 'AAACACTGACTG','CCCCACTGACTG','TTTCACTGACTG','GGGCACTGACTG','NNNCACTGACTG','ACTCACTGACTG','ACTCACTGACTG', 'TCTNACTGACTG' ]
        for original in originals:   
            c = sc.compress(original)
            sc.persist(c, guid=original )

        sc.compress_relative_to_consensus(guid = 'AAACACTGACTG', cutoff_proportion = 0.8)
        initial_consensi_keys = set(sc.consensi.keys())
        for original in originals:
            self.assertEqual(original, sc.uncompress(sc.seqProfile[original]))

        # add more changing the consensus by adding at T
        more_seqs = [ 'QQACACTGACTG','TAACACTGACTG','TCCCACTGACTG','TTTCACTGACTG','TGGCACTGACTG','TNNCACTGACTG','TCTCACTGACTG','TCTCACTGACTG', 'TCTNACTGACTG' ]
        for more_seq in more_seqs:   
            c = sc.compress(more_seq)
            sc.persist(c, guid=more_seq )

        sc.compress_relative_to_consensus(guid = 'AAACACTGACTG', cutoff_proportion = 0.8)
  
        for original in originals:
            self.assertEqual(original, sc.uncompress(sc.seqProfile[original]))
        for more_seq in more_seqs:
            self.assertEqual(more_seq, sc.uncompress(sc.seqProfile[more_seq]))

        # there should be one consensus
        later_consensi_keys = set(sc.consensi.keys())
        self.assertNotEqual(initial_consensi_keys, later_consensi_keys)
        self.assertEqual(len(later_consensi_keys), 1)
 
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
        