#!/usr/bin/env python3

# python code to compare fasta sequences
import unittest
import os
import glob
import sys
from Bio import SeqIO
import datetime
import pickle
import hashlib
import collections
import math
import multiprocessing
import uuid
import json
import psutil
from gzip import GzipFile
import random
import itertools


class seqComparer():
    def __init__(self,
                    reference,
                    maxNs,
                    snpCeiling,
                    debugMode=False,
                    excludePositions=set(),
                    snpCompressionCeiling = 250
                ):

        """ instantiates the sequence comparer, an object whcih manages in-memory reference compressed sequences.
        
        It does not manage persistence, nor does it automatically load sequences.
        
        reference is a string consisting of the reference sequence.
        This is required because, as a data compression technique,
        only differences from the reference are stored.
               
        excludePositions contains a set of bases which should not be considered at all in the sequence comparisons.
        Any bases which are always N should be added to this set.
        Not doing so will substantially degrade the algorithm's performance.
        
        If debugMode==True, the server will only load 500 samples.
        
        If the number of Ns are more than maxNs, no data from the sequence is stored.
        If the number of Ns exceeds nCompressionCutoff, Ns are stored not as single base positions but as ranges.  This markedly reduces memory
        usage if the Ns are in long blocks, but slows the comparison rate from about 5,000 per second to about 150/second.
        
        Results > snpCeiling are not returned or stored.
        
        If snpCompressionCeiling is not None, then will consider samples up to snpCompressionCeiling
        when performing deltas compression relative to close neighbours.
        
        David Wyllie, University of Oxford, June 2018
        
        - to run unit tests, do
        python -m unittest seqComparer
        """
        
        # we support three kinds of sequences.
        # sequence in strings;
        # reference based compression relative to reference 'compressed_sequence';
        # reference based compression relative to a consensus 'patch_and_consensus'.
        # we detected the latter two by their keys.
        
        self.compressed_sequence_keys = set(['invalid','A','C','G','T', 'N'])
        self.patch_and_consensus_keys=set(['consensus_md5','patch'])
        self.patch_keys = set(['add','subtract'])
        self.consensus_keys = set(['A','C','G','T', 'N'])
        
        # store snpCeilings.
        self.snpCeiling = snpCeiling
        self.snpCompressionCeiling = snpCompressionCeiling
              
        # sequences with more than maxNs Ns will be considered invalid and their details (apart from their invalidity) will not be stored.
        self.maxNs=maxNs

        # check composition of the reference.
        self.reference=str(reference)           # if passed a Bio.Seq object, coerce to string.
        letters=collections.Counter(self.reference)   
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
      
        # initiate object used to track memory
        self.memtracker = psutil.Process(os.getpid())

    def persist(self, object, guid):
        """ keeps a reference compressed object into RAM.
            Note: the sequences are stored on disc/db relative to the reference.
            Compression relative to each other is carried out post-hoc in ram
            """
        self.seqProfile[guid]=object
    def load(self, guid):
        """ recovers (loads) a variable containing a reference compressed object into RAM.
            Note: the sequences are stored on disc/db relative to the reference.
            Compression relative to each other is carried out post-hoc in ram
            """
        return self.seqProfile[guid]
      
    def _refresh(self):
        self._seq1=None
        self._seq2=None
        self.seq1md5=None
        self.seq2md5=None
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
    def memoryInfo(self):
        """ returns information about memory usage by the machine """ 
        return(self.memtracker.memory_info())
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
        compressed_sequence = self._computeComparator(compressed_sequence)    # decompress if it is a patch_consensus 
        if compressed_sequence['invalid']==1:
            raise ValueError("Cannot uncompress an invalid sequence, as it is not stored")
        
        seq = list(self.reference)
        for item in ['A','C','T','G','N']:
            for x in compressed_sequence[item]:
                seq[x]=item
        return(''.join(seq))
    
    def compress(self, sequence):
        """ reads a string sequence and extracts position - genome information from it.
        returns a dictionary consisting of zero-indexed positions of non-reference bases.
        
        """
        if not len(sequence)==len(self.reference):
            raise TypeError("sequence must of the same length as reference")
        
        # we consider - characters to be the same as N
        sequence=sequence.replace('-','N')
        
        # we only record differences relative to to refSeq.
        # anything the same as the refSeq is not recorded.
        diffDict={ 'A':set([]),'C':set([]),'T':set([]),'G':set([]),'N':set([])}        

        for i in self.included:     # for the bases we need to compress

            if not sequence[i]==self.reference[i]:
                diffDict[sequence[i]].add(i)
                 
        # convert lists to sets (faster to do this all at once)

        for key in ['A','C','G','T']:
            diffDict[key]=set(diffDict[key])
            
        if len(diffDict['N'])>self.maxNs:
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
        """ stores a reference compressed sequence (no patch) in self._seq1 """
        self._seq1=self._computeComparator(sequence)
           
    def setComparator2(self,sequence):
        """ stores a reference compressed sequence (no patch) in self._seq2 """
        self._seq2=self._computeComparator(sequence)
        
        
    def _setStats(self, sortedset1, sortedset2):
        """ compares sortedset1, which contains a series of ranges {(0,1) (10,11)}
        with         sortedset2, which also contains  ranges       {(1,2)}
        
        returns *
        * the number of elements in sortedset1  4
        * the number of elements in sortedset2  2
        * the number of elements in the union of sortedset1 and sortedset2 5
        * sortedSet1 {0,1}
        * sortedSet2 {10,11}
        * the union of sorted set1 and sortedset2   {0,1,2,10,11)}
        
        """
        
        if type(sortedset1)==set and type(sortedset2)==set:
            # then we can just use a standard set union operation.
            retVal=sortedset1 | sortedset2
            return(len(sortedset1), len(sortedset2), len(retVal), sortedset1, sortedset2, retVal)
        
        retVal1=set()        
        if type(sortedset1)==set:
                retVal1=sortedset1
        else:
            raise TypeError("Cannot merge something of class {0}".type(sortedset1))
        
        retVal2=set()
        if type(sortedset2)==set:
                retVal2= sortedset2
        else:
            raise TypeError("Cannot merge something of class {0}".type(sortedset2))

        retVal=retVal2 | retVal1

        return(len(retVal1), len(retVal2), len(retVal), retVal1, retVal2, retVal)   

    def countDifferences_byKey(self, keyPair, cutoff=None):
        """ compares the in memory refCompressed sequences at
        self.seqProfile[key1] and self.seqProfile[key2]

        Returns the number of SNPs between self._seq1 and self._seq2, and, if it is less than cutoff,
        the number of Ns in the two sequences and the union of their positions.
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
        self.setComparator1(self.seqProfile[key1])
        self.setComparator2(self.seqProfile[key2])
        nDiff=self.countDifferences(cutoff=cutoff)
        
        if nDiff is None:
            return((key1, key2, nDiff, None, None, None, None, None, None))
        elif nDiff<=cutoff:
            (n1, n2, nboth, N1pos, N2pos, Nbothpos) = self._setStats(self._seq1['N'],self._seq2['N'])
            return((key1, key2, nDiff, n1,n2,nboth, N1pos, N2pos, Nbothpos))
        else:
            return((key1, key2, nDiff, None, None, None, None, None, None))

    
    def countDifferences(self,cutoff=None):
        """ compares self._seq1 with self._seq2;
        these are set with self.setComparator1 and 2 respectively.
        Returns the number of SNPs between self._seq1 and self._seq2.
        
        Transparently decompresses any sequences stored as deltas relative to a consensus
        scan rate about 25000 per second."""
        #  if cutoff is not specified, we use snpCeiling
        if cutoff is None:
            cutoff = self.snpCeiling
     
        nDiff=0
        if self._seq1['invalid']==1 or self._seq2['invalid']==1:
            return(None)
         
        # compute positions which differ;
        differing_positions = set()
        for nucleotide in ['C','G','A','T']:
       
            # we do not consider differences relative to the reference if the other nucleotide is an N
            nonN_seq1=self._seq1[nucleotide]-self._seq2['N']
            nonN_seq2=self._seq2[nucleotide]-self._seq1['N']
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
          
        retVal = {'add':add_positions, 'subtract':subtract_positions}
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
            if item in patch['add']:
                add_these = patch['add'][item]
            else:
                add_these = set()
            if item in patch['subtract']:
                subtract_these = patch['subtract'][item]
            else:
                subtract_these = set()
                
            compressed_sequence[item]=  (consensus[item]|add_these)-subtract_these
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
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([]), 'invalid':0})
class test_seqComparer_4(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'

        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)

        retVal=sc.compress(sequence='ACTN')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([3]), 'invalid':0})
class test_seqComparer_5(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        retVal=sc.compress(sequence='ACT-')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([3]), 'invalid':0})         
class test_seqComparer_6(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'

        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)

        retVal=sc.compress(sequence='TCT-')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([0]), 'N': set([3]), 'invalid':0})
class test_seqComparer_7(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'

        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        retVal=sc.compress(sequence='ATT-')
        self.assertEqual(retVal,{ 'G': set([]), 'A': set([]), 'C': set([]), 'T': set([1]), 'N': set([3]), 'invalid':0})

class test_seqComparer_6b(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'

        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        originals = [ 'AAAA','CCCC','TTTT','GGGG','NNNN','ACTG','ACTC', 'TCTN']
        for original in originals:

            compressed_sequence=sc.compress(sequence=original)
          
            roundtrip = sc.uncompress(compressed_sequence)
            self.assertEqual(original, roundtrip)
            
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

class test_seqComparer_saveload3(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, snpCeiling = 20,reference=refSeq)
        compressedObj =sc.compress(sequence='ACTT')
        sc.persist(compressedObj, 'one' )     
        retVal=sc.load(guid='one' )
        self.assertEqual(compressedObj,retVal)        

class test_seqComparer_24(unittest.TestCase):
    """ tests N compression """
    def runTest(self):
        
        refSeq=                     'ACTGTTAATTTTTTTTTGGGGGGGGGGGGAA'
        sc=seqComparer(maxNs = 1e8, snpCeiling = 20,reference=refSeq)

        retVal=sc.compress(sequence='ACTGTTAANNNNNNNNTGGGGGGGGGGGGAA')
        self.assertEqual(retVal,{ 'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([8,9,10,11,12,13,14,15]), 'invalid':0})
        retVal=sc.compress(sequence='NNTGTTAANNNNNNNNTGGGGGGGGGGGGAA')
        self.assertEqual(retVal,{ 'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([0,1,8,9,10,11,12,13,14,15]), 'invalid':0})
       
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
          
  
class test_seqComparer_40(unittest.TestCase):
    """ tests the computation of a hash of a compressed object """
    def runTest(self):

        # generate compressed sequences
        refSeq='ACTG'
        sc=seqComparer( maxNs = 1e8, reference=refSeq, snpCeiling =10)
        compressed_sequence = sc.compress(sequence='TTAA')

        res = sc.compressed_sequence_hash(compressed_sequence)
        self.assertEqual(res, "23b867b142bad108b848b87ad4b79633")
        
        
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

        originals = [ 'AAAA','CCCC','TTTT','GGGG','NNNN','ACTG','ACTC', 'TCTN' ]
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
        more_seqs = [ 'TAACACTGACTG','TCCCACTGACTG','TTTCACTGACTG','TGGCACTGACTG','TNNCACTGACTG','TCTCACTGACTG','TCTCACTGACTG', 'TCTNACTGACTG' ]
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
 