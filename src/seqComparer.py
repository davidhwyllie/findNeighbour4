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
from banyan import SortedSet, OverlappingIntervalsUpdator


class seqComparer():
    def __init__(self,
                    reference,
                    NCompressionCutoff,
                    maxNs,
                    snpCeiling,
                    persistenceDir=None,
                    persistenceStore=['localmemory'],
                    startAfresh=False,
                    debugMode=False,
                    excludeFile=os.path.join("..","reference","TB-exclude.txt")
                ):

        """ instantiates the sequence comparer.
        
        reference is a string consisting of the reference sequence.
        This is required because, as a data compression technique,
        only differences from the reference are stored.
        
        persistenceStore is one or more of
            'localmemory' (stores on disc as backup, and in local memory);
            'tofile' - to the persistenceDir;
           
        persistenceDir is a directory into which the parsed sequences will be written; required if persistenceStore is either 'localmemory' or 'tofile'.         
        startAfresh=True causes the persistence store to be emptied, which is useful for benchmarking.
        
        excludeFile contains a list of bases which should not be considered at all in the sequence comparisons.  Any bases which are always N should be added to this list.
        Not doing so will substantially degrade the algorithm's performance.
        
        If debugMode==True, the server will only load 500 samples.
        

        If the number of Ns are more than maxNs, no data from the sequence is stored.
        If the number of Ns exceeds nCompressionCutoff, Ns are stored not as single base positions but as ranges.  This markedly reduces memory
        usage if the Ns are in long blocks, but slows the comparison rate from about 5,000 per second to about 150/second.
        
        Results > snpCeiling are not returned or stored.
        
        David Wyllie, University of Oxford, Jan 2017
        
        - to run unit tests, do
        python -m unittest seqComparer
        """
        # store snpCeiling.
        self.snpCeiling = snpCeiling

        # check the nature of persistence store; if it is a single item, store that as a list. 
        if type(persistenceStore) is str:
            persistenceStore=[persistenceStore]     # make it a list
        self.persistenceStore=persistenceStore      # what is used to store sequences
        
        # sequences with more than maxNs Ns will be considered invalid and their details (apart from their invalidity) will not be stored.
        self.maxNs=maxNs

        self.NCompressionCutoff=NCompressionCutoff

        # check composition of the reference.
        self.reference=str(reference)           # if passed a Bio.Seq object, coerce to string.
        letters=collections.Counter(self.reference)   
        if not set(letters.keys())==set(['A','C','G','T']):
            raise TypeError("Reference sequence supplied contains characters other than ACTG: {0}".format(letters))
        
        ## verify that the storage systems stated are ready for use.
        # check the existence of the persistenceDir and make it if it does not exist.
        self.persistenceDir=persistenceDir
        if 'localmemory' in self.persistenceStore or 'tofile' in self.persistenceStore:     
            if persistenceDir is not None:        
                self.persistenceDir=os.path.join(persistenceDir)
                if not os.path.exists(self.persistenceDir):
                        os.makedirs(self.persistenceDir)
                
        # load the excluded bases
        self.excluded=set()
        if excludeFile is not None:
            with open(excludeFile,'rt') as f:
                rows=f.readlines()
            for row in rows:
                self.excluded.add(int(row))
            print("Excluded {0} positions.".format(len(self.excluded)))
        
        # define what is included
        self.included=set(range(len(self.reference)))-self.excluded
        
        # initialise pairwise sequences for comparison.
        self._refresh()

        # empty out stores if appropriate (e.g. for benchmarking)
        for this_persistenceStore in self.persistenceStore:
            if startAfresh==True:
                if this_persistenceStore=='localmemory':
                        self.seqProfile={}
                        self.emptyPersistenceDir()
                elif this_persistenceStore=='tofile':
                        self.emptyPersistenceDir()
                else:
                        raise TypeError("Do not know how to use {0} as a persistence mechanism.".format(this_persistenceStore))
                  
        # load the signatures into memory if directed to do so.
        self.seqProfile={}
        if 'localmemory' in self.persistenceStore:
            if self.persistenceDir is not None:
                print("Loading profiles of existing samples")
                for (i,filepath) in enumerate(glob.glob(os.path.join(self.persistenceDir,'*.pickle'))):       # DirProfile if just profiles
                    guid=os.path.basename(filepath)

                    if (i % 500 ==0):
                        print("Starting up; Loaded {0}; please wait".format(i))

                    if debugMode==True and i>500:
                        break
                    with open(filepath,'rb') as f:
                        # remove the .pickle suffix
                        guid=os.path.basename(filepath).replace('.pickle','')
                        self.seqProfile[guid]=pickle.load(f)
        
        # initiate object used to track memory
        self.memtracker = psutil.Process(os.getpid())

    def _refresh(self):
        self.seq1=None
        self.seq2=None
        self.seq1md5=None
        self.seq2md5=None
    def emptyPersistenceDir(self):
        """ deletes all *.pickle files in the persistence directory """
        if not self.persistenceDir is None:
            i=0
            for (i,pickleFile) in enumerate(glob.glob(os.path.join(self.persistenceDir,'*.pickle'))):
                os.unlink(pickleFile)
            print("Emptied persistence directory, removing {0} files.".format(i))
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
    def iscachedtofile(self,guid):
        """ returns true or false depending whether we have a compressed copy on disc (name=guid) in this machine """
        expectedFilename=os.path.join(self.persistenceDir, "{0}.pickle".format(guid))
        return(os.path.exists(expectedFilename))
    def guidscachedtofile(self):
        """ returns all guids with sequence profiles currently in this machine """
        guids=set()
        cachedfiles=glob.glob(os.path.join(self.persistenceDir,'*.pickle'))
        for cachedfile in cachedfiles:
            guid=os.path.basename(cachedfile).replace('.pickle','')
            guids.add(guid)

        return(guids)          

    def _guid(self):
        """ returns a guid """
        return(str(uuid.uuid1()))
    def memoryInfo(self):
        """ returns information about memory usage by the machine """ 
        return(self.memtracker.memory_info())
    def _delta(self,x):
        """ returns the difference between two numbers in a tuple x """
        return(x[1]-x[0])
    def _ranges(self, i):
        """ converts a set or list of integers, i,  into a list of tuples defining the start and stop positions
        of continuous blocks of integers.
        
        e.g. [0,1,2,3,4,6,7]  --> [(0,4),(6,7)]
        
        This is derived from the approach here:
        http://stackoverflow.com/questions/4628333/converting-a-list-of-integers-into-range-in-python
        
        """
        i=sorted(i)
        tups=[t for t in enumerate(i)]

        for a,b in itertools.groupby( tups, self._delta):
            b=list(b)
            yield(b[0][1], b[-1][1])

    def compress(self, sequence):
        """ reads a string sequence and extracts position - genome information from it.
        returns a dictionary consisting of zero-indexed positions of non-reference bases. """
        if not len(sequence)==len(self.reference):
            raise TypeError("sequence must of the same length as reference")
        
        # we consider - characters to be the same as N
        sequence=sequence.replace('-','N')
        
        # we only record differences relative to to refSeq.
        # anything the same as the refSeq is not recorded.
        diffDict={'A':set([]),'C':set([]),'T':set([]),'G':set([]),'N':set([])}        

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
            # for Ns, we code 'runs' of Ns as ranges;
            nRuns=list(self._ranges(diffDict['N']))
            
            # If there are large number of Ns, we can optionally index this into
            # a Banyan sorted set (~ in-memory B-tree) for high-speed N exclusion during comparisons;
            # the optimal cutoff point at which it becomes worthwhile doing this needs to be tuned
            if len(diffDict['N'])>self.NCompressionCutoff:
                diffDict['N'] = SortedSet(list(self._ranges(diffDict['N'])), key_type = (int,int), updator=OverlappingIntervalsUpdator)


        return(diffDict)         
    def persist_tofile(self, refCompressedSequence, guid):
        """ writes the reference encoded complete difference object refCompressedSequence to file in self.persistenceDir

        """
        if self.persistenceDir is None:
            raise NotImplementedError("Cannot persist picked object if persistenceDir is not set.")
        
        outputFile=os.path.join(self.persistenceDir, "{0}.pickle".format(guid))       
        if not os.path.exists(outputFile):
            with open(outputFile,'wb') as f:
                pickle.dump(refCompressedSequence, file=f, protocol=2)

        return(guid)
    def persist(self, refCompressedSequence, method, guid=None, **kwargs):
        """ persists the refCompressedSequence to the location chosen, as stated in 'method'.
        Options are: 'tofile','localmemory','toredis_sets','toredis_pickle'. multiple options are possible"""
        ## check that sequence is indeed a dictionary 

        if refCompressedSequence is None:
            raise Typeerror("refCompressed sequence cannot be None")
        if not type(refCompressedSequence)==dict:
            raise TypeError(".persist method needs to be passed a dictionary as generated by .compress() method.")


        if type(method)==str:
            # then we convert to a list;
            method=[method]
            
        if guid is None:
            # assign
            guid=self._guid()

        for this_method in method:
            if this_method=='tofile':
                self.persist_tofile(refCompressedSequence=refCompressedSequence, guid=guid)
            
            elif this_method=='localmemory':
                # store it on disc and in memory
                self.seqProfile[guid]=refCompressedSequence
                self.persist_tofile(refCompressedSequence=refCompressedSequence, guid=guid)
            
            else:
                raise TypeError("Do not know how to persist to '{0}'".format(this_method))       
        return(guid)
    def load_fromfile(self, guid):
        """ loads a picked object from file """
        if self.persistenceDir is not None:
                    filepath=os.path.join(self.persistenceDir,"{0}.pickle".format(guid))
                    with open(filepath,'rb') as f:
                        return(pickle.load(f))      # will raise an error if it does not exist; could trap, but what would we raise?
        else:
            raise TypeError("Cannot recovery from file when presistenceDir is not set")
    def load(self,guid, method='localmemory'):
        """ recovers the full reference- compressed object denoted by the guid on the sequence"""
        if method=='tofile':
            return(self.load_fromfile(guid))
        elif method=='localmemory':
            return(self.seqProfile[guid])
    def setComparator1(self,sequence):
        self.seq1=self.compress(sequence)
    def setComparator2(self,sequence):
        self.seq2=self.compress(sequence)
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
        if type(sortedset1)==SortedSet:
            for (start,stop) in sortedset1:
                for x in range(start,stop+1):
                    retVal1.add(x)
        elif type(sortedset1)==set:
                retVal1=sortedset1
        else:
            raise TypeError("Cannot merge something of class {0}".type(sortedset1))
        
        retVal2=set()
        if type(sortedset2)==SortedSet:
            for (start,stop) in sortedset2:
                for x in range(start,stop+1):
                    retVal2.add(x)
                    
        elif type(sortedset2)==set:
                retVal2= sortedset2
        else:
            raise TypeError("Cannot merge something of class {0}".type(sortedset2))

        retVal=retVal2 | retVal1

        return(len(retVal1), len(retVal2), len(retVal), retVal1, retVal2, retVal)   
    def _sortedSetSubtract(self, set1, sortedset2):
        """ during reference mapping, large numbers of 'N's can be called - i.e. when there is no sequence corresponding to the reference.
        If we store all Ns in the reference, then we may have to store large numbers of positions.
        Typically, Ns occurs in blocks, with a start and end point.
        
        If the number of Ns is small, it may not be worth worrying about this.
        But if they are large, they can be represented efficiently (~ 54 bits) as a range, as opposed to taking up 24 bits per integer.
        We support use of the Banyan SortedSet class to represent a series of range of this type.
        
        When presented with sequence1 with (say) Ts at positions 100, 200 in set1,
        we wish to know whether positions 100, 200 are Ns in the comparator sequence.
        The Ns in the comparator are in sortedset2.
        
        Efficiently removing Ns from set1 iff they are in the ranges in sortedset2 is therefore necessary.
        Doing so rapidly is one of the key determinants of the speed of the whole seqComparer operation.
        
        This function returns set1-sortedset2
    """
        if type(sortedset2)==set:
            # then we can just use standard set subtraction; compression using SortedSets is not operational.
            return(set1-sortedset2)
        elif type(sortedset2)==SortedSet:
            subtracted=set1.copy()
            current_n_range=[]
            for element in set1:
                if len(sortedset2.overlap_point(element))>0:
                    subtracted.remove(element)
            return(subtracted)
        else:
            raise TypeError("Cannot perform sorted set subtraction on an object of class {0}".format(sortedset2))
    def countDifferences(self, cutoff=None, method='one'):            # appears to be faster than method 2
        if method=='one':
            nDiff=self.countDifferences_one(cutoff=cutoff)
        else:
            raise TypeError("Don't know how to perform difference operation using method {0}.".format(method))
        return(nDiff)
    
    def countDifferences_byKey(self,keyPair, cutoff=None):
        """ compares the in memory refCompressed sequences at
        self.seqProfile[key1] and self.seqProfile[key2]
        uses Algorithm 'one', below.
        Returns the number of SNPs between self.seq1 and self.seq2, and, if it is less than cutoff,
        the number of Ns in the two sequences and the union of their positions.
        Typical computational time is less than 0.5 msec."""


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
        self.seq1=self.seqProfile[key1]
        self.seq2=self.seqProfile[key2]
        nDiff=self.countDifferences_one(cutoff=cutoff)
        
        if nDiff is None:
            return((key1, key2, nDiff, None, None, None, None, None, None))
        elif nDiff<=cutoff:
            (n1, n2, nboth, N1pos, N2pos, Nbothpos) = self._setStats(self.seq1['N'],self.seq2['N'])
            return((key1, key2, nDiff, n1,n2,nboth, N1pos, N2pos, Nbothpos))
        else:
            return((key1, key2, nDiff, None, None, None, None, None, None))

    
    def countDifferences_one(self,cutoff=None):
        """ compares self.seq1 with self.seq2;
        these are set with self.setComparator1 and 2 respectively.
        Returns the number of SNPs between self.seq1 and self.seq2.
        

        rate about 7500 per second."""
        # if cutoff is not specified, we use snpCeiling
        if cutoff is None:
            cutoff = self.snpCeiling

            
        nDiff=0
        if self.seq1['invalid']==1 or self.seq2['invalid']==1:
            return(None)
        for nucleotide in ['C','G','A','T']:
       
            # we do not consider differences relative to the reference if the other nucleotide is an N
            # we invoke the SortedSetSubtract method to help us here.
            nonN_seq1=self._sortedSetSubtract(self.seq1[nucleotide],self.seq2['N'])
            nonN_seq2=self._sortedSetSubtract(self.seq2[nucleotide],self.seq1['N'])

            nDiff=nDiff+len(nonN_seq1 ^ nonN_seq2)
            if nDiff>cutoff:
                return(None)
        return(nDiff)    

class test_seqComparer_init1(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, persistenceStore='localmemory', startAfresh=False)
class test_seqComparer_init2(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, persistenceStore='localmemory', startAfresh=True)
class test_seqComparer_init5(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, persistenceStore='tofile', startAfresh=False)      
class test_seqComparer_init6(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, persistenceStore='tofile', startAfresh=True)              
class test_seqComparer_1(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True)
        self.assertEqual(sc.reference,refSeq)     
class test_seqComparer_2(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'

        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True)

        with self.assertRaises(TypeError):
            retVal=sc.compress(sequence='AC')
class test_seqComparer_3(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True, )
        retVal=sc.compress(sequence='ACTG')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([]), 'invalid':0})
class test_seqComparer_4(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'

        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True)

        retVal=sc.compress(sequence='ACTN')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([3]), 'invalid':0})
class test_seqComparer_5(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True)
        retVal=sc.compress(sequence='ACT-')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([3]), 'invalid':0})         
class test_seqComparer_6(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'

        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True)

        retVal=sc.compress(sequence='TCT-')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([0]), 'N': set([3]), 'invalid':0})
class test_seqComparer_7(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'

        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True)
        retVal=sc.compress(sequence='ATT-')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([1]), 'N': set([3]), 'invalid':0})
class test_seqComparer_3a(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True, NCompressionCutoff=0)
        retVal=sc.compress(sequence='ACTG')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': SortedSet([]), 'invalid':0})
class test_seqComparer_4a(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True, NCompressionCutoff=0)
 
        retVal=sc.compress(sequence='ACTN')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': SortedSet([(3,3)]), 'invalid':0})
class test_seqComparer_5a(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True, NCompressionCutoff=0)
 
        retVal=sc.compress(sequence='ACT-')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': SortedSet([(3,3)]), 'invalid':0})         
class test_seqComparer_6a(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True, NCompressionCutoff=0)
 
        retVal=sc.compress(sequence='TCT-')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([0]), 'N': SortedSet([(3,3)]), 'invalid':0})
class test_seqComparer_7a(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True, NCompressionCutoff=0)
 
        retVal=sc.compress(sequence='ATT-')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([1]), 'N': SortedSet([(3,3)]), 'invalid':0})
class test_seqComparer_8(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff=1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True)

        sc.setComparator1(sequence='ACTG')
        sc.setComparator2(sequence='ACTG')
class test_seqComparer_9(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True)
        sc.setComparator1(sequence='ACTG')
        sc.setComparator2(sequence='ACTG')
        self.assertEqual(sc.countDifferences(method='one'),0)
class test_seqComparer_10(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True)
        sc.setComparator1(sequence='TTTG')
        sc.setComparator2(sequence='ACTG')
        self.assertEqual(sc.countDifferences(method='one'),2)
class test_seqComparer_11(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True)
        sc.setComparator1(sequence='TTTG')
        sc.setComparator2(sequence='NNTG')
        self.assertEqual(sc.countDifferences(method='one'),0)
class test_seqComparer_12(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True)

        sc.setComparator2(sequence='TTTG')
        sc.setComparator1(sequence='NNTG')
        self.assertEqual(sc.countDifferences(method='one'),0)
class test_seqComparer_13(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True)
        sc.setComparator2(sequence='TTTG')
        sc.setComparator1(sequence='--TG')
        self.assertEqual(sc.countDifferences(method='one'),0)
class test_seqComparer_14(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True)
        sc.setComparator2(sequence='TTAA')
        sc.setComparator1(sequence='--AG')
        self.assertEqual(sc.countDifferences(method='one'),1)
class test_seqComparer_15(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True)
        sc.setComparator1(sequence='TTAA')
        sc.setComparator2(sequence='--AG')
        self.assertEqual(sc.countDifferences(method='one'),1)

class test_seqComparer_saveload3(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True, persistenceStore='localmemory', persistenceDir=os.path.join('..','unittest_tmp'))
        compressedObj =sc.compress(sequence='ACTT')
        sc.persist(compressedObj, guid='one', method='localmemory')     # no guid supplied
        retVal=sc.load(guid='one', method='localmemory')
        self.assertEqual(compressedObj,retVal)        
class test_seqComparer_saveload4(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True, persistenceStore='tofile', persistenceDir=os.path.join('..','unittest_tmp'))
        compressedObj =sc.compress(sequence='ACTT')
        sc.persist(compressedObj, guid='one', method='tofile')     # no guid supplied
        retVal=sc.load(guid='one', method='tofile')
        self.assertEqual(compressedObj,retVal)
class test_seqComparer_23(unittest.TestCase):
    """ tests N compression """
    def runTest(self):
        
        refSeq='ACTGTTAATTTTTTTTTGGGGGGGGGGGGAA'
        sc=seqComparer(maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True, NCompressionCutoff=0)

        retVal=sc.compress(sequence='ACTGTTAATTTTTTTTTGGGGGGGGGGGGAA')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': SortedSet([]), 'invalid':0})
class test_seqComparer_24(unittest.TestCase):
    """ tests N compression """
    def runTest(self):
        
        refSeq=                     'ACTGTTAATTTTTTTTTGGGGGGGGGGGGAA'
        sc=seqComparer(maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True, NCompressionCutoff=10000)

        retVal=sc.compress(sequence='ACTGTTAANNNNNNNNTGGGGGGGGGGGGAA')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([8,9,10,11,12,13,14,15]), 'invalid':0})
        retVal=sc.compress(sequence='NNTGTTAANNNNNNNNTGGGGGGGGGGGGAA')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set([0,1,8,9,10,11,12,13,14,15]), 'invalid':0})
class test_seqComparer_23a(unittest.TestCase):
    """ tests N compression """
    def runTest(self):
        
        refSeq='ACTGTTAATTTTTTTTTGGGGGGGGGGGGAA'
        sc=seqComparer(maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True, NCompressionCutoff=100000)

        retVal=sc.compress(sequence='ACTGTTAATTTTTTTTTGGGGGGGGGGGGAA')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': set(), 'invalid':0})
class test_seqComparer_24b(unittest.TestCase):
    """ tests N compression """
    def runTest(self):
        
        refSeq=                     'ACTGTTAATTTTTTTTTGGGGGGGGGGGGAA'
        sc=seqComparer(maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True, NCompressionCutoff=0)

        retVal=sc.compress(sequence='ACTGTTAANNNNNNNNTGGGGGGGGGGGGAA')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': SortedSet([(8,15)]), 'invalid':0})
        retVal=sc.compress(sequence='NNTGTTAANNNNNNNNTGGGGGGGGGGGGAA')
        self.assertEqual(retVal,{'G': set([]), 'A': set([]), 'C': set([]), 'T': set([]), 'N': SortedSet([(0,1),(8,15)]), 'invalid':0})

class test_seqComparer_25(unittest.TestCase):
    """ tests N compression """
    def runTest(self):
        
        refSeq=                     'ACTGTTAATTTTTTTTTGGGGGGGGGGGGAA'
        sc=seqComparer(maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True, NCompressionCutoff=0)

        
        Ns=[0,1,2,3,4,5,6,9,10]
        
        retVal=list(sc._ranges(Ns))
        self.assertEqual(retVal, [(0,6),(9,10)])

class test_seqComparer_26(unittest.TestCase):
    """ tests N compression """
    def runTest(self):
        
        refSeq=                     'ACTGTTAATTTTTTTTTGGGGGGGGGGGGAA'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True)
        
        Ns=[]
        
        retVal=list(sc._ranges(Ns))
        self.assertEqual(retVal, [])
        
class test_seqComparer_27(unittest.TestCase):
    """ tests N compression """
    def runTest(self):
        
        refSeq=                     'ACTGTTAATTTTTTTTTGGGGGGGGGGGGAA'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, snpCeiling = 20,reference=refSeq, startAfresh=True)
        
        Ns=[1]
        
        retVal=list(sc._ranges(Ns))
        self.assertEqual(retVal, [(1,1)])
        
        
class test_seqComparer_28(unittest.TestCase):
    """ tests _sortedSubtract """
    def runTest(self):
        
        refSeq=                             'ACTGTTAATTTTTTTTTGGGGGGGGGGGGAA'
        sc=seqComparer (maxNs = 1e8, snpCeiling = 20,reference=refSeq, NCompressionCutoff=0, startAfresh=True)
        compressedObj1=sc.compress(sequence='GGGGTTAANNNNNNNNNGGGGGAAAAGGGAA')
        compressedObj2=sc.compress(sequence='ACTGTTAATTTTTTTTTNNNNNNNNNNNNNN')
        self.assertEqual(compressedObj1, {'invalid': 0, 'C': set(), 'G': {0, 1, 2}, 'T': set(), 'A': {24, 25, 22, 23}, 'N': SortedSet([(8, 16)])})
        self.assertEqual(compressedObj2, {'invalid': 0, 'C': set(), 'G': set(), 'T': set(), 'A': set(), 'N': SortedSet([(17, 30)])})
        nonN_seq1=sc._sortedSetSubtract(compressedObj1['A'],compressedObj2['N'])
        self.assertEqual(nonN_seq1,set())
        #self.assertEqual(intersect1,{22,23,24,25})
        nonN_seq2=sc._sortedSetSubtract(compressedObj1['G'],compressedObj2['N'])
        self.assertEqual(nonN_seq2,{0,1,2})
        #self.assertEqual(intersect2,set())
        
class test_seqComparer_29(unittest.TestCase):
    """ tests _setStats """
    def runTest(self):
        
        refSeq=                             'ACTGTTAATTTTTTTTTGGGGGGGGGGGGAA'
        sc=seqComparer(maxNs = 1e8, snpCeiling = 20,reference=refSeq, NCompressionCutoff=0, startAfresh=True)
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
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, reference=refSeq, startAfresh=True, snpCeiling= 1)
        sc.setComparator1(sequence='ACTG')
        sc.setComparator2(sequence='ACTG')
class test_seqComparer_31(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, reference=refSeq, startAfresh=True, snpCeiling =1 )
        sc.setComparator1(sequence='ACTG')
        sc.setComparator2(sequence='ACTG')
        self.assertEqual(sc.countDifferences(method='one'),0)
class test_seqComparer_32(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, reference=refSeq, startAfresh=True, snpCeiling =1)
        sc.setComparator1(sequence='TTTG')
        sc.setComparator2(sequence='ACTG')
        self.assertEqual(sc.countDifferences(method='one'),None)
class test_seqComparer_33(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, reference=refSeq, startAfresh=True, snpCeiling =1)
        sc.setComparator1(sequence='TTTG')
        sc.setComparator2(sequence='NNTG')
        self.assertEqual(sc.countDifferences(method='one'),0)
class test_seqComparer_34(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, reference=refSeq, startAfresh=True, snpCeiling =1)
        sc.setComparator2(sequence='TTTG')
        sc.setComparator1(sequence='NNTG')
        self.assertEqual(sc.countDifferences(method='one'),0)
class test_seqComparer_13(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, reference=refSeq, startAfresh=True, snpCeiling =1)
        sc.setComparator2(sequence='TTTG')
        sc.setComparator1(sequence='--TG')
        self.assertEqual(sc.countDifferences(method='one'),0)
class test_seqComparer_35(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, reference=refSeq, startAfresh=True, snpCeiling =1)
        sc.setComparator2(sequence='TTAA')
        sc.setComparator1(sequence='--AG')
        self.assertEqual(sc.countDifferences(method='one'),1)
class test_seqComparer_36(unittest.TestCase):
    def runTest(self):
        refSeq='ACTG'
        sc=seqComparer(NCompressionCutoff = 1e8, maxNs = 1e8, reference=refSeq, startAfresh=True, snpCeiling =1)
        sc.setComparator1(sequence='TTAA')
        sc.setComparator2(sequence='--AG')
        self.assertEqual(sc.countDifferences(method='one'),1)
      
