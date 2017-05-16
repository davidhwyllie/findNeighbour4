#!/usr/bin/env python3
""" code for profiling the core of EW2, independent of the server component """

from seqComparer import seqComparer
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
import uuid
import json
import psutil
       
if __name__=="__main__":
        
    # read reference
    inputref=os.path.join("..","reference","reference.fasta")
    with open(inputref,'rt') as f:
        for r in SeqIO.parse(f,'fasta'):
            reference=str(r.seq)
            
    # create sequence comparer object using reference, and a supply a directory to store the persistent objects;
    # instantiation will load the profiles, but not the full sequences, into memory.
    # defaults are supplied by the constructor about what proportion of bases are to be sampled.
    print("Startup")
    sc=seqComparer(reference=reference, \
                       persistenceStore='localmemory', \
                       persistenceDir="/home/dwyllie/data/relatednesstest/TB_FASTA_TESTDIR", \
                       startAfresh=False)
    
    # note timings
    startParse=datetime.datetime.now()
    nTested=0
    keyList=[]
    # we are going to read the sequences from a test file provided by Trien.
    testpath="/home/dwyllie/data/relatednesstest/TB_FASTA/*_v3.fasta" 
    nRead=0
    fastaFiles=glob.glob(testpath)
    for fastaFile in fastaFiles:
        for seq_record in SeqIO.parse(fastaFile, 'fasta'):
            guid=os.path.basename(fastaFile)[0:36]
            print(guid)
            
            seq=str(seq_record.seq)
            nTested+=1
            if nTested % 50==0:
                print("Loaded {0}".format(nTested), )
            if not sc.iscachedinram(guid):
                co=sc.compress(seq)
                sc.persist(refCompressedSequence=co, guid=guid, method='localmemory')     # store the parsed object in ram and on disc

            keyList.append(guid)
            print(guid)
        if nTested>500:
            break
          
    # note  the end of parse time.      
    endParse=datetime.datetime.now()
    #print("LOADED to ",this_persistenceStore, startParse, endParse, nTested, "Time per sample=",(endParse-startParse)/nTested)
    try:
        perEvent=(endParse-startParse)/nTested
    except ZeroDivisionError:
        print("No sequences found")
        exit()
        
    print("Loading speed per sample was (seconds) ",perEvent.total_seconds())

    ## run timings for comparisons
    print("Timing comparison ")
    startCmps=datetime.datetime.now()
    nCmps=0
    
    # compare all vs all:
    with open("/home/dwyllie/datapairlist.txt",'wt') as f:
        f.write("{0}\t{1}\t{2}\n".format('guid1','guid2','dist', 'guid1Npos','guid2Npos','sharedNpos'))
        for key1 in keyList:
            for key2 in keyList:
                if key1<key2:
                    nCmps+=1
                    (guid1,guid2,dist, guid1N, guid2N, sharedN, guid1N, guid2N, sharedNSet)=sc.countDifferences_byKey(keyPair=(key1,key2))
                    if dist is not None:
                        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(guid1,guid2,dist, guid1N, guid2N, sharedN))
        endCmps=datetime.datetime.now()
        perEvent=(endCmps-startCmps)/nCmps
        print("Time per comparison (seconds) =",perEvent.total_seconds())

