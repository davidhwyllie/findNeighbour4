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

from seqComparer import seqComparer


# test code
if __name__ == "__main__":
    
    reffile = os.path.join('..','reference','NC_000962.fasta')
    with open(reffile) as f:
        for record in SeqIO.parse(f, format='fasta'):
            refSeq = record.seq
            
    # create sc object
    sc=seqComparer( maxNs = 1e8,
                       reference=refSeq,
                       startAfresh=True,
                       snpCeiling =20,
                       persistenceDir=os.path.join('..','unittest_tmp'))
    
    # the below data contains test data from https://ora.ox.ac.uk/objects/uuid:82ce6500-fa71-496a-8ba5-ba822b6cbb50
    inputdir = os.path.join("E:/", "dwyllie", "TBTESTDATA")
    globpath = os.path.join(inputdir, '*.fasta')
    nRead = 0
    guids = []
    for fastafile in glob.glob(globpath):
        guid = os.path.basename(fastafile)[0:36]

        with open(fastafile) as f:
            for record in SeqIO.parse(f, format='fasta'):
                this_seq = str(record.seq)
                all_n = this_seq.count('N')+this_seq.count('-')
                if all_n < 400000:
                    print(guid, all_n)
                    guids.append(guid)
                    c = sc.compress(this_seq)
                    sc.persist(c, guid=guid, method='localmemory')
            nRead+=1
            if nRead == 500:
                break
    # estimate memory used
    p1 = len(pickle.dumps(sc.seqProfile))
    p2 = len(pickle.dumps(sc.consensi))
    print('Precompression','',p1, p2, p1+p2)

    covered = set()
    for guid in guids:
        if not guid in covered:
            res = sc.compress_relative_to_consensus(guid)
            covered = covered | set(res)
            if len(res)>1:      # we found matches
                p1 = len(pickle.dumps(sc.seqProfile))
                p2 = len(pickle.dumps(sc.consensi))
                print('Compressing',guid, p1, p2, p1+p2)        
    