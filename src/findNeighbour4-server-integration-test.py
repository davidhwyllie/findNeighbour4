""" integration test for findNeighbour4
assumes a findNeighbour4 server is running, with the connection string stated in ../demos/AC587/config/config_nocl.json.

An example command doing this would be (starting from /src)

pipenv run python3 findNeighbour4-server.py ../demos/AC587/config/config_nocl.json

The test loads the server with data from the AC587 test data
and compares the SNP distances and links with those from a seqComparer
instance.

"""

import os
import glob
import datetime
import pandas as pd
from fn3client import fn3Client
from seqComparer import seqComparer
from preComparer import preComparer
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide
import unittest
from urllib.parse import urlparse as urlparser
from urllib.parse import urljoin as urljoiner
import uuid
import time

if __name__ == '__main__':
        
        # define directory where the fastas are
        fastadir = os.path.join('..','demos','AC587','fasta')

        # instantiate findNeighbour client
        fn3c = fn3Client()      # expects operation on local host; pass baseurl if somewhere else.

        # instantiate seqComparer
        inputfile = "../COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, 'rt') as f:
            for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                refSeq = str(record.seq)
        excludePositions=set()
        with open("../reference/TB-exclude-adaptive.txt",'rt') as f:
            rows=f.readlines()
            for row in rows:
                excludePositions.add(int(row))

        sc = seqComparer(   reference= refSeq, 
                    snpCeiling=2e6, 
                    maxNs=130000,
                    return_none_if_high_snp_distance=False,
                    excludePositions=excludePositions)
        pc = preComparer(**{"selection_cutoff":20,
                            "over_selection_cutoff_ignore_factor":5,
                            "mixed_reporting_cutoff":0,
                            "N_mean":40000,
                            "N_sd":9000,
                            "highN_z_reporting_cutoff":2,
                            "alpha":1e-14,
                            "probN_inflation_factor":3,
                            "n_positions_examined":3812800})
        # reset server
        fn3c.reset()

        # add fasta files to both server and seqComparer instance.  
        guids = set()
        for i,fastafile in enumerate(glob.glob(os.path.join(fastadir, 'controls','*.mfasta.gz'))):
            guid = os.path.basename(fastafile).replace('.mfasta.gz','')
            seq = fn3c.read_fasta_file(fastafile)['seq']
            print("Controls",datetime.datetime.now(), i, guid)
            fn3c.insert(guid=guid,seq=seq)
            
            obj = sc.compress(seq)
            sc.persist(obj, guid)
            pc.persist(obj, guid)
            guids.add(guid)
         
        for i,fastafile in enumerate(sorted(glob.glob(os.path.join(fastadir, 'test', '*.mfasta.gz')))):
            guid = os.path.basename(fastafile).replace('.mfasta.gz','')
            seq = fn3c.read_fasta_file(fastafile)['seq']
            if not guid in guids:   # not already entered
                print("Test",datetime.datetime.now(), i, guid)
                fn3c.insert(guid=guid,seq=seq)
             
                obj = sc.compress(seq)
                if obj['invalid']==0:
                    guids.add(guid)
                sc.persist(obj, guid)
                
                pc.persist(obj, guid)
                if obj['invalid']==0:
                    guids.add(guid)

         
        nFailed =0
        distances = {}
        for guid in guids:
            distances[guid] = fn3c.guid2neighbours(guid, threshold=20, quality_cutoff=0.85, timeout =None)  
        for res in sc.distmat(half=True):       # problem is with distance measurement, not with storage
            neighbours = distances[res[0]]
            pc_result = pc.compare(res[0],res[1])
            if res[0] in guids:
                target_guids = [x[0] for x in neighbours]
                if res[2] <= 20:
                    # guid2 should be present
                    test_success = res[1] in target_guids
                else:
                    # guid2 should be absent
                    test_success = not res[1] in target_guids 
                if not test_success:
                    nFailed +=1
                    print(res[0],res[1],res[2], test_success, pc_result)
                    print(neighbours)
             
        print("Finished.  Failures = {0}".format(nFailed))

