""" illustrates use of findNeighbour3
assumes a findNeighbour3 server is running, with the connection string stated in ../config/config.json.

An example command doing this would be (starting from /src)
python3 findNeighbour3-server.py ../demos/AA041/config/config.json
"""

import os
import glob
import datetime
import pandas as pd
from fn4client import fn4Client

if __name__ == '__main__':
    # define directory where the fastas are
    fastadir = os.path.join('..','demos','AA041','fasta')
    outputdir = os.path.join('..','demos','AA041','output')

    # instantiate client
    fn4c = fn4Client("http://127.0.0.1:5041")      # expects operation on local host; config file runs server on port 5031 for cw version and 5041 for python version.

    # names of the clustering algorithms
    clusters=fn4c.clustering()

    existing_guids = set(fn4c.guids())
    clustering_created = False
    print("There are {0} existing guids".format(len(existing_guids)))
    # add control fasta files.  The system evaluates the %N in terms of the population existing
    # we load 50 randomly selected guids as controls

    for i,fastafile in enumerate(glob.glob(os.path.join(fastadir, 'control','*.mfasta.gz'))):
        guid = "ctrl_"+os.path.basename(fastafile).replace('.mfasta.gz','')
        seq = fn4c.read_fasta_file(fastafile)['seq']
        if not guid in existing_guids:
            fn4c.insert(guid=guid,seq=seq)
            result= 'inserted'
        else:
            result = 'exists, skipped re-insert'
        print("Controls",datetime.datetime.now(), i, guid, result)

    for i,fastafile in enumerate(sorted(glob.glob(os.path.join(fastadir, 'test', '*.mfasta.gz')))):
        
        guid = os.path.basename(fastafile).replace('.mfasta.gz','')
        read_failed = False
        try:
            seq = fn4c.read_fasta_file(fastafile)['seq']
        except IOError:
            read_failed = True
        
        if not read_failed:    
            if not guid in existing_guids:
                fn4c.insert(guid=guid,seq=seq)
                result= 'inserted'
        print("Test",datetime.datetime.now(), i, guid, result)
              

    print("Finished")
