""" illustrates use of findNeighbour3
assumes a findNeighbour3 server is running, with the connection string stated in ../config/config.json.

An example command doing this would be (starting from /src)

python3 findNeighbour3-server.py ../demos/AC587/config/config.json
"""

import os
import glob
import datetime
import pandas as pd
from fn4client import fn4Client

if __name__ == '__main__':
    # define directory where the fastas are
    fastadir = os.path.join('demos','AC587','fasta')
    outputdir = os.path.join('demos','AC587','output')

    # instantiate client
    fn4c = fn4Client("http://localhost:5032")      # expects operation on local host; pass baseurl if somewhere else.

    # names of the clustering algorithms
    clusters=fn4c.clustering()

    # add control fasta files.  The system evaluates the %N in terms of the population existing
    # we load 40 randomly selected guids as controls

    for i,fastafile in enumerate(glob.glob(os.path.join(fastadir, 'controls','*.mfasta.gz'))):
        guid = os.path.basename(fastafile).replace('.mfasta.gz','')
        seq = fn4c.read_fasta_file(fastafile)['seq']
        print("Controls",datetime.datetime.now(), i, guid)
        fn4c.insert(guid=guid,seq=seq)
     
     
    for i,fastafile in enumerate(sorted(glob.glob(os.path.join(fastadir, 'test', '*.mfasta.gz')))):
        guid = os.path.basename(fastafile).replace('.mfasta.gz','')
        seq = fn4c.read_fasta_file(fastafile)['seq']
        print("Test",datetime.datetime.now(), i, guid)
        fn4c.insert(guid=guid,seq=seq)
        
    print("Finished")
