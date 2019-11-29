""" illustrates use of findNeighbour3
assumes a findNeighbour3 server is running, with the connection string stated in ../config/config.json.

An example command doing this would be (starting from /src)

python3 findNeighbour3-server.py ../demos/AC587/config/config.json
"""

import os
import glob
import datetime
import pandas as pd
from fn3client import fn3Client

# define directory where the fastas are
fastadir = os.path.join('..','demos','AC587','fasta')
outputdir = os.path.join('..','demos','AC587','output')

# instantiate client
fn3c = fn3Client("http://localhost:5033")      # expects operation on local host; pass baseurl if somewhere else.

# names of the clustering algorithms
clusters=fn3c.clustering()

# add control fasta files.  The system evaluates the %N in terms of the population existing
# we load 40 randomly selected guids as controls

for i,fastafile in enumerate(glob.glob(os.path.join(fastadir, 'controls','*.mfasta.gz'))):
    guid = os.path.basename(fastafile).replace('.mfasta.gz','')
    seq = fn3c.read_fasta_file(fastafile)['seq']
    print("Controls",datetime.datetime.now(), i, guid)
    fn3c.insert(guid=guid,seq=seq)
 
 
for i,fastafile in enumerate(sorted(glob.glob(os.path.join(fastadir, 'test', '*.mfasta.gz')))):
    guid = os.path.basename(fastafile).replace('.mfasta.gz','')
    seq = fn3c.read_fasta_file(fastafile)['seq']
    print("Test",datetime.datetime.now(), i, guid)
    fn3c.insert(guid=guid,seq=seq)
    
print("Finished")
