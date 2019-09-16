""" illustrates use of findNeighbour3
assumes a findNeighbour3 server is running, with the connection string stated in ../config/config.json.

An example command doing this would be (starting from /src)
python3 findNeighbour3-server.py ../demos/AA041/config/config.json
"""

import os
import glob
import datetime
import pandas as pd
from fn3client import fn3Client

def ll2s(x):
    """ converts a list of lists, e.g. [['guid1',2],['guid2',0]] into a set {'guid1','guid2'} """
    neighbour_set = set()
    for neighbour in x:
        neighbour_set.add(neighbour[0])
    return neighbour_set

# define directory where the fastas are
fastadir = os.path.join('/srv','data','mixfiles','mfasta')

# instantiate client
fn3c = fn3Client("http://localhost:5026")      # expects operation on local host; pass baseurl if somewhere else.

existing_guids = set(fn3c.guids())
clustering_created = False
print("There are {0} existing guids".format(len(existing_guids)))

for i,fastafile in enumerate(sorted(glob.glob(os.path.join(fastadir,  '*.mfasta.gz')))):
    t1 = datetime.datetime.now()
    
    guid = os.path.basename(fastafile).replace('.mfasta.gz','')
    if not guid in existing_guids:
        read_failed = False
        try:
            seq = fn3c.read_fasta_file(fastafile)['seq']
        except IOError:
            read_failed = True
        
        if not read_failed: 
            fn3c.insert(guid=guid,seq=seq)
            result= 'inserted'
        else:
            result = 'read file failed'


    else:
        result = 'exists'

    t2 = datetime.datetime.now()
    i1 = t2-t1
    s = i1.total_seconds()
    
    print(datetime.datetime.now(), i, guid, result, s, "seconds")


