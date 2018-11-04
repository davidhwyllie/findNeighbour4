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
fastadir = os.path.join('..','demos','AA041','fasta')
outputdir = os.path.join('..','demos','AA041','output')

# instantiate client
fn3c = fn3Client()      # expects operation on local host; pass baseurl if somewhere else.

# names of the clustering algorithms
clusters=fn3c.clustering()

existing_guids = set(fn3c.guids())
clustering_created = False
print("There are {0} existing guids".format(len(existing_guids)))
# add control fasta files.  The system evaluates the %N in terms of the population existing
# we load 50 randomly selected guids as controls

for i,fastafile in enumerate(glob.glob(os.path.join(fastadir, 'control','*.fasta'))):
    guid = "ctrl_"+os.path.basename(fastafile).replace('.fasta','')
    seq = fn3c.read_fasta_file(fastafile)['seq']
    if not guid in existing_guids:
        fn3c.insert(guid=guid,seq=seq)
        result= 'inserted'
    else:
        result = 'exists, skipped re-insert'
    print("Controls",datetime.datetime.now(), i, guid, result)

for i,fastafile in enumerate(sorted(glob.glob(os.path.join(fastadir, 'test', '*.fasta')))):
    
    guid = os.path.basename(fastafile).replace('.fasta','')
    seq = fn3c.read_fasta_file(fastafile)['seq']
    if not guid in existing_guids:
        fn3c.insert(guid=guid,seq=seq)
        result= 'inserted'
        for clustering_algorithm in clusters['algorithms']:
        
            df = fn3c.guids2clusters(clustering_algorithm)
            df['step'] = i
            df['clustering_algorithm']=clustering_algorithm
            if not clustering_created:
                clustering_df = df
                clustering_created = True
            else:
                clustering_df = pd.concat([clustering_df,df], ignore_index=True, sort=False)

    else:
        result = 'exists, skipped re-insert'
    print("Test   ",datetime.datetime.now(), i, guid, result)

clustering_df.to_excel(os.path.join(outputdir, "clustering.xlsx"))