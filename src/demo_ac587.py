""" illustrates use of findNeighbour3
assumes a findNeighbour3 server is running, with the connection string stated in ../config/config.json.

An example command doing this would be (starting from /src)

python findNeighbour3-server.py ../demos/AC587/config/config.json
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
fastadir = os.path.join('..','demos','AC587','fasta')
outputdir = os.path.join('..','demos','AC587','output')

# instantiate client
fn3c = fn3Client()      # expects operation on local host; pass baseurl if somewhere else.

# names of the clustering algorithms
clusters=fn3c.clustering()

# add control fasta files.  The system evaluates the %N in terms of the population existing
# we load 40 randomly selected guids as controls

for i,fastafile in enumerate(glob.glob(os.path.join(fastadir, 'controls','*.fasta'))):
    guid = "ctrl_"+os.path.basename(fastafile).replace('.fasta','')
    seq = fn3c.read_fasta_file(fastafile)['seq']
    print("Controls",datetime.datetime.now(), i, guid)
    fn3c.insert(guid=guid,seq=seq)
 
 
for i,fastafile in enumerate(sorted(glob.glob(os.path.join(fastadir, 'test', '*.fasta')))):
    guid = os.path.basename(fastafile).replace('.fasta','')
    seq = fn3c.read_fasta_file(fastafile)['seq']
    print("Test   ",datetime.datetime.now(), i, guid)
    fn3c.insert(guid=guid,seq=seq)
    
    for clustering_algorithm in clusters['algorithms']:
    
        df = fn3c.guids2clusters(clustering_algorithm)
        df['step'] = i
        df['clustering_algorithm']=clustering_algorithm
        if i==0:
            clustering_df = df
        else:
            clustering_df = pd.concat([clustering_df,df], ignore_index=True, sort=False)

        # note the MSA for all clusters this sample is in
        df = fn3c.guids2clusters(clustering_algorithm)
        print(guid)
        print(df.query('guid==@guid'))

    #    df.to_excel(os.path.join(outputdir,"{0}.xlsx".format(guid)))

    # find neighbours
    #neighbour_set = ll2s(fn3c.guid2neighbours(guid, threshold=30))
    #print("found {0} neighbours".format(len(neighbour_set)))
    
    # compute pairwise links with neighbours [slow - we'll run post hoc]
    #print("Assessing mixtures")
    #df = fn3c.assess_mixed(guid, neighbour_set)
    #print("Mixture assessment over")
    #if not df is None:
    #    df['nneighbours']= len(neighbour_set)        
    #    df.to_excel(os.path.join(outputdir,"{0}.xlsx".format(guid)))

clustering_df.to_excel(os.path.join(outputdir, "clustering.xlsx"))