""" runs  findNeighbour3 benchmark
assumes a findNeighbour3 server is running, with the connection string stated in the config file passed as the single argument.

"""

import os
import glob
import datetime
import pandas as pd
import argparse
import json
from fn3client import fn3Client

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class= argparse.RawTextHelpFormatter,
        description="""Runs benchmarks against findNeighbour3-server""")
    parser.add_argument('path_to_server_config_file', type=str, action='store', nargs='?',help='the path to the configuration file',default='')
    args = parser.parse_args()
    if len(args.path_to_server_config_file)>0:
        configFile = args.path_to_server_config_file
    else:
        print("Must supply config file of server as the single argument")
        exit(1)
    with open(configFile,'r') as f:
        CONFIG=f.read()
    CONFIG=json.loads(CONFIG)
    description= CONFIG['SERVERNAME']
			    
    def ll2s(x):
        """ converts a list of lists, e.g. [['guid1',2],['guid2',0]] into a set {'guid1','guid2'} """
        neighbour_set = set()
        for neighbour in x:
            neighbour_set.add(neighbour[0])
        return neighbour_set

    # define directory where the fastas are
    fastadir = os.path.join('..','demos','benchmark','fasta')
    outputdir = os.path.join('..','demos','benchmark','output')
     
    # instantiate client
    baseurl = "http://{0}:{1}".format(CONFIG['IP'],CONFIG['REST_PORT'])
    print("Contacting fn3 server on ",baseurl)
    fn3c = fn3Client(baseurl=baseurl)      # expects operation on local host; pass baseurl if somewhere else.

    # names of the clustering algorithms
    clusters=fn3c.clustering()

    existing_guids = set(fn3c.guids())
    clustering_created = False

    print("There are {0} existing guids".format(len(existing_guids)))

    for i,fastafile in enumerate(sorted(glob.glob(os.path.join(fastadir,  '*.fasta')))):

        guid = os.path.basename(fastafile).replace('.fasta','')
        read_failed = False
        try:
            seq = fn3c.read_fasta_file(fastafile)['seq']
        except IOError:
            read_failed = True
        
        if not read_failed:    
            if not guid in existing_guids:
                fn3c.insert(guid=guid,seq=seq)
                result= 'inserted'
            else:
                result = 'exists, skipped re-insert'
            print("Test   ",datetime.datetime.now(), i, guid, result)

    # recover the change_id
    print("Recovering clustering data")
    clustering_created = False
    for clustering_algorithm in clusters['algorithms']:
        change_id = fn3c.change_id(clustering_algorithm)['change_id']
        print("There were {0} changes with pipeline {1}".format(change_id, clustering_algorithm))
        df = fn3c.guids2clusters(clustering_algorithm, after_change_id=0) 
        df['clustering_algorithm']=clustering_algorithm
        if not clustering_created:
           clustering_df = df
           clustering_created = True
        else:
           clustering_df = pd.concat([clustering_df,df], ignore_index=True, sort=False)
    if clustering_created:
        clustering_df.to_csv(os.path.join(outputdir, "clustering_{0}.csv".format(description)))

    print("Recovering memory usage history")
    df = fn3c.server_memory_usage(nrows=int(100000000))
    df.to_csv(os.path.join(outputdir, "memoryusage_{0}.csv".format(description)))


    print("Complete")

