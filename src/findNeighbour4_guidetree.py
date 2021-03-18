""" clustering for findNeighbour4
assumes a findNeighbour4 server is running, with the connection string stated in ../demos/AC587/config/config_cl.json.

An example command doing this would be (starting from /src)

pipenv run python3 findNeighbour4_server.py ../demos/AC587/config/config_cl.json

The test performs clustering.
"""

# import libraries
import os
import sys
import pathlib
import json
import logging

import logging.handlers
import warnings
import datetime
import glob
import sys
import pandas as pd
import numpy as np
import sentry_sdk
import argparse
import progressbar
import time
from random import sample as random_sample
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sentry_sdk import capture_message, capture_exception
from sentry_sdk.integrations.flask import FlaskIntegration

# logging
from logging.config import dictConfig

# startup
from mongoStore import fn3persistence
from hybridComparer import hybridComparer
from read_config import ReadConfig

if __name__ == '__main__':

    snp_cutoff = 3                                 # how many snp from a sample should we include neighbours
    sample_quality_cutoff = 0.9                    # only build tree from samples with >= 90% ACTG
    max_neighbourhood_mixpore = None               # set to integer, e.g. 50 to do mixpore check on all samples
    
    # command line usage.  Pass the location of a config file as a single argument.
    parser = argparse.ArgumentParser(
        formatter_class= argparse.RawTextHelpFormatter,
        description="""Runs findNeighbour4_guidetree, a findNeighbour4 component.
                                     

Example usage: 
============== 

## does not require findNeighbour4_server to be running

Minimal example:
python findNeighbour4_guidetree.py ../config/myConfigFile.json  

if a config file is not provided, it will terminate. 

"""
)


    parser.add_argument('path_to_config_file', type=str, action='store', nargs='?',
                        help='the path to the configuration file', default=''  )
    args = parser.parse_args()
    
    # an example config file is default_test_config.json

    ############################ LOAD CONFIG ######################################
    print("findNeighbour4 guidetree .. reading configuration file.")

    if len(args.path_to_config_file)>0:
            configFile = args.path_to_config_file
            debugmode = False
            logging.info(configFile)
    else:
            configFile = os.path.join('..','config','default_test_config.json')
            debugmode = True
            warnings.warn("No config file name supplied ; using a configuration ('default_test_config.json') suitable only for testing, not for production. ")

    rc = ReadConfig()
    CONFIG = rc.read_config(configFile)


    ########################### SET UP LOGGING #####################################  
    # create a log file if it does not exist.
    print("Starting logging")
    logdir = os.path.dirname(CONFIG['LOGFILE'])
    pathlib.Path(os.path.dirname(CONFIG['LOGFILE'])).mkdir(parents=True, exist_ok=True)

    # set up logger
    logger = logging.getLogger()
    loglevel=logging.INFO
    if 'LOGLEVEL' in CONFIG.keys():
            if CONFIG['LOGLEVEL']=='WARN':
                    loglevel=logging.WARN
            elif CONFIG['LOGLEVEL']=='DEBUG':
                    loglevel=logging.DEBUG

    # configure logging object 
    logger.setLevel(loglevel)       
    logfile = os.path.join(logdir, "clustering-{0}".format(os.path.basename(CONFIG['LOGFILE'])))
    print("Logging to {0} with rotation".format(logfile))
    file_handler = logging.handlers.RotatingFileHandler(logfile, mode = 'a', maxBytes = 1e7, backupCount = 7)
    formatter = logging.Formatter( "%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s ")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # determine whether a FN_SENTRY_URLenvironment variable is present,
    # if so, the value of this will take precedence over any values in the config file.
    # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
    if os.environ.get("FN_SENTRY_URL") is not None:
        CONFIG["SENTRY_URL"] = os.environ.get("FN_SENTRY_URL")
        print("Set Sentry connection string from environment variable")
    else:
        print("Using Sentry connection string from configuration file.")
        
    # launch sentry if API key provided
    if 'SENTRY_URL' in CONFIG.keys():
            logger.info("Launching communication with Sentry bug-tracking service")
            sentry_sdk.init(CONFIG['SENTRY_URL'], integrations=[FlaskIntegration()])

    ########################### prepare to start clustering ####################################
    # construct data access object
    
    logger.info("Connecting to backend data store")
    try:
            PERSIST=fn3persistence(dbname = CONFIG['SERVERNAME'],
                connString=CONFIG['FNPERSISTENCE_CONNSTRING'],
                debug=0
                                   )  # if in debug mode wipes all data.  This is not what is wanted here, even if we are using unittesting database

    except Exception as e:
            logger.exception("Error raised on creating persistence object")
            raise

     ################################# clustering #############################################
    # open PERSIST and hybridComparer object used by all samples
    # this is only used for data access and msa.
    # inserts are not allowed
    logger.info("Building hybridComparer object")
    hc = hybridComparer(reference=CONFIG['reference'],
        maxNs=CONFIG['MAXN_STORAGE'],
        snpCeiling=  CONFIG['SNPCEILING'],
        excludePositions=CONFIG['excluded'],
        preComparer_parameters={},
        PERSIST=PERSIST,
        disable_insertion = True)
    
    

    print("Loading all samples & annotations.")
    sample_annotations = pd.DataFrame.from_dict(PERSIST.guid2items(None,None), orient='index')
    all_samples = sample_annotations.index.to_list()
    high_quality_samples = sample_annotations.loc[sample_annotations['DNAQuality:propACTG']>sample_quality_cutoff].copy()
    
    sampling_population = high_quality_samples.index.to_list()

    all_samples = set(sampling_population)
    sampling_population = set(sampling_population)          # samples now (note: more may be added during computations by other processes, so we pick a set to work on at the start of the computations)

    print("From all {0} eligible samples currently present, studying {1} high quality (ACTG > 0.9) samples ".format(len( sample_annotations.index.to_list()), len(sampling_population)))
    representative_subsample = set()
    remaining_samples = sampling_population.copy()
    representative_size = {}

    ## DEBUG 
    # remaining_samples = set(random_sample(list(remaining_samples), 5) )

    mixed_samples = set()
    singletons = set()
    iteration = 0
    total_delta = 0

    total_n_samples  = len(remaining_samples)
    bar = progressbar.ProgressBar(max_value=len(remaining_samples))
     
    while len(remaining_samples)>0:

        iteration +=1
        test_sample = random_sample(list(remaining_samples),1)[0]        # randomly sample one
        test_sample_neighbour_info = PERSIST.guid2neighbours(test_sample, snp_cutoff, returned_format=3)        
        
        # filter these neighbours.   Only include neighbours which are part of all_samples.
        # because samples get added to the server all the time, this is not guaranteed.
        test_sample_neighbours = []
        for item in test_sample_neighbour_info['neighbours']:
            if item in all_samples:
                test_sample_neighbours.append(item)

        # maybe put in place a quality filter (total Ns)
        delta = 1           # one less at a minimum (test_sample is never reconsidered);

        if len(test_sample_neighbours) == 0:
            singletons.add(test_sample)

        else: 
            representative_subsample.add(test_sample)
            representative_size[test_sample] = len(test_sample_neighbours)

            # remove sample and neighbours from further consideration
            n1 = len(remaining_samples)
            for item in test_sample_neighbours:
                remaining_samples.discard(item)         # discard as item may be low quality or otherwise pre-excluded
            n2 = len(remaining_samples)
            delta = 1+ n1- n2
            total_delta = total_delta + delta   # total number of samples

        # remove the test sample, if it has not already been removed
        remaining_samples.discard(test_sample)
        if len(representative_subsample)>0:
            ratio_selected = total_delta/float(len(representative_subsample))
        else: 
            ratio_selected = None
          #print("{0} | Iteration: {1} | Remaining samples {2} |  Change now/total {9}/{10} | Sample = {3} |  Mixture statistic {4} | Neighbours {5} | Selected representatives {6} | Mixed samples {7} | Singletons {8} | Approx Subsample 1 in {11:.0f}".format(
          #datetime.datetime.now().isoformat(),
          #iteration,
          #len(remaining_samples),
          #test_sample,
          #'ND',
          #len(test_sample_neighbours),
          #len(representative_subsample),
          #len(mixed_samples),
          # len(singletons),
          #delta, 
          # total_delta,
          #ratio_selected
          
        #)
        bar.update(total_n_samples-len(remaining_samples))


    print("Generating fasta output from {0} samples; ignoring {1} singletons.".format(len(representative_subsample), len(singletons)))
    # export masked fasta (note: holds all selected sequences in ram)
    fasta_outputfile =  'tree_selected.fasta' 
    seqs = []
    for test_sample in representative_subsample:
        seq = hc.uncompress_guid(test_sample)
        seq = seq.replace('N','-')  ## fastTree convention
        sr = SeqRecord(     Seq(seq), 
                                id= test_sample,
                                description="| within {0} snp of {1} others".format(snp_cutoff, representative_size[test_sample]) 
                                )        
        seqs.append(sr)

    with open(fasta_outputfile, 'w') as f:
        SeqIO.write(seqs, f, "fasta")
    print("Selected samples are written to {0}".format(fasta_outputfile))
    print("Finished")
    

