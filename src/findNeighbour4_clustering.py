""" clustering for findNeighbour4
assumes a findNeighbour4 server is running, with the connection string stated in ../demos/AC587/config/config_cl.json.

An example command doing this would be (starting from /src)

pipenv run python3 findNeighbour4_server.py ../demos/AC587/config/config_cl.json

The test performs clustering.
"""

# import libraries
import os
import sys
import requests
import json
import logging
import logging.handlers
import warnings
import datetime
import glob
import sys
import hashlib
import queue
import threading
import gc
import io
import pymongo
import pandas as pd
import numpy as np
import copy
import pathlib
import markdown
import codecs
import sentry_sdk
import matplotlib
import dateutil.parser
import argparse
import progressbar
import time
import progressbar

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_nucleotide
from sentry_sdk import capture_message, capture_exception
from sentry_sdk.integrations.flask import FlaskIntegration

# logging
from logging.config import dictConfig

# utilities for file handling and measuring file size
import psutil

# startup
from mongoStore import fn3persistence
from ma_linkage import MixtureAwareLinkage,MixPOREMixtureChecker, MixtureAwareLinkageResult
from msa import MSAStore
from hybridComparer import hybridComparer
from read_config import ReadConfig
from clusternomenclature import ClusterNameAssigner, ClusterNomenclature

if __name__ == '__main__':

    # command line usage.  Pass the location of a config file as a single argument.
    parser = argparse.ArgumentParser(
        formatter_class= argparse.RawTextHelpFormatter,
        description="""Runs findNeighbour4_clustering, a findNeighbour4 component.
                                     

Example usage: 
============== 

## does not require findNeighbour4_server to be running
Minimal example:
python findNeighbour4_clustering.py ../config/myConfigFile.json  

# Relabel clusters
nohup pipenv run python3 findNeighbour4_clustering.py ../phe_dev/config_phe_dev.json --rebuild_clusters_debug --label_clusters_using=../reference/guid2cluster.xlsx &

if a config file is not provided, it will run (as does findNeighbour4_server) is debug mode: it will run once, and then terminate.  This is useful for unit testing.  If a config file is specified, the clustering will  run until terminated.  

Checks for new sequences are conducted once per minute.



""")
    parser.add_argument('path_to_config_file', type=str, action='store', nargs='?',
                        help='the path to the configuration file', default=''  )
    parser.add_argument('--rebuild_clusters_debug', help='delete existing, and rebuild, clusters.  Only for use in a debug setting', action='store_true')
    parser.add_argument('--label_clusters_using', help='label clusters based on existing membership, as in the two column excel file referred to.  The program will stop after the initial run if this option is used, and will need to be restarted without this option to function normally', action='store')
    args = parser.parse_args()
    
    # an example config file is default_test_config.json

    ############################ LOAD CONFIG ######################################
    print("findNeighbour4 clustering .. reading configuration file.")

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

    # launch sentry if API key provided
    if 'SENTRY_URL' in CONFIG.keys():
            logger.info("Launching communication with Sentry bug-tracking service")
            sentry_sdk.init(CONFIG['SENTRY_URL'], integrations=[FlaskIntegration()])

    ########################### read file containing labels, if it exists ##########
    relabel = False  
    if args.label_clusters_using is not None:
        if os.path.exists(args.label_clusters_using):
            print("File of cluster labels {0} exists; reading it".format(args.label_clusters_using))
            label_df = pd.read_excel(args.label_clusters_using)
            relabel= True
            existing_labels = set()
            # set column headers
            label_df.columns = ['original_label','guid']
            # drop anything which empty labels
            label_df.dropna(inplace=True)
            label_df['label'] = [x[:5] for x in label_df['original_label']]
            previous_guid2cluster_label = {}
            for ix in label_df.index:
                existing_labels.add(label_df.at[ix,'label'])
                previous_guid2cluster_label[label_df.at[ix,'guid']] = [label_df.at[ix,'label']]
            existing_labels = list(existing_labels)
            logging.info("Relabelling samples based on existing labels. In existing data, {0} guids have {1} labels".format(len(previous_guid2cluster_label),len(existing_labels)))


    ########################### prepare to start clustering ####################################
    # construct the required global variables
    
    logger.info("Connecting to backend data store")
    try:
            PERSIST=fn3persistence(dbname = CONFIG['SERVERNAME'],
                connString=CONFIG['FNPERSISTENCE_CONNSTRING'],
                debug=0
                                   )  # if in debug mode wipes all data.  This is not what is wanted here, even if we are using unittesting database

    except Exception as e:
            logger.exception("Error raised on creating persistence object")
            raise

    if args.rebuild_clusters_debug:
        logger.warning("Wiping existing clustering data as --rebuild_clusters_debug is set")
        PERSIST._delete_existing_clustering_data()
        logger.warning("Wiped existing clustering data")
    else:
        logger.info("Working with existing data ... ")

    ################################# clustering #############################################
    # open PERSIST and hybridComparer object used by all samples
    # this is only used for data access and msa.
    # inserts are not allowed
    logger.info("Building hybridComparer object")
    hc = hybridComparer(reference=CONFIG['reference'],
        maxNs=CONFIG['MAXN_STORAGE'],
        snpCeiling=  CONFIG['SNPCEILING'],
        excludePositions=CONFIG['excluded'],
        preComparer_parameters=CONFIG['PRECOMPARER_PARAMETERS'],
        PERSIST=PERSIST,
        disable_insertion = True)

    # get a clustering object's settings
    logger.info("Creating clustering objects ...")
    clusterers = {}
    clusternomenclature = {}
    clusternameassigner= {}
    for clustering_name in CONFIG['CLUSTERING'].keys():
        clustering_setting = CONFIG['CLUSTERING'][clustering_name]
        
        mpmc = MixPOREMixtureChecker(hc, **clustering_setting)  # uses hybridComparer to load samples and compute msas
  
        # check update adds remaining guids
        logger.info("Creating clustering object {0}".format(clustering_name))
        clusterers[clustering_name] = MixtureAwareLinkage(PERSIST=PERSIST, 
                    MIXCHECK = mpmc,
                    mixed_sample_management = clustering_setting['mixed_sample_management'], 
                    snv_threshold=clustering_setting['snv_threshold'],
                    serialisation=None,
                    parameters= clustering_setting,
                    name = clustering_name)
        clusterers[clustering_name].remove_legacy()     # remove any old versions
        logger.info("Created clustering object {0}".format(clustering_name))
        # if applicable, make a cluster nomenclature object;
        if 'cluster_nomenclature_method' in clusterers[clustering_name].parameters:
            if not relabel:     # otherwise, we use the labels provided to us
                existing_labels = clusterers[clustering_name].existing_labels()
            else:
                logging.info("Using pre-supplied existing labels")
            # we are instructed to do cluster naming
            clusternomenclature[clustering_name] = ClusterNomenclature(
                cluster_nomenclature_method = clusterers[clustering_name].parameters['cluster_nomenclature_method'],
                existing_labels = existing_labels)
            clusternameassigner[clustering_name] = ClusterNameAssigner(clusternomenclature[clustering_name])
            logger.info("Created name assigner {0} with method {1}".format(clustering_name, clusterers[clustering_name].parameters['cluster_nomenclature_method']))
        
    # now iterate - on a loop
    while True:
        whitelist= set()
        nbuilt=0
        for clustering_name in CONFIG['CLUSTERING'].keys():         # DEBUG
            clustering_setting = CONFIG['CLUSTERING'][clustering_name]
            clusterers[clustering_name].update()    
            clusterers[clustering_name].cluster()

            # if applicable, generate and store cluster labels
            if 'cluster_nomenclature_method' in clusterers[clustering_name].parameters:
                # we are instructed to do cluster naming
                cluster2guid = clusterers[clustering_name].cluster2names
                if not relabel:     # then we use the existing labels as our source, otherwise we use the labels already computed
                    previous_guid2cluster_label = clusterers[clustering_name].guid2cluster_labels()
                else:
                    logging.info("Using pre-supplied guid to cluster lookup")
                clusterid2clusterlabel = clusternameassigner[clustering_name].assign_new_clusternames(
                                clusterid2guid = clusterers[clustering_name].cluster2names,
                                previous_guid2cluster_label = previous_guid2cluster_label
                            )
                clusterers[clustering_name].apply_cluster_labels(clusterid2clusterlabel)
                
            # store the output
            clusterers[clustering_name].persist(what='graph')
            clusterers[clustering_name].persist(what='output')
            clusterers[clustering_name].remove_legacy()     # remove old versions

            # read the result 
            malr = MixtureAwareLinkageResult(PERSIST=PERSIST, name=clustering_name)   

            ms = MSAStore(PERSIST=PERSIST, in_ram_persistence_time=60)      # persist 60 seconds
            # estimate expected:
            estimated_unk = hc.estimate_expected_unk(sample_size=100, unk_type= clusterers[clustering_name].parameters['uncertain_base_type'])
            if estimated_unk is not None:
                estimated_p1 = estimated_unk / (len(hc.reference)-len(hc.excluded))
            else:   
                estimated_p1 = None     
            # recover existing msas
            stored_msa = ms.existing_tokens()

            # build multisequence alighments
            logger.info("Precomputing clusters for {0}".format(clustering_name))
            cluster_contents=malr.cluster2guid.values()
            bar = progressbar.ProgressBar(max_value=len(cluster_contents))
            
            for i,guids in enumerate(cluster_contents):
    
                bar.update(i+1)
                # token identifies cluster contents, nature of analysis, and outgroup
                token = ms.get_token(malr.parameters['uncertain_base_type'],False, guids)           
                if len(guids)>2:
                    whitelist.add(token)        # we need to retain this msa, if it exists
                    
                    if not token in stored_msa:     # if we haven't already computed it 
                        msa_result = hc.multi_sequence_alignment(guids, expected_p1=estimated_p1, uncertain_base_type= malr.parameters['uncertain_base_type'])
                        ms.persist(token, msa_result)   
                        nbuilt+=1

                    
            bar.finish()


        # cleanup anything we don't need, including old clustering versions and msas
        ms.unpersist(whitelist=whitelist)
        logger.info("Cleanup complete.  Stored data on {0} MSAs; Built {1} new clusters".format(len(whitelist), nbuilt))
    
        if debugmode:
            exit(0)
        if relabel:
            exit(0)
        logger.info("Waiting 60 seconds")
        time.sleep(60)
