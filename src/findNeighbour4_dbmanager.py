 #!/usr/bin/env python
""" 
findNeighbour3 is a server providing relatedness information for bacterial genomes via a Restful API.

This application optimises the database which findNeighbour3 uses, in the background.
Its use is not essential, but running this software will
* reduce the number of documents
* consequently, reduce the index size
* make the database more scalable.

36The main operation this software provides is to repack the database from a format in which there is one document per cell in the snp matrix
to one in which there is one document per row (or, if the matrix is very big, to a small number or rows)

It is designed to be run from the command line by a scheduler, such as cron.
It will repack up to 10 samples per run
"""
 
# import libraries
import os
import sys
import logging
import logging.handlers
import warnings
import sys
import pymongo
import pathlib
import sentry_sdk
import json
import time
import random
import datetime 
import hashlib
import argparse

# logging
from logging.config import dictConfig

# reference based compression, storage and clustering modules
from mongoStore import fn3persistence

# only used for unit testing
import unittest
            

# startup
if __name__ == '__main__':

        # command line usage.  Pass the location of a config file as a single argument.
        # an example config file is default_test_config.json
            # command line usage.  Pass the location of a config file as a single argument.
        parser = argparse.ArgumentParser(
        formatter_class= argparse.RawTextHelpFormatter,
        description="""Runs findNeighbour4_dbmanager, part of the findNeighbour4 server system.

        findNeighbour4_server will run without this, but it it recommended that one or more (see below) processes of findNeighbour4_dbmanager should be 
        running as well.  The findNeighbour4_dbmanager is entirely independent of findNeighbour4_server, and does not require that the findNeighbour4_server
        be running.

        The prupose of the database manager is as follows:
        * It rotates log files for all findNeighbour4 applications
        * It periodically releases database space back to the operating system
        * It rearranges the database continuously, to ensure that neighbours are stored in  the most efficient way.  For details on how it does this, see the mongoStore.repack() method.                                     

Example usage: 
============== 
# show command line options 
python findNeighbour4_dbmanager.py --help  

# run using settings in myConfigFile.json.   
python findNeighbour4_dbmanager.py ../config/myConfigFile.json          # this config file should be the same as that used by the server.  recompress everything needing it.    

python findNeighbour4_server.py ../config/myConfigFile.json \ 
                        --recompress_subset 0123456789abcdef            # effects same as above

# two processes, which can be run simultaneously; they will not recompress the same samples
python findNeighbour4_server.py ../config/myConfigFile.json \ 
                        --recompress_subset 01234567            # one process only compresses samples if the hash of the sample_id begins with 0,1,..7
python findNeighbour4_server.py ../config/myConfigFile.json \ 
                        --recompress_subset 89abcdef            # one process only compresses samples if the hash of the sample_id begins with 8,9, a, ... f

""")
        parser.add_argument('path_to_config_file', type=str, action='store', nargs='?',
                                help='the path to the configuration file', default='')
        parser.add_argument('--recompress_subset', type=str, nargs=1, action='store', default='abcdef0123456789', 
                                help='only compress samples whose hash begins with one of these letters.  Set this if you are using multiple dbmanagers, to make sure they do different work.')
        args = parser.parse_args()
        
        # an example config file is default_test_config.json

                
        ############################ check input ######################################
        print("findNeighbour4_dbmanager server .. reading configuration file.")

        if len(args.path_to_config_file) > 0 :
                configFile = args.path_to_config_file
        else:
                configFile = os.path.join('..','config','default_test_config.json')
                warnings.warn("No config file name supplied ; using a configuration ('default_test_config.json') suitable only for testing, not for production. ")
   
        # the manager can only compress a fraction of the samples' records.  we provide functionality to select 1/16ths of the samples for this instance to compress.
        # this is useful if you want to run multiple dbmanagers.
        # the systems works by taking the hash of the sample name, and then examining the first character of the hexadecimal representation of the hash.
        # the hash yields ~ random first letters, but is deterministic, so different processes hashing the same sample name will get the same answers.
        recompress_subset = 'abcdef0123456789'                # this is the default : this instance analyses all samples names
        if isinstance(args.recompress_subset, list):          # we were passed the parameters
                if not len(args.recompress_subset) == 1:      # check validity
                        raise TypeError("Must pass a single parameter to recompress_subset, if you are passing it.")
                else:
                        recompress_subset = args.recompress_subset[0]   # list element to string
        recompress_subset_chars = sorted(list(recompress_subset))    # split up the string into its components
        recompress_subset_label = 'recompress_all'            # make a label for output in the log
        if not recompress_subset  == 'abcdef0123456789' :
                recompress_subset_label = "recompress_subset_{0}_to_{1}".format(recompress_subset_chars[0], recompress_subset_chars[-1])
       
        ############################ load config ######################################
        # open the config file
        try:
                with open(configFile,'r') as f:
                         CONFIG=f.read()

        except FileNotFoundError:
                raise FileNotFoundError("Passed one parameter, which should be a CONFIG file name; tried to open a config file at {0} but it does not exist ".format(sys.argv[1]))

        if isinstance(CONFIG, str):
                CONFIG=json.loads(CONFIG)   # assume JSON string; convert.

        # check CONFIG is a dictionary  
        if not isinstance(CONFIG, dict):
                raise KeyError("CONFIG must be either a dictionary or a JSON string encoding a dictionary.  It is: {0}".format(CONFIG))
        
        # check that the keys of config are as expected.

        required_keys=set(['IP', 'REST_PORT', 'DEBUGMODE', 'LOGFILE', 'MAXN_PROP_DEFAULT'])
        missing=required_keys-set(CONFIG.keys())
        if not missing == set([]):
                raise KeyError("Required keys were not found in CONFIG. Missing are {0}".format(missing))

        ########################### set up logging #####################################  
        # create a log file if it does not exist.
        print("Starting logging")
        logdir = os.path.dirname(CONFIG['LOGFILE'])
        pathlib.Path(os.path.dirname(logdir)).mkdir(parents=True, exist_ok=True)

        # set up logger
        loglevel=logging.INFO
        if 'LOGLEVEL' in CONFIG.keys():
                if CONFIG['LOGLEVEL']=='WARN':
                        loglevel=logging.WARN
                elif CONFIG['LOGLEVEL']=='DEBUG':
                        loglevel=logging.DEBUG

        # configure logging object
        logger = logging.getLogger()
        logger.setLevel(loglevel)
        logfile = os.path.join(logdir, "dbmanager-{0}".format(os.path.basename(CONFIG['LOGFILE'])))
        print("Logging to {0} with rotation".format(logfile))
        file_handler = logging.handlers.RotatingFileHandler(logfile, mode = 'a', maxBytes = 1e7, backupCount = 7)
 
        formatter = logging.Formatter( "%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s ")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        logger.addHandler(logging.StreamHandler())
 

        ######################### launch sentry if API key provided ##################################
         # determine whether a FN_SENTRY_URLenvironment variable is present,
        # if so, the value of this will take precedence over any values in the config file.
        # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
        if os.environ.get("FN_SENTRY_URL") is not None:
                CONFIG["SENTRY_URL"] = os.environ.get("FN_SENTRY_URL")
                print("Set Sentry connection string from environment variable")
        else:
                print("Using Sentry connection string from configuration file.")

        if 'SENTRY_URL' in CONFIG.keys():
                logger.info("Launching logger")
                sentry_sdk.init(CONFIG['SENTRY_URL'])


        # determine whether a FN_SENTRY_URL environment variable is present,
        # if so, the value of this will take precedence over any values in the config file.
        # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
        if os.environ.get("FN_SENTRY_URL") is not None:
                CONFIG["SENTRY_URL"] = os.environ.get("FN_SENTRY_URL")
                print("Set Sentry connection string from environment variable")
        else:
                print("Using Sentry connection string from configuration file.")      

        ######################### launch mongodb connection ##################################
        # determine whether a FNPERSISTENCE_CONNSTRING environment variable is present,
        # if so, the value of this will take precedence over any values in the config file.
        # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuration file.
        if os.environ.get("FNPERSISTENCE_CONNSTRING") is not None:
                CONFIG["FNPERSISTENCE_CONNSTRING"] = os.environ.get("FNPERSISTENCE_CONNSTRING")
                print("Set connection string to mongodb from environment variable")
        else:
                print("Using connection string to mongodb from configuration file.")


        ######################################################## ######################


        ########################  START Operations ###################################
        logger.info("Connect to database")
        try:
             PERSIST=fn3persistence(dbname = CONFIG['SERVERNAME'],
                                    connString=CONFIG['FNPERSISTENCE_CONNSTRING'],
                                    debug=CONFIG['DEBUGMODE'],
                                    server_monitoring_min_interval_msec = 0)
        except Exception as e:
             logger.exception("Error raised on creating persistence object")
             if e.__module__ == "pymongo.errors":
                 logger.info("Error raised pertains to pyMongo connectivity")
                 raise    

        date_last_log_rotated =  datetime.datetime.now()-datetime.timedelta(hours=25)           # force log rotation on startup  
         
        while True:
             if datetime.datetime.now() > date_last_log_rotated+datetime.timedelta(hours=24):
                     date_last_log_rotated =datetime.datetime.now()
                     logger.info("Rotated mongodb log; deleting old server monitor entries")
                     try:
                            PERSIST.rotate_log()
                     except pymongo.errors.OperationFailure as e:
                             capture_exception(e)
                             logger.warning("Failed to rotate log")                     # can occur if multiple threads try to do it simultenously

                     PERSIST.delete_server_monitoring_entries(before_seconds= (3600 * 24 * 7))        # 7 days
                     pass

             nModified = 0
             # does this guid have any singleton guid2neighbour records which need compressing
             logger.info("(instance = {0}) | Finding compressable records.".format(recompress_subset_label))            
             to_update =PERSIST.singletons(method='approximate', return_top=100)                     # compress in small batches of 100, after which we will reassess what is worst.  

             logger.info("(instance = {1}) | Found records which may need compression; examining {0} with higher singleton estimates.".format(len(to_update.index), recompress_subset_label))

             for guid in to_update.index:
                hashed_guid = hashlib.md5(guid.encode('utf-8')).hexdigest()
                if hashed_guid[0] in recompress_subset_chars:
                        audit_stats = PERSIST.guid2neighbour_repack(guid, always_optimise=False, min_optimisable_size=1)         
                        logger.info("(instance = {1}) | Repacked {0} : {2}".format(guid,recompress_subset_label, audit_stats))
                        if audit_stats['actions_taken']['n_written'] > 0:
                                nModified+=1
                
                        # log database size
                        if nModified % 25 == 0:
                                db_summary = PERSIST.summarise_stored_items()
                                PERSIST.server_monitoring_store(what='dbManager', message="Repacking", guid='-', content=db_summary)

             if nModified == 0:
                # everything has been packed
                logger.info("(instance = {0}) | Nothing found to repack within the specimen number range for this dbmonitor.  Waiting 1 min .. ".format(recompress_subset_label))
                time.sleep(60)    # recheck later
             # otherwise go get some more to compress
         


