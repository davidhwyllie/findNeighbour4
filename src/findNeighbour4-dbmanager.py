 #!/usr/bin/env python
""" 
findNeighbour3 is a server providing relatedness information for bacterial genomes via a Restful API.

This application optimises the database which findNeighbour3 uses, in the background.
Its use is not essential, but running this software will
* reduce the number of documents
* consequently, reduce the index size
* make the database more scalable.

The main operation this software provides is to repack the database from a format in which there is one document per cell in the snp matrix
to one in which there is one document per row (or, if the matrix is very big, to a small number or rows)

It is designed to be run from the command line by a scheduler, such as cron.
It will repack up to 10 samples per run
"""
 
# import libraries
import os
import sys
import logging
import warnings
import sys
import pymongo
import pathlib
import sentry_sdk
import json
import time
import random

# logging
from logging.config import dictConfig

# reference based compression, storage and clustering modules
from mongoStore import fn3persistence

# only used for unit testing
import unittest
			
def repack(guids):
	""" generates a smaller and faster representation in the persistence store
	for the guids in the list. optional"""
	if guids is None:
		guids = PERSIST.guids()  # all the guids
	for this_guid in guids:
		logger.debug("Repacking {0}".format(this_guid))
		PERSIST.guid2neighbour_repack(this_guid)
		

# startup
if __name__ == '__main__':

        # command line usage.  Pass the location of a config file as a single argument.
        # an example config file is default_test_config.json
               
        ############################ LOAD CONFIG ######################################
        print("findNeighbour3-dbmanager server .. reading configuration file.")

        max_batch_size = 100

        if len(sys.argv) == 2:
                configFile = sys.argv[1]
        else:
                configFile = os.path.join('..','config','default_test_config.json')
                warnings.warn("No config file name supplied ; using a configuration ('default_test_config.json') suitable only for testing, not for production. ")
   
        # open the config file
        try:
                with open(configFile,'r') as f:
                         CONFIG=f.read()

        except FileNotFoundError:
                raise FileNotFoundError("Passed one parameter, which should be a CONFIG file name; tried to open a config file at {0} but it does not exist ".format(sys.argv[1]))

        if isinstance(CONFIG, str):
                CONFIG=json.loads(CONFIG)	# assume JSON string; convert.

        # check CONFIG is a dictionary	
        if not isinstance(CONFIG, dict):
                raise KeyError("CONFIG must be either a dictionary or a JSON string encoding a dictionary.  It is: {0}".format(CONFIG))
        
        # check that the keys of config are as expected.

        required_keys=set(['IP', 'REST_PORT', 'DEBUGMODE', 'LOGFILE', 'MAXN_PROP_DEFAULT'])
        missing=required_keys-set(CONFIG.keys())
        if not missing == set([]):
                raise KeyError("Required keys were not found in CONFIG. Missing are {0}".format(missing))

        ########################### SET UP LOGGING #####################################  
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
        file_handler = logging.FileHandler(logfile)
        formatter = logging.Formatter( "%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s ")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        logger.addHandler(logging.StreamHandler())
 

        # launch sentry if API key provided
        if 'SENTRY_URL' in CONFIG.keys():
                logger.info("Launching logger")
                sentry_sdk.init(CONFIG['SENTRY_URL'])

        # set min logging interval if not supplied
        if not 'SERVER_MONITORING_MIN_INTERVAL_MSEC' in CONFIG.keys():
               CONFIG['SERVER_MONITORING_MIN_INTERVAL_MSEC']=0

        # determine whether a FNPERSISTENCE_CONNSTRING environment variable is present,
        # if so, the value of this will take precedence over any values in the config file.
        # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuration file.
        if os.environ.get("FNPERSISTENCE_CONNSTRING") is not None:
                CONFIG["FNPERSISTENCE_CONNSTRING"] = os.environ.get("FNPERSISTENCE_CONNSTRING")
                print("Set connection string to mongodb from environment variable")
        else:
                print("Using connection string to mongodb from configuration file.")

        # determine whether a FN_SENTRY_URL environment variable is present,
        # if so, the value of this will take precedence over any values in the config file.
        # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
        if os.environ.get("FN_SENTRY_URL") is not None:
                CONFIG["SENTRY_URL"] = os.environ.get("FN_SENTRY_URL")
                print("Set Sentry connection string from environment variable")
        else:
                print("Using Sentry connection string from configuration file.")      
        #########################  CONFIGURE HELPER APPLICATIONS ######################


        ########################  START Operations ###################################
        logger.info("Collecting samples for repacking")

        logging.info("Connecting to backend data store")
        try:
             PERSIST=fn3persistence(dbname = CONFIG['SERVERNAME'],
									connString=CONFIG['FNPERSISTENCE_CONNSTRING'],
									debug=CONFIG['DEBUGMODE'],
									server_monitoring_min_interval_msec = CONFIG['SERVER_MONITORING_MIN_INTERVAL_MSEC'])
        except Exception as e:
             logger.exception("Error raised on creating persistence object")
             if e.__module__ == "pymongo.errors":
                 logger.info("Error raised pertains to pyMongo connectivity")
                 raise        
        while True:
             nModified = 0
             to_update = set()
             # does this guid have any singleton guid2neighbour records which need compressing
             print("Gathering guids for update .. ")
             for res in PERSIST.db.guid2neighbour.find({'rstat':'s'}):                   
                 to_update.add(res['guid'])
                 if len(to_update)>max_batch_size:
                     break
             to_update = list(to_update)
             random.shuffle(to_update)
             for guid in to_update:
                     logger.info("Repacking {0}".format(guid))
                     repack([guid])
                     nModified += 1

                     # log database size
                     db_summary = PERSIST.summarise_stored_items()
                     PERSIST.server_monitoring_store(what='dbManager', message="Repack one", guid=guid, content=db_summary)

             if nModified == 0:
                     # everything has been packed
                     logger.info("Nothing found to repack.  Waiting 60s .. ")
                     time.sleep(60)	# recheck in 1 minute
		 


