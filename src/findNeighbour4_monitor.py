#!/usr/bin/env python3
""" produces depictions of server activity, based on findNeighbour monitoring data """

# import libraries
import os
import sys
import logging
import logging.handlers
import warnings
import pymongo
import pandas as pd
import numpy as np
import pathlib
import sentry_sdk
import json
import time
import random
import dateutil.parser
import datetime
import unittest
from bokeh.embed import file_html
from bokeh.plotting import show
from bokeh.resources import CDN

# fn3 storage module
from mongoStore import fn3persistence
from depictStatus import MakeHumanReadable, DepictServerStatus

# startup
if __name__ == '__main__':

        # command line usage.  Pass the location of a config file as a single argument.
        # an example config file is default_test_config.json
               
        ############################ LOAD CONFIG ######################################
        print("findNeighbour4_monitor .. reading configuration file.")

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
        pathlib.Path(os.path.dirname(CONFIG['LOGFILE'])).mkdir(parents=True, exist_ok=True)

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
        logfile = os.path.join(logdir, "monitor-{0}".format(os.path.basename(CONFIG['LOGFILE'])))
        print("Logging to {0} with rotation".format(logfile))
        file_handler = logging.handlers.RotatingFileHandler(logfile, mode = 'a', maxBytes = 1e7, backupCount = 7)
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
                logger.info("Set connection string to mongodb from environment variable")
        else:
                logger.info("Using connection string to mongodb from configuration file.")
        
        
        # determine whether a FN_SENTRY_URL environment variable is present,
        # if so, the value of this will take precedence over any values in the config file.
        # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
        if os.environ.get("FN_SENTRY_URL") is not None:
                CONFIG["SENTRY_URL"] = os.environ.get("FN_SENTRY_URL")
                logger.info("Set Sentry connection string from environment variable")
        else:
                logger.info("Using Sentry connection string from configuration file.")
                
                
        #########################  CONFIGURE HELPER APPLICATIONS ######################


        ########################  START Operations ###################################
        logger.info("Preparing to produce visualisations")

        print("Connecting to backend data store at {0}".format(CONFIG['SERVERNAME']))
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
        dss1 = DepictServerStatus(logfile= CONFIG['LOGFILE'],
                                    server_url=CONFIG['IP'],
                                    server_port=CONFIG['REST_PORT'],
                                    server_description=CONFIG['DESCRIPTION'])
        while True:
            insert_data = PERSIST.recent_server_monitoring(selection_field="context|info|message", selection_string="About to insert", max_reported=500)
            recent_data = PERSIST.recent_server_monitoring(selection_field="content|activity|whatprocess", selection_string="server", max_reported=100)
            page_content = dss1.make_report(insert_data, recent_data)
            for item in page_content.keys():       
                html = file_html(page_content[item], CDN, item)
                PERSIST.monitor_store(item, html) 
 
                # with open("test.html",'wt') as f:
                #    f.write(html)
                  
            time.sleep(120)	# rerun in 2 minutes



