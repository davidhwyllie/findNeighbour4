#!/usr/bin/env python3
""" produces depictions of server activity, based on findNeighbour4 monitoring data 

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.  See see <https://www.gnu.org/licenses/>.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

"""

# import libraries
import os
import sys
import argparse
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

# config
from findn.common_utils import ConfigManager

# fn3 storage module
from findn.mongoStore import fn3persistence
from findn.depictStatus import MakeHumanReadable, DepictServerStatus

# startup
if __name__ == '__main__':
        parser = argparse.ArgumentParser(
                formatter_class= argparse.RawTextHelpFormatter,
                description="""Runs findNeighbour4_monitor, a service which every 2 minutes makes an report on server status. This report is available from the server."""
                )                
        parser.add_argument('path_to_config_file', type=str, action='store', nargs='?',
                                help='the path to the configuration file', default=None)
        parser.add_argument('--run_once_only', action='store_true', 
                                help='just produce one depiction; do not keep running.  Default False.  Mainly useful for unit testing.')
        args = parser.parse_args()


        ############################ LOAD CONFIG ######################################

        config_file = args.path_to_config_file
        if config_file is None:
                config_file =  os.path.join('config','default_test_config.json')
                debugmode = 1
                warnings.warn("No config file name supplied ; using a configuration ('default_test_config.json') suitable only for testing, not for production. ")
        else:
                debugmode = 0
        cfm = ConfigManager(config_file)  
        CONFIG = cfm.read_config()
                    

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
                logger.info("Launching sentry client")

        #########################  CONFIGURE HELPER APPLICATIONS ######################


        ########################  START Operations ###################################
        logger.info("Preparing to produce visualisations")

        logger.info("Connecting to backend data store at {0}".format(CONFIG['SERVERNAME']))
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

        processes = ['server','dbmanager','clustering']
        logfiles = {}
        for process in processes:
                logfiles[process]= os.path.join(logdir, "{0}-{1}".format(process, os.path.basename(CONFIG['LOGFILE'])))


        dss1 = DepictServerStatus(logfiles= logfiles,
                                        server_url=CONFIG['IP'],
                                        server_port=CONFIG['REST_PORT'],
                                        server_description=CONFIG['DESCRIPTION'])
        while True:
                insert_data = PERSIST.recent_server_monitoring(selection_field="context|info|message", selection_string="About to insert", max_reported=500)
                recent_data = PERSIST.recent_server_monitoring(selection_field="content|activity|whatprocess", selection_string="server", max_reported=100)
                db_data = PERSIST.recent_server_monitoring(selection_field="content|activity|whatprocess", selection_string="dbManager", max_reported=1000)

                page_content = dss1.make_report(insert_data, recent_data, db_data)
                for item in page_content.keys():       
                        html = file_html(page_content[item], CDN, item)
                        PERSIST.monitor_store(item, html) 

                # with open("test.html",'wt') as f:
                #    f.write(html)
                        
                if args.run_once_only:
                        logger.info("Exiting as run_once_only specified")
                        exit(0)

                        
                if debugmode == 1:
                        logger.info("Exiting as debugmode")
                        exit(0)

                time.sleep(120)	# rerun in 2 minutes



