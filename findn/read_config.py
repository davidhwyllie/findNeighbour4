#!/usr/bin/env python
""" classes reading the findNeighbour config file 

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
import os
import unittest
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class ReadConfig():
    """ reads, and where approporiate modified from environmental variables containing secret 
	configuration parameters etc, a findNeighbour config string.
    also reads from reference files, where appropriate.

    CONFIG files are json files with fomrat similar to the below:


		An example CONFIG is below:
		
		{			
		"DESCRIPTION":"A test server operating in ../unittest_tmp, only suitable for testing",
		"IP":"127.0.0.1",
		"INPUTREF":"../reference/TB-ref.fasta",
		"EXCLUDEFILE":"../reference/TB-exclude.txt",
		"DEBUGMODE":0,
		"SERVERNAME":"TBSNP",
		"FNPERSISTENCE_CONNSTRING":"mongodb://127.0.0.1",
		"MAXN_STORAGE":100000,
		"MAXN_PROP_DEFAULT":0.70,
		"PRECOMPARER_PARAMETERS":{},
		"LOGFILE":"../unittest_tmp/logfile.log",
		"LOGLEVEL":"INFO",
		"SNPCEILING": 20,
		"SERVER_MONITORING_MIN_INTERVAL_MSEC":0,
		"SENTRY_URL":"https://c******************@sentry.io/1******",
		"CLUSTERING":{'SNV12_ignore' :{'snv_threshold':12, 'mixed_sample_management':'ignore', 'mixture_criterion':'pvalue_1', 'cutoff':0.001},
					  'SNV12_include':{'snv_threshold':12, 'mixed_sample_management':'include', 'mixture_criterion':'pvalue_1', 'cutoff':0.001}
					 },
		"LISTEN_TO":"127.0.0.1"
		}
		
		where the database connection binds to
		FNPERSISTENCE_CONNSTRING
		Note: if a FNPERSISTENCE_CONNSTRING environment variable is present, then the value of this will take precedence over any values in the config file.
		This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuration file.
		
		related to what monitoring the server uses
		SERVER_MONITORING_MIN_INTERVAL_MSEC (optional)
		
		related to error handling
		SENTRY_URL (optional)
		Note: if a FN_SENTRY URL environment variable is present, then the value of this will take precedence over any values in the config file.
		This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
		PERSIST is a storage object needs to be supplied.  The fn3Persistence class in mongoStore is one suitable object.
		PERSIST=fn3persistence(connString=CONFIG['FNPERSISTENCE_CONNSTRING'])

    """

    def __init__(self):
        """ initializes the object """
        pass

    def read_config(self, configFile):
        """ read a config file 

        input:      configFile, a path to a configuration file
        returns:    the config file as a dictionary"""

        try:
            with open(configFile,'r') as f:
                CONFIG=f.read()

        except FileNotFoundError:
            raise FileNotFoundError("Passed a positional parameter, which should be a CONFIG file name; tried to open a config file at {0} but it does not exist ".format(configFile))

        if isinstance(CONFIG, str):
            CONFIG=json.loads(CONFIG)	# assume JSON string; convert.

        # check CONFIG is a dictionar
        if not isinstance(CONFIG, dict):
            raise KeyError("CONFIG must be either a dictionary or a JSON string encoding a dictionary.  It is: {0}".format(CONFIG))
	    
        # check that the keys of config are as expected.
        required_keys=set(['IP', 'REST_PORT', 'DEBUGMODE', 'LOGFILE', 'MAXN_PROP_DEFAULT'])
        missing=required_keys-set(CONFIG.keys())
        if not missing == set([]):
            raise KeyError("Required keys were not found in CONFIG. Missing are {0}".format(missing))

        # determine whether a FNPERSISTENCE_CONNSTRING environment variable is present,
        # if so, the value of this will take precedence over any values in the config file.
        # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
        if os.environ.get("FNPERSISTENCE_CONNSTRING") is not None:
            CONFIG["FNPERSISTENCE_CONNSTRING"] = os.environ.get("FNPERSISTENCE_CONNSTRING")
            print("Set mongodb connection string  from environment variable")
        else:
            print("Using mongodb connection string from configuration file.")

        # determine whether a FN_SENTRY_URLenvironment variable is present,
        # if so, the value of this will take precedence over any values in the config file.
        # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
        if os.environ.get("FN_SENTRY_URL") is not None:
            CONFIG["SENTRY_URL"] = os.environ.get("FN_SENTRY_URL")
            print("Set Sentry connection string from environment variable")
        else:
            print("Using Sentry connection string from configuration file.")

        excluded=set()
        if CONFIG['EXCLUDEFILE'] is not None:
            with open(CONFIG['EXCLUDEFILE'],'rt') as f:
	            rows=f.readlines()
            for row in rows:
	            excluded.add(int(row))
        CONFIG['excluded']=excluded

        # load reference
        with open(CONFIG['INPUTREF'],'rt') as f:
	        for r in SeqIO.parse(f,'fasta'):
		        CONFIG['reference']=str(r.seq)

        return CONFIG

class Test_RC_1(unittest.TestCase):
    """ tests the ReadConfig class"""
    def runTest(self):

        rc = ReadConfig()
        res = rc.read_config("../config/default_test_config.json")

        self.assertTrue(isinstance(res,dict))
 
  
