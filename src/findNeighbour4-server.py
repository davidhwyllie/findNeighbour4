#!/usr/bin/env python
""" 
A server providing relatedness information for bacterial genomes via a Restful API.

Implemented in pure Python3, it uses in-memory data storage backed by MongoDb.
It loads configuration from a config file, which must be set in production.

If no config file is provided, it will run in  'testing' mode with the  parameters
in default_test_config.json.  This expects a mongodb database to be running on
the default port on local host.  As a rough guide to the amount of space required in mongodb,
about 0.5MB of database is used per sequence, or about 2,000 sequences per GB.

If no config file is provided, it will run in  'testing' mode with the  parameters
in default_test_config.json.  This expects a mongodb database to be running on
the default port on local host.  As a rough guide to the amount of space required in mongodb,
about 0.5MB of database is used per sequence, or about 2,000 sequences per GB.

All internal modules, and the restful API, are covered by unit testing.
Unit testing can be achieved by:

# starting a test RESTFUL server
python3 findNeighbour4-server.py

# And then (e.g. in a different terminal) launching unit tests with
python3 -m unittest findNeighbour4-server
"""

# import libraries
import os
import sys
import requests
import json
import logging
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
import networkx as nx
import progressbar
from sentry_sdk import capture_message, capture_exception
from sentry_sdk.integrations.flask import FlaskIntegration

# flask
from flask import Flask, make_response, jsonify, Markup
from flask import request, abort, send_file
from flask_cors import CORS		# cross-origin requests are not permitted except for one resource, for testing

# logging
from logging.config import dictConfig

# utilities for file handling and measuring file size
import psutil

# reference based compression, storage and clustering modules
from NucleicAcid import NucleicAcid
from mongoStore import fn3persistence
from hybridComparer import hybridComparer	
from clustering import snv_clustering
from guidLookup import guidSearcher  		# fast lookup of first part of guids
from ma_linkage import MixtureAwareLinkageResult

# network visualisation
from visualiseNetwork import snvNetwork

# server status visualisation
from depictStatus import MakeHumanReadable

# only used for unit testing
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide
import unittest
from urllib.parse import urlparse as urlparser
from urllib.parse import urljoin as urljoiner
import uuid
import time

class findNeighbour4():
	""" a server based application for maintaining a record of bacterial relatedness using SNP distances.
	
		The high level arrangement is that
		- This class interacts with  sequences, partly held in memory
		  [handled by the hybridComparer class] 
		- cached data is accessed by an fn3persistence object
		- methods in findNeighbour4() return native python3 objects.
		
		- a web server, currently flask, handles the inputs and outputs of this class
		- in particular, native python3 objects returned by this class are serialised by the Flask web server code.
		"""
		
	def __init__(self,CONFIG, PERSIST):
		""" Using values in CONFIG, starts a server with CONFIG['NAME'] on port CONFIG['PORT'].

		CONFIG contains Configuration parameters relevant to the reference based compression system which lies
		at the core of the server.
			INPUTREF:       the path to fasta format reference file.
			EXCLUDEFILE:    a file containing the zero-indexed positions in the supplied sequences which should be ignored in all cases.
							Typically, this is because the software generating the mapped fasta file has elected not to call these regions,
							in any samples, e.g. because of difficulty mapping to these regions.
							Such regions can occupy up 5- 20% of the genome and it is important for efficient working of this software
							that these regions are supplied for exclusion on sequence loading.  Not doing so will slow loading, and markedly increase
							memory requirements, but will not alter the results produced.
			DEBUGMODE:      Controls operation of the server:

							DEBUGMODE =                                                              0       1        2
							Run server in production mode (errors logged, not returned to client)    Y       N        N
							Run server in debug mode (errors reported to client)                     N       Y        Y
							Create Database if it does not exist                                     Y       Y        Y
							Delete all data on startup                                               N       N        Y
							Enable /restart endpoint, which restarts empty server (for testing)      N       N        Y

			SERVERNAME:     the name of the server. used as the name of mongodb database which is bound to the server.
			PRECOMPARER_PARAMETERS:		parameters passed to precomparer.  ** DOCS NEEDED ON THIS**.  It is important that these values are either set to 'fail always'  -route all precompared samples to the seqcomparer - or are the result of calibration.  preComparer_calibration will do the calibration process automatically.
			FNPERSISTENCE_CONNSTRING: a valid mongodb connection string. if shard keys are set, the 'guid' field is suitable key.
							Note: if a FNPERSISTENCE_CONNSTRING environment variable is present, then the value of this will take precedence over any values in the config file.
							This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
		
			MAXN_STORAGE:   The maximum number of Ns in the sequence <excluding those defined in > EXCLUDEFILE which should be indexed.
							Other files, e.g. those with all Ns, will be tagged as 'invalid'.  Although a record of their presence in the database
							is kept, they are not compared with other sequences.
			MAXN_PROP_DEFAULT: if the proportion not N in the sequence exceeds this, the sample is analysed, otherwise considered invalid.
			LOGFILE:        the log file used
			LOGLEVEL:		default logging level used by the server.  Valid values are DEBUG INFO WARNING ERROR CRITICAL
			SNPCEILING: 	links between guids > this are not stored in the database
			
			REPACK_FREQUENCY: see /docs/repack_frequency.md
			CLUSTERING:		a dictionary of parameters used for clustering.  In the below example, there are two different
							clustering settings defined, one named 'SNV12_ignore' and the other 'SNV12_include.
							{'SNV12_ignore' :{'snv_threshold':12, 'mixed_sample_management':'ignore', 'mixture_criterion':'p_value1', 'cutoff':0.001},
							 'SNV12_include':{'snv_threshold':12, 'mixed_sample_management':'include', 'mixture_criterion':'p_value1', 'cutoff':0.001}
							}
							Each setting is defined by four parameters:
							snv_threshold: clusters are formed if samples are <= snv_threshold from each other
							mixed_sample_management: this defines what happens if mixed samples are detected.
								Suppose there are three samples, A,B and M.  M is a mixture of A and B.
								A and B are > snv_threshold apart, but their distance to M is zero.
								If mixed_sample_management is
								'ignore', one cluster {A,B,M} is returned
								'include', two clusters {A,M} and {B,M}
								'exclude', three clusters are returns {A},{B},{C}
							mixture_criterion: sensible values include 'p_value1','p_value2','p_value3' but other output from  seqComparer._msa() is also possible.
								 these p-values arise from three different tests for mixtures.  Please see seqComparer._msa() for details.
							cutoff: samples are regarded as mixed if the mixture_criterion is less than or equal to this value.
			SENTRY_URL:  optional.  If provided, will launch link Sentry to the flask application using the API key provided.  See https://sentry.io for a description of this service. 

								Note: if a FN_SENTRY_URL environment variable is present, then the value of this will take precedence over any values in the config file.
								This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
					LISTEN_TO:   optional.  If missing, will bind to localhost (only) on 127.0.0.1.  If present, will listen to requests from the IP stated.  if '0.0.0.0', the server will respond to all external requests.
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

		Some of these settings are read when the server is first-run, stored in a database, and the server will not
		change the settings on re-start even if the config file is changed.  Examples are:
		SNPCEILING
		MAXN_PROP_DEFAULT
		EXCLUDEFILE
		INPUTREF
		CLUSTERING
		PRECOMPARER_PARAMETERS
		
		These settings cannot be changed because they alter the way that the data is stored; if you want to change
		the settings, the data will have to be re-loaded. 
		
		However, most other settings can be changed and will take effect on server restart.  These include:
		server location
		IP
		SERVERNAME
		REST_PORT
		LISTEN_TO (optional)
		
		internal logging	
		LOGFILE
		LOGLEVEL
		
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
		
		# store the persistence object as part of this object
		self.PERSIST=PERSIST
		
		# check input
		if isinstance(CONFIG, str):
			self.CONFIG=json.loads(CONFIG)	# assume JSON string; convert.
		elif isinstance(CONFIG, dict):
			self.CONFIG=CONFIG
		else:
			raise TypeError("CONFIG must be a json string or dictionary, but it is a {0}".format(type(CONFIG)))
		
		# check it is a dictionary	
		if not isinstance(self.CONFIG, dict):
			raise KeyError("CONFIG must be either a dictionary or a JSON string encoding a dictionary.  It is: {0}".format(CONFIG))
		
		# check that the keys of config are as expected.
		required_keys=set(['IP','INPUTREF','EXCLUDEFILE','DEBUGMODE','SERVERNAME',
						   'FNPERSISTENCE_CONNSTRING', 'MAXN_STORAGE',
						    "SNPCEILING", 'MAXN_PROP_DEFAULT', 'REST_PORT',
						   'LOGFILE','LOGLEVEL', 'CLUSTERING', "PRECOMPARER_PARAMETERS"])
		missing=required_keys-set(self.CONFIG.keys())
		if not missing == set([]):
			raise KeyError("Required keys were not found in CONFIG. Missing are {0}".format(missing))

		if len(self.CONFIG['PRECOMPARER_PARAMETERS'])>0:
			# these are supplied
			observed_keys = set(list(self.CONFIG['PRECOMPARER_PARAMETERS'].keys()))
			expected_keys = set(["selection_cutoff",
		    					"over_selection_cutoff_ignore_factor",
                    			"mixed_reporting_cutoff",
                    			"N_mean",
                    			"N_sd",
                    			"highN_z_reporting_cutoff",
                    			"alpha",
                    			"probN_inflation_factor",           
                    			"n_positions_examined"])

			missing = expected_keys - observed_keys
			if not missing == set([]):
				raise KeyError("Precomparer parameters were supplied, but the required keys were not found . Missing are {0}".format(missing))

		# the following keys are not stored in any database backend, as a server could be moved, i.e.
		# running on the same data but with different IP etc
		
		do_not_persist_keys=set(['IP',"SERVERNAME",'FNPERSISTENCE_CONNSTRING',
								 'LOGFILE','LOGLEVEL','REST_PORT',
								 'SENTRY_URL', 'SERVER_MONITORING_MIN_INTERVAL_MSEC'])
				
		# determine whether this is a first-run situation.
		if self.PERSIST.first_run():
			self.first_run(do_not_persist_keys)

		# load global settings from those stored at the first run.
		cfg = self.PERSIST.config_read('config')
		
		# set easy to read properties from the config
		self.reference = cfg['reference']
		self.excludePositions = set(cfg['excludePositions'])
		self.debugMode = cfg['DEBUGMODE']
		self.maxNs = cfg['MAXN_STORAGE']
		self.snpCeiling = cfg['SNPCEILING']
		self.maxn_prop_default = cfg['MAXN_PROP_DEFAULT']
		self.clustering_settings = cfg['CLUSTERING']
		
		## start setup
		self.write_semaphore = threading.BoundedSemaphore(1)        # used to permit only one process to INSERT at a time.
		
		# initialise nucleic acid analysis object
		self.objExaminer=NucleicAcid()
		
		# formatting utility
		self.mhr = MakeHumanReadable()
		
		# load in-memory sequences
		self.gs = guidSearcher()
		self._load_in_memory_data()
		
		print("findNeighbour4 is ready.")
	
	def _load_in_memory_data(self):
		""" loads in memory data into the hybridComparer object from database storage """
		
		# set up clustering
		# clustering is performed by findNeighbour4-cluster, and results are stored in the mongodb
		# the clustering object used here is just a reader for the clustering status.
		app.logger.info("findNeighbour4 is loading clustering data.")
		
		self.clustering={}		# a dictionary of clustering objects, one per SNV cutoff/mixture management setting
		for clustering_name in self.clustering_settings.keys():
			self.clustering[clustering_name] = MixtureAwareLinkageResult(
						PERSIST=self.PERSIST, 
						name=clustering_name,
						serialisation = None)
			app.logger.info("set up clustering access object {0}".format(clustering_name))
			
		# initialise hybridComparer, which manages in-memory reference compressed data
		# preComparer_parameters will be read from disc
		self.hc = hybridComparer(reference=self.reference,
							maxNs=self.maxNs,
							snpCeiling= self.snpCeiling,
							excludePositions=self.excludePositions,
							preComparer_parameters={},
                    		PERSIST=self.PERSIST)
		
		app.logger.info("In-RAM data store set up.")
		
		# determine how many guids there in the database
		guids = self.PERSIST.refcompressedsequence_guids()
		if len(guids)==0:
			self.server_monitoring_store(message='There is nothing to load')
			app.logger.info("Nothing to load")
		else:	
			self.server_monitoring_store(message='Starting load of sequences into memory from database')
			app.logger.info("Loading sequences from database")
			
			nLoaded = 0
			nRecompressed = 0
			bar = progressbar.ProgressBar(max_value = len(guids))
			for guid in guids:
				nLoaded+=1
				self.gs.add(guid)
				self.hc.repopulate(guid=guid)
				bar.update(nLoaded)

				if nLoaded % 100 == 0:
					self.server_monitoring_store(message='Server restarting; loaded {0} from database ..'.format(nLoaded))


			bar.finish()
			app.logger.info("findNeighbour4 has finished loaded {0} sequences from database".format(len(guids)))

		self.server_monitoring_store(message='Garbage collection.')
		app.logger.info("Garbage collecting")		
		gc.collect()		# free up ram		
		self.server_monitoring_store(message='Load from database complete.')

	def reset(self):
		""" restarts the server, deleting any existing data """
		if not self.debugMode == 2:
			return		 # no action taken by calls to this unless debugMode ==2
		else:
			print("Deleting existing data and restarting")
			self.PERSIST._delete_existing_data()
			time.sleep(2) # let the database recover
			self._load_in_memory_data()

	def server_monitoring_store(self, message="No message supplied", guid=None):
		""" reports server memory information to store """
		hc_summary = {}
		try:
			hc_summary =  self.hc.summarise_stored_items()
		except AttributeError:		# no hc object, occurs during startup
			pass

		db_summary =  self.PERSIST.summarise_stored_items()
		mem_summary = self.PERSIST.memory_usage()
		self.PERSIST.server_monitoring_store(message=message, what='server', guid= guid, content={**hc_summary, **db_summary, **mem_summary})

	def first_run(self, do_not_persist_keys):
		""" actions taken on first-run only.
		Include caching results from CONFIGFILE to database, unless they are in do_not_persist_keys"""
		
		app.logger.info("First run situation: parsing inputs, storing to database. ")

		# create a config dictionary
		config_settings= {}
		
		# store start time 
		config_settings['createTime']= datetime.datetime.now()
		
		# store description
		config_settings['description']=self.CONFIG['DESCRIPTION']

		# store clustering settings
		self.clustering_settings=self.CONFIG['CLUSTERING']
		config_settings['clustering_settings']= self.clustering_settings

		# store precomparer settings
		if len(self.CONFIG['PRECOMPARER_PARAMETERS'])>0:
			self.PERSIST.config_store('preComparer', self.CONFIG['PRECOMPARER_PARAMETERS'])

		# load the excluded bases
		excluded=set()
		if self.CONFIG['EXCLUDEFILE'] is not None:
			with open(self.CONFIG['EXCLUDEFILE'],'rt') as f:
				rows=f.readlines()
			for row in rows:
				excluded.add(int(row))

		app.logger.info("Noted {0} positions to exclude.".format(len(excluded)))
		config_settings['excludePositions'] = list(sorted(excluded))
		
		# load reference
		with open(self.CONFIG['INPUTREF'],'rt') as f:
			for r in SeqIO.parse(f,'fasta'):
				config_settings['reference']=str(r.seq)

		# persist other config settings.
		for item in self.CONFIG.keys():
			if not item in do_not_persist_keys:
				config_settings[item]=self.CONFIG[item]
				
		res = self.PERSIST.config_store('config',config_settings)
		self.server_monitoring_store(message='First run complete.')

		app.logger.info("First run actions complete.")
	

	def insert(self,guid,dna):
		""" insert DNA called guid into the server,
		persisting it in both RAM and on disc, and updating any clustering.
		"""
		
		# clean, and provide summary statistics for the sequence
		app.logger.info("Preparing to insert: {0}".format(guid))

		if not self.hc.iscachedinram(guid):                   # if the guid is not already there
			self.server_monitoring_store(message='About to insert',guid=guid)
			app.logger.info("Guid is not present: {0}".format(guid))
			
			# insert sequence into the sequence store.
			self.objExaminer.examine(dna)  					  # examine the sequence				
			cleaned_dna=self.objExaminer.nucleicAcidString.decode()
			refcompressedsequence =self.hc.compress(cleaned_dna)          # compress it and store it in RAM
			self.server_monitoring_store(message='Compression complete',guid=guid)
			self.write_semaphore.acquire()				    # addition should be an atomic operation
			
			try:
				loginfo = self.hc.persist(refcompressedsequence, 
						guid,
						{'DNAQuality':self.objExaminer.composition})	
			
			except Exception as e:
				self.write_semaphore.release()                  # release the write semaphore
				app.logger.exception("Error raised on persisting {0}".format(guid))
				app.logger.exception(e)
				capture_exception(e)
				# Rollback anything which could leave system in an inconsistent state
				# remove the guid from RAM is the only step necessary
				self.hc.remove(guid)	
				app.logger.info("Guid {0}  removed from preComparer. {0}".format(guid))
				abort(503,e)		# the mongo server may be refusing connections, or busy.  This is observed occasionally in real-world use
				
			# release semaphore
			self.write_semaphore.release()                  # release the write semaphore
			app.logger.info("Insert succeeded {0}".format(guid))
			for info in loginfo:
				app.logger.info(info)		# performance info
			self.server_monitoring_store(message='Stored compressed sequence',guid=guid)
			return "Guid {0} inserted.".format(guid)		
		else:
			app.logger.info("Already present, no insert needed: {0}".format(guid))

			return "Guid {0} is already present".format(guid)
			
	def exist_sample(self,guid):
		""" determine whether the sample exists in RAM"""
		
		## this call measures presence on disc
		return self.PERSIST.guid_exists(guid)

	def server_time(self):
		""" returns the current server time """
		return {"server_name":self.CONFIG['SERVERNAME'], "server_time":datetime.datetime.now().isoformat()}

	def server_name(self):
		""" returns information about the server """
		return {"server_name":self.CONFIG['SERVERNAME'],
				"server_description":self.CONFIG['DESCRIPTION']
				}
	def server_config(self):
		""" returns the config file with which the server was launched
		
		This may be highly undesirable, and is only available in DEBUG mode.
		as it reveals the internal server architecture  including
		backend databases and perhaps connection strings with passwords.
		"""
		
		if self.debugMode==2:
			return self.CONFIG
		else:
			return None
	def server_nucleotides_excluded(self):
		""" returns the nucleotides excluded by the server """
		# TODO: make this solely dependent on the config file so we don't need memory access
		return {"exclusion_id":self.hc.excluded_hash(), "excluded_nt":list(self.hc.excluded)}
	
	def server_memory_usage(self, max_reported=None):
		""" reports recent server memory activity """
		if max_reported is None:
			max_reported =100		# a default
		return self.PERSIST.recent_server_monitoring(max_reported= max_reported)
	
	def neighbours_within_filter(self, guid, snpDistance, cutoff=0.85, returned_format=1):
		""" returns a list of guids, and their distances, by a sample quality cutoff	
			returns links either as
			format 1 [[otherGuid, distance]]
			or as
			format 2 [[otherGuid, distance, N_just1, N_just2, N_either]]
			or as
			format 3 [otherGuid, otherGuid2, otherGuid3]
			or as
			format 4 [{'guid':otherGuid, 'snv':distance}, {'guid':otherGuid2, 'snv':distance2}]
		"""

		# check the query is of good quality
		inScore = self.PERSIST.guid_quality_check(guid,float(cutoff))
		if inScore == None:
			raise KeyError("{0} not found".format(guid))	# that's an error, maybe should raise KeyError
		elif inScore == False:
			return []		# bad sequence; just to report no links

		# if it is of good quality, then we look for links
		idList=list()

		# gets the similar sequences from the database;
		retVal = self.PERSIST.guid2neighbours(guid=guid, cutoff=snpDistance, returned_format=returned_format)
		
		# run a quality check on the things our sample is like.
		# extract the matching guids, independent of the format requested.
		sampleList=retVal['neighbours']
		idList=[]
		for sa in sampleList:
			if isinstance(sa, list):
				idList.append(sa[0])		# add the guid
			elif isinstance(sa, str):
				idList.append(sa)
			elif isinstance(sa, dict):
				idList.append(sa['guid'])
			else:
				raise TypeError("Unknown format returned {0} {1}".format(type(sa),sampleList))
		
		guid2qual=self.PERSIST.guid2quality(idList)
					  
		# Filter to get good matching guids
		goodGuids=set()
		cutoff=float(cutoff)
		for guid in guid2qual.keys():
			if guid2qual[guid]>=cutoff:
				goodGuids.add(guid)
		
		# note one could put a filter to identify samples based on Ns here: these statistics are return in the sampleList
		
		# assemble output by filtering sampleList
		finalOutput = list()
		for sa in sampleList:
			if isinstance(sa, list):
				guid = sa[0]
			elif isinstance(sa, str):
				guid = sa
			elif isinstance(sa, dict):
				guid = sa['guid']
		
			if guid in goodGuids:
				finalOutput.append(sa)
				
		return finalOutput
	
	def get_all_guids(self):
		return self.PERSIST.guids()
	
	def guids_with_quality_over(self,cutoff=0.66):
		rs=self.PERSIST.guid2propACTG_filtered(float(cutoff))
		if rs==None:
			return []
		else:
			return list(rs.keys())
		
	def get_all_guids_examination_time(self):
		res = self.PERSIST.guid2ExaminationDateTime()
		# isoformat all the keys, as times are not json serialisable
		retDict = res
		for key in retDict:
			retDict[key]=retDict[key].isoformat()
		return(retDict)
	
	def get_all_annotations(self):
		return self.PERSIST.guid_annotations()

	def get_one_annotation(self, guid):
		return self.PERSIST.guid_annotation(guid)
		
	def sequence(self, guid):
		""" gets masked sequence for the guid, in format sequence|fasta """
		if not self.hc.iscachedinram(guid):
			return None
		try:		
			seq = self.hc.uncompress_guid(guid)
			return {'guid':guid, 'invalid':0,'comment':'Masked sequence, as stored','masked_dna':seq}
		except ValueError:
				return {'guid':guid, 'invalid':1,'comment':'No sequence is available, as invalid sequences are not stored'}
			
# default parameters for unit testing only.
RESTBASEURL   = "http://127.0.0.1:5020"
ISDEBUG = True
LISTEN_TO = '127.0.0.1'		# only local addresses

# initialise Flask 
app = Flask(__name__)
CORS(app)	# allow CORS
app.logger.setLevel(logging.DEBUG)

			
def isjson(content):
		""" returns true if content parses as json, otherwise false. used by unit testing. """
		try:
			x = json.loads(content.decode('utf-8'))
			return True
 
		except json.decoder.JSONDecodeError:
			return False

def tojson(content):
	""" json dumps, formatting dates as isoformat """
	def converter(o):
		if isinstance(o, datetime.datetime):
			return o.isoformat()
		else:
			return json.JSONEncoder.default(o)
	return(json.dumps(content, default=converter))

# --------------------------------------------------------------------------------------------------
@app.errorhandler(404)
def not_found(error):
	json_err = jsonify({'error': 'Not found (custom error handler for mis-routing)'})
	return make_response(json_err, 404)
# --------------------------------------------------------------------------------------------------
 
@app.teardown_appcontext
def shutdown_session(exception=None):
	fn.PERSIST.closedown()		# close database connection

def do_GET(relpath):
	""" makes a GET request  to relpath.
		Used for unit testing.   """
	
	url = urljoiner(RESTBASEURL, relpath)
	print("GETing from: {0}".format(url))

	session = requests.Session()
	session.trust_env = False

	# print out diagnostics
	print("About to GET from url {0}".format(url))
	response = session.get(url=url, timeout=None)

	print("Result:")
	print("code: {0}".format(response.status_code))
	print("reason: {0}".format(response.reason))
	try:     
		print("text: {0}".format(response.text[:100]))
	except UnicodeEncodeError:
		# which is what happens if you try to display a gz file as text, which it isn't
		print("Response cannot be coerced to unicode ? a gz file.  The response content had {0} bytes.".format(len(response.text)))
		print("headers: {0}".format(response.headers))

	session.close()
	return(response)

def do_POST(relpath, payload):
	""" makes a POST request  to relpath.
		Used for unit testing.
		payload should be a dictionary"""
	
	url = urljoiner(RESTBASEURL, relpath)

	# print out diagnostics
	print("POSTING to url {0}".format(url))
	if not isinstance(payload, dict):
		raise TypeError("not a dict {0}".format(payload))
	response = requests.post(url=url, data=payload)

	print("Result:")
	print("code: {0}".format(response.status_code))
	print("reason: {0}".format(response.reason))
	print("content: {0}".format(response.content))
		
	return(response)

def render_markdown(md_file):
	""" render markdown as html
	"""
	with codecs.open(md_file, mode="r", encoding="utf-8") as f:
		text = f.read()
		html = markdown.markdown(text, extensions = ['tables'])
	return html

@app.route('/', methods=['GET'])
def routes():
	""" returns server info page
	"""
	routes_file = os.path.join("..","doc","rest-routes.md")
	return make_response(render_markdown(routes_file))

@app.route('/ui/info', methods=['GET'])
def server_info():
	""" returns server info page
	"""
	routes_file = os.path.join("..","doc","serverinfo.md")
	return make_response(render_markdown(routes_file))

@app.route('/api/v2/raise_error/<string:component>/<string:token>', methods=['GET'])
def raise_error(component, token):
	""" * raises an error internally.  Can be used to test error logging.  Disabled unless in debug mode.

	/api/v2/raise_error/*component*/*token*/

	Valid values for component are:
	main - raise error in main code
	persist - raise in PERSIST object
	clustering - raise in clustering
	seqcomparer - raise in seqcomparer.
	"""

	if not fn.debugMode == 2:
		# if we're not in debugMode==2, then this option is not allowed
		abort(404, 'Calls to /raise_error are only allowed with debugMode == 2' )

	if component == 'main':
		raise ZeroDivisionError(token)
	elif component == 'clustering':
		clustering_names = list(fn.clustering_settings.keys())

		if len(clustering_names)==0:
			self.fail("no clustering settings defined; cannot test error generation in clustering")
		else:
			clustering_name = clustering_names[0]
			fn.clustering_settings[clustering_name].raise_error(token)
	elif component == 'seqcomparer':
		fn.sc.raise_error(token)
	elif component == 'persist':
		fn.PERSIST.raise_error(token)
	else:
		raise KeyError("Invalid component called.  Allowed: main;persist;clustering;seqcomparer.")
	
def construct_msa(guids, output_format, what):
	
	""" constructs multiple sequence alignment for guids
		and returns in one of 'fasta' 'json-fasta', 'html', 'json' or 'json-records' format.
		
		what is one of 'N','M','N_or_M'
	
	"""
	res = fn.hc.multi_sequence_alignment(guids, output='df_dict', uncertain_base_type=what)
	df = pd.DataFrame.from_dict(res,orient='index')
	html = df.to_html()
	fasta= ""
	for guid in df.index:
		fasta=fasta + ">{0}\n{1}\n".format(guid, df.loc[guid,'aligned_seq'])
		
	if output_format == 'fasta':
		return make_response(fasta)
	elif output_format == 'json-fasta':
		return make_response(json.dumps({'fasta':fasta}))
	elif output_format == 'html':
		return make_response(html)
	elif output_format == 'json':
		return make_response(json.dumps(res))
	elif output_format == 'json-records':
		if len(df.index)>0:
			df['guid'] = df.index
		return make_response(df.to_json(orient='records'))
	
@app.route('/api/v2/reset', methods=['POST'])
def reset():
	""" deletes any existing data from the server """
	if not fn.debugMode == 2:
		# if we're not in debugMode==2, then this option is not allowed
		abort(404, 'Calls to /reset are only allowed with debugMode == 2' )
	else:
		fn.reset()
		return make_response(json.dumps({'message':'reset completed'}))

@app.route('/api/v2/monitor', methods=['GET'])
@app.route('/api/v2/monitor/<string:report_type>', methods=['GET'])
def monitor(report_type = 'Report' ):
	""" returns an html/bokeh file, generated by findNeighbour4-monitor,
	and stored in a database. If not report_type is specified, uses 'Report',
	which is the name of the default report produced by findNeighbour4-monitor"""
	
	html = fn.PERSIST.monitor_read(report_type)
	if html is None:
		html = "No report called {0} is available.  Check that findNeighbour4-monitor.py is running.".format(report_type)
	return html
@app.route('/api/v2/clustering/<string:clustering_algorithm>/<int:cluster_id>/network',methods=['GET'])
@app.route('/api/v2/clustering/<string:clustering_algorithm>/<int:cluster_id>/minimum_spanning_tree',methods=['GET'])
def cl2network(clustering_algorithm, cluster_id):
	""" produces a cytoscape.js compatible graph from a cluster ,
	either from the network (comprising all edges < snp cutoff)
	or as a minimal spanning tree.
	"""
	# validate input
	fn.clustering[clustering_algorithm].refresh()
	try:
		res = fn.clustering[clustering_algorithm].clusters2guidmeta(after_change_id = None)		
	except KeyError:
		# no clustering algorithm of this type
		return make_response(tojson("no clustering algorithm {0}".format(clustering_algorithm)), 404)
		
	# check guids
	df = pd.DataFrame.from_records(res)

	# check guids
	df = pd.DataFrame.from_records(res)
	
	if len(df.index)==0:
		return make_response(
								tojson(
									{'success':0, 'message':'No samples exist for that cluster'}
								)
							)
	else:
		df = df[df["cluster_id"]==cluster_id]		# only if there are records
		missing_guids = []
		guids = sorted(df['guid'].tolist())
					
		# data validation complete.  construct outputs
		snv_threshold = fn.clustering[clustering_algorithm].snv_threshold
		snvn = snvNetwork(snv_threshold = snv_threshold)
		E=[]
		for guid in guids:
			is_mixed = int(fn.clustering[clustering_algorithm].is_mixed(guid))
			snvn.G.add_node(guid, is_mixed=is_mixed)     
		for guid in guids:
			res = fn.PERSIST.guid2neighbours(guid, cutoff=snv_threshold, returned_format=1)
			for (guid2, snv) in res['neighbours']:
				if guid2 in guids:		# don't link outside the cluster
					E.append((guid,guid2))
					snvn.G.add_edge(guid, guid2, weight=snv, snv=snv)
				
		if request.base_url.endswith('/minimum_spanning_tree'):
			snvn.G = nx.minimum_spanning_tree(snvn.G)
			retVal =snvn.network2cytoscapejs()
			retVal['message']='{0} cluster #{1}. Minimum spanning tree is shown.  Red nodes are mixed.'.format(clustering_algorithm,cluster_id)
		else:
			retVal = snvn.network2cytoscapejs()
			retVal['message']='{0} cluster #{1}. Network of all edges < cutoff shown.  Red nodes are mixed.'.format(clustering_algorithm,cluster_id)
		retVal['success']=1
		return make_response(tojson(retVal))

	
@app.route('/api/v2/multiple_alignment/guids', methods=['POST'])
def msa_guids():
	""" performs a multiple sequence alignment on a series of POSTed guids,
	delivered in a dictionary, e.g.
	{'guids':'guid1;guid2;guid3',
	'output_format':'json'}
	
	Valid values for output_format are:
	json
	json-records
	html
	json-fasta
	fasta
	
	Valid values for what are
	N
	M
	N_or_M
	"""

	# validate input
	request_payload = request.form.to_dict()
	if 'output_format' in request_payload.keys() and 'guids' in request_payload.keys():
		guids = request_payload['guids'].split(';')		# coerce both guid and seq to strings
		output_format= request_payload['output_format']
		if 'what' in request_payload.keys():
			what = request_payload['what']
		else:
			what = 'N'		# default to N
		if not what in ['N','M','N_or_M']:
			abort(404, 'what must be one of N M N_or_M, not {0}'.format(what))
		if not output_format in ['html','json','fasta', 'json-fasta', 'json-records']:
			abort(404, 'output_format must be one of html, json, json-records or fasta not {0}'.format(output_format))
	else:
		abort(501, 'output_format and guids are not present in the POSTed data {0}'.format(data_keys))
	
	# check guids
	missing_guids = []
	for guid in sorted(guids):
		try:
			result = fn.exist_sample(guid)
		except Exception as e:
			capture_exception(e)
			abort(500, e)
		if not result is True:
			missing_guids.append(guid)
	
	if len(missing_guids)>0:
		capture_message("asked to perform multiple sequence alignment with the following missing guids: {0}".format(missing_guids))		
		abort(501, "asked to perform multiple sequence alignment with the following missing guids: {0}".format(missing_guids))
	
	# data validation complete.  construct outputs
	return construct_msa(guids, output_format, what)


@app.route('/api/v2/multiple_alignment_cluster/<string:clustering_algorithm>/<int:cluster_id>/<string:output_format>',methods=['GET'])
def msa_guids_by_cluster(clustering_algorithm, cluster_id, output_format):
	""" performs a multiple sequence alignment on the contents of a cluster
	
	Valid values for format are:
	json
	json-records
	fasta
	html
	"""
	
	# validate input
	fn.clustering[clustering_algorithm].refresh()
	try:
		res = fn.clustering[clustering_algorithm].clusters2guidmeta(after_change_id = None)		
	except KeyError:
		# no clustering algorithm of this type
		return make_response(tojson("no clustering algorithm {0}".format(clustering_algorithm)), 404)
		
	if not output_format in ['html','json','json-records','fasta','json-fasta']:
		abort(501, 'output_format must be one of html, json, json-records fasta or json-fasta not {0}'.format(output_format))

	# check guids
	df = pd.DataFrame.from_records(res)

	if len(df.index)==0:
		return make_response(
								json.dumps(
									{'status':'No samples exist for that cluster'}
								)
							)
	else:
		df = df[df["cluster_id"]==cluster_id]
		missing_guids = []
		guids = []
		for guid in sorted(df['guid'].tolist()):
			try:
				result = fn.exist_sample(guid)
			except Exception as e:
				capture_exception(e)
				abort(500, e)
			if not result is True:
				missing_guids.append(guid)
			else:
				guids.append(guid)
		
		if len(missing_guids)>0:
			abort(501, "asked to perform multiple sequence alignment with the following missing guids: {0}".format(missing_guids))
			
		# data validation complete.  construct outputs
		return construct_msa(guids, output_format, what=fn.clustering[clustering_algorithm].uncertain_base_type)



@app.route('/api/v2/server_config', methods=['GET'])
def server_config():
	""" returns server configuration.

		returns the config file with which the server was launched.
		This may be highly undesirable,
		as it reveals the internal server architecture  including
		backend databases and perhaps connection strings with passwords.

	"""
	res = fn.server_config()
	if res is None:		# not allowed to see it
		return make_response(tojson({'NotAvailable':"Endpoint is only available in debug mode"}), 404)
	else:
		return make_response(tojson(CONFIG))



@app.route('/api/v2/server_memory_usage', defaults={'nrows':100, 'output_format':'json'}, methods=['GET'])
@app.route('/api/v2/server_memory_usage/<int:nrows>', defaults={'output_format':'json'}, methods=['GET'])
@app.route('/api/v2/server_memory_usage/<int:nrows>/<string:output_format>', methods=['GET'])
def server_memory_usage(nrows, output_format):
	""" returns server memory usage information, as list.
	The server notes memory usage at various key points (pre/post insert; pre/post recompression)
	and these are stored. """
	try:
		result = fn.server_memory_usage(max_reported = nrows)

	except Exception as e:
		
		capture_exception(e)
		abort(500, e)
	
	# reformat this into a long, human readable format.
	
	resl = pd.melt(pd.DataFrame.from_records(result), id_vars = ['_id','context|time|time_now', 'context|info|message']).dropna()       # drop any na values
	resl.columns = ['_id','event_time','info_message','event_description','value']
	resl = resl[resl.event_description.astype(str).str.startswith('server')]		# only server
	resl['descriptor1']='Server'
	resl['descriptor2']='RAM'
	resl['detail']= [fn.mhr.convert(x) for x in resl['event_description'].tolist()]
	resl = resl.drop(['event_description'], axis =1)

	if output_format == 'html':
		return(resl.to_html())
	elif output_format == 'json':
		return make_response(resl.to_json(orient='records'))
	else:
		abort(500, "Invalid output_format passed")
		
@app.route('/ui/server_status', defaults={'absdelta':'absolute', 'stats_type':'mstat', 'nrows':1}, methods=['GET'])
@app.route('/ui/server_status/<string:absdelta>/<string:stats_type>/<int:nrows>', methods=['GET'])
def server_storage_status(absdelta, stats_type, nrows):
	""" returns server memory usage information, as list.
	The server notes memory usage at various key points (pre/post insert; pre/post recompression)
	and these are stored."""
	try:
		result = fn.server_memory_usage(max_reported = nrows)
		df = pd.DataFrame.from_records(result, index='_id')  #, coerce_float=True

		# identify target columns
		valid_starts = ['clusters',
						'guid2meta',
						'guid2neighbour',
						'refcompressedseq',
						'server',
						'mstat',
						'scstat']
		
		if not stats_type in valid_starts:
			abort(404, "Valid stats_type values are {0}".format(valid_starts))
		
		if len(df.columns.values)==0:
			return("No column data found from database query")
		
		target_columns = []
		
		target_string = "{0}".format(stats_type)
		for col in df.columns.values:
			if col.find(target_string)>=0:
				if absdelta=='delta' and col.endswith('|delta'):
					target_columns.append(col)
				elif absdelta=='absolute' and not col.endswith('|delta'):
					target_columns.append(col)
		if nrows<1:
			return("More than one row must be requested.")
		if len(target_columns)==0:
			return("No column data found matching this selection. <p>This may be normal if the server has just started up.<p>We tried to select from {2} rows of data, with {3} columns.  We looked for '{4}'.<p>Valid values for the three variables passed in the URL are as follows: <p> stats_type: {0}. <p> absdelta: ['absolute', 'delta']. <p> nrows must be a positive integer. <p> The columns available for selection from the server's monitoring log are: {1}".format(valid_starts,df.columns.values, len(df.index), len(df.columns.values), target_string))
		if len(df.index)==0:
			return("No row data found matching this selection. <p>This may be normal if the server has just started up.<p> We tried to select from {2} rows of data, with {3} columns.  We looked for '{4}'.<p>Valid values for the three variables passed in the URL are as follows: <p> stats_type: {0}. <p> absdelta: ['absolute', 'delta']. <p> nrows must be a positive integer. <p> The columns available for selection from the server's monitoring log are: {1}".format(valid_starts,df.columns.values, len(df.index), len(df.columns.values), target_string))
		
		# convert x-axis to datetime
		for ix in df.index:
			try:
				df.loc[ix,'context|time|time_now']= dateutil.parser.parse(df.loc[ix,'context|time|time_now'])
			except TypeError:
				app.logger.warning("Attempted date conversion on {0} with _id = {1}, but this failed".format(df.loc[ix,'time|time_now'], ix))
				df.loc[ix,'context|time|time_now']=None
	
		# if values are not completed, then use the previous non-null version
		# see https://stackoverflow.com/questions/14399689/matplotlib-drawing-lines-between-points-ignoring-missing-data
		select_cols = target_columns.copy()
		select_cols.append('context|time|time_now')
		dfp = df[select_cols]
		dfp = dfp.dropna()
		if len(dfp.index)==0:
			return("No non-null row data found matching this selection. <p>This may be normal if the server has just started up.<p> We tried to select from {2} rows of data, with {3} columns.  We looked for '{4}'.<p>Valid values for the three variables passed in the URL are as follows: <p> stats_type: {0}. <p> absdelta: ['absolute', 'delta']. <p> nrows must be a positive integer. <p> The columns available for selection from the server's monitoring log are: {1}".format(valid_starts,df.columns.values, len(df.index), len(df.columns.values), target_string))
		
		# construct a dictionary mapping current column names to human readable versions
		mapper={}
		new_target_columns = []
		for item in target_columns:
			mapper[item] = fn.mhr.convert(item)
			new_target_columns.append(mapper[item])
		dfp.rename(mapper, inplace=True, axis='columns')

		# plot
		plts = dfp.plot(kind='line', x='context|time|time_now', subplots=True, y=new_target_columns)
		for plt in plts:
			fig = plt.get_figure()
			fig.set_figheight(len(target_columns)*2)
			fig.set_figwidth(8)
			img = io.BytesIO()
			fig.savefig(img)
			matplotlib.pyplot.close('all')		# have to explicitly close, otherwise memory leaks 
			img.seek(0)
			return send_file(img, mimetype='image/png')

	except Exception as e:
		capture_exception(e)
		abort(500, e)
		
	return make_response(tojson(result))


@app.route('/api/v2/snpceiling', methods=['GET'])
def snpceiling():
	""" returns largest snp distance stored by the server """
	try:
		result = {"snpceiling":fn.snpCeiling}
		
	except Exception as e:
		capture_exception(e)
		abort(500, e)
	return make_response(tojson(result))


@app.route('/api/v2/server_time', methods=['GET'])
def server_time():
	""" returns server time """
	try:
		result = fn.server_time()

	except Exception as e:
		capture_exception(e)
		abort(500, e)
	return make_response(tojson(result))

@app.route('/api/v2/server_name', methods=['GET'])
def server_name():
	""" returns server name """
	try:
		result = fn.server_name()

	except Exception as e:
		capture_exception(e)
		abort(500, e)
	return make_response(tojson(result))

	
@app.route('/api/v2/guids', methods=['GET'])
def get_all_guids(**debug):
	""" returns all guids.  other params, if included, is ignored."""
	try:
		result = list(fn.get_all_guids())
	except Exception as e:
		capture_exception(e)
		abort(500, e)
	return(make_response(tojson(result)))


@app.route('/api/v2/guids_with_quality_over/<float:cutoff>', methods=['GET'])
@app.route('/api/v2/guids_with_quality_over/<int:cutoff>', methods=['GET'])
def guids_with_quality_over(cutoff, **kwargs):
	""" returns all guids with quality score >= cutoff."""
	try:
		result = fn.guids_with_quality_over(cutoff)	
	except Exception as e:
		capture_exception(e)
		abort(500, e)
	return make_response(tojson(result))


@app.route('/api/v2/guids_and_examination_times', methods=['GET'])
def guids_and_examination_times(**kwargs):
	""" returns all guids and their examination (addition) time.
	reference, if passed, is ignored."""
	try:	
		result =fn.get_all_guids_examination_time()	
	except Exception as e:
		capture_exception(e)
		abort(500, e)
	return make_response(tojson(result))



@app.route('/api/v2/guids_beginning_with/<string:startstr>', methods=['GET'])
def get_matching_guids(startstr, max_returned=30):
	""" returns all guids matching startstr.
	A maximum of max_returned matches is returned.
	If > max_returned records match, then an empty list is returned.
	"""
	try:
		result = fn.gs.search(search_string= startstr, max_returned=max_returned)
		app.logger.debug(result)
	except Exception as e:
		capture_exception(e)
		abort(500, e)
	return(make_response(tojson(result)))





@app.route('/api/v2/annotations', methods=['GET'])
def annotations(**kwargs):
	""" returns all guids and associated meta data.
	This query can be slow for very large data sets.
	"""
	try:
		result = fn.get_all_annotations()
		
	except Exception as e:
		capture_exception(e)
		abort(500, e)
		
	return(tojson(result))


@app.route('/api/v2/<string:guid>/exists', methods=['GET'])
def exist_sample(guid, **kwargs):
	""" checks whether a guid exists.
	reference and method are ignored."""
	
	try:
		result = fn.exist_sample(guid)
		
	except Exception as e:
		capture_exception(e)
		abort(500, e)
		
	return make_response(tojson(result))


@app.route('/api/v2/<string:guid>/clusters', methods=['GET'])
def clusters_sample(guid):
	""" returns clusters in which a sample resides """

	clustering_algorithms = sorted(fn.clustering.keys())
	retVal=[]
	for clustering_algorithm in clustering_algorithms:
		fn.clustering[clustering_algorithm].refresh()
		res = fn.clustering[clustering_algorithm].clusters2guidmeta(after_change_id = None)
		for item in res:	
			if item['guid']==guid:
				item['clustering_algorithm']=clustering_algorithm
				retVal.append(item)
	if len(retVal)==0:
		abort(404, "No clustering information for guid {0}".format(guid))
	else:
		return make_response(tojson(retVal))




@app.route('/api/v2/<string:guid>/annotation', methods=['GET'])
def annotations_sample(guid):
	""" returns annotations of one sample """
	
	try:
		result = fn.PERSIST.guid_annotation(guid)
	
	except Exception as e:
		capture_exception(e)
		abort(500, e)
	
	if len(result.keys())==0:
		abort(404, "guid does not exist {0}".format(guid))
	retVal = result[guid]
	retVal['guid'] = guid
	return make_response(tojson(retVal))


@app.route('/api/v2/insert', methods=['POST'])
def insert():
	""" inserts a guids with sequence"""
	try:
		data_keys = set()
		for key in request.form.keys():
			data_keys.add(key)
		payload = {}
		for key in data_keys:
			payload[key]= request.form[key]
	
		if 'seq' in data_keys and 'guid' in data_keys:
			guid = str(payload['guid'])
			seq  = str(payload['seq'])
			result = fn.insert(guid, seq)
		else:
			abort(501, 'seq and guid are not present in the POSTed data {0}'.format(data_keys))
		
	except Exception as e:
		capture_exception(e)
		abort(500, e)
		
	return make_response(tojson(result))

@app.route('/api/v2/mirror', methods=['POST'])
def mirror():
	""" receives data, returns the dictionary it was passed. Takes no other action.
	Used for testing that gateways etc don't remove data."""

	retVal = {}
	for key in request.form.keys():
		retVal[key]=request.form[key]
	return make_response(tojson(retVal))

@app.route('/api/v2/clustering', methods=['GET'])
def algorithms():
	"""  returns the available clustering algorithms """
	res = sorted(fn.clustering.keys())		
	return make_response(tojson({'algorithms':res}))


@app.route('/api/v2/clustering/<string:clustering_algorithm>/what_tested', methods=['GET'])
def what_tested(clustering_algorithm):
	"""  returns what is tested (N, M, N_or_M) for clustering_algorithm.
		 Useful for producing reports of what clustering algorithms are doing
	"""
	try:
		fn.clustering[clustering_algorithm].refresh()
		
	except KeyError:
		# no clustering algorithm of this type
		abort(404, "no clustering algorithm {0} (failed update)".format(clustering_algorithm))
	
	try:
		
		res = fn.clustering[clustering_algorithm].uncertain_base_type
	except KeyError:
		# no clustering algorithm of this type
		abort(404, "no clustering algorithm {0} (failed base recovery)".format(clustering_algorithm))
		
	return make_response(tojson({'what_tested': res, 'clustering_algorithm':clustering_algorithm}))




@app.route('/api/v2/clustering/<string:clustering_algorithm>/change_id', methods=['GET'])
def change_id(clustering_algorithm):
	"""  returns the current change_id number, which is incremented each time a change is made.
		 Useful for recovering changes in clustering after a particular point."""
	try:
		fn.clustering[clustering_algorithm].refresh()
		res = fn.clustering[clustering_algorithm].change_id		
	except KeyError:
		# no clustering algorithm of this type
		abort(404, "no clustering algorithm {0}".format(clustering_algorithm))
		
	return make_response(tojson({'change_id': res, 'clustering_algorithm':clustering_algorithm}))

@app.route('/api/v2/clustering/<string:clustering_algorithm>/guids2clusters', methods=['GET'])
@app.route('/api/v2/clustering/<string:clustering_algorithm>/guids2clusters/after_change_id/<int:change_id>', methods=['GET'])
def g2c(clustering_algorithm, change_id=None):
	"""  returns a guid -> clusterid dictionary for all guids """
	try:

		fn.clustering[clustering_algorithm].refresh()
		
	except KeyError:
		app.logging.info("No algorithm {0}".format(clustering_algorithm))
		abort(404, "no clustering algorithm {0}".format(clustering_algorithm))
	
	res = fn.clustering[clustering_algorithm].clusters2guidmeta(after_change_id = change_id)		
	
	return make_response(tojson(res))


@app.route('/api/v2/clustering/<string:clustering_algorithm>/clusters', methods=['GET'])
@app.route('/api/v2/clustering/<string:clustering_algorithm>/members', methods=['GET'])
@app.route('/api/v2/clustering/<string:clustering_algorithm>/summary', methods=['GET'])
@app.route('/api/v2/clustering/<string:clustering_algorithm>/<int:cluster_id>', methods=['GET'])
def clusters2cnt(clustering_algorithm, cluster_id = None):
	"""  returns a dictionary containing
		 'summary':a clusterid -> count dictionary for all cluster_ids for clustering_algorithm,
		 'members':a list of all guids and clusterids for clustering algorithm
		 
		 * If cluster_id is specified, only returns details for one cluster_id.
		 * If /clusters is requested, returns a dictionary with both 'summary' and 'members' keys.
		 * If /members is requested, returns a dictionary with only 'members' key
		 * If /summary is requested, returns a dictionary with only 'summary' key"""

	try:
		fn.clustering[clustering_algorithm].refresh()
		all_res = fn.clustering[clustering_algorithm].clusters2guidmeta(after_change_id = None)

	except KeyError:
		# no clustering algorithm of this type
		abort(404, "no clustering algorithm {0}".format(clustering_algorithm))

	# if no cluster_id is specified, then we return all data.
	if cluster_id is None:
		res = all_res
	else:
		
		res = []
		for item in all_res:
			
			if item['cluster_id'] == cluster_id:
				res.append(item)
		if len(res) == 0:
			# no cluster exists of that name
			abort(404, "no cluster {1} exists for algorithm {0}".format(clustering_algorithm, cluster_id))
			
	d= pd.DataFrame.from_records(res)
	
	try:
		df = pd.crosstab(d['cluster_id'],d['is_mixed'])
		df = df.add_prefix('is_mixed_')
		df['cluster_id']=df.index
		summary = json.loads(df.to_json(orient='records'))
		detail  = json.loads(d.to_json(orient='records'))
		#print(request.url, request.url.endswith('summary'), request.url.endswith('members'))
		if cluster_id is not None:
			retVal = {"summary":summary, "members":detail}
		elif request.url.endswith('clusters'):
			retVal = {"summary":summary, "members":detail}			
		elif request.url.endswith('summary'):
			retVal = {"summary":summary}
		elif request.url.endswith('members'):
			retVal = {"members":detail}
		else:
			abort(404, "url not recognised: "+request.url)
	except KeyError:  # no data
		if cluster_id is not None:
			retVal = {"summary":[], "members":[]}
		elif request.url.endswith('clusters'):
			retVal = {"summary":[], "members":[]}			
		elif request.url.endswith('summary'):
			retVal = {"summary":[]}
		elif request.url.endswith('members'):
			retVal = {"members":[]}
		else:
			abort(404, "url not recognised: "+request.url)

	return make_response(tojson(retVal))

@app.route('/api/v2/clustering/<string:clustering_algorithm>/cluster_ids', methods=['GET'])
def g2cl(clustering_algorithm):
	"""  returns a guid -> clusterid dictionary for all guids """
	try:
		fn.clustering[clustering_algorithm].refresh()
		res = fn.clustering[clustering_algorithm].clusters2guidmeta(after_change_id = None)		
	except KeyError:
		# no clustering algorithm of this type
		abort(404, "no clustering algorithm {0}".format(clustering_algorithm))
	cluster_ids = set()
	for item in res:
		try:
			cluster_ids.add(item['cluster_id'])
		except KeyError:
			# there's no cluster_id
			pass
	retVal = sorted(list(cluster_ids))
	return make_response(tojson(retVal))


		
		



@app.route('/api/v2/<string:guid>/neighbours_within/<int:threshold>', methods=['GET'])
@app.route('/api/v2/<string:guid>/neighbours_within/<int:threshold>/with_quality_cutoff/<float:cutoff>', methods=['GET'])
@app.route('/api/v2/<string:guid>/neighbours_within/<int:threshold>/with_quality_cutoff/<int:cutoff>', methods=['GET'])
@app.route('/api/v2/<string:guid>/neighbours_within/<int:threshold>/with_quality_cutoff/<float:cutoff>/in_format/<int:returned_format>', methods=['GET'])
@app.route('/api/v2/<string:guid>/neighbours_within/<int:threshold>/with_quality_cutoff/<int:cutoff>/in_format/<int:returned_format>', methods=['GET'])
@app.route('/api/v2/<string:guid>/neighbours_within/<int:threshold>/in_format/<int:returned_format>', methods=['GET'])
def neighbours_within(guid, threshold, **kwargs):
	""" get a guid's neighbours, within a threshold """
	# we support optional cutoff and threshold parameters.
	# we also support 'method' and 'reference' parameters but these are ignored.
	# the default for cutoff and format are 0.85 and 1, respectively.
	if not 'cutoff' in kwargs.keys():
		cutoff = CONFIG['MAXN_PROP_DEFAULT']
	else:
		cutoff = kwargs['cutoff']
		
		
	if not 'returned_format' in kwargs.keys():
		returned_format = 1
	else:
		returned_format = kwargs['returned_format']
		
	# validate input
	if not returned_format in set([1,2,3,4]):
		abort(500, "Invalid format requested, must be 1, 2, 3 or 4.")
	if not ( 0 <= cutoff  and cutoff <= 1):
		abort(500, "Invalid cutoff requested, must be between 0 and 1")
		
	try:
		result = fn.neighbours_within_filter(guid, threshold, cutoff, returned_format)
	except KeyError as e:
		# guid doesn't exist
		abort(404, e)
	except Exception as e:
		capture_exception(e)
		abort(500, e)
	
	return make_response(tojson(result))
	

@app.route('/api/v2/<string:guid>/sequence', methods=['GET'])
def sequence(guid):
	""" returns the masked sequence as a string """	
	result = fn.sequence(guid)
	if result is None:  # no guid exists
		return make_response(tojson('guid {0} does not exist'.format(guid)), 404)
	else:
		return make_response(tojson(result))


@app.route('/api/v2/nucleotides_excluded', methods=['GET'])
def nucleotides_excluded():
	""" returns all nucleotides excluded by the server.
	Useful for clients which need to to ensure that server
	and client masking are identical. """
	
	try:
		result = fn.server_nucleotides_excluded()
		
	except Exception as e:
		capture_exception(e)
		abort(500, e)

	return make_response(tojson(result))


# startup
if __name__ == '__main__':

	# command line usage.  Pass the location of a config file as a single argument.
	parser = argparse.ArgumentParser(
		formatter_class= argparse.RawTextHelpFormatter,
		description="""Runs findNeighbour4-server, a service for bacterial relatedness monitoring.
									 

Example usage: 
============== 
# show command line options 
python findNeighbour4-server.py --help 	

# run with debug settings; only do this for unit testing.
python findNeighbour4-server.py 	

# run using settings in myConfigFile.json.  Memory will be recompressed after loading. 
python findNeighbour4-server.py ../config/myConfigFile.json		

# run using settings in myConfigFile.json; 
# recompress RAM every 20000 samples loaded 
# (only needed if RAM is in short supply and data close to limit) 
# enabling this option will slow up loading 

python findNeighbour4-server.py ../config/myConfigFile.json	\ 
                        --on_startup_recompress-memory_every 20000 

""")
	parser.add_argument('path_to_config_file', type=str, action='store', nargs='?',
						help='the path to the configuration file', default='')
	parser.add_argument('--on_startup_recompress_memory_every', type=int, nargs=1, action='store', default=[None], 
						help='when loading, recompress server memory every so many samples.')
	args = parser.parse_args()
	
	# an example config file is default_test_config.json

	############################ LOAD CONFIG ######################################
	print("findNeighbour4 server .. reading configuration file.")

	if len(args.path_to_config_file)>0:
			configFile = args.path_to_config_file
	else:
			configFile = os.path.join('..','config','default_test_config.json')
			warnings.warn("No config file name supplied ; using a configuration ('default_test_config.json') suitable only for testing, not for production. ")

	# open the config file
	try:
			with open(configFile,'r') as f:
					 CONFIG=f.read()

	except FileNotFoundError:
			raise FileNotFoundError("Passed a positional parameter, which should be a CONFIG file name; tried to open a config file at {0} but it does not exist ".format(sys.argv[1]))

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
	app.logger.setLevel(loglevel)       
	file_handler = logging.FileHandler(CONFIG['LOGFILE'])
	formatter = logging.Formatter( "%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s ")
	file_handler.setFormatter(formatter)
	app.logger.addHandler(file_handler)

	# log a test error on startup
	# app.logger.error("Test error logged on startup, to check logger is working")

	# launch sentry if API key provided
	if 'SENTRY_URL' in CONFIG.keys():
			app.logger.info("Launching communication with Sentry bug-tracking service")
			sentry_sdk.init(CONFIG['SENTRY_URL'], integrations=[FlaskIntegration()])

	########################### prepare to launch server ###############################################################
	# construct the required global variables
	LISTEN_TO = '127.0.0.1'
	if 'LISTEN_TO' in CONFIG.keys():
		LISTEN_TO = CONFIG['LISTEN_TO']

	RESTBASEURL = "http://{0}:{1}".format(CONFIG['IP'], CONFIG['REST_PORT'])

	#########################  CONFIGURE HELPER APPLICATIONS ######################
	## once the flask app is running, errors get logged to app.logger.  However, problems on start up do not.
	## configure mongodb persistence store

	# plotting engine
	matplotlib.use('agg')		#  prevent https://stackoverflow.com/questions/27147300/how-to-clean-images-in-python-django

	if 'SENTRY_URL' in CONFIG.keys():
			app.logger.info("Launching communication with Sentry bug-tracking service")
			sentry_sdk.init(CONFIG['SENTRY_URL'], integrations=[FlaskIntegration()])

	print("Connecting to backend data store")
	try:
			PERSIST=fn3persistence(dbname = CONFIG['SERVERNAME'],			
						connString=CONFIG['FNPERSISTENCE_CONNSTRING'],
						debug=CONFIG['DEBUGMODE'])
	except Exception as e:
			app.logger.exception("Error raised on creating persistence object")
			if e.__module__ == "pymongo.errors":
				  app.logger.info("Error raised pertains to pyMongo connectivity")
			raise

	# instantiate server class
	print("Loading sequences into server, please wait ...")
	try:
		fn = findNeighbour4(CONFIG, PERSIST)
	except Exception as e:
			app.logger.exception("Error raised on instantiating findNeighbour4 object")
			raise


	########################  START THE SERVER ###################################
	if CONFIG['DEBUGMODE']>0:
			flask_debug = True
			app.config['PROPAGATE_EXCEPTIONS'] = True
	else:
			flask_debug = False

	app.logger.info("Launching server listening to IP {0}, debug = {1}, port = {2}".format(LISTEN_TO, flask_debug, CONFIG['REST_PORT']))
	app.run(host=LISTEN_TO, debug=flask_debug, port = CONFIG['REST_PORT'])




