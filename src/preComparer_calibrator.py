#!/usr/bin/env python
""" 
Calibrate preComparer and report on performance of preComparer on real data
"""
 
# import libraries [ not all are needed TODO]
import os
import sys
import requests
import json
import logging
import warnings
import datetime
import glob
import time
import sys
import hashlib
import queue
import threading
import random
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
import bokeh
import pickle
import sqlalchemy 
from scipy.stats import binom_test, median_absolute_deviation
from sentry_sdk import capture_message, capture_exception
from sentry_sdk.integrations.flask import FlaskIntegration

# logging
from logging.config import dictConfig

# utilities for file handling and measuring file size
import psutil

# reference based compression, storage and clustering modules
from NucleicAcid import NucleicAcid
from mongoStore import fn3persistence
from seqComparer import seqComparer		# import from seqComparer_mt for multithreading
from preComparer import preComparer
from clustering import snv_clustering
from guidLookup import guidSearcher  # fast lookup of first part of guids

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

class Calibrator():
	""" analyses existing sequences in a findNeighbour3/4 , derives and tests sequences for
		 """ 
		
	def __init__(self,CONFIG, PERSIST):
		""" Using values in CONFIG, estimates based on distances in the server with CONFIG['NAME'] on port CONFIG['PORT'].
		
		Note that the findNeighbour3 server itself does not have to be operational for this class to work, but the underlying database does have to be present.

		related to error handling:
		SENTRY_URL (optional)
		Note: if a FN_SENTRY URL environment variable is present, then the value of this will take precedence over any values in the config file.
		This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
		PERSIST is a storage object needs to be supplied.  The fn3Persistence class in mongoStore is one suitable object.
		PERSIST=fn3persistence(connString=CONFIG['FNPERSISTENCE_CONNSTRING'])

		"""
		
		# store the persistence object as part of the object
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
						   'LOGFILE','LOGLEVEL'])
		missing=required_keys-set(self.CONFIG.keys())
		if not missing == set([]):
			raise KeyError("Required keys were not found in CONFIG. Missing are {0}".format(missing))

		# load reference
		cfg = self.PERSIST.config_read('config')
		self.reference = cfg['reference']
		self.snpCeiling = cfg['SNPCEILING']
		self.maxNs = cfg['MAXN_STORAGE']

		# set easy to read properties from the config
		self.analysed_reference_length = len( cfg['reference']) - len(set(cfg['excludePositions']))

		# we start without any variation model
		self._reset()
		

	def _reset(self):
		""" clears existing variation model and pca result """ 
		pass

	def guids(self):
		""" returns list of guids currently in the findNeighbour3/4 database"""
		return sorted(self.PERSIST.refcompressedsequence_guids())	# sorting is not essential, but makes for deterministic  testing
		
		
	def derive_settings(self, 
							selection_cutoff = 20, 
							alpha = 1e-14,
							train_on = None, 
							export_raw = None):
		""" 
			input:
				selection_cutoff: SNP distances above this are not reported
				train_on: the number of sequences to analyse
				export_raw: a writable directory into which the data studied will be written.  If none, nothing will be written.
			returns:
				a dictionary containing recommended preComparer settings
			
		"""
		retVal = {'selection_cutoff':selection_cutoff,
					'n_positions_examined':self.analysed_reference_length,
					'alpha':alpha}

		pc = preComparer()
		print("Deriving settings: analysed sequence length",self.analysed_reference_length)

		# determine guids there in the database
		guids = self.guids()
		
	
		if train_on is None:
			train_on = len(guids)

		bar = progressbar.ProgressBar(max_value = train_on)
		nLoaded =0
		for guid in guids:
			nLoaded+=1

			if nLoaded>=train_on:		# trained on enough samples
				break
				pass
			bar.update(nLoaded)
			obj = self.PERSIST.refcompressedsequence_read(guid)	# ref compressed sequence

			if export_raw is not None:
				filename = "{0}.pickle".format(guid)
				filepath = os.path.join(export_raw,filename)
				with open(filepath,'wb') as f:
					pickle.dump(obj,f)

			pc.persist(obj, guid)

		bar.finish()

		composition = pd.DataFrame.from_dict(pc.composition, orient = 'index')
		# allows computation of Z
		valid = composition.query("invalid==0")
		retVal['N_mean'] = np.median(valid['N'])
		retVal['N_sd'] = median_absolute_deviation(valid['N'])

		# next compute over_selection_cutoff_ignore_factor etc.
		# for now, we just supply defaults
		retVal['over_selection_cutoff_ignore_factor'] = 5
		retVal['mixed_reporting_cutoff'] = 0
		retVal['highN_z_reporting_cutoff'] = 2
		retVal['probN_inflation_factor'] = 3

		return retVal

	def test_settings(self, 
							preComparer_settings,
							input_connstring=None,
							train_on = 500,
							output_stem = 'output',
							start_afresh=True):
		""" 
			input:
				selection_cutoff: SNP distances above this are not reported
				train_on: the number of sequences to analyse
				input_connstring: optional; a connection string to a database containing pairwise comparisons of seqComparer (accurate SNP distances) with pcaComparer (approximate distances) created by a previous run of test_settings.  If supplied, test_settings will conduct an analysis on the same samples, in the same order, using the exact distances from the previous analysis.  This speeds up recomputation.
				output_stem: the name and, optionally, path of an sqlite database to which the output will be written.
				start_afresh: discard any previous computations
			returns:
				a connection string to a database containing pairwise comparisons of seqComparer (accurate SNP distances) with pcaComparer (approximate distances)
			
		"""
	
		sql_filename = "{0}.sqlite".format(output_stem)
		if start_afresh:
			try:
				os.unlink(sql_filename)
			except FileNotFoundError:
				pass
		connstring = 'sqlite:///{0}'.format(sql_filename)
		engine = sqlalchemy.create_engine(connstring, echo=False)

		pc = preComparer(**preComparer_settings)
		sc = seqComparer(reference = self.reference, maxNs = self.maxNs , snpCeiling= 1e12)
		print("testing settings")

		# determine guids there in the database
		guids = self.guids()
	
		if train_on is None:
			train_on = len(guids)

		# load samples 
		print("Loading samples")
		bar = progressbar.ProgressBar(max_value = train_on)
		nLoaded =0
		test_guids = set()
		for guid in guids:
			nLoaded+=1

			if nLoaded>=train_on:		# trained on enough samples
				break
				pass
			bar.update(nLoaded)
			obj = self.PERSIST.refcompressedsequence_read(guid)	# ref compressed sequence

			pc.persist(obj, guid)
			sc.persist(obj, guid)
			if obj['invalid']==0:
				test_guids.add(guid)

		bar.finish()

		# compute exact distances
		print("Trying to load exact results")
		try:
			exact_results = pd.read_sql("select * from exact_results;", con=engine)
			print(".. succeeded: using stored exact results ")

		except sqlalchemy.exc.OperationalError:		# doesn't exist

			print(".. none are stored; so computing exact distances")
			bar = progressbar.ProgressBar(max_value = train_on*train_on)
			nLoaded =0
			result = {}
			for guid1 in test_guids:

				for guid2 in test_guids:
					if guid1 < guid2:
						nLoaded+=1
						sc_dist= {'guid1':guid1,'guid2':guid2}
						bar.update(nLoaded)
						t2= datetime.datetime.now()
						[guid1,guid2,dist,n1,n2,nboth,N1pos, N2pos, Nbothpos] = sc.countDifferences_byKey((guid1,guid2))
						sc_dist['dexact']=dist
						t3= datetime.datetime.now()
						i2 = (t3-t2).microseconds
						sc_dist['seqComparer_timing'] = i2
						
						result[nLoaded]=sc_dist
	
			bar.finish()
			per_cmp = pd.DataFrame.from_dict(result, orient='index')
			per_cmp.to_sql('exact_results', con=engine, if_exists = 'append')	
			exact_results = pd.read_sql("select * from exact_results;", con=engine)
			
		# store the settings used
		
		try:
			existing_settings_analysed = pd.read_sql("select * from settings;", con=engine)
			print("Previous settings have been analysed:")
			print(existing_settings_analysed)
			settings_ix= len(existing_settings_analysed.index)
		except sqlalchemy.exc.OperationalError:		# first run, no settings table
			print("No previous settings have been analysed")
			settings_ix = 0
		settings_ix +=1
		print("Storing new settings with index = ",settings_ix)
		new_setting = pd.DataFrame.from_dict({settings_ix: preComparer_settings}, orient='index')
		new_setting.to_sql('settings', con=engine, if_exists = 'append')	

		print("Computing approximate distances")
		print("Settings used are:")
		print(new_setting)

		bar = progressbar.ProgressBar(max_value = train_on*train_on)
		nLoaded =0
		try:
			res, = pd.read_sql("select max(index) from validation;", con=engine)
			print(res)
			per_cmp_ix= 1 + res
		except sqlalchemy.exc.OperationalError:		# first run, no settings table
			per_cmp_ix = 0
	

		result = {}

		for exact_results_ix in exact_results.index:

			guid1 = exact_results.loc[exact_results_ix, 'guid1']
			guid2 = exact_results.loc[exact_results_ix, 'guid2']

			nLoaded+=1
			per_cmp_ix+=1
			bar.update(nLoaded)
			t1= datetime.datetime.now()
			pc_dist = pc.compare(guid1,guid2)
			t2= datetime.datetime.now()
			pc_dist['settings_ix']=settings_ix
			pc_dist['dexact'] = exact_results.loc[exact_results_ix,'dexact']
			pc_dist['seqComparer_timing']=exact_results.loc[exact_results_ix,'seqComparer_timing']

			i1  = (t2-t1).microseconds
			pc_dist['preComparer_timing'] = i1
			pc_dist['ratio_timing'] = pc_dist['seqComparer_timing']/i1
			pc_dist['low_distance'] = (pc_dist['dexact'] <= self.snpCeiling)
			result[per_cmp_ix]=pc_dist
	
		bar.finish()
		per_cmp = pd.DataFrame.from_dict(result, orient='index')
		per_cmp.to_sql('validation', con=engine, if_exists = 'append')		


		### TODO: move out to different reporting function
		# read from database
		summary = pd.read_sql("select settings_ix, low_distance,reported_category, count(*) n from validation  group by settings_ix, low_distance, reported_category;", con=engine)
		print("Results:")
		print(summary)

		failures = pd.read_sql("select settings_ix, reported_category, count(*) n from validation where low_distance=1 and not substr(reported_category,1,9)='No retest' group by settings_ix, reported_category;", con=engine)
		print("The following regimes failed QC")
		print(failures)


		per_category = summary.pivot(index=None, columns='reported_category', values = 'n')

		return connstring

	def retest_settings(self, 
							connstring,
							preComparer_settings,
							output_stem = 'output'):
		""" using output from a previous test_settings() call, examine 
			input:
				connstring: a connection string to a database containing pairwise comparisons of seqComparer (accurate SNP distances) with pcaComparer (approximate distances), as provided by test_settings
				preComparer_settings: a dictionary of precomparer settings, as returned by derive_settings
				output_stem: the name and, optionally, path of an sqlite database to which the output will be written.
			returns:
				a connection string to a database containing pairwise comparisons of seqComparer (accurate SNP distances) with pcaComparer (approximate distances)
			
		"""

		# iterate over sqlite table
		# classification algorithm

		if (self.composition[key1]['Z'] > self.highN_z_reporting_cutoff  or self.composition[key2]['Z'] > self.highN_z_reporting_cutoff ):
		    if penalised_destim > self.over_selection_cutoff_ignore_factor * self.selection_cutoff:
		        res['reported_category'] = 'No retest - HighN'
		        res['no_retest'] = True
		    else:
		        res['reported_category'] = 'High N'
		elif	res['mixed_in_cmp']> self.mixed_reporting_cutoff:
		    if penalised_destim > self.over_selection_cutoff_ignore_factor * self.selection_cutoff:
		        res['reported_category'] = 'No retest - Mixed'
		        res['no_retest'] = True
		    else:
		        res['reported_category'] = 'Mixed'
		        res['no_retest'] = False	
		elif penalised_destim > self.selection_cutoff:
		    res['reported_category'] = 'No retest - No match'
		    res['no_retest'] = True
		else:
		    res['reported_category'] = 'Possible match'
		    res['no_retest'] = False	

		# and report, including a score	
# startup
if __name__ == '__main__':


	# command line usage.  Pass the location of a config file as a single argument.
	parser = argparse.ArgumentParser(
		formatter_class= argparse.RawTextHelpFormatter,
		description="""Runs findNeighbour3-varmod, a service for bacterial relatedness monitoring.
									 

Example usage: 
============== 
# show command line options 
python findNeighbour3-varmod.py --help 	

# run with debug settings; only do this for unit testing.
python findNeighbour3-varmod.py 	

# run using settings in myConfigFile.json.  
python findNeighbour3-varmod.py ../config/myConfigFile.json		

""")
	parser.add_argument('path_to_config_file', type=str, action='store', nargs='?',
						help='the path to the configuration file', default='')
	
	args = parser.parse_args()
	
	# an example config file is default_test_config.json

	############################ LOAD CONFIG ######################################
	print("findNeighbour3 distance estimator .. reading configuration file.")

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

	# determine whether a FN_SENTRY_URL environment variable is present,
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

	########################### prepare to launch server 
	print("Connecting to backend data store")
	try:
			PERSIST=fn3persistence(dbname = 
								CONFIG['SERVERNAME'],
						     	connString=CONFIG['FNPERSISTENCE_CONNSTRING'],
						        debug=CONFIG['DEBUGMODE'],
								server_monitoring_min_interval_msec = 0
								)
	except Exception as e:
			print("Error raised on creating persistence object")
			raise

	if CONFIG['SERVERNAME'] == 'fn3_unittesting':
		# no config file was supplied
		raise ValueError("No config file supplied.")

	# instantiate model builder class
	try:
		calib = Calibrator(CONFIG, PERSIST)
	except Exception as e:
			print("Error raised on instantiating findNeighbour3 distance estimator object")
			raise

	settings = calib.derive_settings(train_on=5000)		# provides default settings

	print("Running first ")
	validation = calib.test_settings(settings, train_on =5000, start_afresh = True)
	print("Running second")
	settings['probN_inflation_factor']=1
	settings['mixed_reporting_cutoff']=100
	validation = calib.test_settings(settings, train_on =5000, input_connstring = validation, start_afresh =False)
	print("Data is in ",validation)
