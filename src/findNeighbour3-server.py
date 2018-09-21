 #!/usr/bin/env python
""" 
A server providing relatedness information for bacterial genomes via a Restful API.

Implemented in pure Python3 3, it uses in-memory data storage backed by MongoDb.
It loads configuration from a config file, which must be set in production.

If no config file is provided, it will run in  'testing' mode with the  parameters
in default_test_config.json.  This expects a mongodb database to be running on
the default port on local host.  As a rough guide to the amount of space required in mongodb,
about 0.5MB of database is used per sequence, or about 2,000 sequences per GB.

All internal modules, and the restful API, are covered by unit testing.
Unit testing can be achieved by:

# starting a test RESTFUL server
python3 findNeighbour3-server-rest.py

# And then (e.g. in a different terminal) launching unit tests with
python3 -m unittest findNeighbour3-server

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
import pymongo
import pandas as pd
import numpy as np
import copy
import pathlib
import sentry_sdk
from sentry_sdk.integrations.flask import FlaskIntegration

# flask
from flask import Flask, make_response, jsonify, Markup
from flask import request, abort

# logging
from logging.config import dictConfig

# utilities for file handling and measuring file size
import psutil

# reference based compression, storage and clustering modules
from NucleicAcid import NucleicAcid
from mongoStore import fn3persistence
from seqComparer import seqComparer
from clustering import snv_clustering

# only used for unit testing
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide
import unittest
from urllib.parse import urlparse as urlparser
from urllib.parse import urljoin as urljoiner
import uuid

class findNeighbour3():
	""" a server based application for maintaining a record of bacterial relatedness using SNP distances.
	
	    The high level arrangement is that
		- This class interacts with in-memory sequences
		  [handled by the seqComparer class] and backends [fn3Persistance class] used by the server
		- methods in findNeighbour3() return native python3 objects.
		
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
            FNPERSISTENCE_CONNSTRING: a valid mongodb connection string. if shard keys are set, the 'guid' field is suitable key.
            MAXN_STORAGE:   The maximum number of Ns in the sequence <excluding those defined in > EXCLUDEFILE which should be indexed.
                            Other files, e.g. those with all Ns, will be tagged as 'invalid'.  Although a record of their presence in the database
                            is kept, they are not compared with other sequences.
            MAXN_PROP_DEFAULT: if the proportion not N in the sequence exceeds this, the sample is analysed, otherwise considered invalid.
            LOGFILE:        the log file used
            LOGLEVEL:		default logging level used by the server.  Valid values are DEBUG INFO WARNING ERROR CRITICAL
            SNPCEILING: 	links between guids > this are not stored in the database
            GC_ON_RECOMPRESS: if 'recompressing' sequences to a local reference, something the server does automatically, perform
                            a full mark-and-sweep gc at this point.  This setting alters memory use and compute time, but not the results obtained.
            RECOMPRESS_FREQUENCY: if recompressable records are detected, recompress every RECOMPRESS_FREQ th detection (e.g. 5).
                            Trades off compute time with mem usage.  This setting alters memory use and compute time, but not the results obtained.
							If zero, recompression is disabled.
            REPACK_FREQUENCY: how the matrix is stored in mongodb.  if REPACK_FREQ=0, there will be one document for every non-empty matrix cell.
			                if REPACK_FREQ>0, then if a guid has REPACK_FREQ-1 neighbours, then a 'repack' operation
							occurs.  This transfers multiple matrix cells into one mongodb document: essentially, part or all of a row
							will be packed into a single document.  This reduces query times, but the repack operation slows inserts.
							Repacking doesn't alter the results at all, and could be performed independently of inserts.
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
		"SNPCOMPRESSIONCEILING":250,
		"MAXN_PROP_DEFAULT":0.70,
		"LOGFILE":"../unittest_tmp/logfile.log",
		"LOGLEVEL":"INFO",
		"SNPCEILING": 20,
		"GC_ON_RECOMPRESS":1,
		"RECOMPRESS_FREQUENCY":5,
		"CLUSTERING":{'SNV12_ignore' :{'snv_threshold':12, 'mixed_sample_management':'ignore', 'mixture_criterion':'pvalue_1', 'cutoff':0.001},
		              'SNV12_include':{'snv_threshold':12, 'mixed_sample_management':'include', 'mixture_criterion':'pvalue_1', 'cutoff':0.001}
					 }
		}

		Some of these settings are read when the server is first-run, stored in a database, and the server will not
		change the settings on re-start even if the config file is changed.  Examples are:
		SNPCEILING
		MAXN_PROP_DEFAULT
		EXCLUDEFILE
		INPUTREF
		CLUSTERING
		These settings cannot be changed because they alter the way that the data is stored; if you want to change
		the settings, the data will have to be re-loaded. 
		
		However, most other settings can be changed and will take effect on server restart.  These include:
		server location
		IP
		SERVERNAME
		REST_PORT
		
		internal logging	
		LOGFILE
		LOGLEVEL
		
		where the database connection binds to
		FNPERSISTENCE_CONNSTRING
		
		related to internal server memory management:
		GC_ON_RECOMPRESS
		RECOMPRESS_FREQUENCY
		SNPCOMPRESSIONCEILING
		
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
						   'SNPCOMPRESSIONCEILING', "SNPCEILING", 'MAXN_PROP_DEFAULT', 'REST_PORT',
						   'LOGFILE','LOGLEVEL','GC_ON_RECOMPRESS','RECOMPRESS_FREQUENCY', 'REPACK_FREQUENCY', 'CLUSTERING'])
		missing=required_keys-set(self.CONFIG.keys())
		if not missing == set([]):
			raise KeyError("Required keys were not found in CONFIG. Missing are {0}".format(missing))

		# the following keys are not stored in any database backend, as a server could be moved, i.e.
		# running on the same data but with different IP etc
		
		do_not_persist_keys=set(['IP','SERVERNAME','FNPERSISTENCE_CONNSTRING',
								 'LOGFILE','LOGLEVEL','REST_PORT',
								 'GC_ON_RECOMPRESS','RECOMPRESS_FREQUENCY', 'REPACK_FREQUENCY'])
				
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
		self.snpCompressionCeiling = cfg['SNPCOMPRESSIONCEILING']
		self.maxn_prop_default = cfg['MAXN_PROP_DEFAULT']
		self.clustering_settings = cfg['CLUSTERING']
		self.recompress_frequency = self.CONFIG['RECOMPRESS_FREQUENCY']
		self.repack_frequency = self.CONFIG['REPACK_FREQUENCY']
		self.gc_on_recompress = self.CONFIG['GC_ON_RECOMPRESS']
			
		## start setup
		self.write_semaphore = threading.BoundedSemaphore(1)        # used to permit only one process to INSERT at a time.
		
		# initialise nucleic acid analysis object
		self.objExaminer=NucleicAcid()
		
		# load in-memory sequences
		self._load_in_memory_data()
		
		print("findNeighbour3 is ready.")
	
	def _load_in_memory_data(self):
		""" loads in memory data into the seqComparer object from database storage """
		# initialise seqComparer, which manages in-memory reference compressed data
		self.sc=seqComparer(reference=self.reference,
							maxNs=self.maxNs,
							snpCeiling= self.snpCeiling,
							debugMode=self.debugMode,
							excludePositions=self.excludePositions,
							snpCompressionCeiling = self.snpCompressionCeiling)
		
		# determine how many guids there in the database
		guids = self.PERSIST.refcompressedsequence_guids()
	
		self.server_monitoring_store(message='Starting load of sequences into memory from database')

		print("Loading {1} sequences from database .. excluding ({0})".format(self.sc.excluded_hash(),len(guids)))
		nLoaded = 0
		for guid in guids:
			nLoaded+=1
			obj = self.PERSIST.refcompressedsequence_read(guid)
			self.sc.persist(obj, guid=guid)
			if nLoaded % 500 ==0:
				print(nLoaded)
				self.server_monitoring_store(message='Loaded {0} from database'.format(nLoaded))

		print("findNeighbour3 has loaded {0} sequences from database.".format(len(guids)))
		self.server_monitoring_store(message='Loaded {0} from database; load complete.'.format(nLoaded))

		# set up clustering
		print("findNeighbour3 is updating clustering.")

		self.clustering={}		# a dictionary of clustering objects, one per SNV cutoff/mixture management setting
		for clustering_name in self.clustering_settings.keys():
			self.server_monitoring_store(message='Loading clustering data into memory for {0}'.format(clustering_name))
			json_repr = self.PERSIST.clusters_read(clustering_name)
			self.clustering[clustering_name] = snv_clustering(saved_result =json_repr)
		
		# ensure that clustering object is up to date.  clustering is an in-memory graph, which is periodically
		# persisted to disc.  It is possible that, if the server crashes/does a disorderly shutdown,
		# the clustering object which is persisted might not include all the guids in the reference compressed
		# database.  This situation is OK, because the clustering object will bring itself up to date when
		# the new guids and their links are loaded into it.
		self.update_clustering()

	def reset(self):
		""" restarts the server, deleting any existing data """
		if not self.debugMode == 2:
			return		 # no action taken by calls to this unless debugMode ==2
		else:
			print("Deleting existing data and restarting")
			self.PERSIST._delete_existing_data()
			self._load_in_memory_data()
			
	def server_monitoring_store(self, message="No message supplied"):
		""" reports server memory information to store """
		self.PERSIST.server_monitoring_store(message=message, content=self.sc.summarise_stored_items())


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

		# create clusters objects
		app.logger.info("Creating clustering objects..")
		
		self.clustering = {}
		expected_clustering_config_keys = set(['snv_threshold',  'mixed_sample_management', 'cutoff', 'mixture_criterion'])
		for clustering_name in self.clustering_settings.keys():
			observed = self.clustering_settings[clustering_name] 
			if not observed.keys() == expected_clustering_config_keys:
				raise KeyError("Got unexpected keys for clustering setting {0}: got {1}, expected {2}".format(clustering_name, observed, expected_clustering_config_keys))
			self.clustering[clustering_name] = snv_clustering(snv_threshold=observed['snv_threshold'] , mixed_sample_management=observed['mixed_sample_management'])
			self.PERSIST.clusters_store(clustering_name, self.clustering[clustering_name].to_dict())
			
		# persist other config settings.
		for item in self.CONFIG.keys():
			if not item in do_not_persist_keys:
				config_settings[item]=self.CONFIG[item]
				
		res = self.PERSIST.config_store('config',config_settings)
		app.logger.info("First run actions complete.")
		
	def repack(self,guids=None):
		""" generates a smaller and faster representation in the persistence store
		for the guids in the list. optional"""
		if guids is None:
			guids = self.PERSIST.guids()  # all the guids
		for this_guid in guids:
			app.logger.debug("Repacking {0}".format(this_guid))
			self.PERSIST.guid2neighbour_repack(this_guid)
	
	def record_server_status(self, message='No message supplied'):
		""" stores server status to database.
		Useful for process monitoring """
		self.server_monitoring_store(message = message)
		
	def insert(self,guid,dna):
		""" insert DNA called guid into the server,
		persisting it in both RAM and on disc, and updating any clustering.
		"""
		
		# clean, and provide summary statistics for the sequence
		app.logger.info("Preparing to insert: {0}".format(guid))
		if not self.sc.iscachedinram(guid):                   # if the guid is not already there
			
			# prepare to insert
			self.objExaminer.examine(dna)  					  # examine the sequence
			cleaned_dna=self.objExaminer.nucleicAcidString.decode()
			refcompressedsequence =self.sc.compress(cleaned_dna)          # compress it and store it in RAM
			self.server_monitoring_store(message='Guid {0} is about to be stored in ram'.format(guid))
			self.sc.persist(refcompressedsequence, guid)			    # insert the DNA sequence into ram.
			self.server_monitoring_store(message='Guid {0} stored in ram'.format(guid))
						
			# construct links with everything existing existing at the time the semaphore was acquired.
			self.write_semaphore.acquire()				    # addition should be an atomic operation

			links={}			
			try:
				# this process reports links less than self.sc.snpCeiling
				app.logger.debug("Finding links: {0}".format(guid))

				to_compress = 0
				for key2 in self.sc.guidscachedinram():
					if not guid==key2:
						(guid1,guid2,dist,n1,n2,nboth, N1pos, N2pos, Nbothpos)=self.sc.countDifferences_byKey(keyPair=(guid,key2),
																											  cutoff = self.snpCompressionCeiling)
					
						link = {'dist':dist,'n1':n1,'n2':n2,'nboth':nboth}
						to_compress +=1
						if dist is not None:
							if link['dist'] <= self.snpCeiling:
								links[guid2]=link			


				## now persist in database.  
				# we have considered what happens if database connectivity fails during the insert operations.
				app.logger.info("Persisting: {0}".format(guid))

				# if the database connectivity fails after this refcompressedseq_store has completed, then 
				# the 'document' will already exist within the mongo file store.
				# in such a case, a FileExistsError is raised.
				# we trap for such errors, logging a warning, but permitting continuing execution since
				# this is expected if refcompressedseq_store succeeds, but subsequent inserts fail.
				try:
					self.PERSIST.refcompressedseq_store(guid, refcompressedsequence)     # store the parsed object to database
				except FileExistsError:
					app.logger.warning("Attempted to refcompressedseq_store {0}, but it already exists.  This is expected only if database connectivity failed during a previous INSERT operation.  Such failures should be noted in earlier logs".format(guid))
				except Exception: 		# something else
					raise			# we don't want to trap other things

				# annotation of guid will update if an existing record exists.  This is OK, and is acceptable if database connectivity failed during previous inserts
				self.PERSIST.guid_annotate(guid=guid, nameSpace='DNAQuality',annotDict=self.objExaminer.composition)						

				# addition of neighbours may cause neighbours to be entered more than once if database connectivity failed during previous inserts.
				# because of the way that extraction of links works, this does not matter, and duplicates will not be reported.
				self.PERSIST.guid2neighbour_add_links(guid=guid, targetguids=links)

			except Exception as e:
				app.logger.exception("Error raised on persisting {0}".format(guid))
				self.write_semaphore.release() 	# ensure release of the semaphore if an error is trapped

				# Rollback anything which could leave system in an inconsistent state
				# remove the guid from RAM is the only step necessary
				self.sc.remove(guid)	
				app.logger.info("Guid successfully removed from ram. {0}".format(guid))
                                    
				if e.__module__ == "pymongo.errors":
					app.logger.info("Error raised pertains to pyMongo connectivity")
					abort(503,e)		# the mongo server may be refusing connections, or busy.  This is observed occasionally in real-world use
				else:
					abort(500,e)		# some other kind of error

				
			# release semaphore
			self.write_semaphore.release()                  # release the write semaphore

			if self.recompress_frequency > 0:				
				if to_compress>= self.recompress_frequency and to_compress % self.recompress_frequency == 0:		# recompress if there are lots of neighbours, every self.recompress_frequency isolates
					app.logger.debug("Recompressing: {0}".format(guid))
					self.server_monitoring_store(message='Guid {0} being recompressed relative to neighbours'.format(guid))
					self.sc.compress_relative_to_consensus(guid)
				if self.gc_on_recompress==1:
					gc.collect()
			app.logger.info("Insert succeeded {0}".format(guid))

			# clean up guid2neighbour; this can readily be done post-hoc, if the process proves to be slow.
			# it is a mongodb reformatting operation which doesn't affect results.

			guids = list(links.keys())
			guids.append(guid)
		
			if self.repack_frequency>0:
				app.logger.info("Repacking around: {0}".format(guid))
				if len(guids) % self.repack_frequency ==0:		# repack if there are repack_frequency-1 neighbours
					self.repack(guids)
			
			# cluster
			app.logger.info("Clustering around: {0}".format(guid))
			self.update_clustering()

			return "Guid {0} inserted.".format(guid)		
		else:
			return "Guid {0} is already present".format(guid)
			app.logger.info("Already present, no insert needed: {0}".format(guid))
	
	def update_clustering(self, store=True):
		""" performs clustering on any samples within the persistence store which are not already clustered
		    If Store=True, writes the clustered object to mongo."""
		
		# update clustering and re-cluster
		for clustering_name in self.clustering_settings.keys():
			
			# ensure that clustering object is up to date.  clustering is an in-memory graph, which is periodically
			# persisted to disc.  It is possible that, if the server crashes/does a disorderly shutdown,
			# the clustering object which is persisted might not include all the guids in the reference compressed
			# database.  This situation is OK, because the clustering object will bring itself up to date when
			# the new guids and their links are loaded into it.
			
			guids = self.PERSIST.refcompressedsequence_guids()			# all guids processed and refernece compressed
			in_clustering_guids = self.clustering[clustering_name].guids()  # all clustered guids
			to_add_guids = guids - in_clustering_guids					# what we need to add
			remaining_to_add_guids = copy.copy(to_add_guids)				# we iterate until there's nothing left to add
			logging.info("Clustering graph {0} contains {2} guids out of {1}; updating.".format(clustering_name, len(guids), len(in_clustering_guids)))
			while len(remaining_to_add_guids)>0:
				to_add_guid = remaining_to_add_guids.pop()				# get the guid
				links = self.PERSIST.guid2neighbours(to_add_guid, returned_format=3)['neighbours']	# and its links	
				self.clustering[clustering_name].add_sample(to_add_guid, links)		# add it to the clustering db
				
			# check any clusters to which to_add_guids have been added for mixtures.
			nMixed = 0
			guids_to_check = set()
			clusters_to_check = set()
			
			for guid in to_add_guids:
				guids_to_check.add(guid)								# then we need to check it and
				for cluster in self.clustering[clustering_name].guid2clusters(guid):
					clusters_to_check.add(cluster)						# everything else in the same cluster as it
			
			cl2guids = 	self.clustering[clustering_name].clusters2guid()	# dictionary allowing cluster -> guid lookup
			for cluster in clusters_to_check:
				guids_for_msa = cl2guids[cluster]							# do msa on the cluster
				logging.debug("** Checking cluster {0}; performing MSA on {1} samples".format(cluster,len(guids_for_msa)))
				msa = self.sc.multi_sequence_alignment(guids_for_msa, output='df')		#  a pandas dataframe; p_value tests mixed
				logging.debug("** Multi sequence alignment is complete")

				if not msa is None:		# no alignment was made
					mixture_criterion = self.clustering_settings[clustering_name]['mixture_criterion']
					mixture_cutoff = self.clustering_settings[clustering_name]['cutoff']
					##################################################################################################
					## NOTE: query_criterion is EVAL'd by pandas.  This potentially a route to attack the server, but
					# to do so requires that you can edit either the CONFIG file (pre-startup) or the Mongodb in which
					# the config is stored post first-run
					##################################################################################################
					query_criterion = "{0} <= {1}".format(mixture_criterion,mixture_cutoff)
					
					msa_mixed = msa.query(query_criterion)
					for mixed_guid in msa_mixed.index:
						if not self.clustering[clustering_name].is_mixed(mixed_guid):   # if it's not known to be mixed
							# set it as mixed
							links = self.PERSIST.guid2neighbours(mixed_guid, returned_format=3)['neighbours']
							self.clustering[clustering_name].set_mixed(mixed_guid, neighbours = links)
						
			in_clustering_guids = self.clustering[clustering_name].guids()
			logging.info("Cluster {0} updated; now contains {1} guids. ".format(clustering_name, len(in_clustering_guids)))
			
			if store==True:
				self.PERSIST.clusters_store(clustering_name, self.clustering[clustering_name].to_dict())
				logging.info("Cluster {0} persisted".format(clustering_name))
			
	def exist_sample(self,guid):
		""" determine whether the sample exists in RAM"""
		
		## this call measures presence on disc
		return self.PERSIST.guid_exists(guid)

	def server_time(self):
		""" returns the current server time """
		return {"server_time":datetime.datetime.now().isoformat()}

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
		return {"exclusion_id":self.sc.excluded_hash(), "excluded_nt":list(self.sc.excluded)}
	
	def server_memory_usage(self, max_reported=None):
		""" reports recent server memory activity """
		if max_reported is None:
			max_reported =100		# a default
		return self.PERSIST.recent_server_monitoring(max_reported= max_reported)
	
	def neighbours_within_filter(self, guid, snpDistance, cutoff=0.85, returned_format=1):
		""" returns a list of guids, and their distances, by a sample quality cutoff	
			returns links either as
			format 1 [otherGuid, distance]
			or as
			format 2 [otherGuid, distance, N_just1, N_just2, N_either]
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
		# wrapper round ewc.neighboursOf()
		# returned value reVal is {'guid':guid, 'neighbours':[ *, *, * ...] }
		# where * is one of the two formats shown in the docstring above.
		
		retVal = self.PERSIST.guid2neighbours(guid=guid, cutoff=snpDistance, returned_format=returned_format)
		
		# run a quality check on the things our sample is like.
		sampleList=retVal['neighbours']
		idList=[]
		for sa in sampleList:
			idList.append(sa[0])		# add the guid
		
		# get the sample qualities from the database
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
		for item in sampleList:
			if returned_format == 1:
				item=[item[0],item[1]]		# just the guid and the distance;
			# otherwise return the whole of item	
			if item[0] in goodGuids:
				finalOutput.append(item)
				
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
	
	def query_get_detail(self, sname1, sname2):
		""" gets detail on the comparison of a pair of samples.  Computes this on the fly """
		ret = self.sc.query_get_detail(sname1,sname2)
		return(ret)

	def sequence(self, guid):
		""" gets masked sequence for the guid, in format sequence|fasta """
		if not self.sc.iscachedinram(guid):
			return None
		try:		
			seq = self.sc.uncompress(self.sc.seqProfile[guid])
			return {'guid':guid, 'invalid':0,'comment':'Masked sequence, as stored','masked_dna':seq}
		except ValueError:
				return {'guid':guid, 'invalid':1,'comment':'No sequence is available, as invalid sequences are not stored'}
			
# default parameters for unit testing only.
RESTBASEURL   = "http://127.0.0.1:5020"
ISDEBUG = True
LISTEN_TO = '127.0.0.1'		# only local addresses

# initialise Flask 
app = Flask(__name__)
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
    fn3.PERSIST.closedown()		# close database connection

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

@app.route('/', methods=['GET'])
def server_info():
	""" returns server info page
	"""
	res = """findNeighbour3 web server operating.<p>Endpoints are in rest-routes.md, in the docs<p>"""
	return make_response(res)

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

	if not fn3.debugMode == 2:
		# if we're not in debugMode==2, then this option is not allowed
		abort(404, 'Calls to /reset are only allowed with debugMode == 2' )

	if component == 'main':
		raise ZeroDivisionError(token)
	elif component == 'clustering':
		clustering_names = list(fn3.clustering_settings.keys())

		if len(clustering_names)==0:
			self.fail("no clustering settings defined; cannot test error generation in clustering")
		else:
			clustering_name = clustering_names[0]
			fn3.clustering_settings[clustering_name].raise_error(token)
	elif component == 'seqcomparer':
		fn3.sc.raise_error(token)
	elif component == 'persist':
		fn3.PERSIST.raise_error(token)
	else:
		raise KeyError("Invalid component called.  Allowed: main;persist;clustering;seqcomparer.")

@unittest.skip("skipped; known issue with error handling within flask")		
class test_raise(unittest.TestCase):
	""" tests route /api/v2/reset
	
	Note: this test currently fails, and has been disabled.
	It appears that (at least as currently configured) errors raised during
	Flask execution are not logged to app.logger.
	
	This is unexpected; the errors raised are printed to STDERR and
	are also logged using Sentry, if configured.
	
	The logger is working, and explicit calls to app.logger.exception() within try/except blocks
	do log.
	
	This remains an unresolved issue.
	"""
	def runTest(self):
		
		# get the server's config - requires that we're running in debug mode
		relpath = "/api/v2/server_config"
		try:
			res = do_GET(relpath)
		except requests.exceptions.HTTPError:
			self.fail("Could not read config. This unit test requires a server in debug mode")
			
		self.assertTrue(isjson(content = res.content))

		config_dict = json.loads(res.content.decode('utf-8'))
		try:
			logfile = config_dict['LOGFILE']
		except KeyError:
			self.fail("No LOGFILE element in config dictionary, as obtained from the server.")
			
		for error_at in ['main','persist','clustering','seqcomparer']:
			guid = uuid.uuid4().hex
			print(guid)
			token = "TEST_ERROR_in_{0}_#_{1}".format(error_at, guid)
			relpath = "/api/v2/raise_error/{0}/{1}".format(error_at, token)
			print(relpath)
			res = do_GET(relpath)
	
			if not os.path.exists(logfile):
				self.fail("No logfile {0} exists.  This test only works when the server is on localhost, and debugmode is 2".format(logfile))
			else:
				with open(logfile, 'rt') as f:
					txt = f.read()
					if not token in txt:
						print("NOT LOGGED: ***** {0} <<<<<<".format(txt[-200:]))
						self.fail("Error was not logged {0}".format(error_at))

	
@app.route('/api/v2/assess_mixed', methods=['POST'])
def assess_mixed():
	""" computes estimates of whether *this_guid* is likely to be mixed, relative to the guids in the
	semicolon separated list *related_guids*.
	
	The payload expected is a dictionary like htis:
	{'this_guid':'guid0',
	 'related_guids':'guid1;guid2;guid3',
	 'sample_size':30}
	
	The strategy used is draw at most max_sample_size unique related_guids, and from them
		analyse all (max_sample_size * (max_sample_size-1))/2 unique pairs.
		For each pair, we determine where they differ, and then
		estimate the proportion of mixed bases in those variant sites.

        Pairs of related_guids which do not differ are uninformative and are ignored.

        The output is a pandas dataframe containing mixture estimates for this_guid for each of a series of pairs.

        The p values reported are derived from exact, two-sided binomial tests as implemented in pythons scipy.stats.binom_test().
        
        TEST 1:
        This tests the hypothesis that the number of Ns in the *alignment*
        is GREATER than those expected from the expected_N in the population of whole sequences.
 
        Does so by comparing the observed number of Ns in the alignment (alignN),
        given the alignment length (4 in the above case) and an expectation of the proportion of bases which will be N.
        The expected number of Ns is estimated by randomly sampling sample_size guids from those stored in the server and
        observing the number of Ns per base.  The estimate_expected_N() function performs this.
        
        This approach determines the median number of Ns in valid sequences, which (if bad samples with large Ns are rare)
        is a relatively unbiased estimate of the median number of Ns in the good quality samples.
        
        If there  are not enough samples in the server to obtain an estimate, p_value is not computed, being
        reported as None.
        
        TEST 2: tests whether the proportion of Ns in the alignment is greater
        than in the bases not in the alignment, for this sequence.

	"""

	# validate input
	request_payload = request.form.to_dict()
	if 'this_guid' in request_payload.keys() and 'related_guids' in request_payload.keys():
		if len(request_payload['related_guids'])==0:
			# no related guids
			return make_response(json.dumps({}))
		
		related_guids = request_payload['related_guids'].split(';')		# coerce both guid and seq to strings
		this_guid= request_payload['this_guid']
	else:
		abort(501, 'this_guid and related_guids are not present in the POSTed data {0}'.format(request_payload.keys()))

	# check guids
	missing_guids = []
	for guid in related_guids:
		try:
			result = fn3.exist_sample(guid)
		except Exception as e:
			abort(500, e)
		if result is False:
			missing_guids.append(guid)
	
	if len(missing_guids)>0:
		abort(501, "asked to perform mixture assessment relative to the following missing guids: {0}".format(missing_guids))
	
	# data validation complete.  construct outputs
	res = fn3.sc.assess_mixed(
		this_guid= this_guid,
		related_guids = related_guids,
		max_sample_size = 10)
	
	if res is None:
		retVal = {}
		
	else:
		retVal = pd.DataFrame.to_dict(res,orient='index')

	return make_response(json.dumps(retVal))

	
class test_assess_mixed_1(unittest.TestCase):
	""" tests route /api/v2/assess_mixed
 		- due for removal
	"""
	def runTest(self):
		relpath = "/api/v2/guids"
		res = do_GET(relpath)
		n_pre = len(json.loads(str(res.text)))		# get all the guids

		inputfile = "../COMPASS_reference/R39/R00000039.fasta"
		with open(inputfile, 'rt') as f:
			for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
					originalseq = list(str(record.seq))
		inserted_guids = []
		
		for i in range(0,50):
			guid_to_insert = "guid_{0}".format(n_pre+i)
			inserted_guids.append(guid_to_insert)
			
			seq = originalseq			
			# make i mutations at position 500,000
			offset = 500000
			for j in range(i):
				mutbase = offset+j
				ref = seq[mutbase]
				if not ref == 'T':
					seq[mutbase] = 'T'
				if not ref == 'A':
					seq[mutbase] = 'A'
			seq = ''.join(seq)
						
			print("Adding TB sequence {2} of {0} bytes with {1} mutations relative to ref.".format(len(seq), i, guid_to_insert))
			self.assertEqual(len(seq), 4411532)		# check it's the right sequence
	
			relpath = "/api/v2/insert"
			res = do_POST(relpath, payload = {'guid':guid_to_insert,'seq':seq})
			self.assertTrue(isjson(content = res.content))
			info = json.loads(res.content.decode('utf-8'))
			self.assertEqual(info, 'Guid {0} inserted.'.format(guid_to_insert))
	
		relpath = "/api/v2/assess_mixed"
		payload = {'this_guid':inserted_guids[0], 'related_guids':';'.join(inserted_guids[1:3])}
		res = do_POST(relpath, payload=payload)
		self.assertTrue(isjson(res.content))
		self.assertEqual(res.status_code, 200)
		d = json.loads(res.content.decode('utf-8'))
		df = pd.DataFrame.from_dict(d,orient='index')
		self.assertEqual(df.index.tolist(), [inserted_guids[0]])

		relpath = "/api/v2/assess_mixed"
		payload = {'this_guid':inserted_guids[0], 'related_guids':';'.join(inserted_guids[1:20])}
		res = do_POST(relpath, payload=payload)
		self.assertTrue(isjson(res.content))
		self.assertEqual(res.status_code, 200)
		d = json.loads(res.content.decode('utf-8'))
		df = pd.DataFrame.from_dict(d,orient='index')
		
def construct_msa(guids, output_format):
	""" constructs multiple sequence alignment for guids
	    and returns in one of 'fasta' 'html' or 'json' format."""
	res = fn3.sc.multi_sequence_alignment(guids, output='df_dict')
	df = pd.DataFrame.from_dict(res,orient='index')
	html = df.to_html()
	fasta= ""
	for guid in df.index:
		fasta=fasta + ">{0}\n{1}\n".format(guid, df.loc[guid,'aligned_seq'])
		
	if output_format == 'fasta':
		return make_response(fasta)
	elif output_format == 'html':
		return make_response(html)
	elif output_format == 'json':
		return make_response(json.dumps(res))

@app.route('/api/v2/reset', methods=['POST'])
def reset():
	""" deletes any existing data from the server """
	if not fn3.debugMode == 2:
		# if we're not in debugMode==2, then this option is not allowed
		abort(404, 'Calls to /reset are only allowed with debugMode == 2' )
	else:
		fn3.reset()
		return make_response(json.dumps({'message':'reset completed'}))
class test_reset(unittest.TestCase):
	""" tests route /api/v2/reset
	"""
	def runTest(self):
		relpath = "/api/v2/guids"
		res = do_GET(relpath)
		n_pre = len(json.loads(str(res.text)))		# get all the guids
		
		guid_to_insert = "guid_{0}".format(n_pre+1)
		
		inputfile = "../COMPASS_reference/R39/R00000039.fasta"
		with open(inputfile, 'rt') as f:
			for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
					seq = str(record.seq)
		
		relpath = "/api/v2/insert"
		res = do_POST(relpath, payload = {'guid':guid_to_insert,'seq':seq})
		self.assertTrue(isjson(content = res.content))
		info = json.loads(res.content.decode('utf-8'))
		self.assertEqual(info, 'Guid {0} inserted.'.format(guid_to_insert))
		
		relpath = "/api/v2/guids"
		res = do_GET(relpath)
		n_post = len(json.loads(str(res.text)))		# get all the guids
		
		relpath = "/api/v2/reset"
		res = do_POST(relpath, payload={})
		
		relpath = "/api/v2/guids"
		res = do_GET(relpath)
		n_post_reset = len(json.loads(str(res.text)))		# get all the guids
			
		self.assertTrue(n_post>0)
		self.assertTrue(n_post_reset==0)
@app.route('/api/v2/multiple_alignment/guids', methods=['POST'])
def msa_guids():
	""" performs a multiple sequence alignment on a series of POSTed guids,
	delivered in a dictionary, e.g.
	{'guids':'guid1;guid2;guid3',
	'output_format':'json'}
	
	Valid values for format are:
	json
	html
	"""

	# validate input
	request_payload = request.form.to_dict()
	if 'output_format' in request_payload.keys() and 'guids' in request_payload.keys():
		guids = request_payload['guids'].split(';')		# coerce both guid and seq to strings
		output_format= request_payload['output_format']
		if not output_format in ['html','json','fasta']:
			abort(501, 'output_format must be one of html, json, or fasta not {0}'.format(format))
	else:
		abort(501, 'output_format and guids are not present in the POSTed data {0}'.format(data_keys))
	
	# check guids
	missing_guids = []
	for guid in guids:
		try:
			result = fn3.exist_sample(guid)
		except Exception as e:
			abort(500, e)
		if not result is True:
			missing_guids.append(guid)
	
	if len(missing_guids)>0:
		abort(501, "asked to perform multiple sequence alignment with the following missing guids: {0}".format(missing_guids))
	
	# data validation complete.  construct outputs
	return construct_msa(guids, output_format)


class test_msa_2(unittest.TestCase):
	""" tests route /api/v2/multiple_alignment/guids, with additional samples.
	"""
	def runTest(self):
		relpath = "/api/v2/guids"
		res = do_GET(relpath)
		n_pre = len(json.loads(str(res.text)))		# get all the guids

		inputfile = "../COMPASS_reference/R39/R00000039.fasta"
		with open(inputfile, 'rt') as f:
			for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
					originalseq = list(str(record.seq))
		inserted_guids = []			
		for i in range(0,3):
			guid_to_insert = "guid_{0}".format(n_pre+i)
			inserted_guids.append(guid_to_insert)
			
			seq = originalseq			
			# make i mutations at position 500,000
			offset = 500000
			for j in range(i):
				mutbase = offset+j
				ref = seq[mutbase]
				if not ref == 'T':
					seq[mutbase] = 'T'
				if not ref == 'A':
					seq[mutbase] = 'A'
			seq = ''.join(seq)
						
			print("Adding TB sequence {2} of {0} bytes with {1} mutations relative to ref.".format(len(seq), i, guid_to_insert))
			self.assertEqual(len(seq), 4411532)		# check it's the right sequence
	
			relpath = "/api/v2/insert"
			res = do_POST(relpath, payload = {'guid':guid_to_insert,'seq':seq})
			self.assertTrue(isjson(content = res.content))
			info = json.loads(res.content.decode('utf-8'))
			self.assertEqual(info, 'Guid {0} inserted.'.format(guid_to_insert))
	
		relpath = "/api/v2/multiple_alignment/guids"
		payload = {'guids':';'.join(inserted_guids),'output_format':'html'}
		res = do_POST(relpath, payload=payload)
		self.assertFalse(isjson(res.content))
		self.assertEqual(res.status_code, 200)
		self.assertTrue(b"</table>" in res.content)
		
		
		payload = {'guids':';'.join(inserted_guids),'output_format':'json'}
		res = do_POST(relpath, payload=payload)
		self.assertTrue(isjson(res.content))
		self.assertEqual(res.status_code, 200)
		self.assertFalse(b"</table>" in res.content)
		d = json.loads(res.content.decode('utf-8'))
		not_present = set(inserted_guids) - set(d.keys())
		self.assertEqual(not_present, set())

		payload = {'guids':';'.join(inserted_guids),'output_format':'fasta'}
		res = do_POST(relpath, payload=payload)
		self.assertFalse(isjson(res.content))
		self.assertEqual(res.status_code, 200)
	
		relpath = "/api/v2/clustering/SNV12_ignore/guids2clusters"
		res = do_GET(relpath)
		self.assertEqual(res.status_code, 200)
		retVal = json.loads(str(res.text))
		self.assertTrue(isinstance(retVal, list))
		res = json.loads(res.content.decode('utf-8'))
		cluster_id=None
		for item in res:
			if item['guid'] in inserted_guids:
				cluster_id = item['cluster_id']
		self.assertTrue(cluster_id is not None)
			
		relpath = "/api/v2/multiple_alignment_cluster/SNV12_ignore/{0}/json".format(cluster_id)
		res = do_GET(relpath)
		self.assertTrue(isjson(res.content))
		self.assertEqual(res.status_code, 200)
		d = json.loads(res.content.decode('utf-8'))
		self.assertEqual(set(inserted_guids)-set(d.keys()),set([]))

		relpath = "/api/v2/multiple_alignment_cluster/SNV12_ignore/{0}/fasta".format(cluster_id)
		res = do_GET(relpath)
		self.assertFalse(isjson(res.content))
		self.assertEqual(res.status_code, 200)

@app.route('/api/v2/multiple_alignment_cluster/<string:clustering_algorithm>/<int:cluster_id>/<string:output_format>',methods=['GET'])
def msa_guids_by_cluster(clustering_algorithm, cluster_id, output_format):
	""" performs a multiple sequence alignment on the contents of a cluster
	
	Valid values for format are:
	json
	html
	"""
	# validate input
	try:
		res = fn3.clustering[clustering_algorithm].clusters2guidmeta(after_change_id = None)		
	except KeyError:
		# no clustering algorithm of this type
		return make_response(tojson("no clustering algorithm {0}".format(clustering_algorithm)), 404)
		
	if not output_format in ['html','json','fasta']:
		abort(501, 'output_format must be one of html, json, or fasta not {0}'.format(format))

	# check guids
	df = pd.DataFrame.from_records(res)
	if len(df.index)==0:
		retVal = {}
	else:
		missing_guids = []
		guids = []
		for guid in df['guid'].tolist():
			try:
				result = fn3.exist_sample(guid)
			except Exception as e:
				abort(500, e)
			if not result is True:
				missing_guids.append(guid)
			else:
				guids.append(guid)
		
		if len(missing_guids)>0:
			abort(501, "asked to perform multiple sequence alignment with the following missing guids: {0}".format(missing_guids))
			
		# data validation complete.  construct outputs
		return construct_msa(guids, output_format)


class test_msa_1(unittest.TestCase):
	""" tests route /api/v2/multiple_alignment/guids, with additional samples.
	"""
	def runTest(self):
		relpath = "/api/v2/guids"
		res = do_GET(relpath)
		n_pre = len(json.loads(str(res.text)))		# get all the guids

		inputfile = "../COMPASS_reference/R39/R00000039.fasta"
		with open(inputfile, 'rt') as f:
			for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
					originalseq = list(str(record.seq))
		inserted_guids = []			
		for i in range(0,3):
			guid_to_insert = "guid_{0}".format(n_pre+i)
			inserted_guids.append(guid_to_insert)
			
			seq = originalseq			
			# make i mutations at position 500,000
			offset = 500000
			for j in range(i):
				mutbase = offset+j
				ref = seq[mutbase]
				if not ref == 'T':
					seq[mutbase] = 'T'
				if not ref == 'A':
					seq[mutbase] = 'A'
			seq = ''.join(seq)
						
			print("Adding TB sequence {2} of {0} bytes with {1} mutations relative to ref.".format(len(seq), i, guid_to_insert))
			self.assertEqual(len(seq), 4411532)		# check it's the right sequence
	
			relpath = "/api/v2/insert"
			res = do_POST(relpath, payload = {'guid':guid_to_insert,'seq':seq})
			self.assertEqual(res.status_code, 200)
		
			self.assertTrue(isjson(content = res.content))
			info = json.loads(res.content.decode('utf-8'))
			self.assertEqual(info, 'Guid {0} inserted.'.format(guid_to_insert))
	
		relpath = "/api/v2/multiple_alignment/guids"
		payload = {'guids':';'.join(inserted_guids),'output_format':'html'}
		res = do_POST(relpath, payload=payload)
		self.assertFalse(isjson(res.content))
		self.assertEqual(res.status_code, 200)
		self.assertTrue(b"</table>" in res.content)
		
		
		payload = {'guids':';'.join(inserted_guids),'output_format':'json'}
		res = do_POST(relpath, payload=payload)
		self.assertTrue(isjson(res.content))
		self.assertEqual(res.status_code, 200)
		self.assertFalse(b"</table>" in res.content)
		d = json.loads(res.content.decode('utf-8'))
		self.assertEqual(set(d.keys()), set(inserted_guids))

		payload = {'guids':';'.join(inserted_guids),'output_format':'fasta'}
		res = do_POST(relpath, payload=payload)
		self.assertFalse(isjson(res.content))
		self.assertEqual(res.status_code, 200)


@app.route('/api/v2/server_config', methods=['GET'])
def server_config():
    """ returns server configuration.

        returns the config file with which the server was launched.
        This may be highly undesirable,
        as it reveals the internal server architecture  including
        backend databases and perhaps connection strings with passwords.

    """
    res = fn3.server_config()
    if res is None:		# not allowed to see it
        return make_response(tojson({'NotAvailable':"Endpoint is only available in debug mode"}), 404)
    else:
        return make_response(tojson(CONFIG))

class test_server_config(unittest.TestCase):
    """ tests route v2/server_config"""
    def runTest(self):
        relpath = "/api/v2/server_config"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))

        config_dict = json.loads(res.content.decode('utf-8'))
       
        self.assertTrue('GC_ON_RECOMPRESS' in config_dict.keys())
        self.assertEqual(res.status_code, 200)


@app.route('/api/v2/server_memory_usage', defaults={'nrows':100}, methods=['GET'])
@app.route('/api/v2/server_memory_usage/<int:nrows>', methods=['GET'])
def server_memory_usage(nrows):
	""" returns server memory usage information, as list.
	The server notes memory usage at various key points (pre/post insert; pre/post recompression)
	and these are stored. """
	try:
		result = fn3.server_memory_usage(max_reported = nrows)

	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
		
	return make_response(tojson(result))

class test_server_memory_usage(unittest.TestCase):
    """ tests route /api/v2/server_memory_usage"""
    def runTest(self):
        relpath = "/api/v2/server_memory_usage"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content = res.content))

        res = json.loads(res.content.decode('utf-8'))
        self.assertTrue(isinstance(res,list))



@app.route('/api/v2/server_time', methods=['GET'])
def server_time():
	""" returns server time """
	try:
		result = fn3.server_time()
		
	except Exception as e:
		abort(500, e)
	return make_response(tojson(result))

class test_server_time(unittest.TestCase):
    """ tests route /api/v2/server_time"""
    def runTest(self):
        relpath = "/api/v2/server_time"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        config_dict = json.loads(res.content.decode('utf-8'))
        self.assertTrue('server_time' in config_dict.keys())
        self.assertEqual(res.status_code, 200)

	
@app.route('/api/v2/guids', methods=['GET'])
def get_all_guids(**debug):
	""" returns all guids.  reference, if included, is ignored."""
	try:
		result = list(fn3.get_all_guids())
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
	return(make_response(tojson(result)))

class test_get_all_guids_1(unittest.TestCase):
    """ tests route /api/v2/guids"""
    def runTest(self):
        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        guidlist = json.loads(str(res.content.decode('utf-8')))
        self.assertTrue(isinstance(guidlist, list))
        self.assertEqual(res.status_code, 200)
        ## TODO: insert guids, check it doesn't fail.

@app.route('/api/v2/guids_with_quality_over/<float:cutoff>', methods=['GET'])
@app.route('/api/v2/guids_with_quality_over/<int:cutoff>', methods=['GET'])
def guids_with_quality_over(cutoff, **kwargs):
	""" returns all guids with quality score >= cutoff."""
	try:
		result = fn3.guids_with_quality_over(cutoff)	
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
	return make_response(tojson(result))

class test_guids_with_quality_over_1(unittest.TestCase):
    """ tests route /api/v2/guids_with_quality_over"""
    def runTest(self):
        relpath = "/api/v2/guids_with_quality_over/0.7"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        guidlist = json.loads(res.content.decode('utf-8'))
        self.assertTrue(isinstance(guidlist, list))
        self.assertEqual(res.status_code, 200)
        

@app.route('/api/v2/guids_and_examination_times', methods=['GET'])
def guids_and_examination_times(**kwargs):
	""" returns all guids and their examination (addition) time.
	reference, if passed, is ignored."""
	try:	
		result =fn3.get_all_guids_examination_time()	
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
	return make_response(tojson(result))


class test_get_all_guids_examination_time_1(unittest.TestCase):
    """ tests route /api/v2/guids_and_examination_times"""
    def runTest(self):
        relpath = "/api/v2/guids_and_examination_times"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        guidlist = json.loads(res.content.decode('utf-8'))
        
        self.assertTrue(isinstance(guidlist, dict))
        self.assertEqual(res.status_code, 200)

        #  test that it actually works
        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))		# get all the guids

        guid_to_insert = "guid_{0}".format(n_pre+1)

        inputfile = "../COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, 'rt') as f:
            for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                    seq = str(record.seq)

        print("Adding TB reference sequence of {0} bytes".format(len(seq)))
        self.assertEqual(len(seq), 4411532)		# check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload = {'guid':guid_to_insert,'seq':seq})
        self.assertEqual(res.status_code, 200)

        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(info, 'Guid {0} inserted.'.format(guid_to_insert))

        relpath = "/api/v2/guids_and_examination_times"
        res = do_GET(relpath)
        et= len(json.loads(res.content.decode('utf-8')))


@app.route('/api/v2/annotations', methods=['GET'])
def annotations(**kwargs):
	""" returns all guids and associated meta data.
	This query can be slow for very large data sets.
	"""
	try:
		result = fn3.get_all_annotations()
		
	except Exception as e:
		abort(500, e)
		
	return(tojson(result))

class test_annotations_1(unittest.TestCase):
    """ tests route /api/v2/annotations """
    def runTest(self):
        relpath = "/api/v2/annotations"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content = res.content))
        inputDict = json.loads(res.content.decode('utf-8'))
        self.assertTrue(isinstance(inputDict, dict)) 
        guiddf = pd.DataFrame.from_dict(inputDict,orient='index')		#, orient='index'
        self.assertTrue(isinstance(guiddf, pd.DataFrame)) 

@app.route('/api/v2/<string:guid>/exists', methods=['GET'])
def exist_sample(guid, **kwargs):
	""" checks whether a guid exists.
	reference and method are ignored."""
	
	try:
		result = fn3.exist_sample(guid)
		
	except Exception as e:
		abort(500, e)
		
	return make_response(tojson(result))

class test_exist_sample(unittest.TestCase):
    """ tests route /api/v2/guid/exists """
    def runTest(self):
        relpath = "/api/v2/non_existent_guid/exists"
        res = do_GET(relpath)
       
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), bool)
        self.assertEqual(info, False)


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
			result = fn3.insert(guid, seq)
		else:
			abort(501, 'seq and guid are not present in the POSTed data {0}'.format(data_keys))
		
	except Exception as e:
		print("Exception raised", e)
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
	res = sorted(fn3.clustering.keys())		
	return make_response(tojson({'algorithms':res}))

class test_algorithms(unittest.TestCase):
	"""  tests return of a change_id number """
	def runTest(self):
		relpath = "/api/v2/clustering"
		res = do_GET(relpath)
		self.assertEqual(res.status_code, 200)
		retDict = json.loads(str(res.text))
		self.assertEqual(retDict, {'algorithms': ['SNV12_ignore', 'SNV12_include']})


@app.route('/api/v2/clustering/<string:clustering_algorithm>/change_id', methods=['GET'])
def change_id(clustering_algorithm):
	"""  returns the current change_id number, which is incremented each time a change is made.
	     Useful for recovering changes in clustering after a particular point."""
	try:
		res = fn3.clustering[clustering_algorithm].change_id		
	except KeyError:
		# no clustering algorithm of this type
		abort(404, "no clustering algorithm {0}".format(clustering_algorithm))
		
	return make_response(tojson({'change_id': res, 'clustering_algorithm':clustering_algorithm}))

@app.route('/api/v2/clustering/<string:clustering_algorithm>/guids2clusters', methods=['GET'])
def g2c(clustering_algorithm):
	"""  returns a guid -> clusterid dictionary for all guids """
	try:
		res = fn3.clustering[clustering_algorithm].clusters2guidmeta(after_change_id = None)		
	except KeyError:
		# no clustering algorithm of this type
		abort(404, "no clustering algorithm {0}".format(clustering_algorithm))
		
	return make_response(tojson(res))

class test_g2c(unittest.TestCase):
	"""  tests return of a change_id number """
	def runTest(self):
		relpath = "/api/v2/clustering/SNV12_ignore/guids2clusters"
		res = do_GET(relpath)
		self.assertEqual(res.status_code, 200)
		retVal = json.loads(str(res.text))
		self.assertTrue(isinstance(retVal, list))

@app.route('/api/v2/clustering/<string:clustering_algorithm>/guids2clusters/after_change_id/<int:change_id>', methods=['GET'])
def g2ca(clustering_algorithm, change_id):
	"""  returns a guid -> clusterid dictionary, with changes occurring after change_id, a counter which is incremented each time a change is made.
	     Useful for recovering changes in clustering after a particular point."""
	try:
		res = fn3.clustering[clustering_algorithm].clusters2guidmeta(after_change_id = change_id)		
	except KeyError:
		# no clustering algorithm of this type
		abort(404, "no clustering algorithm {0}".format(clustering_algorithm))
		
	return make_response(tojson(res))

class test_g2ca(unittest.TestCase):
	"""  tests return of a change_id number """
	def runTest(self):
		relpath = "/api/v2/clustering/SNV12_ignore/guids2clusters/after_change_id/1"
		res = do_GET(relpath)
		self.assertEqual(res.status_code, 200)
		retVal = json.loads(str(res.text))
		self.assertTrue(isinstance(retVal, list))
		
class test_change_id(unittest.TestCase):
	"""  tests return of a change_id number """
	def runTest(self):
		relpath = "/api/v2/clustering/SNV12_ignore/change_id"
		res = do_GET(relpath)
		self.assertEqual(res.status_code, 200)
		retDict = json.loads(str(res.text))
		self.assertEqual(set(retDict.keys()), set(['change_id','clustering_algorithm']))
		self.assertEqual(retDict['clustering_algorithm'],'SNV12_ignore')

		relpath = "/api/v2/clustering/SNV12_ignore/change_id"
		res = do_GET(relpath)
		self.assertEqual(res.status_code, 200)
		retDict = json.loads(str(res.text))
		self.assertEqual(set(retDict.keys()), set(['change_id','clustering_algorithm']))
		self.assertEqual(retDict['clustering_algorithm'],'SNV12_ignore')

		relpath = "/api/v2/clustering/not_exists/change_id"
		res = do_GET(relpath)
		self.assertEqual(res.status_code, 404)
		

class test_insert_1(unittest.TestCase):
    """ tests route /api/v2/insert """
    def runTest(self):
        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))		# get all the guids

        guid_to_insert = "guid_{0}".format(n_pre+1)

        inputfile = "../COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, 'rt') as f:
            for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                    seq = str(record.seq)

        print("Adding TB reference sequence of {0} bytes".format(len(seq)))
        self.assertEqual(len(seq), 4411532)		# check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload = {'guid':guid_to_insert,'seq':seq})
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(info, 'Guid {0} inserted.'.format(guid_to_insert))

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_post = len(json.loads(res.content.decode('utf-8')))
        self.assertEqual(n_pre+1, n_post)
                

        # check if it exists
        relpath = "/api/v2/{0}/exists".format(guid_to_insert)
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), bool)
        self.assertEqual(res.status_code, 200)
        self.assertEqual(info, True)

class test_insert_10(unittest.TestCase):
	""" tests route /api/v2/insert, with additional samples.
		Also provides a set of very similar samples, testing recompression code."""
	def runTest(self):
		relpath = "/api/v2/guids"
		res = do_GET(relpath)
		n_pre = len(json.loads(str(res.text)))		# get all the guids

		inputfile = "../COMPASS_reference/R39/R00000039.fasta"
		with open(inputfile, 'rt') as f:
			for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
					originalseq = list(str(record.seq))
					
		for i in range(1,10):
			guid_to_insert = "guid_{0}".format(n_pre+i)

			seq = originalseq			
			# make i mutations at position 500,000
			offset = 500000
			for j in range(i):
				mutbase = offset+j
				ref = seq[mutbase]
				if not ref == 'T':
					seq[mutbase] = 'T'
				if not ref == 'A':
					seq[mutbase] = 'A'
			seq = ''.join(seq)
						
			print("Adding TB sequence {2} of {0} bytes with {1} mutations relative to ref.".format(len(seq), i, guid_to_insert))
			self.assertEqual(len(seq), 4411532)		# check it's the right sequence
	
			relpath = "/api/v2/insert"
			res = do_POST(relpath, payload = {'guid':guid_to_insert,'seq':seq})
			self.assertTrue(isjson(content = res.content))
			info = json.loads(res.content.decode('utf-8'))
			self.assertEqual(info, 'Guid {0} inserted.'.format(guid_to_insert))
	
			relpath = "/api/v2/guids"
			res = do_GET(relpath)
			n_post = len(json.loads(res.content.decode('utf-8')))
			self.assertEqual(n_pre+i, n_post)
					
			# check if it exists
			relpath = "/api/v2/{0}/exists".format(guid_to_insert)
			res = do_GET(relpath)
			self.assertTrue(isjson(content = res.content))
			info = json.loads(res.content.decode('utf-8'))
			self.assertEqual(type(info), bool)
			self.assertEqual(res.status_code, 200)
			self.assertEqual(info, True)	

class test_insert_60(unittest.TestCase):
	""" tests route /api/v2/insert, with additional samples.
		Also provides a set of very similar samples, testing recompression code."""
	def runTest(self):
		relpath = "/api/v2/guids"
		res = do_GET(relpath)
		n_pre = len(json.loads(str(res.text)))		# get all the guids

		inputfile = "../COMPASS_reference/R39/R00000039.fasta"
		with open(inputfile, 'rt') as f:
			for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
					originalseq = list(str(record.seq))
		guids_inserted = list()			
		for i in range(1,40):
			
			seq = originalseq
			if i % 5 ==0:
				is_mixed = True
				guid_to_insert = "mixed_{0}".format(n_pre+i)
			else:
				is_mixed = False
				guid_to_insert = "nomix_{0}".format(n_pre+i)	
			# make i mutations at position 500,000
			
			offset = 500000
			for j in range(i):
				mutbase = offset+j
				ref = seq[mutbase]
				if is_mixed == False:
					if not ref == 'T':
						seq[mutbase] = 'T'
					if not ref == 'A':
						seq[mutbase] = 'A'
				if is_mixed == True:
						seq[mutbase] = 'N'					
			seq = ''.join(seq)
			guids_inserted.append(guid_to_insert)			
			if is_mixed:
					print("Adding TB sequence {2} of {0} bytes with {1} mutations relative to ref.".format(len(seq), i, guid_to_insert))
			else:
					print("Adding mixed TB sequence {2} of {0} bytes with {1} Ns relative to ref.".format(len(seq), i, guid_to_insert))
				
				
			self.assertEqual(len(seq), 4411532)		# check it's the right sequence
	
			relpath = "/api/v2/insert"
			res = do_POST(relpath, payload = {'guid':guid_to_insert,'seq':seq})
			self.assertEqual(res.status_code, 200)
		
			self.assertTrue(isjson(content = res.content))
			info = json.loads(res.content.decode('utf-8'))
			self.assertEqual(info, 'Guid {0} inserted.'.format(guid_to_insert))
	
			relpath = "/api/v2/guids"
			res = do_GET(relpath)
			n_post = len(json.loads(res.content.decode('utf-8')))
			self.assertEqual(n_pre+i, n_post)
					
			# check if it exists
			relpath = "/api/v2/{0}/exists".format(guid_to_insert)
			res = do_GET(relpath)
			self.assertTrue(isjson(content = res.content))
			info = json.loads(res.content.decode('utf-8'))
			self.assertEqual(type(info), bool)
			self.assertEqual(res.status_code, 200)
			self.assertEqual(info, True)	

		# check: is everything there?
		for guid in guids_inserted:
			relpath = "/api/v2/{0}/exists".format(guid)
			res = do_GET(relpath)
			self.assertTrue(isjson(content = res.content))
			info = json.loads(res.content.decode('utf-8'))
			self.assertEqual(type(info), bool)
			self.assertEqual(res.status_code, 200)
			self.assertEqual(info, True)	

		# is everything clustered?
		relpath = "/api/v2/clustering/SNV12_ignore/guids2clusters"
		res = do_GET(relpath)
		self.assertEqual(res.status_code, 200)
		retVal = json.loads(str(res.text))
		self.assertTrue(isinstance(retVal, list))
		
		print("running mixed checks:")
		for item in retVal:
			if 'mixed_' in item['guid']:
				self.assertTrue(item['is_mixed'])
				print(item['guid'], item['is_mixed'])
		
class test_mirror(unittest.TestCase):
    """ tests route /api/v2/mirror """
    def runTest(self):
        
        relpath = "/api/v2/mirror"
        payload = {'guid':'1', 'seq':"ACTG"}
        res = do_POST(relpath, payload = payload)
        res_dict = json.loads(res.content.decode('utf-8'))
        self.assertEqual(payload, res_dict)


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
	if not returned_format in set([1,2,]):
		abort(500, "Invalid format requested, must be 1 or 2")
	if not ( 0 <= cutoff  and cutoff <= 1):
		abort(500, "Invalid cutoff requested, must be between 0 and 1")
		
	try:
		result = fn3.neighbours_within_filter(guid, threshold, cutoff, returned_format)
	except KeyError as e:
		# guid doesn't exist
		abort(404, e)
	except Exception as e:
		abort(500, e)
	
	return make_response(tojson(result))
	
class test_neighbours_within_1(unittest.TestCase):
    """ tests route /api/v2/guid/neighbours_within/ """
    def runTest(self):
        relpath = "/api/v2/non_existent_guid/neighbours_within/12"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), dict)
        self.assertEqual(res.status_code, 404)

class test_neighbours_within_2(unittest.TestCase):
    """ tests route /api/v2/guid/neighbours_within/ """
    def runTest(self):
        relpath = "/api/v2/non_existent_guid/neighbours_within/12/with_quality_cutoff/0.5"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), dict)
        self.assertEqual(res.status_code, 404)

class test_neighbours_within_3(unittest.TestCase):
    """ tests route /api/v2/guid/neighbours_within/ """
    def runTest(self):
        relpath = "/api/v2/non_existent_guid/neighbours_within/12/with_quality_cutoff/0.5/in_format/1"
        res = do_GET(relpath)
        print(res)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), dict)
        self.assertEqual(res.status_code, 404)

class test_neighbours_within_4(unittest.TestCase):
    """ tests route /api/v2/guid/neighbours_within/ """
    def runTest(self):
        relpath = "/api/v2/non_existent_guid/neighbours_within/12/with_quality_cutoff/0.5/in_format/2"
        res = do_GET(relpath)
        print(res)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), dict)
        self.assertEqual(res.status_code, 404)

class test_neighbours_within_5(unittest.TestCase):
    """ tests route /api/v2/guid/neighbours_within/ """
    def runTest(self):
        relpath = "/api/v2/non_existent_guid/neighbours_within/12/in_format/2"
        res = do_GET(relpath)
        print(res)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), dict)
        self.assertEqual(res.status_code, 404)
 
class test_neighbours_within_6(unittest.TestCase):
   """ tests all the /api/v2/guid/neighbours_within methods using test data """
   def runTest(self):
        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))

        inputfile = "../COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, 'rt') as f:
            for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                    seq = str(record.seq)
                    
        # generate variants
        variants = {}
        for i in range(4):
                 guid_to_insert = "guid_insert_{0}".format(n_pre+i+1)
                 vseq=list(seq)
                 vseq[100*i]='A'
                 vseq=''.join(vseq)
                 variants[guid_to_insert] = vseq

        for guid_to_insert in variants.keys():

                print("Adding mutated TB reference sequence called {0}".format(guid_to_insert))        
                relpath = "/api/v2/insert"
                
                res = do_POST(relpath, payload = {'guid':guid_to_insert,'seq':variants[guid_to_insert]})
                self.assertTrue(isjson(content = res.content))
                info = json.loads(res.content.decode('utf-8'))
                self.assertTrue('inserted' in info)

                # check if it exists
                relpath = "/api/v2/{0}/exists".format(guid_to_insert)
                res = do_GET(relpath)
                self.assertTrue(isjson(content = res.content))
                info = json.loads(res.content.decode('utf-8'))
                self.assertEqual(type(info), bool)
                self.assertEqual(res.status_code, 200)
                self.assertEqual(info, True)
        
        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_post = len(json.loads(res.content.decode('utf-8')))
        self.assertEqual(n_pre+4, n_post)

        test_guid = min(variants.keys())
        print("Searching for ",test_guid)
        
        search_paths = ["/api/v2/{0}/neighbours_within/1",
                        "/api/v2/{0}/neighbours_within/1/with_quality_cutoff/0.5",
                        "/api/v2/{0}/neighbours_within/1/with_quality_cutoff/0.5/in_format/1"
                        "/api/v2/{0}/neighbours_within/1/with_quality_cutoff/0.5/in_format/2"
                        "/api/v2/{0}/neighbours_within/1/in_format/1"
                        "/api/v2/{0}/neighbours_within/1/in_format/2"
                        ]
        
        for search_path in search_paths:
                url = search_path.format(test_guid)
                res = do_GET(relpath)
                self.assertTrue(isjson(content = res.content))
                info = json.loads(res.content.decode('utf-8'))
                self.assertEqual(type(info), list)
                guids_found = set()
                for item in info:
                        guids_found.add(item)
                recovered = guids_found.intersection(variants.keys())
                self.assertEqual(len(recovered),4)
                self.assertEqual(res.status_code, 200)

@app.route('/api/v2/<string:guid>/sequence', methods=['GET'])
def sequence(guid):
	""" returns the masked sequence as a string """	
	result = fn3.sequence(guid)
	if result is None:  # no guid exists
		return make_response(tojson('guid {0} does not exist'.format(guid)), 404)
	else:
		return make_response(tojson(result))

class test_sequence_1(unittest.TestCase):
    """ tests route /api/v2/*guid*/sequence"""
    def runTest(self):
        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))		# get all the guids

        guid_to_insert = "guid_{0}".format(n_pre+1)

        inputfile = "../COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, 'rt') as f:
            for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                    seq = str(record.seq)

        print("Adding TB reference sequence of {0} bytes".format(len(seq)))
        self.assertEqual(len(seq), 4411532)		# check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload = {'guid':guid_to_insert,'seq':seq})
        self.assertEqual(res.status_code, 200)

        relpath = "/api/v2/{0}/sequence".format(guid_to_insert)
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)

        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(info['guid'], guid_to_insert)
        self.assertEqual(info['invalid'], 0)
        self.assertEqual(info['masked_dna'].count('N'), 557291)

class test_sequence_2(unittest.TestCase):
    """ tests route /api/v2/*guid*/sequence"""
    def runTest(self):
 
        relpath = "/api/v2/{0}/sequence".format('no_guid_exists')
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 404)

class test_sequence_3(unittest.TestCase):
    """ tests route /api/v2/*guid*/sequence"""
    def runTest(self):
        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        print(res)
        n_pre = len(json.loads(res.content.decode('utf-8')))		# get all the guids

        guid_to_insert = "guid_{0}".format(n_pre+1)

        inputfile = "../COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, 'rt') as f:
            for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                    seq = str(record.seq)
        seq = 'N'*4411532
        print("Adding TB reference sequence of {0} bytes with {1} Ns".format(len(seq), seq.count('N')))
        self.assertEqual(len(seq), 4411532)		# check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload = {'guid':guid_to_insert,'seq':seq})
        self.assertEqual(res.status_code, 200)

        relpath = "/api/v2/{0}/sequence".format(guid_to_insert)
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)

        info = json.loads(res.content.decode('utf-8'))
        print(info)
        self.assertEqual(info['guid'], guid_to_insert)
        self.assertEqual(info['invalid'], 1)

class test_sequence_4(unittest.TestCase):
    """ tests route /api/v2/*guid*/sequence"""
    def runTest(self):
        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        print(res)
        n_pre = len(json.loads(res.content.decode('utf-8')))		# get all the guids

        inputfile = "../COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, 'rt') as f:
            for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                    seq2 = str(record.seq)

        guid_to_insert1 = "guid_{0}".format(n_pre+1)
        guid_to_insert2 = "guid_{0}".format(n_pre+2)


        seq1 = 'N'*4411532
        print("Adding TB reference sequence of {0} bytes with {1} Ns".format(len(seq1), seq1.count('N')))
        self.assertEqual(len(seq1), 4411532)		# check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload = {'guid':guid_to_insert1,'seq':seq1})
        self.assertEqual(res.status_code, 200)

        print("Adding TB reference sequence of {0} bytes with {1} Ns".format(len(seq2), seq2.count('N')))
        self.assertEqual(len(seq2), 4411532)		# check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload = {'guid':guid_to_insert2,'seq':seq2})
        self.assertEqual(res.status_code, 200)
		
        relpath = "/api/v2/{0}/sequence".format(guid_to_insert1)
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)

        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(info['guid'], guid_to_insert1)
        self.assertEqual(info['invalid'], 1)
        relpath = "/api/v2/{0}/sequence".format(guid_to_insert2)
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)

        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(info['guid'], guid_to_insert2)
        self.assertEqual(info['invalid'], 0)

@app.route('/api/v2/nucleotides_excluded', methods=['GET'])
def nucleotides_excluded():
	""" returns all nucleotides excluded by the server.
	Useful for clients which need to to ensure that server
	and client masking are identical. """
	
	try:
		result = fn3.server_nucleotides_excluded()
		
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)

	return make_response(tojson(result))

class test_nucleotides_excluded(unittest.TestCase):
    """ tests route /api/v2/nucleotides_excluded"""
    def runTest(self):
        relpath = "api/v2/nucleotides_excluded"
        res = do_GET(relpath)
        resDict = json.loads(res.text)
        self.assertTrue(isinstance(resDict, dict))
        self.assertEqual(set(resDict.keys()), set(['exclusion_id', 'excluded_nt']))
        self.assertEqual(res.status_code, 200)
 

# startup
if __name__ == '__main__':

        # command line usage.  Pass the location of a config file as a single argument.
        # an example config file is default_test_config.json
               
        ############################ LOAD CONFIG ######################################
        print("findNeighbour3 server .. reading configuration file.")

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
        app.logger.setLevel(loglevel)       
        file_handler = logging.FileHandler(CONFIG['LOGFILE'])
        formatter = logging.Formatter( "%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s ")
        file_handler.setFormatter(formatter)
        app.logger.addHandler(file_handler)
 
        # log a test error on startup
        # app.logger.error("Test error logged on startup, to check logger is working")

        # launch sentry if API key provided
        if 'SENTRY_URL' in CONFIG.keys():
                app.logger.info("Launching logger")
                sentry_sdk.init(CONFIG['SENTRY_URL'], integrations=[FlaskIntegration()])

        ########################### prepare to launch server ###############################################################
        # construct the required global variables
        LISTEN_TO = '127.0.0.1'
        RESTBASEURL = "http://{0}:{1}".format(CONFIG['IP'], CONFIG['REST_PORT'])

        #########################  CONFIGURE HELPER APPLICATIONS ######################
        ## once the flask app is running, errors get logged to app.logger.  However, problems on start up do not.
        ## configure mongodb persistence store
        print("Connecting to backend data store")
        try:
                PERSIST=fn3persistence(dbname = CONFIG['SERVERNAME'],connString=CONFIG['FNPERSISTENCE_CONNSTRING'], debug=CONFIG['DEBUGMODE'])
        except Exception as e:
                app.logger.exception("Error raised on creating persistence object")
                if e.__module__ == "pymongo.errors":
                      app.logger.info("Error raised pertains to pyMongo connectivity")
                raise

        # instantiate server class
        print("Loading sequences into server, please wait ...")
        try:
        	fn3 = findNeighbour3(CONFIG, PERSIST)
        except Exception as e:
                app.logger.exception("Error raised on instantiating findNeighbour3 object")
                raise


        ########################  START THE SERVER ###################################
        if CONFIG['DEBUGMODE']>0:
                flask_debug = True
                app.config['PROPAGATE_EXCEPTIONS'] = True
        else:
                flask_debug = False

        app.logger.info("Launching server on {0}, debug = {1}, port = {2}".format(LISTEN_TO, flask_debug, CONFIG['REST_PORT']))
        #try:
        app.run(host=LISTEN_TO, debug=flask_debug, port = CONFIG['REST_PORT'])

        #except Exception as e:
        #       app.logger.exception("Failed to launch flask web server")
        #        raise



