#!/usr/bin/env python

### SERVER
from xmlrpc.server import SimpleXMLRPCServer
import uuid
import datetime
import unittest
import os
import glob
import sys
import datetime
import pickle
import hashlib
import collections
import uuid
import subprocess
import threading
import time
import json
import queue
import warnings
import logging
import psutil

# import resource

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Table, Column, Integer, String, Float, DateTime, MetaData, select, func, or_
from sqlalchemy import create_engine, exc as excb
from sqlalchemy.orm import sessionmaker
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship, backref, exc

from FN import *
from ewSetCore import ewSetCore
from seqComparer import seqComparer
from ewsnpstore import db_ewss, ElephantWalkDataSNPaccessor


class MyXMLRPCServer(SimpleXMLRPCServer):
	"""...""" 
	stopped = False

	def serve_forever(self):
			while not self.stopped:
					self.handle_request()

	@staticmethod
	def force_stop():
			MyXMLRPCServer.stopped = True
			urllib.urlopen("http://{0}/".format(args.URL)).read()

			
class ElephantWalk():
	""" a server based application for maintaining a record of bacterial relatedness using SNP distances """
	PERSIST=None
	def __init__(self,CONFIG):
		""" Using values in CONFIG, starts a server with CONFIG['NAME'] on port CONFIG['PORT'].
		
		CONFIG contains Configuration parameters relevant to the reference based compression system which lies
		at the core of the server.
		
		for explanations as to the meanings of these values, please see the documentation in  ewSetCore, to
		which the CONFIG dictionary gets passed.
		
		An example CONFIG is below:
		
		{			
		"DESCRIPTION":"A test server operating in ../unittest_tmp, only suitable for testing",
		"PORT":8184,
		"IP":"127.0.0.1",
		"INPUTREF":"../reference/TB-ref.fasta",
		"PERSISTENCEDIR":"../unittest_tmp/",
		"EXCLUDEFILE":"../reference/TB-exclude.txt",
		"SNPDIR":"../unittest_tmp",
		"DEBUGMODE":0,
		"SERVERNAME":"TBSNP",
		"EDGEDB_CONNSTRING":"sqlite:///<<DEFAULT>>/{0}.db",
		"FNPERSISTENCE_CONNSTRING":"sqlite:///<<DEFAULT>>/findNeighbour.db",
		"MAXN_STORAGE":100000,
		"MAXN_PROP_DEFAULT":0.70,
		"NCOMPRESSIONCUTOFF":100000,
		"LOGFILE":"../unittest_tmp/logfile.log",
		"LOGLEVEL":"INFO",
		"SNPCEILING": 20,
		"MULTIPROCESSING":0
		}

	
		"""
		
		## check input
		if isinstance(CONFIG, str):
			CONFIG=json.loads(CONFIG)	# assume JSON string; convert.
		
		# check it is a dictionary	
		if not isinstance(CONFIG, dict):
			raise KeyError("CONFIG must be either a dictionary or a JSON string encoding a dictionary.  It is: {0}".format(CONFIG))
		
		# check that the keys of config are as expected.
		required_keys=set(['IP','PORT','INPUTREF','PERSISTENCEDIR','EXCLUDEFILE','SNPDIR','DEBUGMODE','SERVERNAME',
						   'EDGEDB_CONNSTRING','FNPERSISTENCE_CONNSTRING', 'MAXN_STORAGE',
						   'NCOMPRESSIONCUTOFF', 'SNPCOMPRESSIONCEILING', 'MAXN_PROP_DEFAULT', 'MULTIPROCESSING', 'REST_PORT',
						   'LOGFILE','LOGLEVEL'])
		missing=required_keys-set(CONFIG.keys())
		if not missing == set([]):
			raise KeyError("Required keys were not found in CONFIG. Missing are {0}".format(missing))
		
		# start process
		self.write_semaphore = threading.BoundedSemaphore(1)
		self.CONFIG=CONFIG
		self.objExaminer=NucleicAcid()
		self.ewsc=ewSetCore(CONFIG=self.CONFIG)
		print ("EW2 is Ready; ({0})".format(self.ewsc.sc.excluded_hash()))
		
	def insert(self,sname,dna):
		""" insert DNA called sname into the server
		
		TODO:
		At present, inadequate consideration has been given to what happens if
		(as happens, but very rarely) the database layer fails.
		
		Therefore, calls to the data base layer should be wrapped in try/catch and
		suitable rollback / error raising arranged.
		"""
		
		# clean, and provide summary statistics for the sequence
		logging.info("Inserting: {0}".format(sname))
		
		already_exists = self.ewsc.exist_sample(sname)
		if not already_exists:
			self.objExaminer.examine(dna)  
			
			# store the summary statistics
			ElephantWalk.PERSIST.annotateFromDict(sequenceGuid=sname, nameSpace='DNAQuality',annotDict=self.objExaminer.composition)
			
			# write 
			cleaned_dna=self.objExaminer.nucleicAcidString.decode()
			self.write_semaphore.acquire()				    # addition should be an atomic operation
			self.ewsc.insert(sname, cleaned_dna)			# insert the DNA sequence into the server.
			self.write_semaphore.release()                  # release the write semaphore
			return json.dumps(["OK"])									# a response object with a 200 code would be better
		else:
			return json.dumps(["Already present"])
	
	def exist_sample(self,sname):
		""" determine whether the sample exists """
		return self.ewsc.exist_sample(sname)

	def server_time(self):
		""" returns the current server time """
		return json.dumps({"server_time":datetime.datetime.now().isoformat()})

	def server_config(self):
		""" returns the config file with which the server was launched
		
		This may be highly undesirable,
		as it reveals the internal server architecture  including
		backend databases and perhaps connection strings with passwords.
		"""
		return json.dumps(self.CONFIG)

	def server_nucleotides_excluded(self):
		""" returns the nucleotides excluded by the server """
		l = sorted(list(self.ewsc.sc.excluded))
		return json.dumps({"exclusion_id":self.ewsc.sc.excluded_hash(), "excluded_nt":l})
	
	def server_memory_usage(self):
		""" returns memory usage by current python process
		
		Uses the resource module.
		Please see:
		https://docs.python.org/3/library/resource.html
		for more information.  """
		
		mem= {'maximum_resident_set_size':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
			  'note':'Values are as returned by the python3 resource module for the server process only.  Please see https://docs.python.org/3/library/resource.html for more information',
			  'memory_units':'bytes'}
		return json.dumps(mem)
	
	def query_get_value_snp_filter(self, sname, snpDistance, cutoff=0.85, returned_format=1):
		""" returns a list of guids, and their distances, by a sample quality cutoff
		
		    returns links either as
			format 1 [otherGuid, distance]
            or as
			format 2 [otherGuid, distance, N_just1, N_just2, N_either]
        """

		# check the query is of good quality
		inScore = ElephantWalk.PERSIST.testIndividualGuidQuality(sname,float(cutoff))
		if inScore == None:
			return json.dumps(['Err','missing sample'])
		elif inScore == False:
			return json.dumps(['Bad','bad sequence'])

		# if it is of good quality, then we look for links
		idList=list()


		# gets the similar sequences from the database;
		# wrapper round ewc.neighboursOf()
		# returned value reVal is {'guid':guid, 'neighbours':[ *, *, * ...] }
		# where * is one of the two formats shown in the docstring above.
		
		retVal = self.ewsc.query_get_values(guid=sname, cutoff=snpDistance, returned_format=returned_format)
		
		# run a quality check on the things our sample is like.
		sampleList=retVal['neighbours']
		idList=[]
		for sa in sampleList:
			idList.append(sa[0])		# add the guid
		
		# get the sample qualities from the database
		guid2qual=ElephantWalk.PERSIST.guid2quality(idList)
					  
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
				
		return json.dumps(finalOutput)
	

	def get_all_guids(self):
		return ElephantWalk.PERSIST.asJson_sequenceGuids()
	
	def get_all_filtered_guids(self,cutoff=0.66):
		rs=ElephantWalk.PERSIST.asJson_propACTG_filteredSequenceGuids(float(cutoff))
		if rs==None:
			return json.dumps([])
		else:
			return rs
		
	def get_all_guids_examination_time(self):
		return ElephantWalk.PERSIST.asJson_guid2ExaminationDateTime()
	
	def get_all_annotations(self):
		return ElephantWalk.PERSIST.asJson_allAnnotations()
	
	def query_get_detail(self, sname1, sname2):
		""" gets detail on the comparison of a pair of samples.  Computes this on the fly """
		ret = self.ewsc.query_get_detail(sname1,sname2)
		return(json.dumps(ret))
	
	def get_all_values(self, cutoff=12):
		ret = self.ewsc.query_get_values_snp(cutoff=cutoff)
		return(json.dumps(ret))
	
if __name__=='__main__':
	
	# command line usage.  Pass the location of a config file as a single argument.
	# an example config file is default_config.json

	############################ LOAD CONFIG ######################################
	logging.basicConfig(level=logging.INFO)
	if len(sys.argv) == 2:			 
		try:
			with open(sys.argv[1],'r') as f:
				CONFIG_STRING=f.read()
				

		except FileNotFoundError:
			raise FileNotFoundError("Passed one parameter, which should be a CONFIG file name; tried to open a config file at {0} but it does not exist ".format(sys.argv[1]))
	else:
		# use default which may be inappropriate in production
		warnings.warn("No config file name supplied as a single argument; using '../config/default_config.json' This configuration suitable only for testing, not for production.")
		try:
			with open(os.path.join('..','config','default_test_config.json'),'r') as f:
				CONFIG_STRING=f.read()
				
		except FileNotFoundError:
			raise FileNotFoundError("No default_test_config.json file found.")

	CONFIG=json.loads(CONFIG_STRING)

	########################### SET UP LOGGING #####################################
	# defaults to INFO.  WARN and DEBUG also supported.
	loglevel=logging.INFO

	if 'LOGLEVEL' in CONFIG.keys():
		if CONFIG['LOGLEVEL']=='WARN':
			loglevel=logging.WARN
		elif CONFIG['LOGLEVEL']=='DEBUG':
			loglevel=logging.DEBUG

	if 'LOGFILE' in CONFIG.keys():
			logfile=os.path.abspath(CONFIG['LOGFILE'])
			print("Logging to {0}".format(logfile))
			logging.basicConfig(filename=logfile, format='%(asctime)s|%(levelname)s|%(message)s', level=loglevel)
	else:
		warnings.warn("No LOGFILE entry in CONFIG, so no logging to file in place.")
		logging.basicConfig(format='%(asctime)s|%(levelname)s|%(message)s', level=loglevel)
		
	
	#########################  CONFIGURE HELPER APPLICATIONS ######################		

	ElephantWalk.PERSIST=FNPersistence(db=db, engineName=CONFIG['FNPERSISTENCE_CONNSTRING'], rootdir=CONFIG['PERSISTENCEDIR'])

	
	########################  START THE SERVER ###################################
	logging.info("Attempting to launch server on {0} and port {1}".format(CONFIG['IP'],CONFIG['PORT']))
	server = MyXMLRPCServer((CONFIG['IP'],CONFIG['PORT']))
	server.register_introspection_functions()
	server.register_instance(ElephantWalk(CONFIG))
	server.serve_forever()
