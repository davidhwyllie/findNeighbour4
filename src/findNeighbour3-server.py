 #!/usr/bin/env python
""" 
provides restful interface to ElephantWalk2 functions
 
The endpoint provided is designed as an internal API.
As implemented, it is not protected by authentication.

requires python3.

loads configuration from a config file.
If no config file is provided, it will run in 'testing' mode with the
parameters in default_test_config.json.

The config file must include three parameters:
REST_PORT - the port on which this server should run.

Unit testing can be achieved by

# starting a test RESTFUL server
python3 findNeighbour3-server-rest.py

# And then (e.g. in a different terminal) launching unit tests with
python3 -m unittest findNeighbour3-server

# all should pass
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

# flask
from flask import Flask, make_response, jsonify, Markup
from flask import request, abort

# logging
from logging.config import dictConfig

# utilities for file handling and measuring file size
import psutil

# measure server memory usage; linux specifc
if not os.name == 'nt':
	import resource

# reference based compression modules
from NucleicAcid import NucleicAcid
from mongoStore import fn3persistence
from seqComparer import seqComparer

# pandas/numpy
import pandas as pd
import numpy as np

# only used for unit testing
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide
import unittest
from urllib.parse import urlparse as urlparser
from urllib.parse import urljoin as urljoiner


class ElephantWalk():
	""" a server based application for maintaining a record of bacterial relatedness using SNP distances.
	
	    The high level arrangement is that
		- This class interacts with in-memory sequences
		  [handled by the seqComparer class] and backends [fn3Persistance class] used by the server
		- methods in ElephantWalk() return native python objects.
		
		- a web server, currently flask, handles the inputs and outputs of this class
		- in particular, native python objects returned by this class are serialised by the Flask web server code.
		"""
		
	def __init__(self,CONFIG, PERSIST):
		""" Using values in CONFIG, starts a server with CONFIG['NAME'] on port CONFIG['PORT'].
		
		CONFIG contains Configuration parameters relevant to the reference based compression system which lies
		at the core of the server.
		
		for7 explanations as to the meanings of these values, please see the documentation in  ewSetCore, to
		which the CONFIG dictionary gets passed.
		
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
		"NCOMPRESSIONCUTOFF":100000,
		"LOGFILE":"../unittest_tmp/logfile.log",
		"LOGLEVEL":"INFO",
		"SNPCEILING": 20,
		"MULTIPROCESSING":0
		}

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
						   'NCOMPRESSIONCUTOFF', 'SNPCOMPRESSIONCEILING', "SNPCEILING", 'MAXN_PROP_DEFAULT', 'MULTIPROCESSING', 'REST_PORT',
						   'LOGFILE','LOGLEVEL'])
		missing=required_keys-set(self.CONFIG.keys())
		if not missing == set([]):
			raise KeyError("Required keys were not found in CONFIG. Missing are {0}".format(missing))

		# the following keys are not stored in any database backend, as a server could be moved, i.e.
		# running on the same data but with different IP etc
		do_not_persist_keys=set(['IP','SERVERNAME','FNPERSISTENCE_CONNSTRING','LOGFILE','LOGLEVEL','REST_PORT'])
				
		# determine whether this is a first-run situation.
		if self.PERSIST.first_run:
			self.first_run(do_not_persist_keys)

		# load global settings from those stored at the first run.  
		cfg = self.PERSIST.config_read('config')
		self.reference = cfg['reference']
		self.excludePositions = set(cfg['excludePositions'])
		self.debugMode = cfg['DEBUGMODE']
		self.maxNs = cfg['MAXN_STORAGE']
		self.snpceiling = cfg['SNPCEILING']
		self.ncompressioncutoff = cfg['NCOMPRESSIONCUTOFF']
		self.snpcompressionCeiling = cfg['SNPCOMPRESSIONCEILING']
		self.maxn_prop_default = cfg['MAXN_PROP_DEFAULT']
		
		# start process
		self.write_semaphore = threading.BoundedSemaphore(1)        # used to permit only one process to INSERT at a time.
		self.objExaminer=NucleicAcid()
		self.sc=seqComparer(reference=self.reference,
							maxNs=self.maxNs,
							snpCeiling= self.snpceiling,
							debugMode=self.debugMode,
							excludePositions=self.excludePositions,
							snpCompressionCeiling = self.snpcompressionCeiling)
		
		# now load compressed sequences into ram.
		# note this does not apply 'double delta' recompression
		# a more sophisticated algorithm will be needed to do that
		guids = self.PERSIST.refcompressedsequence_guids()
		print("EW3 core is loading {1} sequences from database .. excluding ({0})".format(self.sc.excluded_hash(),len(guids)))
		nLoaded = 0
		for guid in guids:
			nLoaded+=1
			obj = self.PERSIST.refcompressedsequence_load(guid)
			self.sc.persist(obj)
			if nLoaded % 500 ==0:
				print(nLoaded)
		print("EW3 core is ready; loaded {0} sequences from database".format(len(guids)))
	
	def first_run(self, do_not_persist_keys):
		""" actions taken on first-run only.
		Include caching results from CONFIGFILE to database, unless they are in do_not_persist_keys"""
		
		logging.info("First run situation: parsing inputs, storing to database. ")

		# create a config dictionary
		config_settings= {}
		
		# store start time 
		config_settings['createTime']= datetime.datetime.now()
		
		# store description
		config_settings['description']=self.CONFIG['DESCRIPTION']
				
		# load the excluded bases
		excluded=set()
		if self.CONFIG['EXCLUDEFILE'] is not None:
			with open(self.CONFIG['EXCLUDEFILE'],'rt') as f:
				rows=f.readlines()
			for row in rows:
				excluded.add(int(row))

		logging.info("Noted {0} positions to exclude.".format(len(excluded)))
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
		
	def repack(self,guids=None):
		""" generates a smaller and faster representation in the persistence store
		for the guids in the list. optional"""
		if guids is None:
			guids = self.PERSIST.guids()  # all the guids
		for this_guid in guids:
			logging.info("Repacking {0}".format(this_guid))
			self.PERSIST.guid2neighbour_repack(this_guid)
			
	def insert(self,guid,dna):
		""" insert DNA called guid into the server
		
		TODO:
		At present, inadequate consideration has been given to what happens if
		(as happens, but very rarely) the database layer fails.
		
		Therefore, calls to the data base layer should be wrapped in try/catch and
		suitable rollback / error raising arranged.
		
		# reconsider this now Mongo is used.
		
		"""
		
		# clean, and provide summary statistics for the sequence
		logging.info("Inserting: {0}".format(guid))
		if not self.sc.iscachedinram(guid):                   # if the guid is not already there
			
			# prepare to insert
			self.objExaminer.examine(dna)  					  # examine the sequence
			cleaned_dna=self.objExaminer.nucleicAcidString.decode()
			refcompressedsequence =self.sc.compress(cleaned_dna)          # compress it and store it in RAM
			self.sc.persist(refcompressedsequence, guid)			    # insert the DNA sequence into ram.
					
			# construct links with everything existing existing at the time the semaphore was acquired.
			self.write_semaphore.acquire()				    # addition should be an atomic operation
			
			links={}
			for key2 in self.sc.guidscachedinram():
				if not guid==key2:
					(guid1,guid2,dist,n1,n2,nboth, N1pos, N2pos, Nbothpos)=self.sc.countDifferences_byKey(keyPair=(guid,key2))
					link = {'dist':dist,'n1':n1,'n2':n2,'nboth':nboth}
					if dist is not None:
						links[guid2]=link

			# write
			## should trap here to return sensible error message if database connectivity is lost.
			self.PERSIST.refcompressedseq_store(guid, refcompressedsequence)     # store the parsed object on disc
			self.PERSIST.guid_annotate(guid=guid, nameSpace='DNAQuality',annotDict=self.objExaminer.composition)						
			self.PERSIST.guid2neighbour_add_links(guid=guid, targetguids=links)

			# release semaphore
			self.write_semaphore.release()                  # release the write semaphore
			
			# clean up guid2neighbour; this can readily be done post-hoc, if the process is slow.  it doesn't affect results.
			guids = list(links.keys())
			guids.append(guid)
			self.repack(guids)
			return "Guid {0} inserted.".format(guid)		# a 200 will be added by flask
		else:
			return "Guid {0} is already present".format(guid)
	
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
		
		if self.debugMode==True:
			return self.CONFIG
		else:
			return None
	def server_nucleotides_excluded(self):
		""" returns the nucleotides excluded by the server """
		return {"exclusion_id":self.sc.excluded_hash(), "excluded_nt":list(self.sc.excluded)}
	
	def server_memory_usage(self):
		""" returns memory usage by current python process
		
		Uses the resource module.
		Please see:
		https://docs.python.org/3/library/resource.html
		for more information.
		
		These calls are linux specific"""
		if os.name == 'nt':
			mem= {'maximum_resident_set_size':-1,
			  'note':'Server is running windows; memory assessment is not available.',
			  'memory_units':'bytes'}			
		else:
			mem= {'maximum_resident_set_size':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
			  'note':'Values are as returned by the python3 resource module for the server process only.  Please see https://docs.python.org/3/library/resource.html for more information',
			  'memory_units':'bytes'}
		return mem
	
	def query_get_value_snp_filter(self, guid, snpDistance, cutoff=0.85, returned_format=1):
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
	
	def get_all_filtered_guids(self,cutoff=0.66):
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


# default parameters for unit testing only.
RESTBASEURL   = "http://127.0.0.1:5000"
ISDEBUG = True
LISTEN_TO = '127.0.0.1'		# only local addresses

# initialise Flask 
app = Flask(__name__)
app.logger.setLevel(logging.DEBUG)

			

def isjson(content):
        """ returns true if content parses as json, otherwise false. used by unit testing. """
        try:
            x = json.loads(content.decode('utf-8'))
            print("JSON DECODE SUCCEEDED YIELDING {1}: {0}".format(x,type(x)))
            return True
 
        except json.decoder.JSONDecodeError:
            print("JSON DECODE FAILED : {0}".format(content.decode('utf-8')))
            return False

def tojson(content):
	""" json dumps, formatting dates as isoformat """
	def converter(o):
		if isinstance(o, datetime.datetime):
			return o.isoformat()
		elif isinstance(o, pd.DataFrame):
			return o.to_json()  #(orient='index', date_format='iso')
		else:
			return o.__str__()
	return(json.dumps(content, default=converter))

# --------------------------------------------------------------------------------------------------
@app.errorhandler(404)
def not_found(error):
    json_err = jsonify({'error': 'Not found (custom error handler for mis-routing)'})
    return make_response(json_err, 404)
# --------------------------------------------------------------------------------------------------
 
@app.teardown_appcontext
def shutdown_session(exception=None):
    ew.PERSIST.closedown()		# close database connection
	
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
	session = requests.Session()
	session.trust_env = False

	# print out diagnostics
	print("POSTING to url {0}".format(url))
	response = session.post(url=url, data=payload, timeout=None)

	print("Result:")
	print("code: {0}".format(response.status_code))
	print("reason: {0}".format(response.reason))
	
	session.close()
	return(response)

@app.route('/', methods=['GET'])
def server_info():
	""" returns server info page
	"""
	res = """# findNeighbour3 web server operating.<p>Endpoints are <p>"""
	for route in ['not available']:
		res = res+"{0}<p>".format(route)
	return make_response(res)

@app.route('/api/v2/server_config', methods=['GET'])
def server_config():
    """ returns server configuration.

        returns the config file with which the server was launched.
        This may be highly undesirable,
        as it reveals the internal server architecture  including
        backend databases and perhaps connection strings with passwords.

    """
    res = ew.server_config()
    if res is None:		# not allowed to see it
        abort(404, "Endpoint only available in debug mode")
    else:
        return make_response(tojson(CONFIG))

class test_server_config(unittest.TestCase):
    """ tests route v2/server_config"""
    def runTest(self):
        relpath = "/api/v2/server_config"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))

        config_dict = json.loads(res.content.decode('utf-8'))
        self.assertTrue('NCOMPRESSIONCUTOFF' in config_dict.keys())
        self.assertEqual(res.status_code, 200)


@app.route('/api/v2/server_memory_usage', methods=['GET'])
def server_memory_usage():
	""" returns server memory usage """
	try:
		result = ew.server_memory_usage()

	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
	return make_response(tojson(result))

class test_server_memory_usage(unittest.TestCase):
    """ tests route /api/v2/server_memory_usage"""
    def runTest(self):
        relpath = "/api/v2/server_memory_usage"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))

        config_dict = json.loads(res.content.decode('utf-8'))
        self.assertTrue('note' in config_dict.keys())
        self.assertEqual(res.status_code, 200)


@app.route('/api/v2/server_time', methods=['GET'])
def server_time():
	""" returns server time """
	try:
		result = ew.server_time()
		
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
def get_all_guids(**kwargs):
	""" returns all guids.  reference, if included, is ignored."""
	try:
		result = list(ew.get_all_guids())
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
def get_all_filtered_guids(cutoff, **kwargs):
	""" returns all guids with quality score >= cutoff."""
	try:
		result = ew.get_all_filtered_guids(cutoff)	
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
	return make_response(tojson(result))

class test_get_all_filtered_guids_1(unittest.TestCase):
    """ tests route /api/v2/guids_with_quality_over"""
    def runTest(self):
        relpath = "/api/v2/guids_with_quality_over/0.7"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        guidlist = json.loads(res.content.decode('utf-8'))
        self.assertTrue(isinstance(guidlist, list))
        self.assertEqual(res.status_code, 200)
        # TODO: insert guids, check it doesn't fail.

@app.route('/api/v2/guids_and_examination_times', methods=['GET'])
def get_guids_examtime(**kwargs):
	""" returns all guids and their examination (addition) time.
	reference, if passed, is ignored."""
	try:	
		result =ew.get_all_guids_examination_time()	
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

        # TODO: test that it actually works


@app.route('/api/v2/annotations', methods=['GET'])
def get_guids_annotations(**kwargs):
	""" returns all guids and associated meta data.
	This query can be slow for very large data sets.
	The reference is ignored."""
	try:
		result = ew.get_all_annotations()
		
	except Exception as e:
		abort(500, e)
		
	return(tojson(result))

class test_get_guids_annotations_1(unittest.TestCase):
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
		result = ew.exist_sample(guid)
		
	except Exception as e:
		print(e)
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
		
		# note: additional testing is performed in test_insert_*


@app.route('/api/v2/insert', methods=['POST'])
def insert():
	""" inserts a guids with sequence, which it expects gzipped."""
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
			result = ew.insert(guid, seq)
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
	try:
		data_keys = set()
		for key in request.form.keys():
			data_keys.add(key)
		payload = {}
		for key in data_keys:
			payload[key]= request.form[key]
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)		
			
	return make_response(tojson(payload))

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

class test_insert_plus(unittest.TestCase):
	""" tests route /api/v2/insert """
	def runTest(self):
		relpath = "/api/v2/guids"
		res = do_GET(relpath)
		n_pre = len(json.loads(str(res.text)))		# get all the guids

		for i in range(1,10):
			guid_to_insert = "guid_{0}".format(n_pre+i)

			inputfile = "../COMPASS_reference/R39/R00000039.fasta"
			with open(inputfile, 'rt') as f:
				for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
						seq = list(str(record.seq))
						
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

class test_mirror(unittest.TestCase):
    """ tests route /api/v2/mirror """
    def runTest(self):
        
        relpath = "/api/v2/mirror"
        payload = {'guid':'1','seq':"ACTG"}
        res = do_POST(relpath, payload = payload)
        res_dict = json.loads(res.content.decode('utf-8'))
        self.assertEqual(payload, res_dict)
        self.assertTrue(isinstance(res_dict, dict))
        print(res.text)

@app.route('/api/v2/<string:guid>/neighbours_within/<int:threshold>', methods=['GET'])
@app.route('/api/v2/<string:guid>/neighbours_within/<int:threshold>/with_quality_cutoff/<float:cutoff>', methods=['GET'])
@app.route('/api/v2/<string:guid>/neighbours_within/<int:threshold>/with_quality_cutoff/<int:cutoff>', methods=['GET'])
@app.route('/api/v2/<string:guid>/neighbours_within/<int:threshold>/with_quality_cutoff/<float:cutoff>/in_format/<int:returned_format>', methods=['GET'])
@app.route('/api/v2/<string:guid>/neighbours_within/<int:threshold>/with_quality_cutoff/<int:cutoff>/in_format/<int:returned_format>', methods=['GET'])
@app.route('/api/v2/<string:guid>/neighbours_within/<int:threshold>/in_format/<int:returned_format>', methods=['GET'])
def query_get_value_snp(guid, threshold, **kwargs):
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
		result = ew.query_get_value_snp_filter(guid, threshold, cutoff, returned_format)
	except KeyError as e:
		# guid doesn't exist
		abort(404, e)
	except Exception as e:
		abort(500, e)
	
	return make_response(tojson(result))
	
class test_query_get_value_snp_1(unittest.TestCase):
    """ tests route /api/v2/guid/neighbours_within/ """
    def runTest(self):
        relpath = "/api/v2/non_existent_guid/neighbours_within/12"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), dict)
        self.assertEqual(res.status_code, 404)

class test_query_get_value_snp_2(unittest.TestCase):
    """ tests route /api/v2/guid/neighbours_within/ """
    def runTest(self):
        relpath = "/api/v2/non_existent_guid/neighbours_within/12/with_quality_cutoff/0.5"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), dict)
        self.assertEqual(res.status_code, 404)

class test_query_get_value_snp_3(unittest.TestCase):
    """ tests route /api/v2/guid/neighbours_within/ """
    def runTest(self):
        relpath = "/api/v2/non_existent_guid/neighbours_within/12/with_quality_cutoff/0.5/in_format/1"
        res = do_GET(relpath)
        print(res)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), dict)
        self.assertEqual(res.status_code, 404)

class test_query_get_value_snp_4(unittest.TestCase):
    """ tests route /api/v2/guid/neighbours_within/ """
    def runTest(self):
        relpath = "/api/v2/non_existent_guid/neighbours_within/12/with_quality_cutoff/0.5/in_format/2"
        res = do_GET(relpath)
        print(res)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), dict)
        self.assertEqual(res.status_code, 404)

class test_query_get_value_snp_5(unittest.TestCase):
    """ tests route /api/v2/guid/neighbours_within/ """
    def runTest(self):
        relpath = "/api/v2/non_existent_guid/neighbours_within/12/in_format/2"
        res = do_GET(relpath)
        print(res)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), dict)
        self.assertEqual(res.status_code, 404)
 
class test_query_get_value_snp_6(unittest.TestCase):
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

                
# @app.route('/api/v2/<string:guid1>/<string:guid2>/detailed_comparison', methods=['GET'])
# def get_detail(guid1, guid2):
# 	""" detailed comparison of two guids """
# 	try:
# 		result = ew.query_get_detail(guid1, guid2)
# 		
# 	except Exception as e:
# 		print("Exception raised", e)
# 		abort(500, e)
# 		
# 	return(result)
# 
# class test_get_detail(unittest.TestCase):
#     """ tests route /detailed_comparison """
#     def runTest(self):
#         relpath = "/api/v2/guid1/guid2/detailed_comparison"
#         res = do_GET(relpath)
# 
#         self.assertTrue(isjson(content = res.content))
#         info = json.loads(res.content.decode('utf-8'))
#         self.assertEqual(type(info), dict)
# 
#         self.assertEqual(info, {"guid1_exists": False, "success": 0, "guid2_exists": False})
#         self.assertEqual(res.status_code, 200)

@app.route('/api/v2/nucleotides_excluded', methods=['GET'])
def get_nucleotides_excluded():
	""" returns all nucleotides excluded by the server.
	Useful for clients which need to to ensure that server
	and client masking are identical. """
	
	try:
		result = ew.server_nucleotides_excluded()
		
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)

	return make_response(tojson(result))

class test_get_nucleotides_excluded(unittest.TestCase):
    """ tests route /api/v2/nucleotides_excluded"""
    def runTest(self):
        relpath = "api/v2/nucleotides_excluded"
        res = do_GET(relpath)
        print(res.text[0:1000])
        resDict = json.loads(res.text)
        self.assertTrue(isinstance(resDict, dict))
        self.assertEqual(set(resDict.keys()), set(['exclusion_id', 'excluded_nt']))
        self.assertEqual(res.status_code, 200)
 

if __name__ == '__main__':

        # command line usage.  Pass the location of a config file as a single argument.
        # an example config file is default_test_config.json
               
        ############################ LOAD CONFIG ######################################
       
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
        # see http://flask.pocoo.org/docs/dev/logging/               
        loglevel=logging.INFO
        if 'LOGLEVEL' in CONFIG.keys():
                if CONFIG['LOGLEVEL']=='WARN':
                        loglevel=logging.WARN
                elif CONFIG['LOGLEVEL']=='DEBUG':
                        loglevel=logging.DEBUG
        
        # configure logging object 
        app.logger.setLevel(loglevel)          
        
        # handles logging both with a stream to stderr and a rotating file
        rfh_handler = logging.handlers.RotatingFileHandler(CONFIG['LOGFILE'], maxBytes=100000, backupCount=5)
        stream_handler = logging.StreamHandler()

        formatter = logging.Formatter( "%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s ")
        rfh_handler.setFormatter(formatter)
        stream_handler.setFormatter(formatter)
        app.logger.addHandler(rfh_handler)
        app.logger.addHandler(stream_handler)
        
        ########################### prepare to launch server ###############################################################
        # construct the required global variables
        LISTEN_TO = '127.0.0.1'
        RESTBASEURL = "http://{0}:{1}".format(CONFIG['IP'], CONFIG['REST_PORT'])

        #########################  CONFIGURE HELPER APPLICATIONS ######################
		## configure mongodb persistence store
        PERSIST=fn3persistence(connString=CONFIG['FNPERSISTENCE_CONNSTRING'], debug=CONFIG['DEBUGMODE'])
        ew = ElephantWalk(CONFIG, PERSIST)
        
        ########################  START THE SERVER ###################################
        if CONFIG['DEBUGMODE']==1:
                app.logger.info("No config file name supplied ; using a configuration ('default_test_config.json') suitable only for testing, not for production. ")

        app.run(host=LISTEN_TO, debug=CONFIG['DEBUGMODE'], port = CONFIG['REST_PORT'])

