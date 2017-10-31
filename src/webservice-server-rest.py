#!/usr/bin/env python
""" 
provides restful interface to ElephantWalk2 functions
 
The endpoint provided is designed as an internal API.
As implemented, it is not protected by authentication.

requires python3.

loads configuration from a config file.
If no config file is provided, it will run in 'testing' mode with the following default parameters:

RESTBASEURL = "http://127.0.0.1:5000"
XMLRPCBASEURL = "http://127.0.0.1:8184"
ISDEBUG = True

The config file must include three parameters:
IP - the IP of this server, and the XMLRPC server which are assumed to be the same;
PORT - the IP of the XMLRPC server;
REST_PORT - the port on which this server should run.

REST_PORT can be added to the CONFIG file used by the findNeighbour2 XMLRPC server
and the same config file used for both servers.


Unit testing can be achieved by
# starting a test XMLRPC server
python3 webservice-server.py &

# starting a test RESTFUL server
python3 webservice-server-rest.py

# And then (e.g. in a different terminal) launching unit tests with
python3 -m unittest webservice-server-rest
# all should pass
"""
 
# import libraries
import os
import sys
import requests
import json
import logging
import xmlrpc.client
import gzip
import warnings
from flask import Flask, make_response, jsonify
from flask import request, abort
from logging.config import dictConfig

# only used for unit testing
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide
import unittest
from urllib.parse import urlparse as urlparser
from urllib.parse import urljoin as urljoiner

# default parameters for unit testing only.
RESTBASEURL   = "http://127.0.0.1:5000"
XMLRPCBASEURL = "http://127.0.0.1:8184"
ISDEBUG = True
LISTEN_TO = '127.0.0.1'		# only local addresses
 
app = Flask(__name__)
app.logger.setLevel(logging.DEBUG)


def isjson(content):
        """ returns true if content parses as json, otherwise false. used by unit testing. """
        try:
            x = json.loads(content.decode('utf-8'))
            print("JSON DECODE SUCCEEDED : {0}".format(content.decode('utf-8')))
            return True
 
        except json.decoder.JSONDecodeError:
            print("JSON DECODE FAILED : {0}".format(content.decode('utf-8')))
            return False


# --------------------------------------------------------------------------------------------------
@app.errorhandler(404)
def not_found(error):
    json_err = jsonify({'error': 'Not found'})
    return make_response(json_err, 404)
# --------------------------------------------------------------------------------------------------
 
def get_client():
	""" instantiates an xmlrpc client """
	try:
		client=xmlrpc.client.ServerProxy(XMLRPCBASEURL)

		try:
			client._()   # Call a fictive method.
			
		except xmlrpc.client.Fault:
			pass
		
		except socket.error:
			abort(500, 'Cannot connect to upstream XMLRPC server')
	
	# untrapped error on instantiating client.
	except Exception as e:
			abort(501,e)

	return(client)

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


@app.route('/v2/server_config', methods=['GET'])
def server_config():
	""" returns server configuration """
	try:
		client=get_client()
		result = client.server_config()
		
	except Exception as e:
		print("Exception raised", e)
		abort(502, e)

	return(result)	
class test_server_config(unittest.TestCase):
    """ tests route v2/server_config"""
    def runTest(self):
        relpath = "/v2/server_config"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))

        config_dict = json.loads(res.content.decode('utf-8'))
        self.assertTrue('MASKER' in config_dict.keys())
        print(res)
        self.assertEqual(res.status_code, 200)


@app.route('/v2/server_memory_usage', methods=['GET'])
def server_memory_usage():
	""" returns server memory usage """
	try:
		client=get_client()
		result = client.server_memory_usage()

	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
	return(result)


class test_server_memory_usage(unittest.TestCase):
    """ tests route /v2/server_memory_usage"""
    def runTest(self):
        relpath = "/v2/server_memory_usage"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))

        config_dict = json.loads(res.content.decode('utf-8'))
        self.assertTrue('note' in config_dict.keys())
        print(res)
        self.assertEqual(res.status_code, 200)


@app.route('/v2/server_time', methods=['GET'])
def server_time():
	""" returns server time """
	try:
		client=get_client()
		result = client.server_time()
		
	except Exception as e:
		abort(500, e)
	return(result)

class test_server_time(unittest.TestCase):
    """ tests route /v2/server_time"""
    def runTest(self):
        relpath = "/v2/server_time"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        config_dict = json.loads(res.content.decode('utf-8'))
        self.assertTrue('server_time' in config_dict.keys())
        self.assertEqual(res.status_code, 200)

	
@app.route('/sample/guids/<string:reference>', methods= ['GET'])
@app.route('/v2/guids', methods=['GET'])
def get_all_guids(**kwargs):
	""" returns all guids.  reference, if included, is ignored."""
	try:
		client=get_client()
		result = client.get_all_guids()	
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
	return(result)


class test_get_all_guids_1(unittest.TestCase):
    """ tests route /v2/guids"""
    def runTest(self):
        relpath = "/v2/guids"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        guidlist = json.loads(str(res.content.decode('utf-8')))
        print(guidlist)
        self.assertEqual(res.status_code, 200)

class test_get_all_guids_2(unittest.TestCase):
    """ tests route /sample/guids"""
    def runTest(self):
        relpath = "/sample/guids/R00039"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        guidlist = json.loads(res.content.decode('utf-8'))
        print(guidlist)
        self.assertEqual(res.status_code, 200)

@app.route('/sample/guids_cutoff/<string:reference>/<float:cutoff>', methods=['GET'])		
@app.route('/v2/guids_with_quality_over/<float:cutoff>', methods=['GET'])
@app.route('/sample/guids_cutoff/<string:reference>/<int:cutoff>', methods=['GET'])		
@app.route('/v2/guids_with_quality_over/<int:cutoff>', methods=['GET'])
def get_all_filtered_guids(cutoff, **kwargs):
	""" returns all guids with quality score >= cutoff.
	reference, if provided, ignored."""
	try:
		client=get_client()
		result = client.get_all_filtered_guids(cutoff)	
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
	return(result)


class test_get_all_filtered_guids_1(unittest.TestCase):
    """ tests route /v2/guids_with_quality_over"""
    def runTest(self):
        relpath = "/v2/guids_with_quality_over/0.7"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        guidlist = json.loads(res.content.decode('utf-8'))
        print(guidlist)
        self.assertEqual(res.status_code, 200)

		
class test_get_all_filtered_guids_2(unittest.TestCase):
    """ tests route /sample/guids_cutoff"""
    def runTest(self):
        relpath = "/sample/guids_cutoff/R00039/0.7"
        res = do_GET(relpath)
        guidlist = json.loads(str(res.text))
        print(guidlist)
        self.assertEqual(res.status_code, 200)

@app.route('/sample/guids_and_time/<string:reference>', methods=['GET'])
@app.route('/v2/guids_and_examination_times', methods=['GET'])
def get_guids_examtime(**kwargs):
	""" returns all guids and their examination (addition) time.
	reference, if passed, is ignored."""
	try:
		client=get_client()
		result = client.get_all_guids_examination_time()	
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
	return(result)


class test_get_all_guids_examination_time_1(unittest.TestCase):
    """ tests route /v2/guids_and_examination_times"""
    def runTest(self):
        relpath = "/v2/guids_and_examination_times"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        guidlist = json.loads(res.content.decode('utf-8'))
        print(guidlist)
        self.assertEqual(res.status_code, 200)
		
class test_get_all_guids_examination_time_2(unittest.TestCase):
    """ tests route /sample/guids_and_time/R00039"""
    def runTest(self):
        relpath = "/sample/guids_and_time/R00039"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
 
        guidlist = json.loads(res.content.decode('utf-8'))
        print(guidlist)
        self.assertEqual(res.status_code, 200)

@app.route('/sample/annotation/<string:reference>', methods=['GET'])
@app.route('/v2/annotations', methods=['GET'])
def get_guids_annotations(**kwargs):
	""" returns all guids and associated meta data.
	This query can be slow for very large data sets.
	The reference is ignored."""
	try:
		client=get_client()
		result = client.get_all_annotations()
		
	except Exception as e:
		abort(500, e)
		
	return(result)


class test_get_guids_annotations_1(unittest.TestCase):
    """ tests route /v2/annotations """
    def runTest(self):
        relpath = "/v2/annotations"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content = res.content))
 
        guidlist = json.loads(res.content.decode('utf-8'))
        print(guidlist)

		
class test_get_guids_annotations_2(unittest.TestCase):
    """ tests route /sample/annotation/R0039 """
    def runTest(self):
        relpath = "/sample/annotation/R0039"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content = res.content))
 
        guidlist = json.loads(res.content.decode('utf-8'))
        print(guidlist)


@app.route('/sample/walks/processed/<string:guid>/<string:reference>/<string:method>', methods=['GET'])
@app.route('/v2/<string:guid>/exists', methods=['GET'])
def exist_sample(guid, **kwargs):
	""" checks whether a guid exists.
	reference and method are ignored."""
	
	try:
		client=get_client()
		result = client.exist_sample(guid)
		
	except Exception as e:
		abort(500, e)
		
	return(json.dumps(result))

class test_exist_sample(unittest.TestCase):
    """ tests route /v2/guid/exists """
    def runTest(self):
        relpath = "/v2/non_existent_guid/exists"
        res = do_GET(relpath)
       
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), bool)
        self.assertEqual(info, False)

@app.route('/sample/walks/snp/<string:guid>/<string:reference>/<int:threshold>/<string:method>', methods = ['GET'])
@app.route('/sample/findneighbour/snp/<string:guid>/<string:reference>/<int:threshold>/<string:method>/<float:cutoff>', methods = ['GET'])
@app.route('/sample/findneighbour/snp/<string:guid>/<string:reference>/<int:threshold>/<string:method>/<int:cutoff>', methods = ['GET'])
@app.route('/sample/neighbours/<string:guid>/<string:reference>/<int:threshold>/<string:method>', methods=['GET'])
@app.route('/v2/<string:guid>/neighbours_within<int:threshold>', methods=['GET'])
@app.route('/v2/<string:guid>/<int:threshold>/neighbours_within/<float:cutoff>', methods=['GET'])
@app.route('/v2/<string:guid>/<int:threshold>/neighbours_within/<float:cutoff>/<int:returned_format>', methods=['GET'])
@app.route('/v2/<string:guid>/<int:threshold>/neighbours_within/<int:cutoff>', methods=['GET'])
@app.route('/v2/<string:guid>/<int:threshold>/neighbours_within/<int:cutoff>/<int:returned_format>', methods=['GET'])
def query_get_value_snp(guid, threshold, **kwargs):
	""" get a guid's neighbours, within a threshold """
	
	# we support optional cutoff and threshold parameters.
	# we also support 'method' and 'reference' parameters but these are ignored.
	# the default for cutoff and format are 0.85 and 1, respectively.
	if not 'cutoff' in kwargs.keys():
		cutoff = 0.85
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
		abort(500, "Invalid format requested, must be between 0 and 1")
		
	try:
		client=get_client()	
		result = client.query_get_value_snp_filter(guid, threshold, cutoff, returned_format)
	
	except Exception as e:
		abort(500, e)
	
	return(result)		# already json

class test_query_get_value_snp_0a(unittest.TestCase):
    """ tests route /sample/findneighbour/snp/nonexistent_guid/R00039/12/elephantwalk/0.85	"""
    def runTest(self):
        relpath = "/sample/findneighbour/snp/nonexistent_guid/R00039/12/elephantwalk/0.85"

        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
       
        self.assertTrue(info[1] == "missing sample")
        self.assertTrue(info[0] == "Err")

class test_query_get_value_snp_0b(unittest.TestCase):
    """ tests route '/sample/neighbours/<string:guid>/<string:reference>/<int:threshold>/<string:method>' """
    def runTest(self):
        relpath = "/sample/neighbours/nonexistent_guid/R00039/12/elephantwalk"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
       
        self.assertTrue(info[1] == "missing sample")
        self.assertTrue(info[0] == "Err")

		
class test_query_get_value_snp_1(unittest.TestCase):
    """ tests route /v2/guid/neighbours_within/ """
    def runTest(self):
        relpath = "/v2/non_existent_guid/neighbours_within/12"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), dict)
        self.assertTrue('error' in info.keys())
        self.assertEqual(res.status_code, 404)

class test_query_get_value_snp_2(unittest.TestCase):
    """ tests route /v2/guid/neighbours_within/ """
    def runTest(self):
        relpath = "/v2/non_existent_guid/neighbours_within/12/0.5"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), dict)
        self.assertTrue('error' in info.keys())
        self.assertEqual(res.status_code, 404)

class test_query_get_value_snp_3(unittest.TestCase):
    """ tests route /v2/guid/neighbours_within/ """
    def runTest(self):
        relpath = "/v2/non_existent_guid/neighbours_within/12/0.5/1"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), dict)
        self.assertTrue('error' in info.keys())
        self.assertEqual(res.status_code, 404)

class test_query_get_value_snp_4(unittest.TestCase):
    """ tests route /v2/guid/neighbours_within/ """
    def runTest(self):
        relpath = "/v2/query_get_value_snp/non_existent_guid/12/0.5/2"
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), dict)
        self.assertTrue('error' in info.keys())
        self.assertEqual(res.status_code, 404)

@app.route('/v2/<string:guid1>/<string:guid2>/detailed_comparison', methods=['GET'])
def get_detail(guid1, guid2):
	""" detailed comparison of two guids """
	try:
		client=get_client()
		result = client.query_get_detail(guid1, guid2)
		
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
		
	return(result)

class test_query_get_detail(unittest.TestCase):
    """ tests route /query_get_detail """
    def runTest(self):
        relpath = "/v2/guid1/guid2/detailed_comparison"
        res = do_GET(relpath)

        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), dict)

        self.assertEqual(info, {"guid1_exists": False, "success": 0, "guid2_exists": False})
        self.assertEqual(res.status_code, 200)



@app.route('/v2/insert', methods=['POST'])
def insert():
	""" inserts a guids with sequence, which it expects gzipped."""
	try:
		data_keys = set()
		for key in request.form.keys():
			data_keys.add(key)
		payload = {}
		for key in data_keys:
			payload[key]= request.form[key]
		client=get_client()
	
		if 'seq' in data_keys and 'guid' in data_keys:
			guid = str(payload['guid'])
			seq  = str(payload['seq'])
			result = client.insert(guid, seq)
		else:
			abort(501, 'seq and guid are not present in the POSTed data {0}'.format(data_keys))
		
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
		
	return(result)

@app.route('/v2/mirror', methods=['POST'])
def mirror():
	""" receives data, returns the dictionary it was passed. Takes no other action.
	Used for testing that gateways etc don't remove data."""
	try:
		print(request.headers)		
		data_keys = set()
		for key in request.form.keys():
			data_keys.add(key)
		payload = {}
		for key in data_keys:
			payload[key]= request.form[key]
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)		
			
	return(json.dumps(payload))

class test_insert(unittest.TestCase):
    """ tests route /v2/insert """
    def runTest(self):
        relpath = "/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))

        guid_to_insert = "guid_{0}".format(n_pre+1)

        inputfile = "../COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, 'rt') as f:
            for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                    seq = str(record.seq)

        print("Adding TB reference sequence of {0} bytes".format(len(seq)))
        self.assertEqual(len(seq), 4411532)		# check it's the right sequence

        relpath = "/v2/insert"
        res = do_POST(relpath, payload = {'guid':guid_to_insert,'seq':seq})
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(info, ['OK'])

        relpath = "/v2/guids"
        res = do_GET(relpath)
        n_post = len(json.loads(res.content.decode('utf-8')))
        self.assertEqual(n_pre+1, n_post)
                

        # check if it exists
        relpath = "/v2/{0}/exists".format(guid_to_insert)
        res = do_GET(relpath)
        self.assertTrue(isjson(content = res.content))
        info = json.loads(res.content.decode('utf-8'))
        self.assertEqual(type(info), bool)
        self.assertEqual(res.status_code, 200)
        self.assertEqual(info, True)


class test_mirror(unittest.TestCase):
    """ tests route /v2/mirror """
    def runTest(self):
        
        relpath = "/v2/mirror"
        payload = {'guid':'1','seq':"ACTG"}
        res = do_POST(relpath, payload = payload)
        res_dict = json.loads(res.content.decode('utf-8'))
        self.assertEqual(payload, res_dict)
        self.assertTrue(isinstance(res_dict, dict))
        print(res.text)


@app.route('/v2/neighbours_within/<int:threshold>', methods=['GET'])
def get_all_values(threshold):
	""" get all pairwise distances, within and including a SNP threshold"""
	try:
		client=get_client()
		result = client.get_all_values(threshold)
		
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
		
	return(str(result))
class test_get_all_values(unittest.TestCase):
    """ tests route /v2/neighbours_within/threshold"""
    def runTest(self):
        relpath = "v2/neighbours_within/12"
        res = do_GET(relpath)
        resList = json.loads(res.text)
        self.assertTrue(isinstance(resList, list))
        self.assertEqual(res.status_code, 200)

@app.route('/v2/nucleotides_excluded', methods=['GET'])
def get_nucleotides_excluded():
	""" returns all nucleotides excluded by the server.
	Useful for clients which need to to ensure that server
	and client masking are identical. """
	
	try:
		client=get_client()
		result = client.server_nucleotides_excluded()
		
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)

	return(str(result))

class test_get_nucleotides_excluded(unittest.TestCase):
    """ tests route /v2/nucleotides_excluded"""
    def runTest(self):
        relpath = "v2/nucleotides_excluded"
        res = do_GET(relpath)
        print(res.text[0:100])
        resDict = json.loads(res.text)
        self.assertTrue(isinstance(resDict, dict))
        self.assertEqual(set(resDict.keys()), set(['exclusion_id', 'excluded_nt']))
        self.assertEqual(res.status_code, 200)
 

if __name__ == '__main__':

        # command line usage.  Pass the location of a config file as a single argument.
        # an example config file is default_test_config.json
               
        ############################ LOAD CONFIG ######################################
        root = logging.getLogger()
        root.setLevel(logging.DEBUG)
        
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
        required_keys=set(['IP','REST_PORT', 'PORT', 'DEBUGMODE', 'LOGFILE'])
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
        
        
        # if we are debug mode, then we remove any log file and note we are in debugmode to the logfile
        if CONFIG['DEBUGMODE']==1:
                try:
                        os.unlink(CONFIG['LOGFILE'])
                except FileNotFoundError:
                        pass            # if we can't delete the file, that's OK
                
        # configure logging object 
        app.logger.setLevel(loglevel)          
        root.setLevel(loglevel)
        
        # handles logging both with a stream to stderr and a rotating file
        rfh_handler = logging.handlers.RotatingFileHandler(CONFIG['LOGFILE'], maxBytes=100000, backupCount=5)
        stream_handler = logging.StreamHandler()

        formatter = logging.Formatter( "%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s ")
        rfh_handler.setFormatter(formatter)
        stream_handler.setFormatter(formatter)
        
        # we add these log handler to the root logger.  That way, all module output (independent of flask) will go there too
        root.addHandler(rfh_handler)
        root.addHandler(stream_handler)
        
        ########################### prepare to launch server ###############################################################
        # construct the required global variables
        LISTEN_TO = '0.0.0.0'
        RESTBASEURL = "http://{0}:{1}".format(CONFIG['IP'], CONFIG['REST_PORT'], CONFIG['PORT'])
        XMLRPCBASEURL = "http://{0}:{2}".format(CONFIG['IP'], CONFIG['REST_PORT'], CONFIG['PORT'])
        
      
        ########################  START THE SERVER ###################################
        app.logger.info("REST Server operating on {0};\n expects communication with XMLRPC server on {1}".format(RESTBASEURL, XMLRPCBASEURL))
        if CONFIG['DEBUGMODE']==1:
                app.logger.info("No config file name supplied ; using a configuration ('default_test_config.json') suitable only for testing, not for production. ")
                app.logger.info("This test configuration assumes the findNeighbour2 XMLRPC server is running at XMLRPCBASEURL (see config file above)")
       
        app.run(host=LISTEN_TO, debug=CONFIG['DEBUGMODE'], port = CONFIG['REST_PORT'])

