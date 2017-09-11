
#!/usr/bin/env python
""" 
provides restful interface to ElephantWalk2 functions
 
The endpoint provides designed as an internal API.
As implemented, it is not protected by authentication.

"""
 
# import libraries
import sys
import requests
import json
import logging
import xmlrpc.client
import gzip

# only used for testing
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide
import unittest


# cope with differential naming in python2/3
try:
    from urllib.parse import urlparse as urlparser
    from urllib.parse import urljoin as urljoiner
except ImportError:   
    from urlparse import urlparse as urlparser
    from urlparse import urljoin as urljoiner
 
from flask import Flask, make_response, jsonify
from flask import request, abort
 
RESTBASEURL = "http://127.0.0.1:5000"
XMLRPCBASEURL = "http://127.0.0.1:8184"          
app = Flask(__name__)
app.logger.setLevel(logging.DEBUG)

 
# --------------------------------------------------------------------------------------------------
@app.errorhandler(404)
def not_found(error):
    json_err = jsonify({'error': 'Not found'})
    return make_response(json_err, 404)
# --------------------------------------------------------------------------------------------------
 
 
# --------------------------------------------------------------------------------------------------

#@app.route('/findneighbour/<string:instance>/<string:reference>/annotation/', methods=['POST', 'GET'])
#@app.route('/findneighbour/<string:instance>/<string:reference>/guids_and_time/', methods=['POST', 'GET'])
 
#@app.route('/findneighbour/<string:instance>/<string:reference>/guids_cutoff/<float:quality>', methods=['POST', 'GET'])
#@app.route('/findneighbour/<string:instance>/<string:reference>/guids_cutoff/<int:quality>', methods=['POST', 'GET'])
 
#@app.route('/findneighbour/<string:instance>/<string:reference>/snp/<string:guid>/<int:SNPCutoff>/<float:qualityCutoff>', methods=['POST', 'GET'])
#@app.route('/findneighbour/<string:instance>/<string:reference>/snp/<string:guid>/<int:SNPCutoff>/<int:qualityCutoff>', methods=['POST', 'GET'])
 
#@app.route('/findneighbour/<string:instance>/<string:reference>/server_memory_usage', methods=['POST', 'GET'])

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
	print("About to contact url {0}".format(url))
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
	print("POSTING to: {0}".format(url))

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

@app.route('/server_config', methods=['GET'])
def server_config():
	""" returns server configuration """
	try:
		client=get_client()
		result = client.server_config()
		
	except Exception as e:
		print("Exception raised", e)
		abort(502, e)
	return(str(result))	
class test_server_config(unittest.TestCase):
    """ tests route /server_config"""
    def runTest(self):
        relpath = "/server_config"
        res = do_GET(relpath)
        config_dict = json.loads(str(res.text))
        self.assertTrue('MASKER' in config_dict.keys())
        print(res)
        self.assertEqual(res.status_code, 200)


@app.route('/server_memory_usage', methods=['GET'])
def server_memory_usage():
	""" returns server memory usage """
	try:
		client=get_client()
		result = client.server_memory_usage()
		
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
	return(str(result))		
class test_server_memory_usage(unittest.TestCase):
    """ tests route /server_memory_usage"""
    def runTest(self):
        relpath = "/server_memory_usage"
        res = do_GET(relpath)
        config_dict = json.loads(str(res.text))
        self.assertTrue('note' in config_dict.keys())
        print(res)
        self.assertEqual(res.status_code, 200)


@app.route('/server_time', methods=['GET'])
def server_time():
	""" returns server time """
	try:
		client=get_client()
		result = client.server_time()
		
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
	return(str(result))		
class test_server_time(unittest.TestCase):
    """ tests route /server_time"""
    def runTest(self):
        relpath = "/server_time"
        res = do_GET(relpath)
        config_dict = json.loads(str(res.text))
        self.assertTrue('server_time' in config_dict.keys())
        print(res)
        self.assertEqual(res.status_code, 200)

@app.route('/get_all_guids', methods=['GET'])
def get_all_guids():
	""" returns all guids """
	try:
		client=get_client()
		result = client.get_all_guids()	
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
	return(str(result))
class test_get_all_guids(unittest.TestCase):
    """ tests route /get_all_guids"""
    def runTest(self):
        relpath = "/get_all_guids"
        res = do_GET(relpath)
        guidlist = json.loads(str(res.text))
        print(guidlist)
        self.assertEqual(res.status_code, 200)


@app.route('/get_all_filtered_guids/<float:cutoff>', methods=['GET'])
def get_all_filtered_guids(cutoff):
	""" returns all guids with quality score >= cutoff """
	try:
		client=get_client()
		result = client.get_all_filtered_guids(cutoff)	
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
	return(str(result))
class test_get_all_filtered_guids(unittest.TestCase):
    """ tests route /get_all_filtered_guids"""
    def runTest(self):
        relpath = "/get_all_filtered_guids/0.7"
        res = do_GET(relpath)
        guidlist = json.loads(str(res.text))
        print(guidlist)
        self.assertEqual(res.status_code, 200)


@app.route('/get_all_guids_examination_time', methods=['GET'])
def get_guids_examtime():
	""" returns all guids and their examination (addition) time """
	try:
		client=get_client()
		result = client.get_all_guids_examination_time()	
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
	return(str(result))
class test_get_all_guids_examination_time(unittest.TestCase):
    """ tests route /get_all_guids_examination_time"""
    def runTest(self):
        relpath = "/get_all_guids_examination_time"
        res = do_GET(relpath)
        guidlist = json.loads(str(res.text))
        print(guidlist)
        self.assertEqual(res.status_code, 200)


@app.route('/get_all_annotations', methods=['GET'])
def get_guids_annotations():
	""" returns all guids and associated meta data.
	This query can be slow for very large data sets."""
	try:
		client=get_client()
		result = client.get_all_annotations()
		
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
		
	return(str(result))
class test_get_guids_annotations(unittest.TestCase):
    """ tests route /get_guids_annotations """
    def runTest(self):
        relpath = "/get_all_guids_examination_time"
        res = do_GET(relpath)
        guidlist = json.loads(str(res.text))
        print(guidlist)
        self.assertEqual(res.status_code, 200)


@app.route('/exist_sample/<string:guid>', methods=['GET'])
def exist_sample(guid):
	""" checks whether a guid exists"""
	try:
		client=get_client()
		result = client.exist_sample(guid)
		
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
		
	return(str(result))
class test_exist_sample(unittest.TestCase):
    """ tests route /exist_sample """
    def runTest(self):
        relpath = "/exist_sample/non_existent_guid"
        res = do_GET(relpath)
       
        self.assertEqual(res.text, 'False')
        self.assertEqual(res.status_code, 200)

@app.route('/query_get_value_snp/<string:guid>/<int:threshold>', methods=['GET'])
def query_get_value_snp(guid, threshold):
	""" get a guid's neighbours, within a threshold"""
	try:
		client=get_client()
		result = client.query_get_value_snp(guid, threshold)
		
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
		
	return(str(result))
class test_query_get_value_snp(unittest.TestCase):
    """ tests route /query_get_value_snp """
    def runTest(self):
        relpath = "/query_get_value_snp/non_existent_guid"
        res = do_GET(relpath)
       
        self.assertEqual(res.text, '{\n  "error": "Not found"\n}\n')
        self.assertEqual(res.status_code, 404)

query_get_value_snp

@app.route('/insert', methods=['POST'])
def insert():
	""" inserts a guids with sequence, which it expects gzipped."""
	try:
		client=get_client()
		
		seq_data = request.form
		if 'seq' in seq_data.keys() and 'guid' in seq_data.keys():
			seq = seq_data['seq']				
			result = client.insert(seq_data['guid'], seq_data['seq'])
		else:
			abort(501, 'seq and guid are not present in the POSTed data {0}'.seq_data.keys())
		
	except Exception as e:
		print("Exception raised", e)
		abort(500, e)
		
	return(str(result))
class test_insert(unittest.TestCase):
    """ tests route /insert """
    def runTest(self):
        relpath = "/get_all_guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))

        inputfile = "../COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, 'rt') as f:
            for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                    seq = str(record.seq)

        print("Adding TB reference sequence of {0} bytes".format(len(seq)))
        self.assertEqual(len(seq), 4411532)		# check it's the right sequence

        relpath = "/insert"
        res = do_POST(relpath, payload = {'guid':'guid','seq':seq})
        
        relpath = "/get_all_guids"
        res = do_GET(relpath)
        n_post = len(json.loads(str(res.text)))
        
        self.assertEqual(n_pre+1, n_post)

        relpath = "/exist_sample/guid"
        res = do_GET(relpath)
        self.assertEqual(res.text, 'True')


# @app.route('/sample_get_snp/ndneighbour/<string:instance>/<string:reference>/server_memory_usage', methods=['POST', 'GET'])
# def sample_get_snp(self,sample,reference,threshold,method):
# 
# 	url=compassconfig.COMPASSCFG['elephantwalk'][reference.lower()]
# 
# 	if method == "elephantwalk":
# 		url=compassconfig.COMPASSCFG[method][reference.lower()]
# 		client=None
# 	
# 
# 		##check connection
# 		try:
# 			client=xmlrpclib.ServerProxy(url)
# 
# 			try:
# 				client._()   # Call a fictive method.
# 			except xmlrpclib.Fault:
# 				pass
# 			except socket.error:
# 				raise Exception("Error in the server side")
# 
# 			##query to walk
# 			print str(sample), threshold
# 			result = str(client.query_get_value_snp(str(sample),threshold))
# 			result = json.loads(result.replace('\'','\"'))
# 
# 			error = str(result[0])
# 			snps = str(result[1]).replace('\'','\"').replace('u','')
# 
# 			if error == 'Err':
# 				return "[\"Err\", \"{0}\"]".format(result[1])	
# 
# 			res = snps
# 
# 		except Exception,e:
# 			return "[\"Err\",\"{0}\"]".format(e)
# 		
# 		return res
# 
# @app.route('/findneighbour/<string:instance>/<string:reference>/server_memory_usage', methods=['POST', 'GET'])
# 	
# 	def GET(self,sample,reference,threshold,method):
# 
# 		url=compassconfig.COMPASSCFG['elephantwalk'][reference.lower()]
# 
# 		if method == "elephantwalk":
# 			url=compassconfig.COMPASSCFG[method][reference.lower()]
# 			client=None
# 		
# 
# 			##check connection
# 			try:
# 				client=xmlrpclib.ServerProxy(url)
# 
# 				try:
# 					client._()   # Call a fictive method.
# 				except xmlrpclib.Fault:
# 					pass
# 				except socket.error:
# 					raise Exception("Error in the server side")
# 
# 				##query to walk
# 				print str(sample), threshold
# 				result = str(client.query_get_value_snp(str(sample),threshold))
# 				result = json.loads(result.replace('\'','\"'))
# 
# 				error = str(result[0])
# 				snps = str(result[1]).replace('\'','\"').replace('u','')
# 
# 				if error == 'Err':
# 					return "[\"Err\", \"{0}\"]".format(result[1])	
# 
# 				res = snps
# 	
# 			except Exception,e:
# 				return "[\"Err\",\"{0}\"]".format(e)
# 			
# 			return res
# 
# 
# def action(*args, **kwargs):
#  
#     app.logger.debug("command args: %s", str(kwargs))
#  
#     # get digest from message and make one yourself for comparison
#     supposed_digest = request.headers['{0}-HMAC-signature'.format(ENDPOINT_ID)].encode('utf-8')
#     #print(request.headers['{0}-HMAC-signature'.format(ENDPOINT_ID)])
#     actual_digest = make_digest(request.headers['{0}-HMAC-Input'.format(ENDPOINT_ID)])
#  
#     # abort with code 401 of they don't match
#     if not supposed_digest == actual_digest:
#         app.logger.debug('Digests do not match: received = {0}; expected = {1}'.format(supposed_digest, actual_digest))
#         abort(401)
#     else:
#         app.logger.debug('Digests match; analysing request')
#        
#     # re-jsonify input part of message
#     HmacInput = json.loads(request.headers['{0}-HMAC-Input'.format(ENDPOINT_ID)])
#  
#     # get time message was sent to compare with now -> problem if server time incorrect!
#     message_time = HmacInput['timestamp']
#     pat = '%Y-%m-%d %H:%M:%S.%f'
#     message_time = int(time.mktime(time.strptime(message_time, pat)))
#     nowish = int(time.time())
#  
#     app.logger.debug('time now: %s, message sent: %s', nowish, message_time)
#  
#     # get further information from input part of message
#     sender_account = HmacInput['accountkey'].encode('utf-8')
#     sent_to_relpath = HmacInput['url']
#  
#     opar = urlparser(request.url)
#  
#     app.logger.critical("sent from account %s, authorised account %s", sender_account, ACCOUNT)
#     app.logger.critical("url in message: %s, url in request: %s", sent_to_relpath, opar.path)
#  
#     # if message was sent more than 2 minutes ago or accounts or urls don't match
#     if (message_time + 120) < nowish or sender_account != ACCOUNT or sent_to_relpath != opar.path:
#         # abort with code 401
#         app.logger.critical("Auth failed")
#         abort(401)
#        
#     else:
#         # forward legitimate requests; we positively identify such
#         fwdurl = None
#         app.logger.critical("analysing %s", opar.path)
#          
#         if 'forestapi' in opar.path:
#             fwdurl =urljoiner(FORESTAPI, opar.path)
#            
#         elif 'fileserver' in opar.path:
#             fwdurl = urljoiner(FILESERVERAPI, opar.path)
#            
#         elif 'findneighbour' in opar.path:
#             app.logger.critical(kwargs.keys())
#             app.logger.critical(kwargs)
#             if 'instance' in kwargs.keys() and 'reference' in kwargs.keys():
#                 # try to construct a key from the 'instance' and 'reference' portions of the URL
#                 ew_key =       "{0}:{1}".format(kwargs['instance'],kwargs['reference'])
#                 ew_urlstart =  "/findneighbour/{0}/{1}".format(kwargs['instance'],kwargs['reference'])
#                 app.logger.critical("constructed {0} {1}".format(ew_key, ew_urlstart))
#          
#                 if ew_key in ELEPHANTWALKS.keys():
#                     ip = ELEPHANTWALKS[ew_key]
#                     fwdurl = opar.path.replace(ew_urlstart, ip)
#                 else:
#                     app.logger.critical("Endpoint requested combination of {0}, but this is not supported.  it may need to be added to the ELEPHANTWALKS dictionary.".format(ew_key))
#                     abort(403)      # combination is forbidden
#                    
#         if fwdurl is not None:
#  
#             app.logger.critical("legitimate request forwarded to %s", fwdurl)
#  
#             try:           
#                 if request.method == 'POST':
#                        resp = requests.post(url=fwdurl, timeout=None)
#                 elif request.method== 'GET':
#                        resp = requests.get(url=fwdurl, timeout=None)
#                 else:
#                     abort(405)      # method not allowed
#             except requests.exceptions.ConnectionError:       # upstream connectivity failure
#                 app.logger.error("URL timed out: %s", fwdurl)
#                 abort(504)      # gateway timeout
#                
#             # we got a response
#             app.logger.critical("response code: %s", resp.status_code)
#             app.logger.critical("response reason: %s", resp.reason)
#             app.logger.info("response text <clipped>: %s", resp.text[:100])
#             app.logger.info("response headers: %s", resp.headers)
#    
#             # return the response from the forwarding url
#             return (resp.text, resp.status_code, tuple())
#    
#         else:
#             abort(401)

 
# --------------------------------------------------------------------------------------------------
 
if __name__ == '__main__':

	app.run(debug=True)
