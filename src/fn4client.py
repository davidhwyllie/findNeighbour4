#!/usr/bin/env python3
""" A python3 client for findNeighbour3-server.

Provides a class which allows access to all routes described in:
../doc/rest-routes.md for a list.
  
Unit testing:
* launch a test server
python findNeighbour3-server.py

* run fn4Client unit tests
python -m unittest fn4client
"""

import glob
import hashlib
import hmac
import base64
import datetime
import collections
import requests
import json
import urllib.parse
import codecs
import unittest
import os
import time
import pymongo
import warnings
import logging
import gzip
import io
import pandas as pd

# used for loading fasta files.
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# used for unittesting
import uuid

class fn4Client():
    """ python3 API to the findNeighbour3 + -server REST endpoint.
    
        All endpoints are supported.
        See ../doc/rest-routes.md for a list.
                      
        """

    def __init__(self,
                baseurl="http://127.0.0.1:5020",
                ):

        """ set up logging """
        logging.getLogger()
        
        # load config
        self.baseurl = baseurl
        
        # run connection check
        res =  self.server_time()
        logging.info("Connection established at server time", res['server_time'])

    def _decode(self, response):
        """ checks that a response object has response code in the 200s.
        If it doesn't, raises an error.
        If it does, decodes the object and returns json."""
        response.raise_for_status()
        return json.loads(response.content.decode('utf-8'))
  
    def _isjson(self, content):
        """ returns true if content parses as json, otherwise false """
        try:
            x = json.loads(content.decode('utf-8'))
            return True

        except json.decoder.JSONDecodeError:
            return False

    def absurl(self, relpath):
        """ constructs an absolute URL from the relative path requested """
        # removed: self.version+
        absurl = urllib.parse.urljoin(self.baseurl,relpath)
        return(absurl)

    def getpost(self, relpath, method='GET', payload=None, timeout=None):
        """ issues GET or POST against url.  returns a response object
	    will raise errors if generated, but does not raise errors if a valid
            response object is returned, even if it has a status_code indicating the call failed
        """
        
        session = requests.Session()
        session.trust_env = False
    
        sc = -1
        url = self.absurl(relpath)
        if method == 'GET':
            response = session.get(url=url)

        elif method == 'POST':
            response = session.post(url=url,
                                        data = payload)
                          
        else:
            raise NotImplementedError("either GET or POST is required as a method.  was passed {0}".format(method))
        
        session.close() 
        
        if response.status_code >= 500: 
            response.raise_for_status()     # raise error if there's an error.  Sub 500 (404 etc) are let through and handled by the client routines
        return response 

    def get(self, relpath,  timeout=None):
        """ issues GET against url.  returns a response object.
        raises errors if fails"""
        response = self.getpost(relpath=relpath, timeout=timeout, method='GET')
        return response

    def post(self, relpath, payload,  timeout=None):
        """ issues  POST against url.  returns a response object. """
        response = self.getpost(relpath= relpath, payload = payload, timeout=timeout, method='POST')
        return response 
    
    def mirror(self, payload, timeout =0.2):
        """ returns payload via server round trip.
        payload much be a dictionary, with key- value pairs where values are strings.
        Useful for testing server """
        r = self.getpost(relpath='/api/v2/mirror', timeout=timeout, payload=payload, method='POST')   
        rd = self._decode(r)
        return rd

    def server_config(self,  timeout =None):
        """ returns server config as a dictionary """
        return self._decode(self.getpost('/api/v2/server_config', timeout=timeout, method='GET'))
    def server_time(self,  timeout =None):
        """ returns server time as an isoformat string"""
        return self._decode(self.getpost('/api/v2/server_time',timeout=timeout, method='GET'))
    def nucleotides_excluded(self,  timeout =0.2):
        """ returns the nucleotides excluded (i.e. the mask) """
        return self._decode(self.getpost('/api/v2/nucleotides_excluded', timeout=timeout, method='GET'))
    def guids(self,  timeout =None):
        """ returns all guids in the server """
        return self._decode(self.getpost('/api/v2/guids', method='GET'))
    def annotations(self,  timeout =None):
        """ returns all guids and their annotations as a pandas dataframe"""
        retVal = self._decode(self.getpost('/api/v2/annotations', timeout=timeout, method='GET'))
        return pd.DataFrame.from_dict(retVal, orient='index')   
    def clustering(self,  timeout =None):
        """ return clustering pipelines available """
        return self._decode(self.getpost('/api/v2/clustering', timeout=timeout, method='GET'))
    def guids_and_examination_times(self,  timeout =1):
        return self._decode(self.getpost('/api/v2/guids_and_examination_times', timeout=timeout, method='GET'))
    def guid_exists(self,  guid, timeout =None):
        """ returns True or False depending on whether the guid exists """
        if not isinstance(guid,str):  
            raise TypeError("guid {0} passed must be a string, not a {1}".format(guid,type(guid)))
        return self._decode(self.getpost('/api/v2/{0}/exists'.format(guid), timeout=timeout, method='GET'))
    def sequence(self,  guid, timeout =None):
        """ returns masked sequence of an existing guid """
        if not isinstance(guid,str):  
            raise TypeError("guid {0} passed must be a string, not a {1}".format(guid,type(guid)))
        res = self.getpost('/api/v2/{0}/sequence'.format(guid), timeout=timeout, method='GET')
        if res.status_code == 404:
            return None
        else:
            return self._decode(res)
    def change_id(self,  clustering_algorithm, timeout =None):
        """ returns the current change_id associated with clustering_algorithm """
        if not isinstance(clustering_algorithm, str):
            raise TypeError("clustering_algorithm must be str not {0}".format(type(clustering_algorithm)))      
        return self._decode(self.getpost('/api/v2/clustering/{0}/change_id'.format(clustering_algorithm), timeout=timeout, method='GET')) 
    def guids2clusters(self,  clustering_algorithm, after_change_id= None, timeout =None):
        """ returns a guid2cluster lookup """
        if not isinstance(clustering_algorithm, str):
            raise TypeError("clustering_algorithm must be str not {0}".format(type(clustering_algorithm)))
        if after_change_id is None:
            res = self._decode(self.getpost('/api/v2/clustering/{0}/guids2clusters'.format(clustering_algorithm), timeout=timeout, method='GET')) 
            return pd.DataFrame.from_records(res)
           
        elif isinstance(after_change_id, int):
            res = self._decode(self.getpost('/api/v2/clustering/{0}/guids2clusters/after_change_id/{1}'.format(clustering_algorithm, after_change_id), timeout=timeout, method='GET')) 
            return pd.DataFrame.from_records(res)
        else:
            raise TypeError("after must be None or an integer, not {0}".format(type(after_change_id)))
    def clusters(self,  clustering_algorithm, timeout =None):
        """ returns a clusters for a given clustering_algorithm """
        if not isinstance(clustering_algorithm, str):
            raise TypeError("clustering_algorithm must be str not {0}".format(type(clustering_algorithm)))
        
        res = self._decode(self.getpost('/api/v2/clustering/{0}/clusters'.format(clustering_algorithm), timeout=timeout, method='GET')) 
        return(pd.DataFrame.from_records(res['members']))
    def cluster_members(self,  clustering_algorithm, timeout =None):
        """ synonym for clusters """
        return self.clusters(clustering_algorithm, timeout=timeout)
    def cluster_summary(self,  clustering_algorithm, timeout =None):
        """ returns a clusters and counts in each cluster for a given clustering_algorithm """
        if not isinstance(clustering_algorithm, str):
            raise TypeError("clustering_algorithm must be str not {0}".format(type(clustering_algorithm)))
        
        res = self._decode(self.getpost('/api/v2/clustering/{0}/summary'.format(clustering_algorithm), timeout=timeout, method='GET')) 
        return(pd.DataFrame.from_records(res['summary']))
    
    def cluster_ids(self,  clustering_algorithm, timeout =None):
        """ returns a cluster_ids for a given clustering_algorithm """
        if not isinstance(clustering_algorithm, str):
            raise TypeError("clustering_algorithm must be str not {0}".format(type(clustering_algorithm)))
        
        res = self._decode(self.getpost('/api/v2/clustering/{0}/cluster_ids'.format(clustering_algorithm), timeout=timeout, method='GET')) 
        return res
    def network(self,  clustering_algorithm, cluster_id, timeout =None):
        """ returns a cytoscape compatible network """
        if not isinstance(clustering_algorithm, str):
            raise TypeError("clustering_algorithm must be str not {0}".format(type(clustering_algorithm)))
        if not isinstance(cluster_id, int):
            raise TypeError("cluster_id must be int not {0}".format(type(cluster_id)))
                
        res = self._decode(self.getpost('/api/v2/clustering/{0}/{1}/network'.format(clustering_algorithm, cluster_id), timeout=timeout, method='GET')) 
        return res
    def server_memory_usage(self,  nrows = 100, timeout =None):
        if not isinstance(nrows, int):
            raise TypeError("nrows must be integer not {0}".format(type(nrows)))      
        retVal= self._decode(self.getpost('/api/v2/server_memory_usage/{0}'.format(nrows), timeout=timeout, method='GET')) 
        return(pd.DataFrame.from_records(retVal))
    def guids_with_quality_over(self,  cutoff = 0, timeout =None):
        if not type(cutoff) in [int, float]:
            raise TypeError("cutoff must be float or int not {0}".format(type(cutoff)))      
        return self._decode(self.getpost('/api/v2/guids_with_quality_over/{0}'.format(cutoff), timeout=timeout, method='GET')) 
    def guid2neighbours(self,  guid, threshold, quality_cutoff=None, timeout =None):
        """ returns a guid2cluster lookup """
        if not isinstance(guid, str):
            raise TypeError("guid must be str not {0}".format(type(guid)))
        if not isinstance(threshold, int):
            raise TypeError("threshold must be int not {0}".format(type(threshold)))
        if quality_cutoff is None:
            quality_cutoff = 0  # no cutoff
        if isinstance(quality_cutoff, float) or isinstance(quality_cutoff, int):
            return self._decode(self.getpost('/api/v2/{0}/neighbours_within/{1}/with_quality_cutoff/{2}'.format(guid,threshold, quality_cutoff), timeout=timeout, method='GET')) 
        else:
            raise TypeError("unhandled: quality_cutoff must be None or float, not {0}".format(type(quality_cutoff)))
    def msa_cluster(self, clustering_algorithm, cluster_id, output_format='json', timeout=None):
        """ does MSA on a cluster """

        res= self.getpost('/api/v2/multiple_alignment_cluster/{0}/{1}/{2}'.format(clustering_algorithm, cluster_id, output_format), timeout=timeout, method='GET') 
        if output_format=='json':
            retDict = self._decode(res)
            return pd.DataFrame.from_dict(retDict, orient='index')          
        else:            
            return res.content
        
    def msa(self, guids, output_format = 'json-records', what='N', timeout=None):
        """ performs msa
        
        valid values for 'what', which determines how the p-values are computed, are
        M
        N
        N_or_M
        """
        guidstring = ';'.join(guids)
        payload = {'guids':guidstring,'output_format':output_format, 'what':what}
        res = self.post('/api/v2/multiple_alignment/guids', payload = payload, timeout= timeout)
        if output_format in ['json-records','json']:
            retList = self._decode(res)
            return pd.DataFrame.from_records(retList)          
        else:            
            return res.content
        
    def reset(self, timeout=30):
        """ resets the server to a state with no data """      
        return self.post('/api/v2/reset', payload = {}, timeout=timeout)
       
    def insert(self, guid, seq, timeout=None):
        """ inserts a sequence seq with guid 
        """
        
        # check input
        if not isinstance(guid,str):  
            raise TypeError("guid {0} passed must be a string, not a {1}".format(guid,type(guid)))
        if not isinstance(seq,str):  
            raise TypeError("sequence passed must be a string, not a {0}".format(type(seq)))
        return self.post('/api/v2/insert', payload = {'guid':guid,'seq':seq}, timeout=timeout)
       
    def read_fasta_file(self, fastafile):
        """ reads the content of a fasta file into memory.
        returns a dictionary {seqid:(first part of defline), seq:(nucleic acid), content:(entire file content)}.
        Supports both .gz and uncompressed files transparently.
        Does not support multi-fasta files.  Will raise an error if such are detected.
        """
        # first determine whether it is a .gz file or not; read into RAM.
        if fastafile.endswith('.gz'):
            # we decompress it on the fly.
            with gzip.open(fastafile,'r') as f:
                content = f.read().decode('utf-8')
        else:
            with open(fastafile, 'rt') as f:
                content = f.read()
        
        # use BioPython3 SeqIO library to read the file.      
        nFiles = 0 
        with io.StringIO(content) as f:
           for record in SeqIO.parse(f,'fasta'):
               nFiles +=1
               if nFiles > 1:       # that's a multifasta, and we don't support that
                    raise ValueError("Multifasta file is present in {0}.  Multifasta files are not supported".format(fastafile))
               else:
                    res = {'seq': str(record.seq), 'seqid':str(record.id), 'content':content }
                    return(res)
        raise IOError("no content parsed from result of length {0}".format(len(content)))


## note: unittests
# these are focused around testing the client, not comprehensively testing the function of the server.
# to do the latter, run the unittests in findNeighbour3-server.

class test_fn4_client_init1(unittest.TestCase):
    def runTest(self):
        """ initialise an fn4 object.  Will unittest against localhost mongodb server on 5020 """
        fn4c =fn4Client()         # expect success

class test_fn4_client_init2(unittest.TestCase):
    def runTest(self):
        """ initialise an fn4 client, against mongodb which doesn't exist """
        with self.assertRaises(Exception):
            fn4c =fn4Client(baseurl = 'http://example.com')         # expect failure as URL doesn't exist
        
class test_fn4_client_mirror(unittest.TestCase):
    def runTest(self):
        """ tests mirror function - which just sends back the POSTed payload """
        fn4c =fn4Client()         # expect failure as URL doesn't exist
        sent_payload = {'guid':'g1','seq':"ACTG"}
           
        retVal = fn4c.mirror(payload = sent_payload)
        self.assertEqual(retVal, sent_payload)
class test_fn4_client_guids(unittest.TestCase):
    def runTest(self):
        """ tests guids"""
        fn4c =fn4Client()         # expect failure as URL doesn't exist    
        retVal = fn4c.guids()
        self.assertEqual(type(retVal),list)
class test_fn4_client_annotations(unittest.TestCase):
    def runTest(self):
        """ tests guids"""
        fn4c =fn4Client()         # expect failure as URL doesn't exist    
        retVal = fn4c.annotations()
        self.assertEqual(type(retVal),pd.DataFrame)
class test_fn4_client_clustering(unittest.TestCase):
    def runTest(self):
        """ tests guids"""
        fn4c =fn4Client()         # expect failure as URL doesn't exist    
        retVal = fn4c.clustering()
        self.assertEqual(type(retVal),dict)
      
class test_fn4_client_server_memory_usage(unittest.TestCase):
    def runTest(self):
        """ tests guids"""
        fn4c =fn4Client()         # expect failure as URL doesn't exist    
        retVal = fn4c.server_memory_usage()
        self.assertEqual(type(retVal),pd.DataFrame)
class test_fn4_client_exists(unittest.TestCase):
    def runTest(self):
        """ tests guid/exists"""
        fn4c =fn4Client()         # expect failure as URL doesn't exist    
        retVal = fn4c.guid_exists("no")
        self.assertEqual(retVal,False)
class test_fn4_client_masked_sequence(unittest.TestCase):
    def runTest(self):
        """ tests guid/exists"""
        fn4c =fn4Client()         # expect failure as URL doesn't exist
        retVal = fn4c.sequence("no")
        self.assertEqual(retVal,None)
class test_fn4_client_config(unittest.TestCase):
    def runTest(self):
        """ tests server_config endpoint """
        fn4c =fn4Client()         # expect failure as URL doesn't exist
        retVal = fn4c.server_config()
        self.assertTrue('INPUTREF' in retVal.keys())
class test_fn4_client_server_time(unittest.TestCase):
    def runTest(self):
        """ tests server_time endpoint """
        fn4c =fn4Client()         # expect failure as URL doesn't exist
        retVal = fn4c.server_time()
        self.assertTrue('server_time' in retVal.keys())
class test_fn4_client_nucleotides_excluded(unittest.TestCase):
    def runTest(self):
        """ tests server_time endpoint """
        fn4c =fn4Client()         # expect failure as URL doesn't exist
        retVal = fn4c.nucleotides_excluded()
        self.assertTrue(type(retVal), list)
class test_fn4_client_guids_and_examination_times(unittest.TestCase):
    def runTest(self):
        """ tests server_time endpoint """
        fn4c =fn4Client()         # expect failure as URL doesn't exist
        retVal = fn4c.guids_and_examination_times()
        self.assertTrue(type(retVal), list)
class test_fn4_client_insert_fasta_1(unittest.TestCase):
    def runTest(self):
        """ initialise a gapi object """
        fn4c =fn4Client()         # expect success
        res = fn4c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t1.fasta" ))
        
        self.assertEqual(res['content'][0], '>')
        self.assertEqual(res['seqid'],'t1')
        self.assertEqual(res['seq'][0:5],'NNNNN')
        self.assertEqual(len(res['seq']), 4411532)
class test_fn4_client_insert_fasta_2(unittest.TestCase):
    def runTest(self):
        
        fn4c =fn4Client()         # expect success
        res = fn4c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t1.fasta" ))
        seq1 = res['seq']
        res = fn4c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t2.fasta" ))
        seq2 = res['seq']

        uuid1 = uuid.uuid1().hex
        uuid2 = uuid.uuid1().hex
        
        fn4c.insert(uuid1, seq1)
        fn4c.insert(uuid2, seq2)
        
        self.assertTrue(uuid1 in fn4c.guids())
        self.assertTrue(uuid2 in fn4c.guids())
        
        res1 = fn4c.guids_with_quality_over(0.25)
        res2 = fn4c.guids_with_quality_over(1.1)
        
        self.assertTrue(len(res1)>0)
        self.assertEqual(len(res2),0)
        
class test_fn4_client_change_id(unittest.TestCase):
    def runTest(self):
        
        fn4c =fn4Client()         # expect success
        res = fn4c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t1.fasta" ))
        seq1 = res['seq']
        res = fn4c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t2.fasta" ))
        seq2 = res['seq']

        uuid1 = uuid.uuid1().hex
        uuid2 = uuid.uuid1().hex
        
        fn4c.insert(uuid1, seq1)
        
        clustering = fn4c.clustering()
        c1 = fn4c.change_id(clustering['algorithms'][0])
        fn4c.insert(uuid2, seq2)
        c2 = fn4c.change_id(clustering['algorithms'][0])
               
        self.assertTrue(uuid1 in fn4c.guids())
        self.assertTrue(uuid2 in fn4c.guids())

        self.assertTrue(c1['change_id']==c2['change_id'])       # only updates on clustering
class test_fn4_client_msa(unittest.TestCase):
    def runTest(self):
        
        fn4c =fn4Client()         # expect success
        res = fn4c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t1.fasta" ))
        seq1 = res['seq']
        res = fn4c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t2.fasta" ))
        seq2 = res['seq']

        uuid1 = uuid.uuid1().hex
        uuid2 = uuid.uuid1().hex
        uuid3 = uuid.uuid1().hex
        fn4c.insert(uuid1, seq1)
        fn4c.insert(uuid2, seq2)
        fn4c.insert(uuid3, seq1)
        res = fn4c.msa([uuid1,uuid2,uuid3])
        
        self.assertTrue(isinstance(res, pd.DataFrame))
  
class test_fn4_client_guids2clusters(unittest.TestCase):
    
    """ tests various endpoints to do with clustering """
    def runTest(self):
        
        fn4c =fn4Client()         # expect success
        res = fn4c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t1.fasta" ))
        seq1 = res['seq']
        res = fn4c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t2.fasta" ))
        seq2 = res['seq']

        uuid1 = uuid.uuid1().hex
        uuid2 = uuid.uuid1().hex
        uuid3 = uuid.uuid1().hex
        fn4c.insert(uuid1, seq1)
        
        clustering = fn4c.clustering()
        c1 = fn4c.change_id(clustering['algorithms'][0])
        fn4c.insert(uuid2, seq2)
        c2 = fn4c.change_id(clustering['algorithms'][0])            
        fn4c.insert(uuid3, seq2)
        c3 = fn4c.change_id(clustering['algorithms'][0])
       
        self.assertTrue(uuid1 in fn4c.guids())
        self.assertTrue(uuid2 in fn4c.guids())
        self.assertTrue(uuid3 in fn4c.guids())


        # recover clustering
        res1 = fn4c.guids2clusters(clustering['algorithms'][0])
        res2 = fn4c.guids2clusters(clustering['algorithms'][0], after_change_id = c2['change_id'])

        self.assertTrue(isinstance(res1, pd.DataFrame))
        self.assertTrue(isinstance(res2, pd.DataFrame))
        
        # find clusters
        cluster_ids= set()
        for ix in res1.index:
            cluster_ids.add(res1.loc[ix, 'cluster_id'])
        
        # check clusters endpoint
        res3 = fn4c.clusters(clustering['algorithms'][0])
        self.assertTrue(isinstance(res3, pd.DataFrame))
        res4 = fn4c.cluster_members(clustering['algorithms'][0])
        self.assertTrue(isinstance(res4, pd.DataFrame))
        res5 = fn4c.cluster_summary(clustering['algorithms'][0])
        self.assertTrue(isinstance(res5, pd.DataFrame))
        res6 = fn4c.cluster_ids(clustering['algorithms'][0])
        self.assertTrue(isinstance(res6, list))
        self.assertTrue(set(res6)==cluster_ids)  # same results both ways
        
        # recover neighbours
        res = fn4c.guid2neighbours(guid= uuid1, threshold = 250)
        self.assertTrue(isinstance(res, dict))

        res = fn4c.guid2neighbours(guid= uuid1, threshold = 250, quality_cutoff = 0.7)
        self.assertTrue(isinstance(res, list))

        

class test_fn4_client_network(unittest.TestCase):
    """ tests network generation """
    def runTest(self):
        
        fn4c =fn4Client()         # expect success
        res = fn4c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t1.fasta" ))
        seq1 = res['seq']
        res = fn4c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t2.fasta" ))
        seq2 = res['seq']

        uuid1 = uuid.uuid1().hex
        uuid2 = uuid.uuid1().hex
        uuid3 = uuid.uuid1().hex
        fn4c.insert(uuid1, seq1)
        fn4c.insert(uuid2, seq2)
        fn4c.insert(uuid3, seq2)

        clustering = fn4c.clustering()
        self.assertTrue('algorithms' in clustering.keys())
        
        for algorithm in clustering['algorithms']:
            cluster_ids = fn4c.cluster_ids(algorithm)
            for cluster_id in cluster_ids:
                network = fn4c.network(algorithm, cluster_id)
                #print(network.keys())
                #print(algorithm, network['nNodes'], network['nEdges'])
                #print(network['elements'])
