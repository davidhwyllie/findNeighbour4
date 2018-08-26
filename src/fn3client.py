#!/usr/bin/env python3
""" A python3 client for findNeighbour3-server.

Provides a class which allows access to the following routes:

Testing the server
------------------
* returns the dictionary posted to the server. Can be used for testing network connectivity.   
/api/v2/mirror  (requires POST)

Describing server configuration
---------------------------------
* Describe server config.  Disabled if not in debug mode.  
/api/v2/server_config
* Return the [last *nrows* of the] server's log of internal memory usage  
/api/v2/server_memory_usage  
/api/v2/server_memory_usage/*nrows*
* Return server time  
/api/v2/server_time
* List nucleotides masked (ignored) by the server in distance computations  
/api/v2/nucleotides_excluded  

Insert into server   
-------------------
/api/v2/insert requires POST; see docs for details  

Search/describe all sequences in the server, each identified by a guid
-----------------------------------------------------------------------
* list all guids (sequence identifiers) in the server  
/api/v2/guids
* list all guids with quality (proportion of Ns in the sequence) over *cutoff*  
/api/v2/guids_with_quality_over/*cutoff*
* list all guids and their examination (i.e. insertion) time  
/api/v2/guids_and_examination_times
* describe annotations (e.g. quality) for all sequences  
/api/v2/annotations

Describe properties/neighbours of a single sequence, identified by a guid
-------------------------------------------------------------------------
* test whether it exists  
/api/v2/*guid*/exists

* specifies threshold, uses default quality cutoff and output format  
/api/v2/*guid*/neighbours_within/*threshold*

* specifying quality cutoff  
uses default output format, as specified in MAXN_PROP_DEFAULT in config file  
/api/v2/*guid*/neighbours_within/*threshold*/with_quality_cutoff/*cutoff*

* specify quality cutoff and output format  
/api/v2/*guid*/neighbours_within/*threshold*/in_format/*returned_format*

Recover masked sequences
------------------------
* recover masked sequences for *guid*  
/api/v2/*guid*/sequence

Mixtures
----------------------------
* compare a given sequence with a set of neighbours, estimating mixtures of recent origin
/api/v2/assess_mixed (requires POST)  

Multiple sequence alignments
----------------------------
* return multiple sequence alignment for an arbitrary set of sequences, either in json or html format.
/api/v2/multiple_alignment/guids   requires POST; see docs for details

* return multiple sequence alignments of members of cluster  
/api/v2/multiple_alignment/*clustering_algorithm*/*cluster_id*/*output_format*

Clustering
----------
* List the clustering settings operating  
/api/v2/clustering  
* Return the change_id, an incrementing integer which rises are changes are made to clusters  
/api/v2/clustering/*clustering_algorithm*/change_id  
* Return a guid -> cluster lookup
/api/v2/clustering/*clustering_algorithm*/guids2clusters
/api/v2/clustering/*clustering_algorithm*/guids2clusters/after_change_id/*change_id*

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
from Bio.Alphabet import generic_nucleotide

# used for unittesting
import uuid

class fn3Client():
    """ python3 API to the findNeighbour3-server REST endpoint.
    
        All endpoints are supported:
        
        Testing the server
        ------------------
        * returns the dictionary posted to the server. Can be used for testing network connectivity.   
        /api/v2/mirror  (requires POST)
        
        Describing server configuration
        ---------------------------------
        * Describe server config.  Disabled if not in debug mode.  
        /api/v2/server_config
        * Return the [last *nrows* of the] server's log of internal memory usage  
        /api/v2/server_memory_usage  
        /api/v2/server_memory_usage/*nrows*
        * Return server time  
        /api/v2/server_time
        * List nucleotides masked (ignored) by the server in distance computations  
        /api/v2/nucleotides_excluded  
        
        Insert into server   
        -------------------
        /api/v2/insert requires POST; see docs for details  
        
        Search/describe all sequences in the server, each identified by a guid
        -----------------------------------------------------------------------
        * list all guids (sequence identifiers) in the server  
        /api/v2/guids
        * list all guids with quality (proportion of Ns in the sequence) over *cutoff*  
        /api/v2/guids_with_quality_over/*cutoff*
        * list all guids and their examination (i.e. insertion) time  
        /api/v2/guids_and_examination_times
        * describe annotations (e.g. quality) for all sequences  
        /api/v2/annotations
        
        Describe properties/neighbours of a single sequence, identified by a guid
        -------------------------------------------------------------------------
        * test whether it exists  
        /api/v2/*guid*/exists
        
        * specifies threshold, uses default quality cutoff and output format  
        /api/v2/*guid*/neighbours_within/*threshold*
        
        * specifying quality cutoff  
        uses default output format, as specified in MAXN_PROP_DEFAULT in config file  
        /api/v2/*guid*/neighbours_within/*threshold*/with_quality_cutoff/*cutoff*
        
        * specify quality cutoff and output format  
        /api/v2/*guid*/neighbours_within/*threshold*/in_format/*returned_format*
        
        Recover masked sequences
        ------------------------
        * recover masked sequences for *guid*  
        /api/v2/*guid*/sequence
        
        Mixtures
        ----------------------------
        * compare a given sequence with a set of neighbours, estimating mixtures of recent origin
        /api/v2/assess_mixed (requires POST)  
        
        Multiple sequence alignments
        ----------------------------
        * return multiple sequence alignment for an arbitrary set of sequences, either in json or html format.
        /api/v2/multiple_alignment/guids   requires POST; see docs for details
        
        * return multiple sequence alignments of members of cluster  
        /api/v2/multiple_alignment/*clustering_algorithm*/*cluster_id*/*output_format*
        
        Clustering
        ----------
        * List the clustering settings operating  
        /api/v2/clustering  
        * Return the change_id, an incrementing integer which rises are changes are made to clusters  
        /api/v2/clustering/*clustering_algorithm*/change_id  
        * Return a guid -> cluster lookup
        /api/v2/clustering/*clustering_algorithm*/guids2clusters
        /api/v2/clustering/*clustering_algorithm*/guids2clusters/after_change_id/*change_id*
                   
        """

    def __init__(self,
                baseurl="http://127.0.0.1:5000",
                ):

        """ set up logging """
        logging.getLogger()
        
        # load config
        self.baseurl = baseurl
        
        # run connection check
        res =  self._decode(self.getpost('/api/v2/server_time', method='GET'))
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
        will raise errors if they are generated
        """
        
        session = requests.Session()
        session.trust_env = False
    
        sc = -1
        url = self.absurl(relpath)
        #print("** URL", url)
        if method == 'GET':
            response = session.get(url=url)
            sc=response.status_code               

        elif method == 'POST':
            #print("**POST PAYLOAD [in getpost]",json.dumps(payload))
            response = session.post(url=url,
                                        data = payload)
            sc = response.status_code             
                           
        else:
            raise NotImplementedError("either GET or POST is required as a method.  was passed {0}".format(method))
        
        session.close() 
        return response 

    def get(self, relpath,  timeout=None):
        """ issues GET against url.  returns a response object.
        raises errors if fails"""
        response = self.getpost(relpath=relpath, raiseError=raiseError, timeout=timeout, method='GET')
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
    def server_memory_usage(self,  nrows = 100, timeout =None):
        if not isinstance(nrows, int):
            raise TypeError("nrows must be integer not {0}".format(type(nrows)))      
        retVal= self._decode(self.getpost('/api/v2/server_memory_usage/{0}'.format(nrows), timeout=timeout, method='GET')) 
        return(pd.DataFrame.from_records(retVal))
    def guids_with_quality_over(self,  cutoff = 0.8, timeout =None):
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
            return self._decode(self.getpost('/api/v2/{0}/neighbours_within/{1}'.format(guid, threshold), timeout=timeout, method='GET')) 
        elif isinstance(quality_cutoff, float):
            return self._decode(self.getpost('/api/v2/{0}/neighbours_within/{1}/with_quality_cutoff/{2}'.format(guid,threshold, quality_cutoff), timeout=timeout, method='GET')) 
        else:
            raise TypeError("quality_cutoff must be None or float, not {0}".format(type(quality_cutoff)))
    def multiple_alignment_cluster(self, clustering_algorithm, cluster_id, output_format='json', timeout=None):
        """ does MSA on a cluster """

        res= self.getpost('/api/v2/multiple_alignment_cluster/{0}/{1}/{2}'.format(clustering_algorithm, cluster_id, output_format), timeout=timeout, method='GET') 
        if output_format=='json':
            retDict = self._decode(res)
            return pd.DataFrame.from_dict(retDict, orient='index')          
        else:            
            return res.content
        
        
        TODO = """
           Insert into server   
       Search/describe all sequences in the server, each identified by a guid
   
        Multiple sequence alignments
        ----------------------------
        * return multiple sequence alignment for an arbitrary set of sequences, either in json or html format.
        /api/v2/multiple_alignment/guids   requires POST; see docs for details
        
           """

    def assess_mixed(self,guid, neighbours, sample_size=10, timeout=None):
        """ assesses whether a guid is mixed relative to neighbours
        """
        neighbourstring = ';'.join(neighbours)
        payload = {'this_guid':guid,'related_guids':neighbourstring,'sample_size':sample_size}
        retVal = self.post('/api/v2/assess_mixed', payload = payload, timeout= timeout)
        res = self._decode(retVal)
        
        if res=={}:
            return None     # no information
        else:
            df= pd.DataFrame.from_dict(res, orient='index')
        return df
    
    def msa(self, guids, output_format = 'json', timeout=None):
        """ performs msa
        """
        guidstring = ';'.join(guids)
        payload = {'guids':guidstring,'output_format':output_format}
        res = self.post('/api/v2/multiple_alignment/guids', payload = payload, timeout= timeout)
        if output_format=='json':
            retDict = self._decode(res)
            return pd.DataFrame.from_dict(retDict, orient='index')          
        else:            
            return res.content
        
    def reset(self):
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
        
        # now use BioPython3 SeqIO library to read the file.      
        nFiles = 0 
        with io.StringIO(content) as f:
           for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):
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

class test_fn3_client_init1(unittest.TestCase):
    def runTest(self):
        """ initialise an fn3 object.  Will unittest against localhost mongodb server on 5000 """
        fn3c =fn3Client()         # expect success

class test_fn3_client_init2(unittest.TestCase):
    def runTest(self):
        """ initialise an fn3 client, against mongodb which doesn't exist """
        with self.assertRaises(Exception):
            fn3c =fn3Client(baseurl = 'http://example.com')         # expect failure as URL doesn't exist
        
class test_fn3_client_mirror(unittest.TestCase):
    def runTest(self):
        """ tests mirror function - which just sends back the POSTed payload """
        fn3c =fn3Client()         # expect failure as URL doesn't exist
        sent_payload = {'guid':'g1','seq':"ACTG"}
           
        retVal = fn3c.mirror(payload = sent_payload)
        self.assertEqual(retVal, sent_payload)
class test_fn3_client_guids(unittest.TestCase):
    def runTest(self):
        """ tests guids"""
        fn3c =fn3Client()         # expect failure as URL doesn't exist    
        retVal = fn3c.guids()
        self.assertEqual(type(retVal),list)
class test_fn3_client_annotations(unittest.TestCase):
    def runTest(self):
        """ tests guids"""
        fn3c =fn3Client()         # expect failure as URL doesn't exist    
        retVal = fn3c.annotations()
        self.assertEqual(type(retVal),pd.DataFrame)
class test_fn3_client_clustering(unittest.TestCase):
    def runTest(self):
        """ tests guids"""
        fn3c =fn3Client()         # expect failure as URL doesn't exist    
        retVal = fn3c.clustering()
        self.assertEqual(type(retVal),dict)
        print(retVal)       
class test_fn3_client_server_memory_usage(unittest.TestCase):
    def runTest(self):
        """ tests guids"""
        fn3c =fn3Client()         # expect failure as URL doesn't exist    
        retVal = fn3c.server_memory_usage()
        self.assertEqual(type(retVal),pd.DataFrame)
class test_fn3_client_exists(unittest.TestCase):
    def runTest(self):
        """ tests guid/exists"""
        fn3c =fn3Client()         # expect failure as URL doesn't exist    
        retVal = fn3c.guid_exists("no")
        self.assertEqual(retVal,False)
class test_fn3_client_masked_sequence(unittest.TestCase):
    def runTest(self):
        """ tests guid/exists"""
        fn3c =fn3Client()         # expect failure as URL doesn't exist
        retVal = fn3c.sequence("no")
        self.assertEqual(retVal,None)
class test_fn3_client_config(unittest.TestCase):
    def runTest(self):
        """ tests server_config endpoint """
        fn3c =fn3Client()         # expect failure as URL doesn't exist
        retVal = fn3c.server_config()
        self.assertTrue('INPUTREF' in retVal.keys())
class test_fn3_client_server_time(unittest.TestCase):
    def runTest(self):
        """ tests server_time endpoint """
        fn3c =fn3Client()         # expect failure as URL doesn't exist
        retVal = fn3c.server_time()
        self.assertTrue('server_time' in retVal.keys())
class test_fn3_client_nucleotides_excluded(unittest.TestCase):
    def runTest(self):
        """ tests server_time endpoint """
        fn3c =fn3Client()         # expect failure as URL doesn't exist
        retVal = fn3c.nucleotides_excluded()
        self.assertTrue(type(retVal), list)
class test_fn3_client_guids_and_examination_times(unittest.TestCase):
    def runTest(self):
        """ tests server_time endpoint """
        fn3c =fn3Client()         # expect failure as URL doesn't exist
        retVal = fn3c.guids_and_examination_times()
        self.assertTrue(type(retVal), list)
class test_fn3_client_insert_fasta_1(unittest.TestCase):
    def runTest(self):
        """ initialise a gapi object """
        fn3c =fn3Client()         # expect success
        res = fn3c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t1.fasta" ))
        
        self.assertEqual(res['content'][0], '>')
        self.assertEqual(res['seqid'],'t1')
        self.assertEqual(res['seq'][0:5],'NNNNN')
        self.assertEqual(len(res['seq']), 4411532)
class test_fn3_client_insert_fasta_2(unittest.TestCase):
    def runTest(self):
        
        fn3c =fn3Client()         # expect success
        res = fn3c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t1.fasta" ))
        seq1 = res['seq']
        res = fn3c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t2.fasta" ))
        seq2 = res['seq']

        uuid1 = uuid.uuid1().hex
        uuid2 = uuid.uuid1().hex
        
        fn3c.insert(uuid1, seq1)
        fn3c.insert(uuid2, seq2)
        
        self.assertTrue(uuid1 in fn3c.guids())
        self.assertTrue(uuid2 in fn3c.guids())
        
        res1 = fn3c.guids_with_quality_over(0.25)
        res2 = fn3c.guids_with_quality_over(1.1)
        
        self.assertTrue(len(res1)>0)
        self.assertEqual(len(res2),0)
        
class test_fn3_client_change_id(unittest.TestCase):
    def runTest(self):
        
        fn3c =fn3Client()         # expect success
        res = fn3c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t1.fasta" ))
        seq1 = res['seq']
        res = fn3c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t2.fasta" ))
        seq2 = res['seq']

        uuid1 = uuid.uuid1().hex
        uuid2 = uuid.uuid1().hex
        
        fn3c.insert(uuid1, seq1)
        
        clustering = fn3c.clustering()
        c1 = fn3c.change_id(clustering['algorithms'][0])
        fn3c.insert(uuid2, seq2)
        c2 = fn3c.change_id(clustering['algorithms'][0])
               
        self.assertTrue(uuid1 in fn3c.guids())
        self.assertTrue(uuid2 in fn3c.guids())

        self.assertTrue(c1['change_id']<c2['change_id'])
class test_fn3_client_assess_mixed(unittest.TestCase):
    def runTest(self):
        
        fn3c =fn3Client()         # expect success
        res = fn3c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t1.fasta" ))
        seq1 = res['seq']
        res = fn3c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t2.fasta" ))
        seq2 = res['seq']

        uuid1 = uuid.uuid1().hex
        uuid2 = uuid.uuid1().hex
        uuid3 = uuid.uuid1().hex
        fn3c.insert(uuid1, seq1)
        fn3c.insert(uuid2, seq2)
        fn3c.insert(uuid3, seq1)
        res = fn3c.assess_mixed(uuid1, [uuid2,uuid3])
        self.assertTrue(isinstance(res, pd.DataFrame))
class test_fn3_client_msa(unittest.TestCase):
    def runTest(self):
        
        fn3c =fn3Client()         # expect success
        res = fn3c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t1.fasta" ))
        seq1 = res['seq']
        res = fn3c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t2.fasta" ))
        seq2 = res['seq']

        uuid1 = uuid.uuid1().hex
        uuid2 = uuid.uuid1().hex
        uuid3 = uuid.uuid1().hex
        fn3c.insert(uuid1, seq1)
        fn3c.insert(uuid2, seq2)
        fn3c.insert(uuid3, seq1)
        res = fn3c.msa([uuid1,uuid2,uuid3])
        self.assertTrue(isinstance(res, pd.DataFrame))
  

class test_fn3_client_guids2clusters(unittest.TestCase):
    def runTest(self):
        
        fn3c =fn3Client()         # expect success
        res = fn3c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t1.fasta" ))
        seq1 = res['seq']
        res = fn3c.read_fasta_file(fastafile = os.path.join("..","testdata", "fasta", "t2.fasta" ))
        seq2 = res['seq']

        uuid1 = uuid.uuid1().hex
        uuid2 = uuid.uuid1().hex
        uuid3 = uuid.uuid1().hex
        fn3c.insert(uuid1, seq1)
        
        clustering = fn3c.clustering()
        c1 = fn3c.change_id(clustering['algorithms'][0])
        fn3c.insert(uuid2, seq2)
        c2 = fn3c.change_id(clustering['algorithms'][0])            
        fn3c.insert(uuid3, seq2)
        c3 = fn3c.change_id(clustering['algorithms'][0])
       
        self.assertTrue(uuid1 in fn3c.guids())
        self.assertTrue(uuid2 in fn3c.guids())
        self.assertTrue(uuid3 in fn3c.guids())

        self.assertTrue(c1['change_id']<c2['change_id'])

        # recover clustering
        res1 = fn3c.guids2clusters(clustering['algorithms'][0])
        res2 = fn3c.guids2clusters(clustering['algorithms'][0], after_change_id = c2['change_id'])

        self.assertTrue(isinstance(res1, pd.DataFrame))
        self.assertTrue(isinstance(res2, pd.DataFrame))
        # find clusters
        cluster_ids= set()
        for ix in res1.index:
            cluster_ids.add(res1.loc[ix, 'cluster_id'])
           
        # recover neighbours
        res = fn3c.guid2neighbours(guid= uuid1, threshold = 250)
        self.assertTrue(isinstance(res, list))

        res = fn3c.guid2neighbours(guid= uuid1, threshold = 250, quality_cutoff = 0.7)
        self.assertTrue(isinstance(res, list))

        # check multiple_alignment_cluster
        msa = fn3c.multiple_alignment_cluster(clustering['algorithms'][0],min(cluster_ids),'json')

        self.assertTrue(isinstance(msa, pd.DataFrame))
        msa = fn3c.multiple_alignment_cluster(clustering['algorithms'][0],min(cluster_ids),'fasta')
