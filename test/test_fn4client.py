""" 
findNeighbour4 is a  server providing relatedness information for bacterial genomes via a Restful API.
See documentation for full details of its functionality.

There are unit tests for its python client fn4_client, which is
just a wrapper around calls to the findNeighbour4_server rest api.

note: unittests are focused around testing the client, not comprehensively testing the function of the server.
 to do the latter, run the unittests in test_server

To run these tests:

# starting a test RESTFUL server
pipenv run python3 findNeighbour4_server.py

# And then (e.g. in a different terminal) launching unit tests with
pipenv run python3 -m unittest test/test_fn4client.py
"""

import unittest
import os
import pandas as pd

import uuid

from fn4client import fn4Client

class test_fn4_client_init1(unittest.TestCase):
    def runTest(self):
        """ initialise an fn4 object.  Will unittest against localhost mongodb server on 5020 """
        fn4c =fn4Client()         # expect success
        self.assertIsInstance(fn4c, fn4Client)

class test_fn4_client_init2(unittest.TestCase):
    def runTest(self):
        """ initialise an fn4 client, against mongodb which doesn't exist """
        with self.assertRaises(Exception):
            fn4Client(baseurl = 'http://example.com')         # fails, does not exist

class test_fn4_client_mirror(unittest.TestCase):
    def runTest(self):
        """ tests mirror function - which just sends back the POSTed payload """
        fn4c =fn4Client()         # default
        sent_payload = {'guid':'g1','seq':"ACTG"}
           
        retVal = fn4c.mirror(payload = sent_payload)
        self.assertEqual(retVal, sent_payload)
class test_fn4_client_guids(unittest.TestCase):
    def runTest(self):
        """ tests guids"""
        fn4c =fn4Client()        
        retVal = fn4c.guids()
        self.assertEqual(type(retVal),list)
class test_fn4_client_annotations(unittest.TestCase):
    def runTest(self):
        """ tests guids"""
        fn4c =fn4Client()         # default    
        retVal = fn4c.annotations()
        self.assertEqual(type(retVal),pd.DataFrame)
class test_fn4_client_clustering(unittest.TestCase):
    def runTest(self):
        """ tests guids"""
        fn4c =fn4Client()         # default    
        retVal = fn4c.clustering()
        self.assertEqual(type(retVal),dict)
      
class test_fn4_client_server_memory_usage(unittest.TestCase):
    def runTest(self):
        """ tests guids"""
        fn4c =fn4Client()         # default    
        retVal = fn4c.server_memory_usage()

        self.assertEqual(type(retVal),pd.DataFrame)
class test_fn4_client_server_database_usage(unittest.TestCase):
    def runTest(self):
        """ tests guids"""
        fn4c =fn4Client()         # default    
        retVal = fn4c.server_database_usage()
        self.assertEqual(type(retVal),dict)
class test_fn4_client_exists(unittest.TestCase):
    def runTest(self):
        """ tests guid/exists"""
        fn4c =fn4Client()         # default    
        retVal = fn4c.guid_exists("no")
        self.assertEqual(retVal,False)
class test_fn4_client_masked_sequence(unittest.TestCase):
    def runTest(self):
        """ tests guid/exists"""
        fn4c =fn4Client()         # default
        retVal = fn4c.sequence("no")
        self.assertEqual(retVal,None)
class test_fn4_client_config(unittest.TestCase):
    def runTest(self):
        """ tests server_config endpoint """
        fn4c =fn4Client()         # default
        retVal = fn4c.server_config()
        self.assertTrue('INPUTREF' in retVal.keys())
class test_fn4_client_server_time(unittest.TestCase):
    def runTest(self):
        """ tests server_time endpoint """
        fn4c =fn4Client()         # default
        retVal = fn4c.server_time()
        self.assertTrue('server_time' in retVal.keys())
class test_fn4_client_nucleotides_excluded(unittest.TestCase):
    def runTest(self):
        """ tests server_time endpoint """
        fn4c =fn4Client()         # default
        retVal = fn4c.nucleotides_excluded()
        self.assertTrue(type(retVal), list)
class test_fn4_client_guids_and_examination_times(unittest.TestCase):
    def runTest(self):
        """ tests server_time endpoint """
        fn4c =fn4Client()         # default
        retVal = fn4c.guids_and_examination_times()
        self.assertTrue(type(retVal), list)
class test_fn4_client_insert_fasta_1(unittest.TestCase):
    def runTest(self):
        """ initialise a gapi object """
        fn4c =fn4Client()         # expect success
        res = fn4c.read_fasta_file(fastafile = os.path.join("testdata", "fasta", "t1.fasta" ))
        
        self.assertEqual(res['content'][0], '>')
        self.assertEqual(res['seqid'],'t1')
        self.assertEqual(res['seq'][0:5],'NNNNN')
        self.assertEqual(len(res['seq']), 4411532)
class test_fn4_client_insert_fasta_2(unittest.TestCase):
    def runTest(self):
        
        fn4c =fn4Client()         # expect success
        res = fn4c.read_fasta_file(fastafile = os.path.join("testdata", "fasta", "t1.fasta" ))
        seq1 = res['seq']
        res = fn4c.read_fasta_file(fastafile = os.path.join("testdata", "fasta", "t2.fasta" ))
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
        res = fn4c.read_fasta_file(fastafile = os.path.join("testdata", "fasta", "t1.fasta" ))
        seq1 = res['seq']
        res = fn4c.read_fasta_file(fastafile = os.path.join("testdata", "fasta", "t2.fasta" ))
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
        res = fn4c.read_fasta_file(fastafile = os.path.join("testdata", "fasta", "t1.fasta" ))
        seq1 = res['seq']
        res = fn4c.read_fasta_file(fastafile = os.path.join("testdata", "fasta", "t2.fasta" ))
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
        res = fn4c.read_fasta_file(fastafile = os.path.join( "testdata", "fasta", "t1.fasta" ))
        seq1 = res['seq']
        res = fn4c.read_fasta_file(fastafile = os.path.join( "testdata", "fasta", "t2.fasta" ))
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
       
        self.assertIsInstance(c1, dict)
        self.assertIsInstance(c2, dict)
        self.assertIsInstance(c3, dict)
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
    
        self.assertTrue(isinstance(res, list))

        res = fn4c.guid2neighbours(guid= uuid1, threshold = 250, quality_cutoff = 0.7)
        self.assertTrue(isinstance(res, list))

        

class test_fn4_client_network(unittest.TestCase):
    """ tests network generation """
    def runTest(self):
        
        fn4c =fn4Client()         # expect success
        res = fn4c.read_fasta_file(fastafile = os.path.join( "testdata", "fasta", "t1.fasta" ))
        seq1 = res['seq']
        res = fn4c.read_fasta_file(fastafile = os.path.join( "testdata", "fasta", "t2.fasta" ))
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
                self.assertIsInstance(network, dict)
                #print(network.keys())
                #print(algorithm, network['nNodes'], network['nEdges'])
                #print(network['elements'])
