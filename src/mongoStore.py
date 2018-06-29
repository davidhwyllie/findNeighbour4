#!/usr/bin/env python
""" FNhelper, a class which provides a storage layer for meta-data (but not SNP distances) about sequences
  linked to this are SQL alchemy table definition classes FN*.
 """
          
import os
import datetime
import hashlib
import uuid
import json
import pandas as pd
import logging
import pymongo

# used for unit testing only
import unittest
from NucleicAcid import NucleicAcid 


## TODO: add classes which persist
# compressed sequences
# snv

class fn3persistence():
        """ System for persisting results from  large numbers of sequences stored in FindNeighbour.
        Uses Mongodb.
        
        in the current schema there are four collections.
        -'refcompressedseq' contains the reference compressed sequences.
        -'guid2meta', contains guid -> metadata
        -'guid2neighbours', contains links between guids, including snv
        -'config', containing configuration information
        
        """
        
        def refcompressedsequence_store(self, guid, object):
            """ stores object into refcompressedseq collection.
            It is assumed object is a dictionary"""
            return self._store('refcompressedseq',guid, object)
        def refcompressedsequence_read(self, guid, object):
            """ loads object from refcompressedseq collection.
            It is assumed object is a dictionary"""
            return self._load('refcompressedseq',key)
        def refcompressedsequence_guids(self):
            """ loads guids from refcompressedseq collection.
            """
            return self._load_ids('refcompressedseq')
        def config_store(self, key, object):
            """ stores object into config collection
            It is assumed object is a dictionary"""

            return self._store('config',key, object)
        def config_read(self, key):
            """ loads object from refcompressedseq.
                It is assumed object is a dictionary"""
            return self._load('config',key)
        def _store(self, collection, key, object):
            """ stores key:object in collection. It is assumed object is a dictionary.  Updates if appropriate."""
            if not isinstance(object, dict):
                raise TypeError(" .. TypeError: object passed must be a dictionary".format(object))
            
            object['_id'] = key

            res = self.db[collection].replace_one({'_id':key}, object, upsert=True)
            if not res.acknowledged is True:
                raise IOError("Mongo {0} did not acknowledge write of data: {1}".format(self.db, object))        
            return res
        def _load(self,collection, key):
            """ loads object from collection[key] """
            return self.db[collection].find_one({'_id':key})
        def _load_ids(self,collection):
            """ loads guids from collection """
            retVal = set()
            for item in self.db[collection].find({}):
                retVal.add(item['_id'])
            return(retVal)
                
        def __init__(self, connString, dbname = 'fn3', debug = 0, config_settings={}):
            """ Creates a connection to a MongoDb database.
            
            connString : the mongoDb connection string
            dbname: the name of the mongoDb database to use.
            if debug = 0, the database is opened or created.
            if debug is 1, any existing collections are deleted.
            config_settings: only used on db creation; optional dictionary to note items in the database's config collection.
            """
            
            self.logger = logging.getLogger()
            
            # client should trap for connection errors etc      
            self.client = pymongo.MongoClient(connString)
            self.dbname = dbname
            self.db = self.client[dbname]
            
            # delete any pre-existing data if we are in debug mode.
            collections = ['refcompressedseq','guid2meta','guid2neighbours','config']
            if debug == 1:
                self.logger.warning("Debug mode operational; deleting all data from collections.")
                for collection in collections:
                    self.db[collection].delete_many({})
        
        def first_run(self):
            """ if there is no config entry, it is a first-run situation """
            if self.db.config.find_one({'_id':'config'}) is None:
                return True
            else:
                return False
        def __del__(self):
            """ closes any session """
            self.closedown() 

        def closedown(self):
            """ closes any session """
            try:
                self.client.close() 
            except:
                pass

                
        def guid_annotate(self, guid, nameSpace, annotDict):
            """ adds multiple annotations of guid from a dictionary;
            all annotations go into a namespace.
            creates the record if it does not exist"""
            
            # check whethere there is an existing metadata object for this
            metadataObj = self.db.guid2meta.find_one({'_id':guid})
            if metadataObj is None:
                # it doesn't exist.  we create a new one.
                metadataObj = {'_id':guid, 'sequence_meta':{nameSpace:annotDict}}

            if not 'sequence_meta' in metadataObj.keys():
                metadataObj['sequence_meta']={}
            if not nameSpace in metadataObj['sequence_meta'].keys():
                metadataObj['sequence_meta'][nameSpace] = {}
                metadataObj['sequence_meta'][nameSpace] = {**metadataObj['sequence_meta'][nameSpace], **annotDict}
                
            res = self.db.guid2meta.replace_one({'_id':guid}, metadataObj, upsert=True)
            if not res.acknowledged is True:
                raise IOError("Mongo {0} did not acknowledge write of data: {1}".format(self.db, self.metadataObj))        
        def guids(self):
            """ returns all registered guids """
            retVal = [x['_id'] for x in self.db.guid2meta.find({}, {'_id':1})]
            return(set(retVal))
        def guid_exists(self, guid):
            """ checks the presence of a single guid """
            res = self.db.guid2meta.find_one({'_id':guid},{'sequence_meta':1})
            if res is None:
                return False
            else:
                return True
        def guid_quality_check(self,guid,cutoff):
         """ Checks whether the quality of one guid exceeds the cutoff.
         
         If the guid does not exist, returns None.
         If the guid does exist and has quality< cutoff, returns False.
         Otherwise, returns True.
         """
         
         # test input
         if not type(cutoff) in [float,int]:
                 raise TypeError ("Cutoff should be either floating point or integer, but it is %s" % type(cutoff))
         if not type(guid)==str:
                 raise TypeError ("The guid passed should be as string, not %s" % str(guid))

         # recover record, compare with quality
         res = self.db.guid2meta.find_one({'_id':guid},{'sequence_meta':1})
         if res is None:        # no entry for this guid
                 return None
         else:
            try:
                dnaq = res['sequence_meta']['DNAQuality']
            except KeyError:
                raise KeyError("DNA quality is not present in the sequence metadata {0}: {1}".format(guid, res))
            
            # check the DNA quality metric expected is present
            if not 'propACTG' in dnaq.keys():
                raise KeyError("propACTG is not present in DNAQuality namespace of guid {0}: {1}".format(guid, dnaq))
            
            # report whether it is larger or smaller than cutoff
            return dnaq['propACTG']>=cutoff
        def guid2item(self, guidList, namespace, tag):
            """ returns the item in namespace:tag for all guids in guidlist.
            To do this, a table scan is performed - indices are not used.
            If guidList is None, all items are returned.
            An error is raised if namespace and tag is not present in each record.   
            """
            
            retDict={}
            if guidList is None:
                results = self.db.guid2meta.find({},{'sequence_meta':1})
            else:
                results = self.db.guid2meta.find({'_id':{"$in":guidList}},{'sequence_meta':1})
            
            if results is None:        # nothing found
                    return None
                
            for res in results:
               try:
                   namespace_content = res['sequence_meta'][namespace]
               except KeyError:
                   raise KeyError("{2} is not present in the sequence metadata {0}: {1}".format(guid, res, namespace))
               
               # check the DNA quality metric expected is present
               if not tag in namespace_content.keys():
                   raise KeyError("{2} is not present in {3} namespace of guid {0}: {1}".format(guid, namespace_content, tag, namespace))
               
               # return property
               retDict[res['_id']] = namespace_content[tag]
            return(retDict)
        def guid2ExaminationDateTime(self, guidList=None):
            """ returns quality scores for all guids in guidlist.  If guidList is None, all results are returned. """
            return self.guid2item(guidList,'DNAQuality','examinationDate')
        def guid2quality(self, guidList=None):
            """ returns quality scores for all guids in guidlist (if guidList is None)"""
            return self.guid2item(guidList,'DNAQuality','propACTG') 
        def guid2propACTG_filtered(self, cutoff=0.85):
            """ recover guids which have good quality, > cutoff.
            These are in the majority, so we run a table scan to find these.
            """
            allresults = self.guid2quality(None)        # get all results
            retDict = {}
            for guid in allresults.keys():
                if allresults[guid]>=cutoff:
                    retDict[guid]=cutoff
            return retDict      # note: slightly different from previous api
        def guid2items(self, guidList, namespaces):
            """ returns all items in namespaces, which is a list, as a pandas dataframe.
            If namespaces is None, all namespaces are returned.
            If guidList is None, all items are returned.
            To do this, a table scan is performed - indices are not used.
            """
            
            retDict={}
            if guidList is None:
                results = self.db.guid2meta.find({},{'sequence_meta':1})
            else:
                results = self.db.guid2meta.find({'_id':{"$in":guidList}},{'sequence_meta':1})
            
            if results is None:        # nothing found
                    return None
    
            for res in results:
               row = {}
               sought_namespaces = set(res['sequence_meta'].keys())
               if namespaces is not None:       # we only want a subset
                   sought_namespaces = sought_namespaces.intersection(namespaces)   # what we want, intersect what we've got
                   
               for sought_namespace in sought_namespaces:
                    for tag in res['sequence_meta'][sought_namespace].keys():
                        col_name = "{0}:{1}".format(sought_namespace, tag)
                        row[col_name] = res['sequence_meta'][sought_namespace][tag]
               retDict[res['_id']]=row
                       
            return(retDict)       
        def guid_annotations(self):
            return self.guid2items(None,None)
        
## persistence unit tests
UNITTEST_MONGOCONN = "mongodb+srv://findNeighbour3_daemon:vuhCWDaUUIfUSh3Kv6Br@findneighbour3-test-lvnyq.mongodb.net/admin"

class Test_SeqMeta_version(unittest.TestCase):
    """ tests version of library.  only tested with > v3.0"""   
    def runTest(self): 
        self.assertTrue(pymongo.__version__>='3.0')
        
class Test_SeqMeta_guid_annotate_1(unittest.TestCase):
    """ tests insert of new data item""" 
    def runTest(self): 
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 1)
        
        # test there is no 'test' item; insert, and confirm insert
        guid = 1
        namespace= 'ns'
        payload = {'one':1, 'two':2}
        res = p.guid_annotate(guid= guid, nameSpace=namespace, annotDict = payload)
        res = p.db.guid2meta.find_one({'_id':1})
        self.assertEqual(res['sequence_meta']['ns'], payload)

class Test_SeqMeta_guid_exists_1(unittest.TestCase):
    """ tests insert of new data item and existence check""" 
    def runTest(self): 
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 1)
        
        # test there is no 'test' item; insert, and confirm insert
        guid = 1
        namespace= 'ns'
        payload = {'one':1, 'two':2}
        res = p.guid_annotate(guid= guid, nameSpace=namespace, annotDict = payload)
        res = p.guid_exists(guid)
        self.assertEqual(res, True)
        res = p.guid_exists(-1)
        self.assertEqual(res, False)
        
class Test_SeqMeta_guid_annotate_2(unittest.TestCase):
    """ tests update of existing data item with same namespace""" 
    def runTest(self): 
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 1)
        
        # test there is no 'test' item; insert, and confirm insert
        guid = 1
        namespace= 'ns'
        payload1 = {'one':1, 'two':2}
        res = p.guid_annotate(guid= guid, nameSpace=namespace, annotDict = payload1)
        res = p.db.guid2meta.find_one({'_id':1})
        self.assertEqual(res['sequence_meta']['ns'], payload1)
        payload2 = {'one':1, 'two':2}
        res = p.guid_annotate(guid= guid, nameSpace=namespace, annotDict = payload2)
        res = p.db.guid2meta.find_one({'_id':1})
        self.assertEqual(res['sequence_meta']['ns'], payload2)
        
class Test_SeqMeta_guid_annotate_3(unittest.TestCase):
    """ tests update of existing data item with different namespace""" 
    def runTest(self): 
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 1)
        
        # test there is no 'test' item; insert, and confirm insert
        guid = 1
        payload1 = {'one':1, 'two':2}
        
        res = p.guid_annotate(guid= guid, nameSpace='ns1', annotDict = payload1)
        res = p.db.guid2meta.find_one({'_id':1})
        self.assertEqual(res['sequence_meta']['ns1'], payload1)
        
        payload2 = {'one':1, 'two':2}
        res = p.guid_annotate(guid= guid, nameSpace='ns2', annotDict = payload2)
        res = p.db.guid2meta.find_one({'_id':1})
        
        payloads = {'ns1':payload1, 'ns2':payload2}
        self.assertEqual(res['sequence_meta'], payloads)
           
class Test_SeqMeta_init(unittest.TestCase):
    """ tests version of library.  only tested with > v3.0""" 
    def runTest(self): 
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 1)
        
        res = p.db.config.find_one({'_id':'config'})
        self.assertTrue(res is not None)
        
        # test there is no 'test' item
        res = p.db.guid2meta.find_one({'_id':'_startup'})
        self.assertEqual(res, None)
        startup = {'_id':'_startup',
                   'payload':"me"}
        res = p.db.guid2meta.insert_one(startup)
        res = p.db.guid2meta.find_one({'_id':'_startup'})
        self.assertEqual(res, startup)
        
        
class Test_SeqMeta_guids(unittest.TestCase):
    """ tests recovery of sequence guids""" 
    def runTest(self): 
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 1)
        
        # test there is no 'test' item
        startup = {'_id':1}
        res = p.db.guid2meta.insert_one(startup)
        startup = {'_id':2}
        res = p.db.guid2meta.insert_one(startup)
        res= p.guids()
        self.assertEqual(res, set([1,2]))


class Test_SeqMeta_Base(unittest.TestCase):
    """ sets up a connection for unit testing""" 
    def setUp(self): 
        self.t = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 1)
        
        res = self.t.db.config.find_one({'_id':'config'})
        self.assertTrue(res is not None)        # check it worked
        
class Test_SeqMeta_guid_quality_check_1(Test_SeqMeta_Base):
        def runTest(self):
               """ tests return of sequences and their qualities """
               # set up nucleic acid object
               na=NucleicAcid()
               na.examine('ACGTACGTNN')         # 20% bad

               self.t.guid_annotate(guid='g1',nameSpace='DNAQuality',annotDict=na.composition)
               na.examine('ACGTACNNNN')         # 40% bad
               self.t.guid_annotate(guid='g2',nameSpace='DNAQuality',annotDict=na.composition)
               na.examine('ACGTNNNNNN')         # 60% bad
               self.t.guid_annotate(guid='g3',nameSpace='DNAQuality',annotDict=na.composition)
               
               r1=self.t.guid_quality_check('g1',0.80)             # valid
               r2=self.t.guid_quality_check('g2',0.80)             # invalid  
               r3=self.t.guid_quality_check('g3',0.80)             # invalid  
               r4=self.t.guid_quality_check('g4',0.80)             # invalid; does not exist.
               
               self.assertEqual(r1, True)
               self.assertEqual(r2, False)
               self.assertEqual(r3, False)               
               self.assertEqual(r4, None)


class Test_SeqMeta_guid2quality1(Test_SeqMeta_Base):
        def runTest(self):
               """ tests return of sequences and their qualities """
               # set up nucleic acid object
               na=NucleicAcid()
               na.examine('ACGTACGTNN')         # 20% bad

               self.t.guid_annotate(guid='g1',nameSpace='DNAQuality',annotDict=na.composition)
               na.examine('ACGTACNNNN')         # 40% bad
               self.t.guid_annotate(guid='g2',nameSpace='DNAQuality',annotDict=na.composition)
               na.examine('ACGTNNNNNN')         # 60% bad
               self.t.guid_annotate(guid='g3',nameSpace='DNAQuality',annotDict=na.composition)
               
               r1=self.t.guid_quality_check('g1',0.80)             # valid
               r2=self.t.guid_quality_check('g2',0.80)             # invalid  
               r3=self.t.guid_quality_check('g3',0.80)             # invalid  
               r4=self.t.guid_quality_check('g4',0.80)             # invalid; does not exist.
               
               self.assertEqual(r1, True)
               self.assertEqual(r2, False)
               self.assertEqual(r3, False)               
               self.assertEqual(r4, None)
               
               resDict=self.t.guid2quality(None)                          # restrict to nothing - return all
               self.assertTrue(resDict is not None)
               self.assertEqual(resDict['g1'],0.80)                
               self.assertEqual(resDict['g2'],0.60)                
               self.assertEqual(resDict['g3'],0.40)

class Test_SeqMeta_guid2quality2(Test_SeqMeta_Base):
        def runTest(self):
               """ tests return of sequences and their qualities """
               # set up nucleic acid object
               na=NucleicAcid()
               na.examine('ACGTACGTNN')         # 20% bad

               self.t.guid_annotate(guid='g1',nameSpace='DNAQuality',annotDict=na.composition)
               na.examine('ACGTACNNNN')         # 40% bad
               self.t.guid_annotate(guid='g2',nameSpace='DNAQuality',annotDict=na.composition)
               na.examine('ACGTNNNNNN')         # 60% bad
               self.t.guid_annotate(guid='g3',nameSpace='DNAQuality',annotDict=na.composition)
               
               r1=self.t.guid_quality_check('g1',0.80)             # valid
               r2=self.t.guid_quality_check('g2',0.80)             # invalid  
               r3=self.t.guid_quality_check('g3',0.80)             # invalid
               
               self.assertEqual(r1, True)
               self.assertEqual(r2, False)
               self.assertEqual(r3, False)                                # check the db insert works
               
               resDict=self.t.guid2quality(['g1','g2','g3'])
               self.assertTrue(resDict is not None)
               self.assertEqual(resDict['g1'],0.80)                
               self.assertEqual(resDict['g2'],0.60)                
               self.assertEqual(resDict['g3'],0.40)
               
class Test_SeqMeta_Base1(unittest.TestCase):
        """ initialise FN persistence and adds data """     
        def setUp(self):
                self.t = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 1)
     
                dna=NucleicAcid()

                # add some sequences
                seqs={'guid1':'ACGT','guid2':'NACT', 'guid3':'TTTT', 'guid4':'NNNN'}
                for guid in seqs.keys():
                        seq=seqs[guid]
                        dna.examine(seq)
                        self.t.guid_annotate(guid=guid, nameSpace='DNAQuality',annotDict=dna.composition)
                        
class Test_SeqMeta_guid2ExaminationDateTime(Test_SeqMeta_Base1):        
         """ recovering guids and examination times; """
         def runTest(self):
                 res = self.t.guid2ExaminationDateTime()
                 expected=4
                 self.assertEqual(len(res.keys()),expected)
                 
class Test_SeqMeta_propACTG_filteredSequenceGuids(Test_SeqMeta_Base1):      
        """  recovered guids filtered by the propACTG criterion """
        def runTest(self):        
                n=0
                for guid in self.t.guid2propACTG_filtered(cutoff=0.85):
                        n+=1
                expected=2
                self.assertEqual(n,expected)

class Test_SeqMeta_allAnnotations(Test_SeqMeta_Base1):
        """ tests recovery of all annoations """
        def runTest(self):
                df = self.t.guid_annotations()     
                self.assertEqual(len(df.keys()),4)
  