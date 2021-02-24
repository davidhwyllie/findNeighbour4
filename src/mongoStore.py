#!/usr/bin/env python
""" fnPersistence, a class which provides a storage layer for meta-data and snv distances in mongodb """
          
import os
import datetime
import hashlib
import uuid
import json
import pandas as pd
import logging
import pymongo
from bson.objectid import ObjectId
import gridfs
import pickle
import psutil
import copy
import io
import numpy as np

# used for unit testing only
import unittest
from NucleicAcid import NucleicAcid 
import time

class NPEncoder(json.JSONEncoder):
    """ encodes Numpy types as jsonisable equivalents """
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NpEncoder, self).default(obj)


class fn3persistence():
        """ System for persisting results from  large numbers of sequences stored in FindNeighbour 3+.
        Uses Mongodb.
        
        in the current schema there are the following collections:
        -'config', containing configuration information
        
        -'refcompressedseq' contains reference compressed sequences. note that this is a gridfs 'filesystem'.  Keys are guids.

        -'clusters' contains a graph of inter-guid links. note that this is a gridfs 'filesystem'.  Keys are names of clustering algorithms.
        
        -'guid2meta', contains guid -> metadata
        
        -'guid2neighbour', contains links between guids, including snv
            Here, individuals documents are identified by mongo assigned unique ids.
            Each document contains three keys:
            {'guid':'a1234', 'rstat':'s', 'neighbours':{}}
        
            Up to max_neighbours_per_document neighbours can be stored per document.
        
        *max_neighbours_per_document* should be less than 5,000, because there is a max. document size in mongodb.
        In debug mode, it is automatically set to 3.
        if max_neighbours_per_document exist in the document, 'rstat' is set to 'f' (full).
        If there is a single item only, rstat is set to 's' (single); if there are multiple items, it is set to 'm'.
        
        Indices exist on (i) guid - allowing you to find all the documents contains guid X's neighbours and
                         (ii) guid/rstat combination- allowing one to find guid X's most recent document, useful for addition.
                         
        This class provides methods to access these four entities.
        
        NOTE:  regarding sharding, the most important collection is guid2neighbour.
        A hashed sharding based on guid should work well when ensuring database scalability.
        
        """
        
        # code handling startup and shutdown.
        def __init__(self,
                     connString,
                     dbname = 'fn3_unittesting',
                     debug=0,
                     config_settings={},
                     max_neighbours_per_document=5000,
                     server_monitoring_min_interval_msec=0):
            """ Creates a connection to a MongoDb database.
            
            connString : the mongoDb connection string
            dbname: the name of the mongoDb database to use.
            if debug = 0 or 1, the database is opened or created.
            if debug = 2, any existing collections are deleted.
            config_settings: only used on db creation; optional dictionary to note items in the database's config collection.
            """
            
            self.logger = logging.getLogger()
            self.logger.setLevel(logging.INFO)
            logging.info("Created connection to mongodb db named {0}".format(dbname))
            
            # client calling mongostore should trap for connection errors etc 
            self.connString = connString     
            self.dbname = dbname
            self._connect()		# will raise ConnectionError if fails

            # can check what exists with connection.database_names()
            self.expected_collections = ['server_monitoring',
                                         'guid2meta','guid2neighbour',
                                         'config',
                                         'refcompressedseq.chunks','refcompressedseq.files',
                                         'clusters.chunks','clusters.files',
                                         'msa.chunks','msa.files']
            self.expected_clustering_collections = [ 'clusters.chunks','clusters.files',
                                         'msa.chunks','msa.files']
            self.max_neighbours_per_document = max_neighbours_per_document
            self.server_monitoring_min_interval_msec = server_monitoring_min_interval_msec
            self.previous_server_monitoring_data = {}
            self.previous_server_monitoring_time = None
            
            # delete any pre-existing data if we are in debug mode.
            if debug == 2:
                self.logger.warning("Debug mode operational [DEBUG={0}]; deleting all data from collections.".format(debug))
                self._delete_existing_clustering_data()
                self._delete_existing_data()

                self.max_neighbours_per_document =2             # used for unittests
            else:
                self.logger.info("Using stored data in mongostore")
                

            # create indices on guid2neighbours; note will do nothing if index already exists
            ix1 = pymongo.IndexModel([("guid",pymongo.ASCENDING),("rstat", pymongo.ASCENDING)], name='by_guid_full')
            ix2 = pymongo.IndexModel([("rstat", pymongo.ASCENDING)], name='by_rstat')
    
            self.db['guid2neighbour'].create_indexes([ix1, ix2])            
 
            # create indices on msa; note will do nothing if index already exists
            ix3 = pymongo.IndexModel([("filename",pymongo.ASCENDING),("uploadDate", pymongo.ASCENDING)], name='filename_date')
            self.db['msa.files'].create_indexes([ix3]) 
           
            # create indices on guid2meta, allowing recovery of valid and invalid specimens rapidly.
            ix4 = pymongo.IndexModel([("sequence_meta.DNAQuality.invalid", pymongo.ASCENDING)], name='guid_valid')

            # note: if additional metadata is added, such as sequence names etc which might be searched for, then we need to add additional indices here.
            
            self.db['guid2meta'].create_indexes([ix4]) 

            # create index on server_monitoring insert times
            ix5 = pymongo.IndexModel([('context|time|time_now', pymongo.ASCENDING)])
            self.db['server_monitoring'].create_indexes([ix5])

        def delete_server_monitoring_entries(self, before_seconds):
            """ deletes server monitoring entries more than before_seconds ago """
            now = datetime.datetime.now()
            earliest_allowed = now - datetime.timedelta(seconds = before_seconds)

            earliest_allowed_str = str(earliest_allowed.isoformat())
            self.db['server_monitoring'].delete_many( {  "context|time|time_now" : { "$lt" : earliest_allowed_str } } )

            
        def summarise_stored_items(self):
            """ counts how many sequences exist of various types """
            retVal = {}
            collections_present = self.db.list_collection_names()
            for this_collection in self.expected_collections:
                if this_collection in collections_present:
                    res = self.db.command('collstats', this_collection)
                    for relevant_metric in ['totalIndexSize','storageSize','count','avgObjSize']:
                         if relevant_metric in res.keys():
                              target_key = "dstats|{0}|{1}".format(this_collection.replace('.','-'), relevant_metric)
                              retVal[target_key] = res[relevant_metric]
            return(retVal)    
 
        def connect(self):
            """ test whether the database is connected, and if not, tries to connect.
            if the connection fails, raises pymongo.errors.ConnectionFailure """
            if not self.is_connected():
                  self._connect()

        def _connect(self):
            """ connect to the database """

            # try to close any existing session, if it exists
            self.closedown()

            # open new client
            self.client = pymongo.MongoClient(self.connString, retryWrites=True)
            self.db = self.client[self.dbname]

            # open gridfs systems
            self.rcs = gridfs.GridFS(self.db, collection='refcompressedseq')       
            self.clusters = gridfs.GridFS(self.db, collection='clusters')       
            self.monitor = gridfs.GridFS(self.db, collection='monitor')       
            self.msa = gridfs.GridFS(self.db, collection='msa') 

            # enable sharding at database level
            #self.client.admin.command('enableSharding', self.dbname)
    
        def is_connected(self):
            """ Tests whether db is connected cf
            http://api.mongodb.com/python/current/api/pymongo/mongo_client.html """
            try:
                # The ismaster command is cheap and does not require auth.
                self.client.admin.command('ismaster')
                # success
                return True
            except pymongo.errors.ConnectionFailure:
                return False

        def raise_error(self,token):
            """ raises a ZeroDivisionError, with token as the message.
            useful for unit tests of error logging """
            raise ZeroDivisionError(token)
    
        def _delete_existing_data(self):
            """ deletes existing data from the databases """
            for collection in self.expected_collections:
                self.db[collection].delete_many({})

        def _delete_existing_clustering_data(self):
            """ deletes any clustering data from the databases """
            for collection in self.expected_clustering_collections:
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

        # generic routines to handle insertion and read from standard mongodb stores
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
        
        def memory_usage(self):
                """ returns memory usage by current python3 process  
                Uses the psutil module, as the resource module is not available in windows.
                """       
                memdict = psutil.virtual_memory()._asdict()
                sm = {'server|mstat|'+k: v for k, v in memdict.items()}
                return(sm)

        # methods for the config collection
        def config_store(self, key, object):
            """ stores object into config collection
            It is assumed object is a dictionary"""

            return self._store('config',key, object)
        
        def config_read(self, key):
            """ loads object from config.
                It is assumed object is a dictionary"""

            return self._load('config',key)


        # methods for the server_monitoring
        def recent_server_monitoring(self, max_reported = 100, selection_field = None, selection_string = None):
            """ returns a list containing recent server monitoring, in reverse order (i.e. tail first).
                The _id field is an integer reflecting the order added.  Lowest numbers are most recent.
                
                Inputs
                max_reported - return this number of lines, at most.
                selection_field - if not None, will only return lines containing selection_string
                                  in the 'selection_field' key of the returned dictionary.
                selection_string -if selection_field is not None, only returns rows if
                                  selection_string is present in the 'selection_field' key of the
                                  monitoring element. If None, this constraint is ignored.
            """

            if not isinstance(max_reported, int):
                raise TypeError("limit must be an integer, but it is a {0}".format(type(max_reported)))
            if not max_reported>=0:
                raise ValueError("limit must be more than or equal to zero")

            if max_reported == 0:
                return []
        
            n= 0
            retVal = []
            
            if selection_field is None:
                    formerly_cursor = self.db['server_monitoring'].find({}).sort('_id', pymongo.DESCENDING)
            else:
                    formerly_cursor = self.db['server_monitoring'].find({selection_field:selection_string}).sort('_id', pymongo.DESCENDING)
                 
            for formerly in formerly_cursor:
                n+=1
                formerly['_id']=n
                retVal.append(formerly)

                if n>=max_reported:
                        break
            return(retVal)
        
        def server_monitoring_store(self, message = 'No message provided', what=None, guid=None, content={}):
            """ stores content, a dictionary, into the server monitoirng log"""
            now = dict(**content)
            if what is not None:
                now['content|activity|whatprocess']= what
            if guid is not None:
                now['content|activity|guid']= guid
            now['context|info|message'] = message
            current_time = datetime.datetime.now()
            now['context|time|time_now']=str(current_time.isoformat())
            now['context|time|time_boot']=datetime.datetime.fromtimestamp(psutil.boot_time()).strftime("%Y-%m-%d %H:%M:%S")  

            # should we write this data?  We have the option not to log all messages, to prevent the store getting very full.
            write_content = False       
            if self.previous_server_monitoring_time is None:
                write_content = True   # yes if this is the first record written.
            else:
                time_since_last_write = current_time - self.previous_server_monitoring_time  # yes if it's after the server_moni
                t= 1000*float(time_since_last_write.seconds)+float(time_since_last_write.microseconds)/1000
                if t >= self.server_monitoring_min_interval_msec:
                        write_content = True
                        
            if write_content:      

                self.db['server_monitoring'].insert_one(now)
                self.previous_server_monitoring_time = current_time
                self.previous_server_monitoring_data = now
                return True
            else:
                return False

        # methods for monitor, which store the contents of an html file
        # in a gridFS store.
        def monitor_store(self, monitoring_id, html):
                """ stores the monitor output string html.  Overwrites any prior object.
                 """
                self.monitor.delete(monitoring_id)
                with io.BytesIO(html.encode('utf-8')) as f:
                        id = self.monitor.put(f, _id=monitoring_id, filename=monitoring_id)
                        return id

        def monitor_read(self, monitoring_id):
                """ loads stored string (e.g. html object) from the monitor collection. """
                try:
                    res = self.monitor.get(monitoring_id)
                except gridfs.errors.NoFile:
                    return None

                if res is None:
                    return None
                else:
                    return res.read().decode('utf-8')
  
        def msa_store(self, msa_token, msa):
                """ stores the msa object msa under token msa_token. """
                
                if not isinstance(msa, dict):
                        raise TypeError("Can only store dictionary objects, not {0}".format(type(dict)))

                res = self.msa.find_one({'_id':msa_token})
                if res is None:                
                    json_repr = json.dumps(msa).encode('utf-8')
                    with io.BytesIO(json_repr) as f:
                            self.msa.put(f, _id=msa_token, filename=msa_token)

                    return msa_token

        def msa_read(self, msa_token):
                """ loads object from msa collection.
                It is assumed object is a dictionary"""
                
                res = self.msa.find_one({'_id':msa_token})
                if res is None:
                    return None
                json_repr = json.loads(res.read().decode('utf-8'))
                return json_repr

        def msa_delete(self, msa_token):
                """ deletes the msa with token msa_token
                """
                
                self.msa.delete(msa_token)

        def msa_stored_ids(self):
                """ returns a list of  msa tokens of all objects stored """
                return [stored_msa._id for stored_msa in self.msa.find({})]
 
        def msa_delete_unless_whitelisted(self, whitelist):
                """ deletes the msa unless the id is in whitelist
                """
                to_delete= set()
                for id in self.msa_stored_ids():
                    if not id in whitelist:
                        to_delete.add(id)
                for msa_token in to_delete:
                    self.msa.delete(msa_token)


        def cluster_store(self, clustering_key, obj):
                """ stores the clustering object obj.  retains previous version.  To clean these up, call cluster_delete_legacy.

                    obj: a dictionary to store
                    clustering_key: the name of the clustering, e.g. TBSNP12-graph

                    Returns: 
                    current cluster version

                    Note; does not replace previous version, but stores a new one. 

                    cf. warning in Mongo docs:
                    Do not use GridFS if you need to update the content of the entire file atomically. 
                    As an alternative you can store multiple versions of each file and specify the current version of the file in the metadata. 
                    You can update the metadata field that indicates “latest” status in an atomic update after uploading the new version of the file, 
                    and later remove previous versions if needed.
                 """
                

                if not isinstance(obj, dict):
                        raise TypeError("Can only store dictionary objects, not {0}".format(type(dict)))
                json_repr = json.dumps(obj, cls=NPEncoder).encode('utf-8')

                with io.BytesIO(json_repr) as f:
                        id = self.clusters.put(f, filename=clustering_key)
                        return id       # this is the current cluster version

        def cluster_read(self, clustering_key):
                """ loads object from clusters collection corresponding to the most recent version of 
                the clustering, saved with filename = 'clustering_key'.
                """
                
                cursor = self.clusters.find({'filename':clustering_key}).sort('uploadDate',-1).limit(1)
                for res in cursor:
                    json_repr = json.loads(res.read().decode('utf-8'))
                    return json_repr
                # nothing there
                return None

        def cluster_read_update(self, clustering_key, current_cluster_version):
                """ loads object from clusters collection corresponding to the most recent version
                    of the clustering, saved with filename = 'clustering_key'.
                    it will read only if the current version is different from current_cluster_version; other wise, it returns None
                    It is assumed object is a dictionary"""
                latest_version = self.cluster_latest_version(clustering_key)
                if latest_version == current_cluster_version:
                    # no update
                    return None
                else:
                    return self.cluster_read(clustering_key)
 
        def cluster_latest_version(self, clustering_key):
                """ returns id of latest version """
                cursor = self.clusters.find({'filename':clustering_key}).sort('uploadDate',-1).limit(1)
                for res in cursor:
                    return res._id
                return None

        def cluster_keys(self, clustering_name=None):
                """ lists  clustering keys beginning with clustering_name.  If clustering_name is none, all clustering keys are returned.
                    
                """
                
                cursor = self.clusters.find({})
                filenames = set()
                retVal=[]
                for res in cursor:     
                    filenames.add(res.filename)
                
                if clustering_name is not None:     # only report keys starting with clustering_name
                    retVal = [x for x in sorted(filenames) if x.startswith(clustering_name)]
                else:
                    retVal = list(sorted(filenames))
                return retVal

        def cluster_versions(self, clustering_key):
                """ lists ids and storage dates corresponding to versions of clustering identifed by clustering_key.
                    the newest version is first.
                """
                
                cursor = self.clusters.find({'filename':clustering_key}).sort('uploadDate',-1)
                retVal=[]
                for res in cursor:
                    
                    retVal.append(res._id)
                return retVal

        def cluster_delete_all(self, clustering_key):
                """ delete all clustering objects, including the latest version, stored under clustering_key
                """
                ids = self.cluster_versions(clustering_key)
                for this_id in ids:
                    self.clusters.delete(this_id)

        def cluster_delete_legacy_by_key(self, clustering_key):
                """ delete all clustering objects, except latest version, stored with key clustering_key
                """
                ids = self.cluster_versions(clustering_key)
                ids = ids[1:]
                for i,this_id in enumerate(ids):
                    logging.info("Removing historical data for {0} {1} / {2}".format(clustering_key, i, len(ids)))
                    self.clusters.delete(this_id)

        def cluster_delete_legacy(self, clustering_name):
                """ delete all clustering objects, except latest version, stored with  clustering_name
                """
                clustering_keys = self.cluster_keys(clustering_name=clustering_name)
                for clustering_key in clustering_keys:
                    self.cluster_delete_legacy_by_key(clustering_key)

        def refcompressedseq_store(self, guid, obj):
                """ stores the pickled object obj with guid guid.
                Issues an error FileExistsError
                if the guid already exists. """
                pickled_obj = pickle.dumps(obj, protocol=2)
                if guid in self.rcs.list():
                        raise FileExistsError("Attempting to overwrite {0}".format(guid))
                id = self.rcs.put(pickled_obj, _id=guid, filename=guid)

                # do a functional test to verify write
                recovered_obj = self.refcompressedsequence_read(guid)
                if not recovered_obj == obj:
                    raise IOError("Integrity check failed on reference compressed item write/read for {0}".format(guid))
                return id

        def refcompressedsequence_read(self, guid):
                """ loads object from refcompressedseq collection.
                It is assumed object is a dictionary"""
             
                res = self.rcs.find_one({'_id':guid})
                if res is None:
                    return None
                
                return pickle.loads(res.read())
        
        def refcompressedsequence_guids(self):
            """ loads guids from refcompressedseq collection.
            """
            
            return(set(self.rcs.list()))

        # methods for guid2meta        
        def guid_annotate(self, guid, nameSpace, annotDict):
            """ adds multiple annotations of guid from a dictionary;
            all annotations go into a namespace.
            creates the record if it does not exist"""
            
            # check whethere there is an existing metadata object for this
            
            metadataObj = self.db.guid2meta.find_one({'_id':guid})
            if metadataObj is None:
                # it doesn't exist.  we create a new one.
                metadataObj = {'_id':guid, 'sequence_meta':{nameSpace:annotDict}}

            if not 'sequence_meta' in metadataObj.keys():       # this is key is mandatory and is always present
                metadataObj['sequence_meta']={}

            # if the namespace does not exist as a subsidiary of sequence_meta, then we create it
            if not nameSpace in metadataObj['sequence_meta'].keys():
                metadataObj['sequence_meta'][nameSpace] = {}

            # we add any annotations to the existing data
            metadataObj['sequence_meta'][nameSpace] = {**metadataObj['sequence_meta'][nameSpace], **annotDict}
                
            res = self.db.guid2meta.replace_one({'_id':guid}, metadataObj, upsert=True)
            if not res.acknowledged is True:
                raise IOError("Mongo {0} did not acknowledge write of data: {1}".format(self.db, self.metadataObj))
        
        def guids(self):
            """ returns all registered guids """
            
            retVal = [x['_id'] for x in self.db.guid2meta.find({}, {'_id':1})]
            return(set(retVal))
 
        def _guids_selected_by_validity(self, validity):
            """ returns  registered guids, selected on their validity

                0 = guid is valid
                1 = guid is invalid

            """
            if not validity in [0,1]:
                raise ValueError("Validity must be 0 or 1, not {0}".format(validity))

            retVal = [x['_id'] for x in self.db.guid2meta.find({"sequence_meta.DNAQuality.invalid":validity}, {'_id':1})]
            return(set(retVal))
        
        def guids_valid(self):
            """ return all registered valid guids.  

                Validity is determined by the contents of the DNAQuality.invalid field, on which there is an index """
            return self._guids_selected_by_validity(0)
        def guids_invalid(self):
            """ return all invalid guids 

                Validity is determined by the contents of the DNAQuality.invalid field, on which there is an index """
            return self._guids_selected_by_validity(1)

        def guid_exists(self, guid):
            """ checks the presence of a single guid """
            
            res = self.db.guid2meta.find_one({'_id':guid},{'sequence_meta':1})
            if res is None:
                return False
            else:
                return True
        
        def guid_valid(self, guid):
            """ checks the validity of a single guid 

            Parameters:
            guid: the sequence identifier

            Returns
            -1   The guid does not exist
            0    The guid exists and the sequence is valid
            1    The guid exists and the sequence is invalid
            -2    The guid exists, but there is no DNAQuality.valid key"""
            
            res = self.db.guid2meta.find_one({'_id':guid},{'sequence_meta':1})
            if res is None:
                return -1
            else:

                try:
                    return int(res['sequence_meta']['DNAQuality']['invalid'])
                except KeyError:
                    return -2
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
            """ return all annotations of all guids """
            
            return self.guid2items(None,None)           # no restriction by namespace or by guid.
        def guid_annotation(self, guid):
            """ return all annotations of one guid """
            
            return self.guid2items([guid],None)           # restriction by guid.

                   
        def guid2neighbour_add_links(self,guid, targetguids):
                """ adds links between guid and targetguids
                
                guid: the 'source' guid for the matches eg 'guid1'
                targetguids: what is guid linked to, eg
                {
                        'guid2':{'dist':12},
                        'guid3':{'dist':2}
                }
                
                This stores links in the guid2neighbour collection;
                each stored document links one guid to one target.
                
                The function guid2neighbour_repack() reduces the number of documents
                required to store the same information.
                
                if annotation is not None, will additionally write an annotation dictionary into 
                """                
                                 
                # find guid2neighbour entry for guid.
                to_insert = []

                for targetguid in targetguids.keys():
                        payload = targetguids[targetguid]
                        payload1 = {guid:payload}
                        payload2 = {targetguid:payload}
                        
                        to_insert.append({'guid':guid, 'rstat':'s', 'neighbours': {targetguid:payload}})
                        to_insert.append({'guid':targetguid, 'rstat':'s', 'neighbours':{guid:payload}})

                # when complete, do update
                if len(to_insert)>0:
                    res = self.db.guid2neighbour.insert_many(to_insert, ordered=False)
                    if not res.acknowledged is True :
                       raise IOError("Mongo {0} did not acknowledge write of data: {1}".format(self.db, to_insert))
                # check there is a metadata object for the guid
                metadataObj = self.db.guid2meta.find_one({'_id':guid})
                if metadataObj is None:
                    # it doesn't exist.  we create a new one.
                    metadataObj = {'_id':guid, 'created':{'created_at':datetime.datetime.now().isoformat()}}
                    res = self.db.guid2meta.insert_one({'_id':guid}, metadataObj)
                    if not res.acknowledged is True :
                       raise IOError("Mongo {0} did not acknowledge write of data: {1}".format(self.db, metadataObj))

       
        def guid2neighbour_repack(self,guid):
                """ alters the mongodb representation of the links of guid.

                This stores links in the guid2neighbour collection;
                each stored document links one guid to one target.
                
                This function reduces the number of documents
                required to store the same information.
                
                Internally, the documents in guid2neighbour are of the form
                {'guid':'guid1', 'rstat':'s', 'neighbours':{'guid2':{'dist':12, ...}}} OR
                {'guid':'guid1', 'rstat':'m', 'neighbours':{'guid2':{'dist':12, ...}, 'guid3':{'dist':5, ...}} OR
                {'guid':'guid1', 'rstat':'f', 'neighbours':{'guid2':{'dist':12, ...}, 'guid3':{'dist':5, ...}} 
                  
                The last example occurs when the maximum number of neighbours permitted per record has been reached.
                
                This operation moves rstat='s' record's neighbours into rstat='m' type records.
                It guarantees that no guid-guid links will be lost, but more than one record of the same link may exist
                during processing.  This does not matter, because guid2neighbours() deduplicates.
                """                
                
                # determine whether there are any rstat 's' entries for this guid.
                # these include only one 'cell' of the distance matrix.
                
                s_ids=[]
                s_ids = [res["_id"] for res in self.db.guid2neighbour.find({'guid':guid, 'rstat':'s'})]

                if len(s_ids)==0:
                        return 1

                # determine whether there are any rstat 'm' entries for this guid.
                # these contain multiple cells on each row/column of the matrix.
                m_ids = [x['_id'] for x in self.db.guid2neighbour.find({'guid':guid, 'rstat':'m'})]
                                
                # iterate over the s_ids
                # setup
                current_m_id= None
                current_m = None
                processed_s_ids = []
                
                # iterate; for each 's' type record
                while len(s_ids)>0:
                        # read the record with this s_id
                        s_id = s_ids.pop()
                        processed_s_ids.append(s_id)
                        s = self.db.guid2neighbour.find_one({'_id':s_id})       # '_id':item
                        if s is not None:
                            # it is possible that another process has moved it into an m record already.                                         
              
                            # make sure we have a record to write into
                            if len(m_ids)==0 and current_m is None:
                                # create a record to write into
                                to_insert = {'guid':guid, 'rstat':'m', 'neighbours': {}}
                                current_m_id = self.db.guid2neighbour.insert_one(to_insert).inserted_id
                                if current_m_id is None:
                                        raise IOError("Failed to create a rstat m record")
                                current_m = self.db.guid2neighbour.find_one({'_id':current_m_id})
                
                            elif len(m_ids)>0 and current_m is None:
                                # we can use an existing record
                                current_m_id = m_ids.pop()
                                current_m = self.db.guid2neighbour.find_one({'_id':current_m_id})
                                
                            if current_m is None:
                                raise IOError("could not read or create record of id {0}".format(current_m_id))
           
                            # add the new neighbours to the existing neighbours    
                            for key in s['neighbours'].keys():
                                 current_m['neighbours'][key] = s['neighbours'][key]
                        
                            # if we've reached the maximum size permitted or there are none left to process
                            if len(current_m['neighbours'].keys()) >= self.max_neighbours_per_document or len(s_ids)==0:
                                if len(current_m['neighbours'].keys()) >= self.max_neighbours_per_document:                        
                                        current_m['rstat']= 'f'    # full
                                res = self.db.guid2neighbour.replace_one({'_id':current_m_id}, current_m)
                                if not res.acknowledged is True:
                                        raise IOError("Mongo {0} did not acknowledge write of data: {1}".format(self.db, current_m)) 
                                current_m = None
                                current_m_id=None
                # delete those processed single records
                self.db.guid2neighbour.delete_many({'_id':{'$in':processed_s_ids}})
                
        def guid2neighbours(self, guid, cutoff =20, returned_format=2):
                """ returns neighbours of guid with cutoff <=cutoff.
                    Returns links either as
                    
                    format 1 [[otherGuid, distance],[otherGuid2, distance2],...]
                    or as
                    format 2 [[otherGuid, distance, N_just1, N_just2, N_either],[],...]
                    or as
                    format 3 [otherGuid1, otherGuid2, otherGuid3]
                    or as
                    format 4 [{'guid':otherguid, 'snv':distance}]
                        
                        Internally, the documents in guid2neighbour are of the form
                        {'guid':'guid1', 'rstat':'s', 'neighbours':{'guid2':{'dist':12, ...}}} OR
                        {'guid':'guid1', 'rstat':'m', 'neighbours':{'guid2':{'dist':12, ...}, 'guid3':{'dist':5, ...}} OR
                        {'guid':'guid1', 'rstat':'f', 'neighbours':{'guid2':{'dist':12, ...}, 'guid3':{'dist':5, ...}} 
                          
                        However, irrespective of their internal representation, this function always returns
                        exactly one item for each link of 'guid'; duplicates are not possible.
                        The last example occurs when the maximum number of neighbours permitted per record has been reached.
                        """                
                
                retVal=[]
                formatting = {1:['dist'], 2:['dist','N_just1','N_just2','N_either'],3:[], 4:['dist']}
                desired_fields = formatting[returned_format]
                results=  self.db.guid2neighbour.find({'guid':guid})
                reported_already = set()
                for result in results:
                        
                        for otherGuid in result['neighbours'].keys():
                                if not otherGuid in reported_already:           # exclude duplicates
                                        if result['neighbours'][otherGuid]['dist']<=cutoff:        # if distance < cutoff
                                                available_fields = result['neighbours'][otherGuid].keys()
                                                reported_fields={}
                                                for desired_field in desired_fields:
                                                        try:
                                                                observed = result['neighbours'][otherGuid][desired_field]
                                                        except KeyError:                # doesn't exist
                                                                observed = None
                                                        reported_fields[desired_field]=observed
                                                         
                                                if returned_format == 1:
                                                    returned_data=[otherGuid, reported_fields['dist']]
                            
                                                elif returned_format == 2:
                                                    raise NotImplementedError("format 2 is no longer supported")
                                                elif returned_format == 3:
                                                    returned_data = otherGuid
                                                
                                                elif returned_format == 4:
                                                    returned_data={'guid':otherGuid, 'snv':reported_fields['dist']}
                            
                                                else:
                                                    raise ValueError("Unable to understand returned_format = {0}".format(returned_format))                          
                                                
                                                reported_already.add(otherGuid)
                                                retVal.append(returned_data)
                                        
                # recover the guids          
                return({'guid':guid, 'neighbours':retVal})

                
## persistence unit tests
UNITTEST_MONGOCONN = "mongodb://localhost"
class Test_Server_Monitoring_1(unittest.TestCase):
        """ adds server monitoring info"""
        def runTest(self):
                p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
                p.server_monitoring_store(message='one')
                
                res = p.recent_server_monitoring(100)

                self.assertEqual(len(res),1)
                self.assertTrue(isinstance(res,list))

class Test_Server_Monitoring_2(unittest.TestCase):
        """ adds server monitoring info"""
        def runTest(self):
                p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
                p.server_monitoring_store(message='one')
                p.server_monitoring_store(message='two')
                p.server_monitoring_store(message='three')
                
                   
                res = p.recent_server_monitoring(0)
                self.assertEqual(len(res),0)
                self.assertTrue(isinstance(res,list))

                res = p.recent_server_monitoring(1)
                self.assertEqual(len(res),1)
                self.assertTrue(isinstance(res,list))

                res = p.recent_server_monitoring(3)
                self.assertEqual(len(res),3)
                self.assertTrue(isinstance(res,list))

                res = p.recent_server_monitoring(5)
                self.assertEqual(len(res),3)
                self.assertTrue(isinstance(res,list))
                
                with self.assertRaises(ValueError):
                        res = p.recent_server_monitoring(-1)                        

                with self.assertRaises(TypeError):
                        res = p.recent_server_monitoring("thing")

class Test_Server_Monitoring_3(unittest.TestCase):
        """ checks whether server_monitoring_min_interval_msec control works"""
        def runTest(self):
                p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2, server_monitoring_min_interval_msec= 2000)
                retVal = p.server_monitoring_store(message='one')  # should insert
                self.assertEqual(retVal, True)
                res = p.recent_server_monitoring(100)
                self.assertEqual(len(res),1)
                self.assertTrue(isinstance(res,list))

                retVal = p.server_monitoring_store(message='two') # should not inserted
                self.assertEqual(retVal, False)
                res = p.recent_server_monitoring(100)
                self.assertEqual(len(res),1)
                self.assertTrue(isinstance(res,list))
                
                time.sleep(2)                   # seconds
                retVal = p.server_monitoring_store(message='three')  # should insert
                self.assertEqual(retVal, True)
                res = p.recent_server_monitoring(100)
                self.assertEqual(len(res),2)
                self.assertTrue(isinstance(res,list))

class Test_Server_Monitoring_4(unittest.TestCase):
        """ checks whether delete_server_monitoring_entries"""
        def runTest(self):
                p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2, server_monitoring_min_interval_msec= 0)
                retVal = p.server_monitoring_store(message='one')  # should insert
                self.assertEqual(retVal, True)
                res = p.recent_server_monitoring(100)
                self.assertEqual(len(res),1)
                self.assertTrue(isinstance(res,list))
                p.delete_server_monitoring_entries(1)
                res = p.recent_server_monitoring(100)
                self.assertEqual(len(res),1)
                self.assertTrue(isinstance(res,list))  #  should not have deleted
                
                time.sleep(2)                   # seconds

                p.delete_server_monitoring_entries(1)
                res = p.recent_server_monitoring(100)
                self.assertEqual(len(res),0)
                self.assertTrue(isinstance(res,list))               
class Test_SeqMeta_guid2neighbour_8(unittest.TestCase):
        """ tests guid2neighboursOf"""
        def runTest(self):
                p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
                p.guid2neighbour_add_links("srcguid",{'guid1':{'dist':12}, 'guid2':{'dist':0}, 'guid3':{'dist':3}, 'guid4':{'dist':4}, 'guid5':{'dist':5}})
                
                res1 = p.guid2neighbours('srcguid',returned_format=1)
                self.assertEqual(5, len(res1['neighbours']))
                with self.assertRaises(NotImplementedError):
                    res2 = p.guid2neighbours('srcguid',returned_format=2)

                res3 = p.guid2neighbours('srcguid',returned_format=3)
                self.assertEqual(5, len(res3['neighbours']))
                res4 = p.guid2neighbours('srcguid',returned_format=4)
                self.assertEqual(5, len(res4['neighbours']))
                
class Test_SeqMeta_guid2neighbour_7(unittest.TestCase):
        """ tests guid2neighboursOf"""
        def runTest(self):
                p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
                p.guid2neighbour_add_links("srcguid",{'guid1':{'dist':12}, 'guid2':{'dist':0}, 'guid3':{'dist':3}, 'guid4':{'dist':4}, 'guid5':{'dist':5}})
                with self.assertRaises(NotImplementedError):
                    res1 = p.guid2neighbours('srcguid', returned_format=2)
                res1 = p.guid2neighbours('srcguid', returned_format=1)
                self.assertEqual(5, len(res1['neighbours']))
                p.guid2neighbour_repack(guid='srcguid')
                res2 = p.guid2neighbours('srcguid', returned_format=1)
                self.assertEqual(5, len(res2['neighbours']))                

class Test_SeqMeta_guid2neighbour_6(unittest.TestCase):
        """ tests repack where repack spans multiple containers"""
        def runTest(self):
                p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
                res = p.guid2neighbour_add_links("srcguid",{'guid1':{'dist':12}, 'guid2':{'dist':0}, 'guid3':{'dist':3}, 'guid4':{'dist':4}, 'guid5':{'dist':5}})
                
                # check the insert worked
                res = p.db.guid2neighbour.count_documents({'guid':'guid1'})
                self.assertEqual(res, 1)
                res = p.db.guid2neighbour.count_documents({'guid':'srcguid'})
                self.assertEqual(res, 5)
                
                p.guid2neighbour_repack(guid='srcguid')
                
                # should compress the two entries for 'srcguid' into two
                res = p.db.guid2neighbour.count_documents({'guid':'guid1'})
                self.assertEqual(res, 1)
                res = p.db.guid2neighbour.count_documents({'guid':'srcguid'})
                self.assertEqual(res, 3)
                res = p.db.guid2neighbour.count_documents({'guid':'srcguid', 'rstat':'f'})
                self.assertEqual(res, 2)
                                
                # the src guid entry should contain two keys, 'guid1' and 'guid2' in its neighbours: section
                results= p.db.guid2neighbour.find({'guid':'srcguid'})
                observed = set()
                for result in results:
                        for item in result['neighbours'].keys():
                                observed.add(item)
                self.assertEqual(observed, set(['guid1','guid2', 'guid3','guid4', 'guid5']))
  
                res = p.guid2neighbour_add_links("guid6",{'srcguid':{'dist':12}})
                res = p.db.guid2neighbour.count_documents({'guid':'srcguid'})
                self.assertEqual(res, 4)
                res = p.db.guid2neighbour.count_documents({'guid':'srcguid', 'rstat':'f'})
                self.assertEqual(res, 2)
                res = p.db.guid2neighbour.count_documents({'guid':'srcguid', 'rstat':'s'})
                self.assertEqual(res, 1)             
                res = p.db.guid2neighbour.count_documents({'guid':'srcguid', 'rstat':'m'})
                self.assertEqual(res, 1)
                
                p.guid2neighbour_repack(guid='srcguid')
                res = p.db.guid2neighbour.count_documents({'guid':'srcguid'})
                self.assertEqual(res, 3)
                res = p.db.guid2neighbour.count_documents({'guid':'srcguid', 'rstat':'f'})
                self.assertEqual(res, 3)
                res = p.db.guid2neighbour.count_documents({'guid':'srcguid', 'rstat':'s'})
                self.assertEqual(res, 0)             
                res = p.db.guid2neighbour.count_documents({'guid':'srcguid', 'rstat':'m'})
                self.assertEqual(res, 0)

                results= p.db.guid2neighbour.find({'guid':'srcguid'})
                observed = set()
                for result in results:
                        for item in result['neighbours'].keys():
                                observed.add(item)
                self.assertEqual(observed, set(['guid1','guid2', 'guid3','guid4', 'guid5', 'guid6']))
                                            
class Test_SeqMeta_guid2neighbour_5(unittest.TestCase):
        """ tests repack """
        def runTest(self):
                p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
                res = p.guid2neighbour_add_links("srcguid",{'guid1':{'dist':12}, 'guid2':{'dist':0}})
                
                # check the insert worked
                res = p.db.guid2neighbour.count_documents({'guid':'guid1'})
                self.assertEqual(res, 1)
                res = p.db.guid2neighbour.count_documents({'guid':'srcguid'})
                self.assertEqual(res, 2)
                
                p.guid2neighbour_repack(guid='srcguid')
                
                # should compress the two entries for 'srcguid' into one.
                res = p.db.guid2neighbour.count_documents({'guid':'guid1'})
                self.assertEqual(res, 1)
                res = p.db.guid2neighbour.count_documents({'guid':'srcguid'})
                self.assertEqual(res, 1)
                
                # the src guid entry should contain two keys, 'guid1' and 'guid2' in its neighbours: section
                res= p.db.guid2neighbour.find_one({'guid':'srcguid'})
                self.assertEqual(set(res['neighbours'].keys()), set(['guid1','guid2']))
                
class Test_SeqMeta_guid2neighbour_4(unittest.TestCase):
        """ tests creation of a new guid2neighbour entry """
        def runTest(self):
                p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
                res = p.guid2neighbour_add_links("srcguid",{'guid1':{'dist':12}, 'guid2':{'dist':0}})
                res = p.db.guid2neighbour.count_documents({'guid':'guid1'})
                self.assertEqual(res, 1)
                res = p.db.guid2neighbour.count_documents({'guid':'srcguid'})
                self.assertEqual(res, 2)
                
class Test_SeqMeta_guid2neighbour_3(unittest.TestCase):
        """ tests creation of a new guid2neighbour entry """
        def runTest(self):
                p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
                res = p.guid2neighbour_add_links("srcguid",{'guid1':{'dist':12}})
                res = p.db.guid2neighbour.count_documents({'guid':'guid1'})
                self.assertEqual(res, 1)
                res = p.db.guid2neighbour.count_documents({'guid':'srcguid'})
                self.assertEqual(res, 1)
                
class Test_SeqMeta_guid2neighbour_2(unittest.TestCase):
        """ tests creation of a new guid2neighbour entry """
        def runTest(self):
                p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
                res = p.guid2neighbour_add_links("srcguid",{})
                res = p.db.guid2neighbour.count_documents({'guid':'guid1'})
                self.assertEqual(res, 0)
                
class Test_SeqMeta_version(unittest.TestCase):
    """ tests version of library.  only tested with > v3.0"""   
    def runTest(self): 
        self.assertTrue(pymongo.__version__>='3.0')
        
class Test_SeqMeta_file_store1(unittest.TestCase):
        """ tests storage of pickle files in database """
        def runTest(self):
                p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
                obj1 = {1,2,3}
                guid ="guid1"
                p.rcs.delete({'filename':guid})              # delete if present
                p.refcompressedseq_store(guid, obj1)
                res = p.rcs.find_one({'filename':guid}).read()
                obj2 = pickle.loads(res)
                self.assertEqual(obj1, obj2)
class Test_SeqMeta_file_store2(unittest.TestCase):
        """ tests storage of pickle files in database """
        def runTest(self):
                p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
                obj1 = {1,2,3}
                guid ="guid1"
                p.rcs.delete({'filename':guid})              # delete if present
                pickled_obj = pickle.dumps(obj1, protocol=2)
                p.refcompressedseq_store(guid, obj1)
                with self.assertRaises(FileExistsError):
                        p.refcompressedseq_store(guid, obj1)
class Test_SeqMeta_file_store3(unittest.TestCase):
        """ tests storage of pickle files in database """
        def runTest(self):
                p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
                obj1 = {1,2,3}
                guid ="guid1"
                p.rcs.delete({'filename':"guid1"})              # delete if present
                p.rcs.delete({'filename':"guid2"})              # delete if present
                res1 = p.refcompressedsequence_guids()
                pickled_obj = pickle.dumps(obj1, protocol=2)
                p.refcompressedseq_store(guid, pickled_obj)
                guid ="guid2"
                p.refcompressedseq_store(guid, pickled_obj)
                res2 = p.refcompressedsequence_guids()
                self.assertEqual(res2-res1,set(["guid1","guid2"]))
    
class Test_SeqMeta_guid_annotate_1(unittest.TestCase):
    """ tests insert of new data item""" 
    def runTest(self): 
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
        
        #  add the dictionary 'payload' to the namespace 'ns'
        guid = 1
        namespace= 'ns'
        payload = {'one':1, 'two':2}
        res = p.guid_annotate(guid= guid, nameSpace=namespace, annotDict = payload)
        res = p.db.guid2meta.find_one({'_id':1})
        self.assertEqual(res['sequence_meta']['ns'], payload)

class Test_SeqMeta_guid_annotate_2(unittest.TestCase):
    """ tests addition to existing data item""" 
    def runTest(self): 
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
        
        guid = 1

        # add the dictionary 'payload' to the namespace 'ns'
        namespace= 'ns'
        payload = {'one':1, 'two':2}
        res = p.guid_annotate(guid= guid, nameSpace=namespace, annotDict = payload)
        res = p.db.guid2meta.find_one({'_id':1})
        self.assertEqual(res['sequence_meta']['ns'], payload)

        # add the dictionary 'payload' to the namespace 'ns'
        namespace= 'ns'
        payload = {'three':3}
        res = p.guid_annotate(guid= guid, nameSpace=namespace, annotDict = payload)
        res = p.db.guid2meta.find_one({'_id':1})
        self.assertEqual(res['sequence_meta']['ns'], {'one':1, 'two':2, 'three':3})

class Test_SeqMeta_guid_annotate_3(unittest.TestCase):
    """ tests addition to existing data item""" 
    def runTest(self): 
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
        
        guid = 1

        # add the dictionary 'payload' to the namespace 'ns'
        namespace= 'ns'
        payload = {'one':1, 'two':2}
        res = p.guid_annotate(guid= guid, nameSpace=namespace, annotDict = payload)
        res = p.db.guid2meta.find_one({'_id':1})
        self.assertEqual(res['sequence_meta']['ns'], payload)

        # add the dictionary 'payload' to the namespace 'ns'
        namespace= 'ns'
        payload = {'two':3}
        res = p.guid_annotate(guid= guid, nameSpace=namespace, annotDict = payload)
        res = p.db.guid2meta.find_one({'_id':1})
        self.assertEqual(res['sequence_meta']['ns'], {'one':1, 'two':3})

class Test_SeqMeta_guid_exists_1(unittest.TestCase):
    """ tests insert of new data item and existence check""" 
    def runTest(self): 
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
        
        # test there is no 'test' item; insert, and confirm insert
        guid = 1
        namespace= 'ns'
        payload = {'one':1, 'two':2}
        res = p.guid_annotate(guid= guid, nameSpace=namespace, annotDict = payload)
        res = p.guid_exists(guid)
        self.assertEqual(res, True)
        res = p.guid_exists(-1)
        self.assertEqual(res, False)

class Test_SeqMeta_guid_valid_1(unittest.TestCase):
    """ tests insert of new data item and validity check""" 
    def runTest(self): 
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
        
        guid = "valid"
        namespace= 'DNAQuality'
        payload = {'invalid':0}
        res = p.guid_annotate(guid= guid, nameSpace=namespace, annotDict = payload)
        guid = "invalid"
        namespace= 'DNAQuality'
        payload = {'invalid':1}
        res = p.guid_annotate(guid= guid, nameSpace=namespace, annotDict = payload)
        guid = "missing"
        namespace= 'DNAQuality'
        payload = {'N':1}
        res = p.guid_annotate(guid= guid, nameSpace=namespace, annotDict = payload)
        
        res = p.guid_valid("valid")
        self.assertEqual(res, 0)
        res = p.guid_valid("invalid")
        self.assertEqual(res, 1)
        res = p.guid_valid("missing")
        self.assertEqual(res, -2)
        res = p.guid_valid("noguid")
        self.assertEqual(res, -1)

class Test_SeqMeta_guid_valid_2(unittest.TestCase):
    """ tests insert of new data item and validity check""" 
    def runTest(self): 
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
        
        guid = "valid1"
        namespace= 'DNAQuality'
        payload = {'invalid':0}
        res = p.guid_annotate(guid= guid, nameSpace=namespace, annotDict = payload)
        guid = "valid2"
        namespace= 'DNAQuality'
        payload = {'invalid':0}
        res = p.guid_annotate(guid= guid, nameSpace=namespace, annotDict = payload)
       
        guid = "invalid"
        namespace= 'DNAQuality'
        payload = {'invalid':1}
        res = p.guid_annotate(guid= guid, nameSpace=namespace, annotDict = payload)
        guid = "missing"
        namespace= 'DNAQuality'
        payload = {'N':1}
        res = p.guid_annotate(guid= guid, nameSpace=namespace, annotDict = payload)
        
        res = p.guids_valid()
        self.assertEqual(res, set(['valid1','valid2']))
        res = p.guids_invalid()
        self.assertEqual(res, set(['invalid']))
        
class Test_SeqMeta_guid_annotate_5(unittest.TestCase):
    """ tests update of existing data item with same namespace""" 
    def runTest(self): 
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
        
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
        
class Test_SeqMeta_guid_annotate_6(unittest.TestCase):
    """ tests update of existing data item with different namespace""" 
    def runTest(self): 
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
        
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
    """ tests database creation""" 
    def runTest(self): 
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
        self.assertTrue(p.first_run() == True)
        res = p.config_read('preComparer')
        self.assertEqual(res, None)

        p.config_store('config',{'item':1})
        self.assertTrue(p.first_run() == False)    
        res = p.config_read('config')
        self.assertEqual(res, {'_id':'config','item':1})

        p.config_store('preComparer',{'item':2})
        res = p.config_read('config')
        self.assertEqual(res, {'_id':'config', 'item':1})
        res = p.config_read('preComparer')
        self.assertEqual(res, {'_id':'preComparer','item':2})

        p.config_store('preComparer',{'item':3})
        res = p.config_read('config')
        self.assertEqual(res, {'_id':'config', 'item':1})
        res = p.config_read('preComparer')
        self.assertEqual(res, {'_id':'preComparer','item':3})


      
      
class Test_SeqMeta_guids(unittest.TestCase):
    """ tests recovery of sequence guids""" 
    def runTest(self): 
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
        
        startup = {'_id':1}
        res = p.db.guid2meta.insert_one(startup)
        startup = {'_id':2}
        res = p.db.guid2meta.insert_one(startup)
        res= p.guids()
        self.assertEqual(res, set([1,2]))


class Test_SeqMeta_Base(unittest.TestCase):
    """ sets up a connection for unit testing""" 
    def setUp(self): 
        self.t = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
        self.assertTrue(self.t.first_run() == True)
           
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
                self.t = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)
     
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

class Test_SeqMeta_oneAnnotation(Test_SeqMeta_Base1):
        """ tests recovery of one annotations """
        def runTest(self):
                df = self.t.guid_annotation('guid3')     
                self.assertEqual(len(df.keys()),1)
                df = self.t.guid_annotation('missing')     
                self.assertEqual(len(df.keys()),0)

class Test_Clusters(unittest.TestCase):
        """ tests saving and recovery of dictionaries to Clusters"""
        def runTest(self):
                p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)              
                
                self.assertIsNone(p.cluster_latest_version('cl1'))

                self.assertEqual(0, len(p.cluster_versions('cl1')))

                payload1 = {'one':1, 'two':2}
                x= p.cluster_store('cl1', payload1)

                payload1b = {'2_one':1, '2_two':2}
                y = p.cluster_store('cl2', payload1b)

                self.assertIsNotNone(p.cluster_latest_version('cl1'))
                clv = p.cluster_latest_version('cl1')
                self.assertEqual(x, clv)
                self.assertEqual(1, len(p.cluster_versions('cl1')))
                self.assertIsNone(p.cluster_read_update('cl1',clv))

                payload2 = p.cluster_read('cl1')   
                self.assertEqual(payload1, payload2)

                payload3 = {'one':10, 'two':20}
                p.cluster_store('cl1', payload3)       # this is now the latest version

                payload4 = {'2-one':10, '2-two':20}
                p.cluster_store('cl2', payload4)       # this is now the latest version

                self.assertEqual(p.cluster_keys(), ['cl1','cl2'])
                self.assertEqual(p.cluster_keys(clustering_name='cl1'), ['cl1'])

                self.assertEqual(2, len(p.cluster_versions('cl1')))
                self.assertNotEqual(clv,p.cluster_latest_version('cl1'))
                self.assertIsNotNone(p.cluster_read_update('cl1',clv))

                payload4 = p.cluster_read('cl1')   
                self.assertEqual(payload4, payload3)
               
                p.cluster_delete_legacy_by_key('cl1')

                self.assertEqual(1, len(p.cluster_versions('cl1')))

                p.cluster_delete_all('cl1')

                self.assertEqual(0, len(p.cluster_versions('cl1')))
                self.assertEqual(2, len(p.cluster_versions('cl2')))
                p.cluster_delete_legacy_by_key('cl2')
                self.assertEqual(1, len(p.cluster_versions('cl2')))
          
class Test_MSA(unittest.TestCase):
        """ tests saving and recovery of dictionaries to MSA"""
        def runTest(self):
                p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)              
                payload1 = {'one':1, 'two':2}
                p.msa_store(msa_token='msa1', msa = payload1)
                payload2 = p.msa_read(msa_token = 'msa1')   
                self.assertEqual(payload1, payload2)
                p.msa_delete(msa_token='msa1')
                payload3 = p.msa_read(msa_token = 'msa1')   
                self.assertIsNone(payload3)

                payload1 = {'one':1, 'two':2}
                p.msa_store(msa_token='msa1', msa = payload1)
                payload2 = {'one':3, 'two':4}
                p.msa_store(msa_token='msa2', msa = payload2)
                self.assertEqual(2, len(p.msa_stored_ids()))
                p.msa_delete_unless_whitelisted(whitelist=['msa1','msa2'])
                self.assertEqual(2, len(p.msa_stored_ids()))
                p.msa_delete_unless_whitelisted(whitelist=['msa1'])
                self.assertEqual(1, len(p.msa_stored_ids()))
class Test_Monitor(unittest.TestCase):
        """ tests saving and recovery of strings to monitor"""
        def runTest(self):
                p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)              
                payload1 = "line1"
                p.monitor_store('r1', payload1)
                payload2 = p.monitor_read('r1')   
                self.assertEqual(payload1, payload2)
                payload3 = p.monitor_read('nil')
                self.assertIsNone(payload3) 
class test_Raise_error(unittest.TestCase):
    """ tests raise_error"""
    def runTest(self):
        # generate compressed sequences
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)              
        with self.assertRaises(ZeroDivisionError):
            p.raise_error("token")

class Test_summarise_stored_items(unittest.TestCase):
        """ adds server monitoring info"""
        def runTest(self):
                p = fn3persistence(connString=UNITTEST_MONGOCONN, debug= 2)                
                res = p.summarise_stored_items()




