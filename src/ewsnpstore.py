# ew interface


# necessary libraries
import os
import unittest
import sqlalchemy
import uuid
import logging
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Table, Column, Integer, String, Float, DateTime, MetaData, select, func, or_, and_
from sqlalchemy import create_engine, exc as excb
from sqlalchemy.orm import sessionmaker, aliased,  relationship, backref, exc
from sqlalchemy import ForeignKey
import datetime
import hashlib
import sys
from Bio import SeqIO
from gzip import GzipFile
import random
import glob
import json
import pandas as pd



## define classes
db_ewss=declarative_base() # classes mapping to persistent database inherit from this

class SnpAccessor():
    """ base class from which SNP sources inherit
    all sources have a

    SNPsource
    description
    matrixIsSparse (True/False) - whether a matrix will be returned
    nodes - a list of sequenceIds
    edges - a generator yielding an edge list
    accessedDateTime - time of access """
    
    def __init__(self):
        self.accessDataTime=None
        self.SNPsource=None
        self.SNPsourceUrl=None  # for Squirrelwalk etc
        self.description=None
        self.matrixIsSparse=None
    def nodes(self):
        for node in []:
            yield node
    def edges(self):
        for edge in []:
            yield edge


class ewReference(db_ewss):
    """ SQLalchemy class definition for the reference """
    __tablename__ = 'ewReference'
    ewRefId=Column(Integer, primary_key=True)
    ewReferenceName=Column(String(36), index=True, unique=True)
    ewGuid=relationship('ewGuid', cascade="all, delete",backref='reference')   


class ewGuid(db_ewss):       # the guids
    """ SQLalchemy class definition for a table including guids """
    __tablename__ = 'ewGuids'
    ewGuidId=Column(Integer, primary_key=True)
    ewRefId=Column(Integer, ForeignKey(ewReference.ewRefId), index=True)
    ewGuid=Column(String(36), index=True, unique=True)        
    examined=Column(Integer)
    examinationDateTime=Column(DateTime, index=True, nullable=True)


class ewEdge(db_ewss):
    """ SQLalchemy class definition for a table including edges between guids
    and Ns in the comparison"""
    __tablename__ = 'ewEdges'
    edgeId=Column(Integer, primary_key=True)
    ewGuidId1=Column(Integer,  index=True)
    ewGuidId2=Column(Integer,  index=True)
    dist=Column(Integer, index=True)
    N_just1=Column(Integer)
    N_just2=Column(Integer)
    N_either=Column(Integer)
 
    
class ElephantWalkDataSNPaccessor(SnpAccessor):
    """ a SNPsource using ElephantWalk data.
    
        This receives the results of querying for neighbours of individual SNPs and
        stores them in a database.  It does not itself interrogate EW; the class
        getEWSNP (module getEWSNP) does this.
        
        The output from getEWSNP.info can be passed to the
        .storeRelations method of this class.
        
        Example:
        ## TODO
        
        (end of example)
        """
        
    def __init__(self, db=None, engineName=None, ewdir=os.getcwd(), maxDistance=20, refName='R00000039'):
        """ constructor
        
        snpCutoff = the snpCutoff passed to elephantwalk server for edge detection
        db = instance of the orm
        engineName = connection string for the database.  examples:
    
            sqlite:///../sqlitedb/v7.db.  The explicit file path has been included.
            sqlite:///<<DEFAULT>>/v7.db.  module will replace <<DEFAULT>> with the relevant location
            sqlite:///:memory:  [not recommended : misses the point of having a persistent store.  Good for testing only.]
            
        Other databases, including MySql and SQL server, are supported.
        Testing has only been carried out with Sqlite to date.
        
        ewdir = a directory used for creating SQLite database, if specified.
        used to replace <<DEFAULT>> in engineName, as occurs with the LsStore class
        
        example usage:
        
        ewdir='/home/dwyllie/ewPersistence'
        dbname='ewEdgeCache'
        test_path="<<DEFAULT>>/%s.db" % dbname          # used for testing LsStore object
        test_connString="sqlite:///%s" % test_path
        print test_connString

        # define a Guid to register
        guid='00578919-1e88-4ea1-99af-e71a96b9ebbf'
        
        # create EW edge storage system
        ewc=ElephantWalkDataSNPaccessor(db=db_ewss, engineName=test_connString, ewdir=ewdir, maxDistance=20)
        
        # delete anything present - optional, would not do this production
        # ewc.remove_all()
        
        # register all the nodes from the SW server [once only - otherwise, you could register single guids]
        retval=ewc.registerAllGuids(guid=guid)
        if not retval=='OK':
            raise Exception("elephantwalk error, expected OK but got %s" % retval)
                          
        # display nodes
        print list(ewc.nodes())        # it is a generator
        
        # display edges
        print list(ewc.edges())
              
        """
        
        # logging
        logging.getLogger()
         
        # add properties to conform to the expectations of snpDataSources
        self.SNPsource='ElephantWalk'
        self.description='ElephantWalk - {0}'.format(refName)
        self.matrixIsSparse=True
        self.maxDistance=maxDistance
        self.ewReference=refName
        
        ## create database
        ## sqlite (if used ) file location
        # replace <<DEFAULT>> in the connection string with the location of the sqlitedb
        engineName=engineName.replace('<<DEFAULT>>', ewdir)
        engine = create_engine(engineName)    # connection string for a database
        Base=db
        Base.metadata.create_all(engine)  # create the table(s)
        Session = sessionmaker(bind=engine)    # class
        self.engineName=engineName
        
        #print("Using sqlite db using the following handle:%s" % connString)
        self.session=Session()

        # if necessary, add a reference.
        ewr=ewReference(
            ewReferenceName=refName)
        
        try:
            self.session.add(ewr)
            self.session.commit()       # will fail if it's already there
        except excb.IntegrityError:
            self.session.rollback()
            
        # should now exist.
        self.ewref=self.session.query(ewReference).filter(ewReference.ewReferenceName==refName).one()       
    def storeRelations(self, links):
        """ stores the relationship between one guid and others.
        links is a dictionary of the following form:
        {"neighbours": [["0b23587d-e324-4daf-a1ed-a623c4aaab1f", 0, 1000, 2000, 1250]], "guid": "dfeffe2a-594e-463e-aaef-13e47447b7d9"}
        as produced and returned by getEWSNP.info().
        """

        # first check input
        if not type(links)==dict:
            raise TypeError("links should have type dict but it is of type {0}".format(type(links)))
        
        if not('neighbours' in links.keys() and 'guid' in links.keys()):
            raise KeyError("links must have keys neighbours and guid")
        
        if not(type(links['neighbours'])==list):
            raise TypeError("links must have a list of neighbours, but it has a {0} viz: {1}".format(type(links['neighbours']),links))
               
        # check whether guid is present.
        try:
            guidRecord=self.session.query(ewGuid).\
                filter(ewGuid.ewGuid==links['guid']).\
                filter(ewGuid.ewRefId==self.ewref.ewRefId).\
                one()
        
            # now check whether the guid has been previously examined
            if guidRecord.examined==1:
                # do nothing more.  This assumes that the edges of this guid have already been added.
                return(0)
            else: # we have not examined it, but are about to, so we update the record.  
                guidRecord.examined=1
                guidRecord.examinationDateTime=datetime.datetime.now()
                self.session.commit()

        except exc.NoResultFound:       # there isn't one, or wasn't one when we last added
            guidRecord=ewGuid(ewGuid=links['guid'],ewRefId=self.ewref.ewRefId, examined=1, examinationDateTime=datetime.datetime.now())
            try:
                self.session.add(guidRecord)
                self.session.commit()
            except orm_exc.IntegrityError:
                # it is already present
                self.session.rollback()
         
        # rewritten to fix issue #42
        ## identify and add all the guids to which it is linked;
        new_guids=set()
        found_guids=set()
        
        # iterate over the links['neighbours'] list, adding any guids which don't currently exist.
        linked_guids=set()
        if len(links['neighbours'])>0:

            for (this_guid,this_dist, njust_1, njust_2, n) in links['neighbours']:

                linked_guids.add(this_guid)

            # a query recovering the guids currently present for this reference
            guidsPresent=self.session.query(ewGuid.ewGuid).\
                    filter(ewGuid.ewGuid.in_(linked_guids)).\
                    filter(ewGuid.ewRefId==self.ewref.ewRefId)
            
            for (this_present_guid, ) in guidsPresent:
                found_guids.add(this_present_guid)

            # identify guids which are not present, but should be
            notPresent=linked_guids.difference(found_guids)

            # add these new nodes, identified in notPresent
            for this_guid in notPresent:
                linkedRecord=ewGuid(ewGuid=this_guid, ewRefId=self.ewref.ewRefId, examined=0, examinationDateTime=None)
                self.session.add(linkedRecord)
            self.session.commit()

            #Recovering ids of all linked nodes
            guid2guidId={}
            for (this_guid, this_guidId) in self.session.query(ewGuid.ewGuid, ewGuid.ewGuidId).\
                    filter(ewGuid.ewGuid.in_(linked_guids)).\
                    filter(ewGuid.ewRefId==self.ewref.ewRefId):
                guid2guidId[this_guid]=this_guidId
           
            # Now add edges
            query_guidId=guidRecord.ewGuidId
                    
            for (this_guid,this_dist, N_just1, N_just2, N_either) in links['neighbours']:

                guidId=guid2guidId[this_guid]
                guidPair=sorted([guidId,query_guidId])      # only add one way round
                try:
                        this_edge=self.session.query(ewEdge).filter(ewEdge.ewGuidId1==guidPair[0]).filter(ewEdge.ewGuidId2==guidPair[1]).one()
                        
                except exc.NoResultFound:       # there isn't one

                        this_edge=ewEdge(ewGuidId1=guidPair[0],ewGuidId2=guidPair[1], dist=this_dist, N_just1=N_just1, N_just2=N_just2, N_either=N_either)

                        self.session.add(this_edge)
            self.session.commit()
        return(0)
    
    def remove_all(self):
        """ removes all nodes and edges """
        self.session.query(ewEdge).delete()
        self.session.query(ewGuid).delete()
        self.session.commit()
        return(0)
    
    def guuid(self):
        """ returns a guid (used for testing purposes) """
        return str(uuid.uuid1())  
        
    def nodes(self):
        """ a list of nodes """
        results=self.session.query(ewGuid)
        for result in results:
            yield result.ewGuid


    def nodeExists(self,guid):
        """ a list of nodes """
        n=self.session.query(ewGuid).filter(ewGuid.ewGuid==guid).one_or_none()
        if n is None:
            return(False)
        else:
            return(True)

    def nodesObservedAfter(self, cutoffDate=None, returnValue='guid'):
        """ a list of nodes """
        #print("naa, cutoffDate={0}".format(cutoffDate))
        if cutoffDate is None:
            cutoffDate=datetime.datetime(1900,1,1,12,0)       # a very early time
            
        results=self.session.query(ewGuid).filter(ewGuid.examinationDateTime>cutoffDate).filter(ewGuid.examined==1)
        for result in results:
            #print(result.ewGuidId, result.examinationDateTime)
            if returnValue=='guid':
                yield result.ewGuid
            else:
                yield result.ewGuidId
                
    def subsetOfNodes(self,examinationStatus):
        """ a list of nodes with a particular examination status
        
        options are 
        0 not examined
        1 examined """
        
        results=self.session.query(ewGuid).filter(ewGuid.examined==examinationStatus)
        nodelist=[]
        for result in results:
            nodelist.append(result.ewGuid)
        return(nodelist)
       
    def edges(self, cutoffDate=None, cutoff=20):
        """ returns a list of edges.  If cutoffDate is not none,

        only returns edges where one or other node was noted after cutoffDate.
        """
 
        ewGuid1 = aliased(ewGuid, name='ewGuid1')
        ewGuid2 = aliased(ewGuid, name='ewGuid2')
        retVal=[]
        res=self.session.query(
                               ewGuid1.ewGuid.label('ewGuid1'),
                               ewGuid2.ewGuid.label('ewGuid2'),
                               ewEdge.dist,
                               ewEdge.N_just1,
                               ewEdge.N_just2,
                               ewEdge.N_either).\
                               join(ewEdge, ewEdge.ewGuidId1 == ewGuid1.ewGuidId).\
                               join(ewGuid2, ewEdge.ewGuidId2 == ewGuid2.ewGuidId)
        for result in res:
            retVal.append(list(result))
        return(retVal)
 
    def neighboursOf(self, guid, cutoff=20, returned_format=1):
        """ returns neighbours of guid with cutoff <=cutoff.
            Returns links either as format 1 [otherGuid, distance]
                              or as format 2 [otherGuid, distance, N_just1, N_just2, N_either] """
        retVal=[]
       
        guidId,=self.session.query(ewGuid.ewGuidId).filter(ewGuid.ewGuid==guid).one_or_none()
        if guidId is None:
            return([])      # nothing matches something which isn't there
        else:
            ## should be easy to write a three table join to  recover this, but
            # it seems hard in sql alchemy.
            ewGuid1 = aliased(ewGuid, name='ewGuid1')
            ewGuid2 = aliased(ewGuid, name='ewGuid2')
            retVal=[]
            res=self.session.query(
                                   ewGuid1.ewGuid.label('ewGuid1'),
                                   ewGuid2.ewGuid.label('ewGuid2'),
                                   ewEdge.dist,
                                   ewEdge.N_just1,
                                   ewEdge.N_just2,
                                   ewEdge.N_either).\
                                   join(ewEdge, ewEdge.ewGuidId1 == ewGuid1.ewGuidId).\
                                   join(ewGuid2, ewEdge.ewGuidId2 == ewGuid2.ewGuidId).\
                                   filter((ewEdge.ewGuidId1==guidId) | (ewEdge.ewGuidId2==guidId)).\
                                   filter(ewEdge.dist<=cutoff)

            for result in res:
                    if result.ewGuid1==guid:
                        otherGuid=result.ewGuid2
                    else:
                        otherGuid=result.ewGuid1
                    
                    if returned_format == 1:
                        returned_data=[otherGuid, result.dist]

                    elif returned_format == 2:
                        returned_data=[otherGuid, result.dist, result.N_just1, result.N_just2, result.N_either]

                    else:
                        raise ValueError("Unable to understand returned_format = {0}".format(returned_format))
                    retVal.append(returned_data)            
            # recover the guids          
        return({'guid':guid, 'neighbours':retVal})
    
    def pairwise(self, guid1, guid2):
        """ returns details on the pairwise comparison of guid1 with guid2 """
        retVal=[]
       
        guid1_found=False
        guid2_found=False
       
        res1=self.session.query(ewGuid.ewGuidId).filter(ewGuid.ewGuid==guid1).one_or_none()
        if res1 is not None:
            guid1_found=True
        res2=self.session.query(ewGuid.ewGuidId).filter(ewGuid.ewGuid==guid2).one_or_none()
        if res2 is not None:
            guid2_found=True
 
        if not ( guid1_found and guid2_found ):
            return(
                {'guid1':guid1,
                 'guid2':guid2,
                 'guid1_found':guid1_found,
                 'guid2_found':guid2_found,
                 'pairwise':'Failed',
                 'results':{}
                 })
        else:
           
            ## unpack the guid identifiers; search for matching rows
            guidId1,=res1
            guidId2,=res2
           
            res=self.session.query(ewEdge).\
                    filter(or_(\
                            and_(ewEdge.ewGuidId1==guidId1,ewEdge.ewGuidId2==guidId2),\
                            and_(ewEdge.ewGuidId2==guidId1,ewEdge.ewGuidId1==guidId2)\
                            )).one_or_none()
            if res is None:     # no entry.  this is because the value is high, and we don't store high values.
                return({'guid1':guid1,'guid2':guid2, 'guid1_found':True, 'pairwise':'High, not stored', 'results':{}})
            else:
                return({'guid1':guid1,'guid2':guid2, 'guid2_found':True,  'pairwise':'Stored', 'results':{'dist':res.dist,
                                                                                                          'N_just_guid1':res.N_just1,
                                                                                                          'N_just_guid2':res.N_just2,
                                                                                                          'N_either':res.N_either}})
                
    
## unit tests  
class Test_EWSnpStore_integration(unittest.TestCase):
    """ integration test of loading getEWSNP output using sqlite"""
    def runTest(self):
       
        # define temporary directory and sqlite filename for trial
        ewdir=os.path.join('..','unittest_tmp')
        dbname='ewEdgeCache'
       
        # delete any existing sqlite db if present
        try:
            os.unlink(os.path.join(ewdir,"{0}.db".format(dbname)))
        except:
            pass
         
        # use sqlite for testing
        test_path="<<DEFAULT>>/%s.db" % dbname         
        test_connString="sqlite:///%s" % test_path
      
        # create EW edge storage object
        ewc=ElephantWalkDataSNPaccessor(db=db_ewss, engineName=test_connString, ewdir=ewdir, maxDistance=20)
      
        # get some data to enter.  We are parsing output from getEWSNP and loading it.
        # testdate is located in testDataPath as individual json files
        testDataPath=os.path.join('..','testdata','EWSnpStore','*.snps.json')
        testFiles=glob.glob(testDataPath)
        guidsAdded=set()        # set holds the 'true' answer to the number of guids
        guidsLinkedTo=set()     # and their links
        i=0
        for (i,testFile) in enumerate(testFiles):
            try:
                with open(testFile,'rt') as f:
                    #print(testFile)
                    links=json.load(f)
 
                    guidsAdded.add(links['guid'])
                    
                    for element in links['neighbours']:
                        #print(element)
                        guid,snp,n1,n2,nb=element
                        guidsLinkedTo.add(guid)
            except IOError:
               self.fail('Could not open test data set {0}'.format(testFile))
               
            ewc.storeRelations(links)
 
        # test 0 make sure the input state worked
        self.assertTrue(i>0)            # we did add to our test.
       
        # now recover from the ewc object and compare with the 'truth' from parsing the json input;
        notExamined=set()
        for (i, result) in enumerate(ewc.subsetOfNodes(examinationStatus=0)):       # not examined
            notExamined.add(result)
           
        recoveredGuids=set()
        for (i, result) in enumerate(ewc.subsetOfNodes(examinationStatus=1)):       # not examined
            recoveredGuids.add(result)
   
        # test #1 check that ewc subsetOfNodes function is giving the right guids
        self.assertEqual(recoveredGuids,guidsAdded)
        self.assertEqual(notExamined,guidsLinkedTo)
 
        # test #2 the number of nodes should be the sum of ewc.subsetOfNodes with examinationStatus in (0,1)
        for (i, result) in enumerate(ewc.nodes(),start=1):
            pass
        self.assertEqual(i, len(notExamined)+len(recoveredGuids))
       
        
        # test #3 there are edges
        j=0
        for (node1, node2, distance, n1,n2,n3) in ewc.edges():
            j+=1
        self.assertTrue(j>0,True)
         

        # #test #4 with an addition date cutoff prior to any additions, all edges should be recovered
        # i=0
        # for (node1, node2, distance,n1,n2,n3) in ewc.edges(cutoffDate=datetime.datetime(1900,1,1,12,0)):
        #     i+=1
        # # i should be the number of edges
        # self.assertEqual(i,j)
        #           
        # #test #5 with an addition date cutoff after any additions, no edges should be recovered
        # k=0
        # for (node1, node2, distance,n1,n2,n3) in ewc.edges(cutoffDate=datetime.datetime(2200,1,1,12,0)):
        #     k+=1
        # # k should be the number of edges, which should remain at zero.
        # self.assertEqual(k,0)
        
class Test_EWSnpStore_jsonLoad(unittest.TestCase):
    """ integration test of loading getEWSNP output using sqlite"""
    def runTest(self):
       
        # define temporary directory and sqlite filename for trial
        ewdir=os.path.join('..','unittest_tmp')
        dbname='ewEdgeCache'
       
        # delete any existing sqlite db if present
        try:
            os.unlink(os.path.join(ewdir,"{0}.db".format(dbname)))
        except:
            pass
         
        # use sqlite for testing
        test_path="<<DEFAULT>>/%s.db" % dbname         
        test_connString="sqlite:///%s" % test_path
      
        # create EW edge storage object
        ewc=ElephantWalkDataSNPaccessor(db=db_ewss, engineName=test_connString, ewdir=ewdir, maxDistance=20)
      
        # get some data to enter.  We are parsing output from getEWSNP and loading it.
        # testdate is located in testDataPath as individual json files
        testFile=os.path.join('..','testdata','EWSnpStore','test.snps.json')
        try:
            with open(testFile,'rt') as f:
                links=json.load(f)
               
                # check that the guid recorded is correct
                self.assertEqual(links['guid'],'91273e18-da11-447d-8c54-6590faf6b754')
               
                for (guid,snp,n1,n2,n3) in links['neighbours']:
                    # check two examples of distances read by the module
                    if guid=='7a008326-e646-47f5-9e89-dc603c2700df':
                        self.assertEqual(snp,2)
                    if guid=='5578949f-75d7-4b6b-9f4e-fc5632efc68e':
                        self.assertEqual(snp,20)
                    if guid=='307bec33-713c-49f9-ad95-65a46f54f9b4':
                        self.assertEqual(snp,4)
 
        except IOError:
            self.fail('Could not open test data set {0}'.format(testFile))
        
        # store the data      
        ewc.storeRelations(links)
       
        # test the output
        guids={}
        for res in ewc.session.query(ewGuid).all():
            guids[res.ewGuidId]=res.ewGuid
        for res in ewc.session.query(ewEdge).all():
            guid1=guids[res.ewGuidId1]
            guid2=guids[res.ewGuidId2]
           
            # test that the pairs above are correct
            if guid1=='91273e18-da11-447d-8c54-6590faf6b754' and guid2=='7a008326-e646-47f5-9e89-dc603c2700df':
                self.assertEqual(res.dist,2)
            if guid2=='91273e18-da11-447d-8c54-6590faf6b754' and guid1=='7a008326-e646-47f5-9e89-dc603c2700df':
                self.assertEqual(res.dist,2)
            # test that the pairs above are correct
            if guid1=='91273e18-da11-447d-8c54-6590faf6b754' and guid2=='5578949f-75d7-4b6b-9f4e-fc5632efc68e':
                self.assertEqual(res.dist,20)
            if guid2=='91273e18-da11-447d-8c54-6590faf6b754' and guid1=='5578949f-75d7-4b6b-9f4e-fc5632efc68e':
                self.assertEqual(res.dist,20)
                        # test that the pairs above are correct
            if guid1=='91273e18-da11-447d-8c54-6590faf6b754' and guid2=='307bec33-713c-49f9-ad95-65a46f54f9b4':
                self.assertEqual(res.dist,4)
            if guid2=='91273e18-da11-447d-8c54-6590faf6b754' and guid1=='307bec33-713c-49f9-ad95-65a46f54f9b4':
                self.assertEqual(res.dist,4)
 
class Test_EWSnpStore_neighbours(unittest.TestCase):
    """ tests neighboursOf """
    def runTest(self):
       
        # define temporary directory and sqlite filename for trial
        ewdir=os.path.join('..','unittest_tmp')
        dbname='ewEdgeCache'
       
        # delete any existing sqlite db if present
        try:
            os.unlink(os.path.join(ewdir,"{0}.db".format(dbname)))
        except:
            pass
         
        # use sqlite for testing
        test_path="<<DEFAULT>>/%s.db" % dbname          
        test_connString="sqlite:///%s" % test_path
      
        # create EW edge storage object
        ewc=ElephantWalkDataSNPaccessor(db=db_ewss, engineName=test_connString, ewdir=ewdir, maxDistance=20)
      
        # get some data to enter.  We are parsing output from getEWSNP and loading it.
        # testdate is located in testDataPath as individual json files
        testFile=os.path.join('..','testdata','EWSnpStore','test.snps.json')
        try:
            with open(testFile,'rt') as f:
                links=json.load(f)
               
                # check that the guid recorded is correct
                self.assertEqual(links['guid'],'91273e18-da11-447d-8c54-6590faf6b754')
               
                for (guid,snp,n1,n2,n3) in links['neighbours']:
                    # check two examples of distances read by the module
                    if guid=='7a008326-e646-47f5-9e89-dc603c2700df':
                        self.assertEqual(snp,2)
                    if guid=='5578949f-75d7-4b6b-9f4e-fc5632efc68e':
                        self.assertEqual(snp,20)
                    if guid=='307bec33-713c-49f9-ad95-65a46f54f9b4':
                        self.assertEqual(snp,4)
 
        except IOError:
            self.fail('Could not open test data set {0}'.format(testFile))
        
        # store the data      
        ewc.storeRelations(links)
       
        # test the output
        guids={}
        for res in ewc.session.query(ewGuid).all():
            guids[res.ewGuidId]=res.ewGuid
        for res in ewc.session.query(ewEdge).all():
            guid1=guids[res.ewGuidId1]
            guid2=guids[res.ewGuidId2]
           
            # test that the pairs above are correct
            if guid1=='91273e18-da11-447d-8c54-6590faf6b754' and guid2=='7a008326-e646-47f5-9e89-dc603c2700df':
                self.assertEqual(res.dist,2)
            if guid2=='91273e18-da11-447d-8c54-6590faf6b754' and guid1=='7a008326-e646-47f5-9e89-dc603c2700df':
                self.assertEqual(res.dist,2)
            # test that the pairs above are correct
            if guid1=='91273e18-da11-447d-8c54-6590faf6b754' and guid2=='5578949f-75d7-4b6b-9f4e-fc5632efc68e':
                self.assertEqual(res.dist,20)
            if guid2=='91273e18-da11-447d-8c54-6590faf6b754' and guid1=='5578949f-75d7-4b6b-9f4e-fc5632efc68e':
                self.assertEqual(res.dist,20)
                        # test that the pairs above are correct
            if guid1=='91273e18-da11-447d-8c54-6590faf6b754' and guid2=='307bec33-713c-49f9-ad95-65a46f54f9b4':
                self.assertEqual(res.dist,4)
            if guid2=='91273e18-da11-447d-8c54-6590faf6b754' and guid1=='307bec33-713c-49f9-ad95-65a46f54f9b4':
                self.assertEqual(res.dist,4)
     
        res20=ewc.neighboursOf('91273e18-da11-447d-8c54-6590faf6b754')
        for res in res20['neighbours']:
            self.assertEqual(len(res),2)
        res0=ewc.neighboursOf('91273e18-da11-447d-8c54-6590faf6b754', cutoff=0)
        res12=ewc.neighboursOf('91273e18-da11-447d-8c54-6590faf6b754', cutoff=4)
        self.assertEqual(len(res20['neighbours']),14)
        self.assertEqual(len(res12['neighbours']),2)
        self.assertEqual(len(res0['neighbours']),0)
            
                
class Test_EWSnpStore_neighbours_format2(unittest.TestCase):
    """ tests neighboursOf """
    def runTest(self):
       
        # define temporary directory and sqlite filename for trial
        ewdir=os.path.join('..','unittest_tmp')
        dbname='ewEdgeCache'
       

        # delete any existing sqlite db if present
        try:
            os.unlink(os.path.join(ewdir,"{0}.db".format(dbname)))
        except:
            pass
        
        # use sqlite for testing
        test_path="<<DEFAULT>>/%s.db" % dbname          
        test_connString="sqlite:///%s" % test_path
      
        # create EW edge storage object
        ewc=ElephantWalkDataSNPaccessor(db=db_ewss, engineName=test_connString, ewdir=ewdir, maxDistance=20)
      
        # get some data to enter.  We are parsing output from getEWSNP and loading it.
        # testdate is located in testDataPath as individual json files
        testFile=os.path.join('..','testdata','EWSnpStore','test.snps.json')
        try:
            with open(testFile,'rt') as f:
                links=json.load(f)
               
                # check that the guid recorded is correct
                self.assertEqual(links['guid'],'91273e18-da11-447d-8c54-6590faf6b754')
               
                for (guid,snp,n1,n2,n3) in links['neighbours']:

                    # check two examples of distances read by the module
                    if guid=='7a008326-e646-47f5-9e89-dc603c2700df':
                        self.assertEqual(snp,2)
                    if guid=='5578949f-75d7-4b6b-9f4e-fc5632efc68e':
                        self.assertEqual(snp,20)
                    if guid=='307bec33-713c-49f9-ad95-65a46f54f9b4':
                        self.assertEqual(snp,4)
 
        except IOError:
            self.fail('Could not open test data set {0}'.format(testFile))
        
        # store the data      
        ewc.storeRelations(links)

        # test the output
        guids={}
        for res in ewc.session.query(ewGuid).all():
            guids[res.ewGuidId]=res.ewGuid
        for res in ewc.session.query(ewEdge).all():
            guid1=guids[res.ewGuidId1]
            guid2=guids[res.ewGuidId2]

           
            # test that the pairs above are correct
            if guid1=='91273e18-da11-447d-8c54-6590faf6b754' and guid2=='7a008326-e646-47f5-9e89-dc603c2700df':
                self.assertEqual(res.dist,2)
            if guid2=='91273e18-da11-447d-8c54-6590faf6b754' and guid1=='7a008326-e646-47f5-9e89-dc603c2700df':
                self.assertEqual(res.dist,2)
            # test that the pairs above are correct
            if guid1=='91273e18-da11-447d-8c54-6590faf6b754' and guid2=='5578949f-75d7-4b6b-9f4e-fc5632efc68e':
                self.assertEqual(res.dist,20)
            if guid2=='91273e18-da11-447d-8c54-6590faf6b754' and guid1=='5578949f-75d7-4b6b-9f4e-fc5632efc68e':
                self.assertEqual(res.dist,20)
                        # test that the pairs above are correct
            if guid1=='91273e18-da11-447d-8c54-6590faf6b754' and guid2=='307bec33-713c-49f9-ad95-65a46f54f9b4':
                self.assertEqual(res.dist,4)
            if guid2=='91273e18-da11-447d-8c54-6590faf6b754' and guid1=='307bec33-713c-49f9-ad95-65a46f54f9b4':
                self.assertEqual(res.dist,4)
     
        res20=ewc.neighboursOf('91273e18-da11-447d-8c54-6590faf6b754',returned_format=2)
        res0=ewc.neighboursOf('91273e18-da11-447d-8c54-6590faf6b754', cutoff=0, returned_format=2)
        res12=ewc.neighboursOf('91273e18-da11-447d-8c54-6590faf6b754', cutoff=4, returned_format=2)
        
        for res in res20['neighbours']:
            self.assertEqual(len(res),5)
        self.assertEqual(len(res20['neighbours']),14)
        self.assertEqual(len(res12['neighbours']),2)
        self.assertEqual(len(res0['neighbours']),0)
            
 
 
class Test_EWSnpStore_pairwise(unittest.TestCase):
    """ tests neighboursOf """
    def runTest(self):
       
        # define temporary directory and sqlite filename for trial
        ewdir=os.path.join('..','unittest_tmp')
        dbname='ewEdgeCache'
       
        # delete any existing sqlite db if present
        try:
            os.unlink(os.path.join(ewdir,"{0}.db".format(dbname)))
        except:
            pass
         
        # use sqlite for testing
        test_path="<<DEFAULT>>/%s.db" % dbname         
        test_connString="sqlite:///%s" % test_path
      
        # create EW edge storage object
        ewc=ElephantWalkDataSNPaccessor(db=db_ewss, engineName=test_connString, ewdir=ewdir, maxDistance=20)
      
        # get some data to enter.  We are parsing output from getEWSNP and loading it.
        # testdate is located in testDataPath as individual json files
        testFile=os.path.join('..','testdata','EWSnpStore','test.snps.json')
        try:
            with open(testFile,'rt') as f:
                links=json.load(f)
               
                # check that the guid recorded is correct
                self.assertEqual(links['guid'],'91273e18-da11-447d-8c54-6590faf6b754')
               
                for (guid,snp,n1,n2,n3) in links['neighbours']:
                    # check two examples of distances read by the module
                    if guid=='7a008326-e646-47f5-9e89-dc603c2700df':
                        self.assertEqual(snp,2)
                    if guid=='5578949f-75d7-4b6b-9f4e-fc5632efc68e':
                        self.assertEqual(snp,20,)
                    if guid=='307bec33-713c-49f9-ad95-65a46f54f9b4':
                        self.assertEqual(snp,4)
 
        except IOError:
            self.fail('Could not open test data set {0}'.format(testFile))
        
        # store the data      
        ewc.storeRelations(links)
       
        # test the output
        guids={}
        for res in ewc.session.query(ewGuid).all():
           guids[res.ewGuidId]=res.ewGuid
        for res in ewc.session.query(ewEdge).all():
            guid1=guids[res.ewGuidId1]
            guid2=guids[res.ewGuidId2]
           
            # test that the pairs above are correct
            if guid1=='91273e18-da11-447d-8c54-6590faf6b754' and guid2=='7a008326-e646-47f5-9e89-dc603c2700df':
                self.assertEqual(res.dist,2)
            if guid2=='91273e18-da11-447d-8c54-6590faf6b754' and guid1=='7a008326-e646-47f5-9e89-dc603c2700df':
                self.assertEqual(res.dist,2)
            # test that the pairs above are correct
            if guid1=='91273e18-da11-447d-8c54-6590faf6b754' and guid2=='5578949f-75d7-4b6b-9f4e-fc5632efc68e':
                self.assertEqual(res.dist,20)
            if guid2=='91273e18-da11-447d-8c54-6590faf6b754' and guid1=='5578949f-75d7-4b6b-9f4e-fc5632efc68e':
                self.assertEqual(res.dist,20)
                        # test that the pairs above are correct
            if guid1=='91273e18-da11-447d-8c54-6590faf6b754' and guid2=='307bec33-713c-49f9-ad95-65a46f54f9b4':
                self.assertEqual(res.dist,4)
            if guid2=='91273e18-da11-447d-8c54-6590faf6b754' and guid1=='307bec33-713c-49f9-ad95-65a46f54f9b4':
                self.assertEqual(res.dist,4)
   
                
        res=ewc.pairwise('91273e18-da11-447d-8c54-6590faf6b754','307bec33-713c-49f9-ad95-65a46f54f9b4')
        self.assertEqual(res['pairwise'],'Stored')
        self.assertEqual(res['results']['dist'],4)
        res=ewc.pairwise('307bec33-713c-49f9-ad95-65a46f54f9b4','91273e18-da11-447d-8c54-6590faf6b754')
        self.assertEqual(res['pairwise'],'Stored')
        self.assertEqual(res['results']['dist'],4)
        res=ewc.pairwise('307bec33-713c-49f9-ad95-65a46f54f9b4','Nothing')
        self.assertEqual(res['pairwise'],'Failed')
        self.assertEqual(res['results'],{})
       


