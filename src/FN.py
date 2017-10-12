#!/usr/bin/env python
""" contains two classes which help ElephantWalk2.
* FNhelper, a class which provides a storage layer for meta-data (but not SNP distances) about sequences
  linked to this are SQL alchemy table definition classes FN*
* NucleicAcid, a class which cleans and validates DNA sequences """


          
import os
import unittest
import csv 
import datetime
import time
import hashlib
import uuid
import inspect
import sqlalchemy
import json
import pandas as pd
import logging
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Table, Column, Integer, String, Float, DateTime, Boolean, MetaData, select, func, LargeBinary
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship, backref
from sqlalchemy import exc  # exceptions


### beginning of NucleicAcid code and tests
class NucleicAcid():
        """ methods to validate that a string is a nucleic acid, and that it is the expected length.
        
        Raises errors if incorrect data is passed.
        Returns composition of the string in a dictionary composition.
        Returns the cleaned (uppercased) string in self.nucleicAcidString"""
        def __init__(self):
                """ constructor """
                logging.getLogger()
                self._reset()
                
        def _reset(self):
                """ sets properties to their default state """
                self.isValid=None
                self.composition={'examinationDate':datetime.datetime.now(), 'MD5':None, 'length':0}
                self.nchars=0
                self.nucleicAcidString=None
                
        def examine(self,nucleicAcidString):
                """ evaluates the string passed """
                self._reset()
                if not isinstance(nucleicAcidString,str):
                        self.isValid=False
                        raise TypeError ("nucleicAcidString must be of type str, but it is %s " % type(nucleicAcidString))
                
                self.nucleicAcidString=bytes(nucleicAcidString.upper().encode(encoding='utf-8'))       # string not unicode; uppercase
                
                # compute the md5
                self.composition['MD5']=hashlib.md5(self.nucleicAcidString).hexdigest()
                
                # assign the string length to the Composition
                self.nchars=len(self.nucleicAcidString)

                # input is invalid if it has zero length
                if self.nchars==0:
                        self.isValid=False
                        raise ValueError ("nucleicAcidString must be of type str, and of non-zero length")
             
                # define valid characters
                totalValid=0
                validChars=['A','C','G','T','N','-']
                for thisChar in validChars:
                        self.composition[thisChar]=self.nucleicAcidString.count(bytes(thisChar.encode(encoding='utf-8')))
                        totalValid=totalValid+self.composition[thisChar]
                self.composition['invalid']=self.nchars-totalValid
                self.composition['length']=self.nchars
                self.composition['ACTG']=self.composition['A']+self.composition['C']+self.composition['G']+self.composition['T']
                self.composition['propACTG']=float(self.composition['ACTG'])/float(self.composition['length'])
                if self.composition['invalid']==0:       # all characters are permitted
                        self.isValid=True
                else:
                        self.isValid=False
                        raise ValueError ("nucleicAcidString must be of type str, and consist only of ACGTN-.  There are %i invalid characters" % self.composition['invalid'])             
class Test_NucleicAcid_Base1(unittest.TestCase):
    """ base class used for tests of NucleicAcid class """
    def setUp(self):
        pass
    def tearDown(self):
        pass   
class Test_NucleicAcid_1(Test_NucleicAcid_Base1):
    """ when initialised, are all outputs appropriately set."""
    def runTest(self):
        na=NucleicAcid()
      
        self.assertEqual(na.isValid,None)
        self.assertEqual(na.composition['MD5'],None)      
class Test_NucleicAcid_2(Test_NucleicAcid_Base1):
    """ checks number of nucleic acids correctly."""
    def runTest(self):
        na=NucleicAcid()
        na.examine('ACGT')
        self.assertEqual(na.isValid,True)
        self.assertEqual(na.composition['MD5'],'f1f8f4bf413b16ad135722aa4591043e')
class Test_NucleicAcid_3(Test_NucleicAcid_Base1):
    """ checks number of nucleic acids """
    def runTest(self):
        na=NucleicAcid()
        na.examine('ACGTN-')
        self.assertEqual(na.isValid,True)
class Test_NucleicAcid_4(Test_NucleicAcid_Base1):
    """ checks composition of nucleic acids """
    def runTest(self):
        na=NucleicAcid()
        self.assertRaises(ValueError, na.examine, nucleicAcidString='ACGTN-X')
class Test_NucleicAcid_5(Test_NucleicAcid_Base1):
    """ checks that what is passed is a string """
    def runTest(self):
        na=NucleicAcid()
        self.assertRaises(TypeError, na.examine, nucleicAcidString=35)
class Test_NucleicAcid_6(Test_NucleicAcid_Base1):
    """ checks that what is passed is a string """
    def runTest(self):
        na=NucleicAcid()
        self.assertRaises(ValueError, na.examine, nucleicAcidString="")
class Test_NucleicAcid_7(Test_NucleicAcid_Base1):
    """ checks number of nucleic acids """
    def runTest(self):
        na=NucleicAcid()
        na.examine('AACCCGT')
        self.assertEqual(na.composition['A'],2)   
        self.assertEqual(na.composition['C'],3)
class Test_NucleicAcid_8(Test_NucleicAcid_Base1):
    """ checks number of nucleic acids """
    def runTest(self):
        na=NucleicAcid()
        na.examine('AAAN')
        self.assertEqual(na.composition['propACTG'],0.75 )  
### end of NucleicAcid code and tests

### beginning of FNpersistence definitions
db=declarative_base() # classes mapping to persistent database inherit from this; global entity
class FN_Sequence(db):
        __tablename__ = 'FN_Sequence'
        sequenceIntId=Column(Integer, primary_key=True)
        sequenceGuid=Column(String(36), index=True)
        intAnnotation=relationship('FN_IntAnnotation', lazy='dynamic',backref='FN_Sequence')
        strAnnotation=relationship('FN_StrAnnotation', lazy='dynamic',backref='FN_Sequence')
        floatAnnotation=relationship('FN_FloatAnnotation', lazy='dynamic',backref='FN_Sequence')
        dateAnnotation=relationship('FN_DateAnnotation', lazy='dynamic',backref='FN_Sequence')
         
class FN_IntAnnotation(db):
        __tablename__ = 'FN_IntAnnotation'
        IntAnnotationId=Column(Integer, primary_key=True)
        sequenceIntId=Column(Integer, ForeignKey(FN_Sequence.sequenceIntId))
        nameSpace=Column(String(32), index=True)
        tag=Column(String(16), index=True)
        value=Column(Integer, index=True)
        
class FN_StrAnnotation(db):
        __tablename__ = 'FN_StrAnnotation'
        StrAnnotationId=Column(Integer, primary_key=True)
        sequenceIntId=Column(Integer, ForeignKey(FN_Sequence.sequenceIntId))
        nameSpace=Column(String(32), index=True)
        tag=Column(String(16), index=True)
        value=Column(String(1024))
        
class FN_FloatAnnotation(db):
        __tablename__ = 'FN_FloatAnnotation'
        FloatAnnotationId=Column(Integer, primary_key=True)
        sequenceIntId=Column(Integer, ForeignKey(FN_Sequence.sequenceIntId))
        nameSpace=Column(String(32), index=True)
        tag=Column(String(16), index=True)
        value=Column(Float, index=True)
        
class FN_DateAnnotation(db):
        __tablename__ = 'FN_DateAnnotation'
        DateAnnotationId=Column(Integer, primary_key=True)
        sequenceIntId=Column(Integer, ForeignKey(FN_Sequence.sequenceIntId))
        nameSpace=Column(String(32), index=True)
        tag=Column(String(16), index=True)
        value=Column(DateTime, index=True)
        
class FNPersistence():
        """ System for persisting results from  large numbers of sequences stored in FindNeighbour.
        Uses SQL alchemy, so can couple to multiple databases.
        
        Example usage is as below:
        
                # example usage
        print("Running an example; ")

        dna=NucleicAcid()
        persist=FNPersistence(db=db, rootdir='/home/local/GEL/dwyllie/dev/FNP/test')
        
        # add some sequences
        seqs={'guid1':'ACGT','guid2':'NACT', 'guid3':'TTTT', 'guid4':'NNNN'}
        for guid in seqs.keys():
                seq=seqs[guid]
                
                ## should try/catch here, as if the sequence is invalid, an error will be raised
                dna.examine(seq)
                
                ## if no errors are raised, you can annotate and add to elephantwalk.
                persist.annotateFromDict(sequenceGuid=guid, nameSpace='DNAQuality',annotDict=dna.composition)
  
        
        print("Recovering guid quality for a list of guids")
        print persist.guid2quality(['guid1','guid2','guid3'])

        print("Recovering guid quality for a list of guids which don't exist")
        print persist.guid2quality(['guid10','guid11','guid12'])
            
        # how to recover guids
        print("Recovering all guids")
        for (guid, intId) in persist.sequenceGuids():
                print(guid, intId)
        
        print("Recovering all guids as json")
        print(persist.asJson_sequenceGuids())
        expectedResult='[{"guid": "guid1"}, {"guid": "guid2"}, {"guid": "guid3"}, {"guid": "guid4"}]'
        
        # recovering guids and examination times;
        print("Recovering guids and examination times")
        for (guid,nameSpace,tag,value) in persist.guid2ExaminationDateTime():
                print(guid,nameSpace,tag,value)
        
        print("Recovering guids and examination times as json")
        print(persist.asJson_guid2ExaminationDateTime())
        
        #  recovered guids filtered by the propACTG criterion
        print("Recovering filtered guids")
        for guid in persist.propACTG_filteredSequenceGuids(cutoff=0.85):
                print(guid)

        print("Recovering filtered guids as json")
        print(persist.asJson_propACTG_filteredSequenceGuids())
        expected=[{"guid": "guid1", "examinationTime": "2016-08-30T12:41:05.935435"}, {"guid": "guid2", "examinationTime": "2016-08-30T12:41:06.031308"}, {"guid": "guid3", "examinationTime": "2016-08-30T12:41:06.096251"}, {"guid": "guid4", "examinationTime": "2016-08-30T12:41:06.164199"}]

        # how to recover all sequenceIntIds, which are small (integer) tokens with 1:1 concordance with a guid.
        print("Recovering sequenceIntIds (integer tokens, instead of guids)")
        print persist.allSequenceIntId()                # not very useful, suggest do not expose this.  These are internal PKs to the tables
        expected=set([1,2,3,4])
        
        print("Showing all annotations")
        # not very nice format.  better to cross tabulate and return as json
        for item in persist.allAnnotations():     # load into pandas
                print(item)
                
        print("Showing all annotations as json")
        print(persist.asJson_allAnnotations())
        
        """   
        def __init__(self, db, engineName='sqlite:///<<DEFAULT>>/findNeighbour.db', rootdir=os.getcwd()):
            """ configures the SQLalchemy engine and file directories ready for use
            
            engineName is a python format database connection string, e.g.
            
            sqlite:///../sqlitedb/v7.db.  The explicit file path has been included.
            sqlite:///<<DEFAULT>>/v7.db.  LsStore will replace <<DEFAULT>> with the relevant location
            sqlite:// for in memory
            
            Other databases, including MySql and SQL server, are supported.
            Testing has only been carried out with Sqlite to date.
        
            rootdir is the location for creating files
            """
            # initiate 
            # ensure directory structure exists.
            # Store directories in a dictionary, which is part of self.
            logging.getLogger()
            
            self._requiredSubDirs=['sqlitedb']                
            self.directories={}
            self.rootdir=rootdir
            
            for thisSubDir in self._requiredSubDirs:
                    requiredDir=os.path.join(rootdir, thisSubDir)
                    self.directories[thisSubDir]=requiredDir
                    self.ensureDir(requiredDir)
                    
            
            # replace <<DEFAULT>> in the connection string with the location of the sqlitedb
            engineName=engineName.replace('<<DEFAULT>>', self.directories['sqlitedb'])
            #print("creating database and engine with %s" % engineName)
            # create or open the database
            self.Base = db             # instance of the ORM
            self.engine = create_engine(engineName)    # connection string for a database
            self.Base.metadata.create_all(self.engine)  # create the table(s)
            Session = sessionmaker(bind=self.engine)    # class
            self.session=Session()        
        def __del__(self):
            """ closes any session """
            try:
                self.session.close() 
            except:
                pass
        def ensureDir(self,path):
            """ ensures that a directory exists.  creates if it doesn't """
            if not os.path.exists(path):
                os.makedirs(path)
    
        def annotateFromDict(self, sequenceGuid, nameSpace, annotDict):
            """ adds multiple annotations from a dictionary."""
            
            # checks whether the sequenceGuid exists;
            sequenceGuidExists=self.session.query(FN_Sequence).filter(FN_Sequence.sequenceGuid==sequenceGuid).count()
            if sequenceGuidExists==0:
                    # create it
                    sequence=FN_Sequence(sequenceGuid=sequenceGuid)
                    self.session.add(sequence)
            # recover the sequenceIntId;
            thisSequenceIntId=self.session.query(FN_Sequence.sequenceIntId).filter(FN_Sequence.sequenceGuid==sequenceGuid).one()[0]                
    
            for tag in annotDict:
                value=annotDict[tag]
                valueType=type(value)
                if valueType==int:
                    # add integer, removing any previous versions;
                    self.session.query(FN_IntAnnotation).filter(FN_IntAnnotation.sequenceIntId==thisSequenceIntId).filter(FN_IntAnnotation.nameSpace==nameSpace).filter(FN_IntAnnotation.tag==tag).delete()
                    a=FN_IntAnnotation(sequenceIntId=thisSequenceIntId, nameSpace=nameSpace, tag=tag, value=value)
                    self.session.add(a)   # add the annotation
                if valueType==float:
                    # add float, deleting any existing value.
                    self.session.query(FN_FloatAnnotation).filter(FN_FloatAnnotation.sequenceIntId==thisSequenceIntId).filter(FN_FloatAnnotation.nameSpace==nameSpace).filter(FN_FloatAnnotation.tag==tag).delete()
                    a=FN_FloatAnnotation(sequenceIntId=thisSequenceIntId, nameSpace=nameSpace, tag=tag, value=value)
                    self.session.add(a)   # add the annotation
                if valueType==str:
                    # add string, deleting any prior results
                    self.session.query(FN_StrAnnotation).filter(FN_StrAnnotation.sequenceIntId==thisSequenceIntId).filter(FN_StrAnnotation.nameSpace==nameSpace).filter(FN_StrAnnotation.tag==tag).delete()
                    a=FN_StrAnnotation(sequenceIntId=thisSequenceIntId, nameSpace=nameSpace, tag=tag, value=value)
                    self.session.add(a)   # add the annotation
                if valueType==datetime.datetime:
                    # add date, deleting any prior results
                    self.session.query(FN_DateAnnotation).filter(FN_DateAnnotation.sequenceIntId==thisSequenceIntId).filter(FN_DateAnnotation.nameSpace==nameSpace).filter(FN_DateAnnotation.tag==tag).delete()
                    a=FN_DateAnnotation(sequenceIntId=thisSequenceIntId, nameSpace=nameSpace, tag=tag, value=value)
                    self.session.add(a)   # add the annotation                
            self.session.commit()
            
        def sequenceGuids(self):
                """ returns all registered guids, and their integer representations """
                ## NEEDS UNIT TEST
                for result in self.session.query(FN_Sequence.sequenceGuid, FN_Sequence.sequenceIntId).all():
                        yield result
        
        def testIndividualGuidQuality(self,guid,cutoff):
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

         exists=self.session.query(FN_Sequence.sequenceGuid).join(FN_FloatAnnotation).filter(FN_FloatAnnotation.nameSpace=='DNAQuality').filter(FN_Sequence.sequenceGuid==guid).count()
         if exists==0:
                 # does not exist
                 return None

         result=self.session.query(FN_Sequence.sequenceGuid).join(FN_FloatAnnotation).filter(FN_FloatAnnotation.nameSpace=='DNAQuality').filter(FN_FloatAnnotation.tag=='propACTG').filter(FN_FloatAnnotation.value>=cutoff).filter(FN_Sequence.sequenceGuid==guid).count()
         if result==1:
                 return True
         if result==0:
                 return False
         
         # otherwise, we return an error
         raise LookupError ("The result of querying the database with %s was neither 1 nor 0; it was %s" % (guid, result))

        def guid2quality(self, guidList):
                """ returns all registered guids and their quality score """
                retDict={}
                if len(guidList)>0:
                        for result in self.session.query(FN_Sequence.sequenceGuid,
                                                         FN_FloatAnnotation.value).\
                                                         join(FN_FloatAnnotation).\
                                                         filter(FN_FloatAnnotation.nameSpace=='DNAQuality').\
                                                         filter(FN_FloatAnnotation.tag=='propACTG').\
                                                         filter(FN_Sequence.sequenceGuid.in_(guidList)).all():
                                retDict[str(result[0])]=result[1]
                return retDict

        def asJson_sequenceGuids(self):
                """ wrapper around sequenceGuids(), but returns result as json """
                data=[]
                for (guid, intId) in self.sequenceGuids():
                        data.append({'guid':guid})
                return(json.dumps(data))
        
        def propACTG_filteredSequenceGuids(self, cutoff=0.85):
                for result in self.session.query(FN_Sequence.sequenceGuid, FN_FloatAnnotation.nameSpace, FN_FloatAnnotation.tag, FN_FloatAnnotation.value).join(FN_FloatAnnotation).filter(FN_FloatAnnotation.nameSpace=='DNAQuality').filter(FN_FloatAnnotation.tag=='propACTG').filter(FN_FloatAnnotation.value>=cutoff).all():
                        yield result
        
        def asJson_propACTG_filteredSequenceGuids(self, cutoff=0.85):
                data=[]
                for (guid,nameSpace,propACTG,quality) in self.propACTG_filteredSequenceGuids(cutoff=cutoff):
                        data.append(guid)
                return(json.dumps(data))
                       
        def allSequenceIntId(self):
                """ returns a set containing all sequenceIntIds """
                ## NEEDS UNIT TESTS
                sequenceIntIds=set()
                for (guid, intId) in self.sequenceGuids():
                        sequenceIntIds.add(intId)
                return sequenceIntIds

        def allSequenceGuid(self):
                """ returns a set containing all sequenceIntIds """
                ## NEEDS UNIT TESTS
                sequenceGuids=set()
                for (guid, intId) in self.sequenceGuids():
                        sequenceGuids.add(guid)
                return sequenceGuids

        def guid2ExaminationDateTime(self):
                """ contains all guids and their exmamination dates """
                for result in self.session.query(FN_Sequence.sequenceGuid, FN_DateAnnotation.nameSpace, FN_DateAnnotation.tag, FN_DateAnnotation.value).join(FN_DateAnnotation).filter(FN_DateAnnotation.tag=='examinationDate').filter(FN_DateAnnotation.nameSpace=='DNAQuality').all():
                        yield result
        
        def asJson_guid2ExaminationDateTime(self):
                """ wrapper round guid2ExaminationDateTime; returns result as json """
                data=[]
                for (guid,nameSpace,tag,value) in self.guid2ExaminationDateTime():
                        data.append({'guid':guid, 'examinationTime':value.isoformat()})
                return(json.dumps(data))

        def allAnnotations(self):
                """ returns all annotations of the sequences """
                for result in self.session.query(FN_Sequence.sequenceGuid, FN_FloatAnnotation.nameSpace, FN_FloatAnnotation.tag, FN_FloatAnnotation.value).join(FN_FloatAnnotation).all():
                        yield result
                for result in self.session.query(FN_Sequence.sequenceGuid, FN_IntAnnotation.nameSpace, FN_IntAnnotation.tag, FN_IntAnnotation.value).join(FN_IntAnnotation).all():
                        yield result
                for result in self.session.query(FN_Sequence.sequenceGuid, FN_StrAnnotation.nameSpace, FN_StrAnnotation.tag, FN_StrAnnotation.value).join(FN_StrAnnotation).all():
                        yield result
                for result in self.session.query(FN_Sequence.sequenceGuid, FN_DateAnnotation.nameSpace, FN_DateAnnotation.tag, FN_DateAnnotation.value).join(FN_DateAnnotation).all():
                        yield result
        
        def asJson_allAnnotations(self):
                """ wrapper around allAnnotations, which returns the annotations as a cross-tabulated json object.
                Note: would be suitable for jqueryGrid etc.
                Not clear how slow it would be for vast amounts of data - not tested to date.
                Requires pandas (python data frame library)
                to crosstabulate """
                
                rows=[]
                for (guid, nameSpace,tag, value) in self.allAnnotations():     # load into pandas
                        rows.append({'guid':guid, 'tag':nameSpace+":"+tag, 'value':value})
                df=pd.DataFrame.from_records(rows)
                xt=pd.crosstab(index=df['guid'], columns=df['tag'], values=df['value'], aggfunc=max)  # cross tab; only one entry per cell, so 'max' is OK as a summary function;
                return(xt.to_json(date_format='iso'))

## persistence unit tests
class Test_FN_version(unittest.TestCase):
    """ tests version of library.  only tested with > v1.0"""   
    def runTest(self): 
        self.assertTrue(sqlalchemy.__version__[0:3]>='1.0')
class Test_FN_createTable(unittest.TestCase):
    """ creates a table directly using sqlAlchemy.  Mainly a test of database connectivity """
    def runTest(self):
        Base = declarative_base()      # for the ORM
        engine = create_engine('sqlite://')    # in memory
        class TestTable(db):                  # sqlAlchemy table
                __tablename__ = 'testTable'
                testId=Column(Integer, primary_key=True)                            
        Base.metadata.create_all(engine)  # create the table(s)
class Test_FN_Base(unittest.TestCase):
        """ Initialises FN persistence. adds no data"""
        def setUp(self):
                # in memory test
                test_connstring="sqlite://" 
                self.t=FNPersistence(db=db, engineName=test_connstring)  
        def tearDown(self):
                self.t.session.flush()
                self.t.session.close()             
class Test_FN_Base1(Test_FN_Base):
        """ initialise FN persistence and adds data """     
        def setUp(self):
                # in memory test
                test_connstring="sqlite://" 
                self.t=FNPersistence(db=db, engineName=test_connstring)  
                dna=NucleicAcid()

                # add some sequences
                seqs={'guid1':'ACGT','guid2':'NACT', 'guid3':'TTTT', 'guid4':'NNNN'}
                for guid in seqs.keys():
                        seq=seqs[guid]
                        
                        ## should try/catch here, as if the sequence is invalid, an error will be raised
                        dna.examine(seq)
                        self.t.annotateFromDict(sequenceGuid=guid, nameSpace='DNAQuality',annotDict=dna.composition)
          
        
        def tearDown(self):
                self.t.session.flush()
                self.t.session.close()
         
class Test_FN_sequenceGuids(Test_FN_Base1):       
        """ tests sequenceGuids """
        def runTest(self):
                n=0
                for (guid, intId) in self.t.sequenceGuids():
                        n+=1
                expected=4
                self.assertEqual(n,expected)

class Test_FN_json_sequenceGuids(Test_FN_Base1):
        """ tests asJson sequenceGuids """
        def runTest(self):
               resdict=json.loads(self.t.asJson_sequenceGuids())
               guids=set()
               for r in resdict:
                guids.add(r['guid'])
               self.assertEqual(guids, set(['guid1','guid2','guid3','guid4']))
        
class Test_FN_guid2ExaminationDateTime(Test_FN_Base1):        
        """ recovering guids and examination times; """
        def runTest(self):
 
                n=0
                for (guid,nameSpace,tag,value) in self.t.guid2ExaminationDateTime():
                            n+=1
                expected=4
                self.assertEqual(n,expected)
 
class Test_FN_propACTG_filteredSequenceGuids(Test_FN_Base1):      
        """  recovered guids filtered by the propACTG criterion """
        def runTest(self):
        
               n=0
               for guid in self.t.propACTG_filteredSequenceGuids(cutoff=0.85):
                       n+=1
               expected=2
               self.assertEqual(n,expected)

class Test_FN_asJson_propACTG_filteredSequenceGuids(Test_FN_Base1):      
        """  recovered guids filtered by the propACTG criterion """
        def runTest(self):
               res=self.t.asJson_propACTG_filteredSequenceGuids(cutoff=0.85)
               self.assertTrue(res in ['["guid1", "guid3"]','["guid3", "guid1"]'])
               
class Test_FN_allSequenceIntId(Test_FN_Base1):      
        """ how to recover all sequenceIntIds, which are small (integer) tokens with 1:1 concordance with a guid. """
        def runTest(self):
 
                res=self.t.allSequenceIntId()                # not very useful, suggest do not expose this.  These are internal PKs to the tables
                expected=set([1,2,3,4])
                self.assertEqual(res,expected)

class Test_FN_allAnnotations(Test_FN_Base1):
        """ tests recovery of all annoations """
        def runTest(self):
                n=0
                for item in self.t.allAnnotations():     
                        n+=1
                expected=48
                self.assertEqual(n,expected)
                
class Test_FN_AddIntTag(Test_FN_Base):
        def runTest(self):
               """ tests addition of an integer """
               # check it is not there
               self.assertEqual(self.t.session.query(FN_IntAnnotation).count(),0)
               
               self.t.annotateFromDict(sequenceGuid='123-456',nameSpace='default',annotDict={'nACTG':10})

               # check it is there
               self.assertEqual(self.t.session.query(FN_IntAnnotation).count(),1)
class Test_FN_AddFloatTag(Test_FN_Base):
        def runTest(self):
               """ tests addition of an float """
               # check it is not there
               self.assertEqual(self.t.session.query(FN_FloatAnnotation).count(),0)
               
               self.t.annotateFromDict(sequenceGuid='123-456',nameSpace='default',annotDict={'floatN':0.75})
               
               # check it is there
               self.assertEqual(self.t.session.query(FN_FloatAnnotation).count(),1)
class Test_FN_AddStrTag(Test_FN_Base):
        def runTest(self):
               """ tests addition of an float """
               # check it is not there
               self.assertEqual(self.t.session.query(FN_StrAnnotation).count(),0)
               
               self.t.annotateFromDict(sequenceGuid='123-456',nameSpace='default',annotDict={'myName':'David'})
               
               # check it is there
               self.assertEqual(self.t.session.query(FN_StrAnnotation).count(),1)
class Test_FN_AddDateTag(Test_FN_Base):
        def runTest(self):
               """ tests addition of an float """
               # check it is not there
               self.assertEqual(self.t.session.query(FN_DateAnnotation).count(),0)
               
               self.t.annotateFromDict(sequenceGuid='123-456',nameSpace='default',annotDict={'myName':datetime.datetime(2016,6,15)})
               
               # check it is there
       
               self.assertEqual(self.t.session.query(FN_DateAnnotation).count(),1)
class Test_FN_AddNAoutput(Test_FN_Base):
        def runTest(self):
               """ tests addition of output from the NucleicAcid class. """
               # check nothing is there
               self.assertEqual(self.t.session.query(FN_DateAnnotation).count(),0)
               self.assertEqual(self.t.session.query(FN_StrAnnotation).count(),0)
               self.assertEqual(self.t.session.query(FN_FloatAnnotation).count(),0)
               
               # set up nucleic acid object
               na=NucleicAcid()
               na.examine('ACGTACGTNN')
               self.t.annotateFromDict(sequenceGuid='123-456',nameSpace='default',annotDict=na.composition)
               
               # check it is there
               self.assertEqual(self.t.session.query(FN_DateAnnotation).count(),1)
               self.assertEqual(self.t.session.query(FN_StrAnnotation).count()>0, True)                    
               self.assertEqual(self.t.session.query(FN_FloatAnnotation).count()>0, True)
class Test_FN_GetSequences(Test_FN_Base):
        def runTest(self):
               """ tests return of sequences """
               
               # set up nucleic acid object
               na=NucleicAcid()
               na.examine('ACGTACGTNN')
               self.t.annotateFromDict(sequenceGuid='1',nameSpace='DNAQuality',annotDict=na.composition)
               self.t.annotateFromDict(sequenceGuid='2',nameSpace='DNAQuality',annotDict=na.composition)
               self.t.annotateFromDict(sequenceGuid='3',nameSpace='DNAQuality',annotDict=na.composition)
                                  
               # check all are recovered
               result=sorted(list(self.t.sequenceGuids()))      # returns tuples sequenceGuid,sequentIntId            
               self.assertEqual(str(result[0][0]),'1')
               self.assertEqual(str(result[1][0]),'2')
               self.assertEqual(str(result[2][0]),'3')
class Test_FN_allSequenceIntIds(Test_FN_Base):
        def runTest(self):
               """ tests return of sequences """
               
               # set up nucleic acid object
               na=NucleicAcid()
               na.examine('ACGTACGTNN')
               self.t.annotateFromDict(sequenceGuid='1',nameSpace='DNAQuality',annotDict=na.composition)
               self.t.annotateFromDict(sequenceGuid='2',nameSpace='DNAQuality',annotDict=na.composition)
               self.t.annotateFromDict(sequenceGuid='3',nameSpace='DNAQuality',annotDict=na.composition)
                                  
               # check all are recovered
               
               result=self.t.allSequenceIntId()  # set of sequentIntId            
               self.assertEqual(result,{1,2,3})             
class Test_FN_guid2quality(Test_FN_Base):
        def runTest(self):
               """ tests return of sequences and their qualities """
               # set up nucleic acid object
               na=NucleicAcid()
               na.examine('ACGTACGTNN')         # 20% bad
               self.t.annotateFromDict(sequenceGuid='g1',nameSpace='DNAQuality',annotDict=na.composition)
               na.examine('ACGTACNNNN')         # 40% bad
               self.t.annotateFromDict(sequenceGuid='g2',nameSpace='DNAQuality',annotDict=na.composition)
               na.examine('ACGTNNNNNN')         # 60% bad
               self.t.annotateFromDict(sequenceGuid='g3',nameSpace='DNAQuality',annotDict=na.composition)
               
               self.assertEqual(len(list(self.t.sequenceGuids())),3)          # check 
           
                # check results
               resDict=self.t.guid2quality([u'g1',u'g2',u'g3'])
               self.assertEqual(resDict['g1'],0.80)                
               self.assertEqual(resDict['g2'],0.60)                
               self.assertEqual(resDict['g3'],0.40)
class Test_FN_testIndividualGuidQuality1(Test_FN_Base):
        def runTest(self):
               """ tests return of sequences and their qualities """
               # set up nucleic acid object
               na=NucleicAcid()
               na.examine('ACGTACGTNN')         # 20% bad
               self.t.annotateFromDict(sequenceGuid='g1',nameSpace='DNAQuality',annotDict=na.composition)
               na.examine('ACGTACNNNN')         # 40% bad
               self.t.annotateFromDict(sequenceGuid='g2',nameSpace='DNAQuality',annotDict=na.composition)
               na.examine('ACGTNNNNNN')         # 60% bad
               self.t.annotateFromDict(sequenceGuid='g3',nameSpace='DNAQuality',annotDict=na.composition)
               
               r1=self.t.testIndividualGuidQuality('g1',0.80)             # valid
               r2=self.t.testIndividualGuidQuality('g2',0.80)             # invalid  
               r3=self.t.testIndividualGuidQuality('g3',0.80)             # invalid  
               r4=self.t.testIndividualGuidQuality('g4',0.80)             # invalid; does not exist.
               
               self.assertEqual(r1, True)
               self.assertEqual(r2, False)
               self.assertEqual(r3, False)               
               self.assertEqual(r4, None)

## unit tests 
if __name__ == '__main__':
        
        # example usage
        print("Running an example; ")

        dna=NucleicAcid()

        # startup only: for setup
        dbfile='/home/local/GEL/dwyllie/dev/FNP/test/findNeighbour.db'
        if os.path.exists(dbfile):
                os.remove(dbfile)
   
        persist=FNPersistence(db=db, rootdir='/home/local/GEL/dwyllie/dev/FNP/test')
        
        # add some sequences
        seqs={'guid1':'ACGT','guid2':'NACT', 'guid3':'TTTT', 'guid4':'NNNN'}
        for guid in seqs.keys():
                seq=seqs[guid]
                
                ## should try/catch here, as if the sequence is invalid, an error will be raised
                dna.examine(seq)
                
                ## if no errors are raised, you can annotate and add to elephantwalk.
                persist.annotateFromDict(sequenceGuid=guid, nameSpace='DNAQuality',annotDict=dna.composition)
  
        
        # how to recover guids
        print("Recovering all guids")
        n=0
        for (guid, intId) in persist.sequenceGuids():
                print(guid, intId)
                n+=1
        expected=4
        
        print("Recovering all guids as json")
        print(persist.asJson_sequenceGuids())
        expectedResult='[{"guid": "guid1"}, {"guid": "guid2"}, {"guid": "guid3"}, {"guid": "guid4"}]'
        
        # recovering guids and examination times;
        print("Recovering guids and examination times")
        n=0
        for (guid,nameSpace,tag,value) in persist.guid2ExaminationDateTime():
                print(guid,nameSpace,tag,value)
                n+=1
        expected=4
        
        print("Recovering guids and examination times as json")
        print(persist.asJson_guid2ExaminationDateTime())
        
        #  recovered guids filtered by the propACTG criterion
        print("Recovering filtered guids")
        n=0
        for guid in persist.propACTG_filteredSequenceGuids(cutoff=0.85):
                n+=1
                print(guid)
        expected=2
        
        print("Recovering filtered guids as json")
        print(persist.asJson_propACTG_filteredSequenceGuids(cutoff=0.85))
        expected=[{"guid": "guid1", "examinationTime": "2016-08-30T12:41:05.935435"}, {"guid": "guid2", "examinationTime": "2016-08-30T12:41:06.031308"}, {"guid": "guid3", "examinationTime": "2016-08-30T12:41:06.096251"}, {"guid": "guid4", "examinationTime": "2016-08-30T12:41:06.164199"}]

        # how to recover all sequenceIntIds, which are small (integer) tokens with 1:1 concordance with a guid.
        print("Recovering sequenceIntIds (integer tokens, instead of guids)")
        print(persist.allSequenceIntId())                # not very useful, suggest do not expose this.  These are internal PKs to the tables
        expected=set([1,2,3,4])
        
        print("Showing all annotations")
        # not very nice format.  better to cross tabulate and return as json
        n=0
        for item in persist.allAnnotations():     # load into pandas
                n+=1
                print(item)
        expected=48
        
        print("Showing all annotations as json")
        print(persist.asJson_allAnnotations())

        
        print("Now running unit tests; ")
        unittest.main()            ## test everything
        print("Done!")

        
