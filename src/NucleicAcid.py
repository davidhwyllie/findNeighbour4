#!/usr/bin/env python
""" NucleicAcid, a class which cleans and validates DNA sequences """
          
import os
import unittest
import hashlib
import logging
import datetime

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
                iupac = ['R','r','Y','y','S','s','W','w','K','k','M','m']       # dinucleotide codes
                validChars=['A','C','G','T','N','-']+iupac
                nMixed = 0
                for thisChar in validChars:
                        self.composition[thisChar]=self.nucleicAcidString.count(bytes(thisChar.encode(encoding='utf-8')))
                        totalValid=totalValid+self.composition[thisChar]
                        if thisChar in iupac:
                            nMixed = nMixed + self.composition[thisChar]
                self.composition['invalid']=self.nchars-totalValid
                self.composition['length']=self.nchars
                self.composition['ACTG']=self.composition['A']+self.composition['C']+self.composition['G']+self.composition['T']
                self.composition['mixed']=nMixed
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
class Test_NucleicAcid_9(Test_NucleicAcid_Base1):
    """ checks number of nucleic acids """
    def runTest(self):
        na=NucleicAcid()
        na.examine('AAAR')
        self.assertEqual(na.composition['mixed'],1 )  
        self.assertEqual(na.composition['R'],1 )  
