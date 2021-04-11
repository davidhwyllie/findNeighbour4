#!/usr/bin/env python
""" NucleicAcid, a class which cleans and validates DNA sequences 

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.  See see <https://www.gnu.org/licenses/>.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.
"""

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
                iupac_di = ['R','r','Y','y','S','s','W','w','K','k','M','m']       # dinucleotide codes
                iupac_tri = ['B','D','H','V','b','d','h','v']                      # trinucleotide codes
                iupac = iupac_di + iupac_tri
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
                        raise ValueError ("nucleicAcidString must be of type str, and consist only of ACGTN-.  There are %i invalid characters".format(self.composition['invalid']))       
