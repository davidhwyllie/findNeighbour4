#!/usr/bin/env python
""" returns a hash on a list.  Used for identifying clusters

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
import unittest

class IdentifySequenceSet():
    """ generates a unique key for a given set of sequence identifiers """
    def __init__(self):
        """ create the object """
        self.permitted_collection_types = set(['cluster', 'connections','msa', 'distmat', 'njtree', 'iqtree', 'raxml'])
        self.has_outgroup2label = {True:'og',False:'no_og'}
        self.permitted_what = set(['M','N','N_or_M','-'])

    def _hashComponents(self,x:list)->str:
        """ returns an sha1 hash on a list, x """
        x = sorted(x)

        to_hash = ";".join([str(item) for item in x])
        return hashlib.sha1(to_hash.encode('utf-8')).hexdigest()

    def make_identifier(self, collection_type, what, has_outgroup, x):
        """ returns id for a representation of a collection of samples.
            collection_type can be one of the following:
                msa, distmat, njtree, iqtree, raxml [see self.permitted_collection_types]
            what: the uncertain base type: N,M, or N_or_M
            has_outgroup: bool, whether there is an outgroup
            x: a collection of sample names, excluding any outgroup  """
               
        if not collection_type in self.permitted_collection_types:
            raise ValueError("collection type {0} is not allowed".format(collection_type))
        if not what in self.permitted_what:
            raise ValueError("what {0} is not allowed".format(what))
        if not isinstance(has_outgroup,bool):
            raise TypeError("has_outgroup must be boolean")

        return "{0}|{1}|{2}|{3}".format(
                collection_type, 
                what,
                self.has_outgroup2label[has_outgroup],
                self._hashComponents(x)
                )
class Test_iss_1(unittest.TestCase):
    """ tests the identifier() method"""
    def runTest(self):
        m = IdentifySequenceSet()     
        self.assertEqual(m.make_identifier('msa','N',True,[1,2,3]),'msa|N|og|cfd68e36493b6db0bc44c859e1e49290e07efb00')
        self.assertEqual(m.make_identifier('msa','N',False,[1,2,3]),'msa|N|no_og|cfd68e36493b6db0bc44c859e1e49290e07efb00')
        

class Test_iss(unittest.TestCase):
    """ tests the _hashComponents() method"""
    def runTest(self):
        m = IdentifySequenceSet()     
        self.assertEqual(m._hashComponents([1,2,3]),'cfd68e36493b6db0bc44c859e1e49290e07efb00')
        self.assertEqual(m._hashComponents(['a','b','c']), 'a58bef80a0b3c42054baeea0edc2308989c7562c')



