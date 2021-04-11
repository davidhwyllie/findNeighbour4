
""" runs unittest for guidLookup

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

import unittest
import datetime
import uuid
from findn.guidLookup import guidSearcher


class test_gm_1(unittest.TestCase):
    def runTest(self):

        gs = guidSearcher()
        gs.add('b1')
        self.assertEqual(gs.guids,['b1'])
        gs.add('b3')
        self.assertEqual(gs.guids,['b1','b3'])
        gs.add('b3')
        self.assertEqual(gs.guids,['b1','b3'])          
        gs.add('b2')
        self.assertEqual(gs.guids,['b1','b2','b3']) 
        gs.add('a1')
        self.assertEqual(gs.guids,['a1','b1','b2','b3'])
        gs.add('c1')
        self.assertEqual(gs.guids,['a1','b1','b2','b3','c1'])

class test_gm_2(unittest.TestCase):
    def runTest(self):

        gs = guidSearcher()
        gs.add('b1')
        gs.add('b3')
        gs.add('b3')
        gs.add('b2')
        gs.add('a1')
        gs.add('c1')
        
        retVal = gs.search('b')
        self.assertEqual(retVal, ['b1','b2','b3'])
        retVal = gs.search('b', max_returned=30)
        self.assertEqual(retVal, ['b1','b2','b3'])
        retVal = gs.search('b', max_returned=1)
        print(retVal)
        self.assertEqual(retVal, ['b1'])
        retVal = gs.search('b', max_returned=1, return_subset=True)
        self.assertEqual(len(retVal),1)
        retVal = gs.search('b', max_returned=1, return_subset=False)
        self.assertEqual(len(retVal),0)

        retVal = gs.search('b', max_returned=2, return_subset=True)
        self.assertEqual(len(retVal),2)
        retVal = gs.search('b', max_returned=2, return_subset=False)
        self.assertEqual(len(retVal),0)
        retVal = gs.search('z')
        self.assertEqual(retVal, [])

class test_gm_3(unittest.TestCase):
    def runTest(self):

        gs = guidSearcher()
        gs.add('b1')
         
        retVal = gs.search('b1')
        self.assertEqual(len(retVal),1)
               
@unittest.skip("benchmark disabled")       
class test_gm_benchmark(unittest.TestCase):
    def runTest(self):
        """ tests time to search a list of 1M guids
        
        median addition time : 0.00026 sec (0.26 ms)
        median search time   : 0.00002 sec (0.02 ms)"""

        gs = guidSearcher()
        
        print("generating list of guids")

        added = []
        for i in range(1000000):
            new_guid = str(uuid.uuid4())
            added.append(new_guid)

        t2=  datetime.datetime.now()           
        print("adding")
        for item in added:
            gs.add(item)
            
        
        t3=  datetime.datetime.now()
        print("searching")
        for item in added:
            gs.search(item)
            
        t4=  datetime.datetime.now()
        print("SEARCH PER SAMPLE", (t4-t3)/len(added))
        print("ADD PER SAMPLE",(t3-t2)/len(added))