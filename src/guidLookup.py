#!/usr/bin/env python
""" 
module to search guids
"""
import bisect

# for unit tests only
import uuid
import unittest
import datetime

class guidSearcher():
    """ maintains a list of sample numbers (such as guids) and provides
    methods to search them rapidly """
    
    def __init__(self):
        """ initialise a guidSearcher.
        
        Example usage:
        
        gs = guidSearcher()
        gs.add('b1')
        gs.add('b3')
        gs.add('b3')
        gs.add('b2')
        gs.add('a1')
        gs.add('c1')
        print(gs.guids)
        
        
        """
        self.guids = []
    
    def add(self,guid):
        """ add a string ('guid') into an ordered list """
        insertion_point = bisect.bisect_left(self.guids, guid)
       
        if not isinstance(guid, str):
            raise TypeError("added item must be a string, not a {0}".format(type(guid)))
        
        # test whether the item exists
        already_exists = False
        try:
            if guid == self.guids[insertion_point]:
                already_exists = True
        except IndexError:
            # doesn't exist
            pass
        if not already_exists:
            self.guids.insert(insertion_point, guid)

    def search(self, search_string, max_returned=30, return_subset=False):
        """
        search_string  the substring in self.guids sought at the beginning of the string
        max_returned   the maximum number of matches returned
        return_subset  whether to return max_returned records if more than max_returned records are found.
        """
        
        search_start = bisect.bisect_left(self.guids, search_string)
        retVal = []
        
        for i in range(max_returned+1):
            search_point = i + search_start
            if search_point > len(self.guids)-1:      # got to the end
                break
            else:
                if self.guids[i].startswith(search_string):
                    retVal.append(self.guids[i])
            
        if len(retVal) == max_returned and return_subset == False:
            # don't return partial lists of matches
            return []
        else:
            return retVal

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
        self.assertEqual(retVal, [])
        retVal = gs.search('b', max_returned=1, return_subset=True)
        self.assertEqual(len(retVal),1)
        retVal = gs.search('z')
        self.assertEqual(retVal, [])

class test_gm_2(unittest.TestCase):
    def runTest(self):

        gs = guidSearcher()
        gs.add('b1')
         
        retVal = gs.search('b1')
        self.assertEqual(retVal, ['b1'])
               
@unittest.skip("benchmark disabled")       
class test_gm_benchmark(unittest.TestCase):
    def runTest(self):
        """ tests time to search a list of 1M guids
        
        median addition time : 0.00026 sec (0.26 ms)
        median search time   : 0.00002 sec (0.02 ms)"""

        gs = guidSearcher()
        
        print("generating list of guids")
        t1=  datetime.datetime.now()
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