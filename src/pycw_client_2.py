"""
run catwalk from python.
"""
import argparse
import threading
import subprocess
import shlex
import time
import json
import requests
import logging
import unittest
import os
import psutil

class CatWalkServerDidNotStartError(Exception):
    def __init__(self, expression, message):
        self.expression= expression
        self.message =message

class CatWalk():
    """ start, stop, and communicate with a CatWalk server"""
    def __init__(self, cw_binary_filepath, instance_name, reference_filepath, mask_filepath, max_distance, cw_url='http://127.0.0.1:5000'):
        """
        Start the catwalk process in the background, if it not running.

        """

        # store parameters
        self.cw_url = cw_url        # at present, the catwalk server always runs on 5000
        self.cw_binary_filepath = cw_binary_filepath
        self.instance_name = instance_name
        self.reference_filepath = reference_filepath
        self.mask_filepath = mask_filepath
        self.max_distance = max_distance

        # start up if not running
        if not self.server_is_running():
            self.start()

        if not self.server_is_running():        # startup failed
            raise CatWalkServerDidNotStartError()

    def server_is_running(self):
        """ returns true if  response is received by the server, otherwise returns False """
        try:
            result=self.info()
            return True
        except requests.exceptions.ConnectionError: 
            return False

    def start(self):
        """ starts a catwalk process in the background """
        cw_binary_filepath = shlex.quote(self.cw_binary_filepath)
        instance_name = shlex.quote(self.instance_name)
        reference_filepath = shlex.quote(self.reference_filepath)
        mask_filepath = shlex.quote(self.mask_filepath)
        max_distance = self.max_distance

        reference_name = "CATWALK-{0}".format(self.max_distance)

        cmd = f"nohup {cw_binary_filepath} --instance_name {instance_name} --reference_name {reference_name} --reference_filepath {reference_filepath} --mask_name {mask_filepath} --mask_filepath {mask_filepath} --max_distance {max_distance} &"
        logging.info("Starting CatWalk server : {0}".format(cmd))
        
        os.system(cmd)
        time.sleep(1)       # wait for it to start
        logging.info("CatWalk server started.")

    def stop(self):
        """ tops the catwalk server, if running.  Note that all processes called 'cw_server' will be terminated """
        for proc in psutil.process_iter():
            
            if 'cw_server' in proc.name():
                #print(proc.name(),proc.pid,"Killed")

                proc.kill()

    def info(self):
        """
        Get status information from catwalk
        """
        target_url = "{0}/info".format(self.cw_url)
        return requests.get(target_url).json()

    def add_sample_from_refcomp(self, name, refcomp):
        """
        Add a reference compressed (dict with ACGTN keys and list of positions as values) sample.

        The json dict must have all keys: ACGTN, even if they're empty
        """
        requests.post("{0}/add_sample_from_refcomp".format(self.cw_url),
                      json={ "name": name,
                             "refcomp": json.dumps(refcomp),
                             "keep": True }
                     )

    def neighbours(self, name):
        """ get neighbours
        """
        r = requests.get("{0}/neighbours/{1}".format(self.cw_url, name))
        j = r.json()
        return [(sample_name, int(distance_str)) for (sample_name, distance_str) in j]


# unit testings
class test_cw(unittest.TestCase):
    """ starts server, and shuts it down """
    def setUp(self):
        """ cw_binary_filepath must point to the catwalk server and mask & reference files to the relevant data files. 
            Shuts down **any catwalk server** running initially"""
        self.cw = CatWalk(cw_binary_filepath="/home/phe.gov.uk/david.wyllie/catwalk/src/cw_server",
                    instance_name="test",
                    reference_filepath="../reference/TB-ref.fasta",
                    mask_filepath="../reference/TB-exclude-adaptive.txt",
                    max_distance=20)

        # stop the server if it is running
        self.cw.stop()
        self.assertFalse(self.cw.server_is_running())

        self.cw.start()
        self.assertTrue(self.cw.server_is_running())

    def teardown(self):
        self.cw.stop()
        pass

class test_cw_1(test_cw):
    """ tests server startup, shutdown, info(), and the server_is_running method.  
        Shuts down **any catwalk server** running initially """
    def runTest(self):

        self.cw.start()
        self.assertTrue(self.cw.server_is_running())

        self.assertIsInstance(self.cw.info(), dict)

        self.cw.stop()
        self.assertFalse(self.cw.server_is_running())

class test_cw_2(test_cw):
    """ tests insert  
         """
    def runTest(self):

        payload1 = {'A':[1000,1001,1002], 'G':[], 'T':[], 'C':[], 'N':[20000,20001,20002]}
        payload2 = {'A':[1003,1004,1005], 'G':[], 'T':[], 'C':[], 'N':[20000,20001,20002]}
        self.cw.add_sample_from_refcomp('guid1',payload1)
        self.cw.add_sample_from_refcomp('guid2',payload1)
        self.assertEqual(self.cw.neighbours('guid1'),[('guid2',0)])      
        self.assertEqual(self.cw.neighbours('guid2'),[('guid1',0)])      

def main():
    """
    You can also start it on the command-line, and then skip start()
    """
    p = argparse.ArgumentParser()
    p.add_argument("cw_binary_filepath")
    p.add_argument("instance_name")
    p.add_argument("reference_filepath")
    p.add_argument("mask_filepath")
    p.add_argument("max_distance")
    args = p.parse_args()

    start(*vars(args).values())

    time.sleep(5)
    while True:
        print(json.dumps(info(), indent=4))
        time.sleep(60)

if __name__ == "__main__":
    main()
