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
import uuid

class CatWalkServerDidNotStartError(Exception):
    """ the catwalk server did not start """
    def __init__(self, expression, message):
        self.expression= expression
        self.message =message

class CatWalkBinaryNotAvailableError(Exception):
    """ no catwalk binary """
    def __init__(self, expression, message):
        self.expression= expression
        self.message =message

class CatWalk():
    """ start, stop, and communicate with a CatWalk server"""
    def __init__(self, cw_binary_filepath, reference_name, reference_filepath, mask_filepath, max_distance, bind_host="localhost", bind_port=5000):
        """
        Start the catwalk process in the background, if it not running.

        """
       
        no_catwalk_exe_message = """
The catWalk client could not find a CatWalk server application.  This is because either:
i) the cw_binary_filepath parameter was None
ii) the above, and there is no CW_BINARY_FILEPATH environment variable.  To set this, you insert a line like
CW_BINARY_FILEPATH=/path/to/cw/executable
in either
.bashrc - if you're not using a virtual environment
.env    - a file in the same directory as the PipFile, if you are using a virtual environment.
          This file is not committed into the repository, so you'll have to create it once in your installation.
"""
        # if cw_binary_filepath is None, we check the env. variable CW_BINARY_FILEPATH and use that if present.
        if cw_binary_filepath is None:
            if 'CW_BINARY_FILEPATH' in os.environ:
                cw_binary_filepath = os.environ['CW_BINARY_FILEPATH']
            else:
                raise CatWalkBinaryNotAvailableError(expression = None, message = no_catwalk_exe_message)
        if not os.path.exists(cw_binary_filepath):
                raise FileNotFoundError(expression=  None, message = "Was provided a cw_binary_filepath, but there is no file there {0}".format(cw_binary_filepath))
        # store parameters
        self.bind_host = bind_host
        self.bind_port = bind_port
        self.cw_url = "http://{0}:{1}".format(bind_host,bind_port)      
        self.cw_binary_filepath = cw_binary_filepath
        self.reference_filepath = reference_filepath
        self.mask_filepath = mask_filepath
        self.max_distance = max_distance
        self.reference_name = reference_name
        self.instance_name = "CatWalk-SNV-{0}-PORT-{1}-{2}".format(self.max_distance, self.bind_port,         str(uuid.uuid1()))
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
        except Exception as e:
            raise e      # something else when wrong

    def start(self):
        """ starts a catwalk process in the background """
        cw_binary_filepath = shlex.quote(self.cw_binary_filepath)
        instance_name = shlex.quote(self.instance_name)
        reference_filepath = shlex.quote(self.reference_filepath)
        mask_filepath = shlex.quote(self.mask_filepath)

        cmd = f"nohup {cw_binary_filepath} --instance_name {instance_name}  --bind_host {self.bind_host} --bind_port {self.bind_port} --reference_name {self.reference_name} --reference_filepath {reference_filepath} --mask_name {mask_filepath} --mask_filepath {mask_filepath} --max_distance {self.max_distance} > cw_server_nohup.out &"
        logging.info("Attempting startup of CatWalk server : {0}".format(cmd))       
        os.system(cmd)

        time.sleep(2)
        if self.info() is None:
            raise CatWalkServerDidNotStartError()
 
    def stop(self):
        """ stops the catwalk server launched by this process, if running.  The process is identified by a uuid, so only one catwalk server will be shut down. """
        for proc in psutil.process_iter():
            if 'cw_server' in proc.name():
                if self.instance_name in proc.cmdline():
                    proc.kill()
    def stop_all(self):
        """ stops all catwalk servers """
        for proc in psutil.process_iter():
            if 'cw_server' in proc.name():
                    proc.kill()
    def info(self):
        """
        Get status information from catwalk
        """
        target_url = "{0}/info".format(self.cw_url)
        r = requests.get(target_url)
        r.raise_for_status()
        return r.json()

    def add_sample_from_refcomp(self, name, refcomp):
        """
        Add a reference compressed (dict with ACGTN keys and list of positions as values) sample.

        The json dict must have all keys: ACGTN, even if they're empty
        """
        print("ABOUT TO POST AT ","{0}/add_sample_from_refcomp".format(self.cw_url))
        r = requests.post("{0}/add_sample_from_refcomp".format(self.cw_url),
                      json={ "name": name,
                             "refcomp": json.dumps(refcomp),
                             "keep": True }
                     )
        r.raise_for_status()
        return r.status_code

    def neighbours(self, name):
        """ get neighbours
        """
        r = requests.get("{0}/neighbours/{1}".format(self.cw_url, name))
        r.raise_for_status()
        j = r.json()
        return [(sample_name, int(distance_str)) for (sample_name, distance_str) in j]

    def sample_names(self):
        """ get neighbours
        """
        r = requests.get("{0}/list_samples".format(self.cw_url))
        r.raise_for_status()
        return r.json()



# unit testings
class test_cw(unittest.TestCase):
    """ starts server, and shuts it down """
    def setUp(self):
        """ cw_binary_filepath must point to the catwalk server and mask & reference files to the relevant data files. 
            Shuts down **any catwalk server** running initially.

            Note: requires CW_BINARY_FILEPATH environment variable to point to the catwalk binary."""
        self.cw = CatWalk(
                    cw_binary_filepath=None,
                    reference_name="H37RV",
                    reference_filepath="../reference/TB-ref.fasta",
                    mask_filepath="../reference/TB-exclude-adaptive.txt",
                    max_distance=20)

        # stop the server if it is running
        self.cw.stop_all()
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
        res = self.cw.add_sample_from_refcomp('guid1',payload1)
        self.assertEqual(res, 200)

        res = self.cw.add_sample_from_refcomp('guid2',payload1)
        self.assertEqual(res, 200)

        res = self.cw.add_sample_from_refcomp('guid2',payload1)     # insert twice
        self.assertEqual(res, 201)

        self.assertEqual(self.cw.neighbours('guid1'),[('guid2',0)])      
        self.assertEqual(self.cw.neighbours('guid2'),[('guid1',0)])      

class test_cw_3(test_cw):
    """ tests list_samples
         """
    def runTest(self):

        payload1 = {'A':[1000,1001,1002], 'G':[], 'T':[], 'C':[], 'N':[20000,20001,20002]}
        payload2 = {'A':[1003,1004,1005], 'G':[], 'T':[], 'C':[], 'N':[20000,20001,20002]}
        res = self.cw.add_sample_from_refcomp('guid1',payload1)
        res = self.cw.add_sample_from_refcomp('guid2',payload1)

        self.assertEqual(set(self.cw.sample_names()), set(['guid1','guid2']))	# order doesn't matter      

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
