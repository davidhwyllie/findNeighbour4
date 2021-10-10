#!/usr/bin/env python
"""  findNeighbour4 is a  server providing relatedness information for bacterial genomes via a Restful API.
See documentation for full details of its functionality.

There are unit tests for the lockmanager component.  To run them:

# And then  launching unit tests with
pipenv run pytest test/test_lockmanager.py              # tests the server running on the default testing port, 5020

"""

# import libraries
import os
import time
import unittest
from findn.common_utils import ConfigManager

class test_lockmanager_1(unittest.TestCase):
    """tests whether the lockmanager runs when no lock is in place"""

    def runTest(self):

        # run the lockmonitor
        config_file = os.path.join("config", "default_test_config.json")

        # empty everything out
        cfm = ConfigManager(config_file)
        cfm.read_config()

        # start test
        cfm = ConfigManager(config_file)
        cfm.read_config(not_debug_mode=True)

        # unlock
        cfm.PERSIST.unlock(1, force=True)

        # there is no lock
        lock_details = cfm.PERSIST.lock_details(1)
        self.assertIsNone(lock_details)

        os.system(
            "pipenv run python3 findNeighbour4_lockmanager.py --max_run_time 10 --run_once_only"
        )

        # check the lock is still not held
        lock_details = cfm.PERSIST.lock_details(1)
        self.assertIsNone(lock_details)


class test_lockmanager_2(unittest.TestCase):
    """tests whether the lockmanager runs and releases a lock running lock"""

    def runTest(self):

        # run the lock monitor
        config_file = os.path.join("config", "default_test_config.json")

        # empty everything out
        cfm = ConfigManager(config_file)
        cfm.read_config()

        cfm = ConfigManager(config_file)
        cfm.read_config(not_debug_mode=True)

        # unlock
        cfm.PERSIST.unlock(1, force=True)

        # there is no lock
        lock_details = cfm.PERSIST.lock_details(1)
        self.assertIsNone(lock_details)

        # take a lock
        cfm.PERSIST.lock(1, "guid1")

        # check the lock is held
        lock_details = cfm.PERSIST.lock_details(1)
        self.assertIsNotNone(lock_details)

        os.system(
            "pipenv run python3 findNeighbour4_lockmanager.py --max_run_time 10 --run_once_only"
        )

        # check the lock is still held
        lock_details = cfm.PERSIST.lock_details(1)
        self.assertIsNotNone(lock_details)

        # unlock
        cfm.PERSIST.unlock(1)

        # there is no lock
        lock_details = cfm.PERSIST.lock_details(1)
        self.assertIsNone(lock_details)

        # take a lock
        cfm.PERSIST.lock(1, "guid1")
        lock_details = cfm.PERSIST.lock_details(1)
        self.assertIsNotNone(lock_details)

        # wait 15 seconds
        time.sleep(15)
        os.system(
            "pipenv run python3 findNeighbour4_lockmanager.py --max_run_time 10 --run_once_only"
        )

        # there is no lock
        lock_details = cfm.PERSIST.lock_details(1)
        self.assertIsNone(lock_details)
