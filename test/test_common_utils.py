""" runs unittest for common_utils

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

 

"""

import unittest
import os
import pathlib
from findn.common_utils import ConfigManager, EnvWriter


class Test_RC_1(unittest.TestCase):
    """tests the ReadConfig class"""

    def runTest(self):

        rc = ConfigManager("config/default_test_config.json")  # first run
        res = rc.read_config()
        self.assertTrue(isinstance(res, dict))

        rc = ConfigManager("config/default_test_config.json")  # second run
        res = rc.read_config()
        self.assertTrue(isinstance(res, dict))

        # check catwalk_backupdir set and exists
        if not os.path.exists(rc.catwalk_backupdir):
            self.fail("ConfigManager.catwalk_backupdir does not exist")

        # create a file in the backupdir
        backupfile = os.path.join(rc.catwalk_backupdir, 'rc_test.tar')
        pathlib.Path(backupfile).touch()
        self.assertTrue(os.path.exists(backupfile))

        rc.delete_existing_data()
        self.assertFalse(os.path.exists(backupfile))


class Test_envwriter_1(unittest.TestCase):
    """tests the EnvWriterclass"""

    def runTest(self):

        with self.assertRaises(FileNotFoundError):
            EnvWriter("no file")  # first run

        env_file = "unitTest_tmp/.env"
        if os.path.exists(env_file):
            os.unlink(env_file)

        # create a new .env file
        env_lines = ["one='one'\n", "two='two'\n"]
        with open(env_file, "wt") as f:
            f.writelines(env_lines)

        ev = EnvWriter(env_file)  # should read it
        self.assertEqual(2, ev.number_of_env_vars())

        ev.set_env_var("three", "'3'")
        self.assertEqual(3, ev.number_of_env_vars())

        ev.del_env_var("four")
        self.assertEqual(3, ev.number_of_env_vars())

        ev.del_env_var("one")
        self.assertEqual(2, ev.number_of_env_vars())

        ev.set_env_var("test", "'test'")
        self.assertEqual(3, ev.number_of_env_vars())

        ev.save_changes()
        ev = EnvWriter(env_file)  # should read it
        self.assertEqual(3, ev.number_of_env_vars())
