""" tests startupshutdown.py

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@ukhsa.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

"""
import unittest
import psutil
import os


class Test_startup_shutdown(unittest.TestCase):
    """tests startup and shutdown"""

    def count_processes(self, match_string):
        """returns a list of processes matching a string, such as a port number"""
        servers = {}
        for proc in psutil.process_iter():

            cmdline_parts = proc.cmdline()
            for i, cmdline_part in enumerate(proc.cmdline()):
                if str(match_string) in cmdline_part:
                    servers[proc.pid] = cmdline_parts
        return len(servers.keys()), servers

    def runTest(self):

        # we can't test startstop with our usual default_test_config, because this automatically
        # starts and stops catwalk servers
        shutdown_all = "./fn4_shutdown.sh config/default_test_config_startstop.json"
        shutdown_partial = "./fn4_shutdown.sh config/default_test_config_startstop.json --leave_catwalk_running"
        startup_all = "./fn4_startup.sh config/default_test_config_startstop.json"

        # A shut anything running down
        os.system(shutdown_all)
        n_catwalk, catwalk_servers = self.count_processes("CatWalk-PORT-5998")
        self.assertEqual(n_catwalk, 0)
        n_fn4, fn4_servers = self.count_processes(
            "config/default_test_config_startstop.json"
        )
        self.assertEqual(n_fn4, 0)

        # B start it up
        os.system(startup_all)
        n_catwalk, catwalk_servers = self.count_processes("CatWalk-PORT-5998")
        self.assertEqual(n_catwalk, 1)
        n_fn4, fn4_servers = self.count_processes(
            "config/default_test_config_startstop.json"
        )
        self.assertGreater(n_fn4, 0)

        # C shut anything running down
        os.system(shutdown_all)
        n_catwalk, catwalk_servers = self.count_processes("CatWalk-PORT-5998")
        self.assertEqual(n_catwalk, 0)
        n_fn4, fn4_servers = self.count_processes(
            "config/default_test_config_startstop.json"
        )
        self.assertEqual(n_fn4, 0)

        # D start it up
        os.system(startup_all)
        n_catwalk, catwalk_servers = self.count_processes("CatWalk-PORT-5998")
        self.assertEqual(n_catwalk, 1)
        n_fn4, fn4_servers = self.count_processes(
            "config/default_test_config_startstop.json"
        )
        self.assertGreater(n_fn4, 0)

        # E shut anything running down, except catwalk
        os.system(shutdown_partial)
        n_catwalk, catwalk_servers = self.count_processes("CatWalk-PORT-5998")
        self.assertEqual(n_catwalk, 1)
        n_fn4, fn4_servers = self.count_processes(
            "config/default_test_config_startstop.json"
        )
        self.assertEqual(n_fn4, 0)

        # F start it up
        os.system(startup_all)
        n_catwalk, catwalk_servers = self.count_processes("CatWalk-PORT-5998")
        self.assertEqual(n_catwalk, 1)
        n_fn4, fn4_servers = self.count_processes(
            "config/default_test_config_startstop.json"
        )
        self.assertGreater(n_fn4, 0)

        # G shut anything running down
        os.system(shutdown_all)
        n_catwalk, catwalk_servers = self.count_processes("CatWalk-PORT-5998")
        self.assertEqual(n_catwalk, 0)
        n_fn4, fn4_servers = self.count_processes(
            "config/default_test_config_startstop.json"
        )
        self.assertEqual(n_fn4, 0)
