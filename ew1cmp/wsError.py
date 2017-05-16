#!/usr/bin/env python3
""" error classes used by web service interface """

import unittest


class WSFailureError(Exception):
    """Exception raised  when the Web Service doesn't work as expected, but we didn't notice by raising an error during execution """
    def __init__(self, message):
        self.message = message

class test_WSFailureError(unittest.TestCase):
    """ tests raising of a custom error """
    def runTest(self):
        with self.assertRaises(WSFailureError):
            raise WSFailureError('message')

