""" runs unittest for common_utils

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
from findn.common_utils import ConfigManager


class Test_RC_1(unittest.TestCase):
    """tests the ReadConfig class"""

    def runTest(self):

        rc = ConfigManager("config/default_test_config.json")  # first run
        res = rc.read_config()
        self.assertTrue(isinstance(res, dict))

        rc = ConfigManager("config/default_test_config.json")  # second run
        res = rc.read_config()
        self.assertTrue(isinstance(res, dict))
