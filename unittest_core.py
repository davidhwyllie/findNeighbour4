""" runs unittest for the findneighbour4 system 

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

Some tests, such as test_server and test_fn4client, require a server running.
The unittest_all.sh script will launch and shut down such.

"""

import unittest

from test import test_persistence
from test import test_depiction
from test import test_pycw_client
from test import test_common_utils
from test import test_guid_lookup
from test import test_hybridcomparer
from test import test_identify_sequence_set
from test import test_seqComparer
from test import test_preComparer
from test import test_msa
from test import test_msaviewer
from test import test_nucleicacid 
from test import test_visualisenetwork
from test import test_malinkage
from test import test_clusternomenclature
from test import test_server
from test import test_fn4client

# initialize the test suite
loader = unittest.TestLoader()
suite  = unittest.TestSuite()


# add tests to the test suite

suite.addTests(loader.loadTestsFromModule(test_persistence))
suite.addTests(loader.loadTestsFromModule(test_depiction))
suite.addTests(loader.loadTestsFromModule(test_pycw_client))
suite.addTests(loader.loadTestsFromModule(test_common_utils))
suite.addTests(loader.loadTestsFromModule(test_guid_lookup))
suite.addTests(loader.loadTestsFromModule(test_hybridcomparer))
suite.addTests(loader.loadTestsFromModule(test_identify_sequence_set))
suite.addTests(loader.loadTestsFromModule(test_seqComparer))
suite.addTests(loader.loadTestsFromModule(test_preComparer))
suite.addTests(loader.loadTestsFromModule(test_msa))
suite.addTests(loader.loadTestsFromModule(test_msaviewer))
suite.addTests(loader.loadTestsFromModule(test_nucleicacid))
suite.addTests(loader.loadTestsFromModule(test_visualisenetwork))
suite.addTests(loader.loadTestsFromModule(test_malinkage))
suite.addTests(loader.loadTestsFromModule(test_clusternomenclature))
suite.addTests(loader.loadTestsFromModule(test_server))
suite.addTests(loader.loadTestsFromModule(test_fn4client))


# initialize a runner, pass it your suite and run it
runner = unittest.TextTestRunner(verbosity=3)
result = runner.run(suite)

## TODO: start up server, run integration tests


## this then completes CI
