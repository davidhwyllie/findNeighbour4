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

"""

import unittest

# import your test modules
from findn import pycw_client
from findn import hybridComparer
from findn import seqComparer
from findn import preComparer
from findn import mongoStore
from findn import msa
from findn import msaviewer
from findn import NucleicAcid
from findn import visualiseNetwork
from findn import read_config
from snpclusters import ma_linkage
from findn import guidLookup
from findn import depictStatus
from findn import identify_sequence_set
from snpclusters import clusternomenclature
from findn import common_utils

# initialize the test suite
loader = unittest.TestLoader()
suite  = unittest.TestSuite()

# add tests to the test suite
suite.addTests(loader.loadTestsFromModule(pycw_client))
suite.addTests(loader.loadTestsFromModule(hybridComparer))
suite.addTests(loader.loadTestsFromModule(seqComparer))
suite.addTests(loader.loadTestsFromModule(preComparer))
suite.addTests(loader.loadTestsFromModule(NucleicAcid))
suite.addTests(loader.loadTestsFromModule(mongoStore))
suite.addTests(loader.loadTestsFromModule(msa))
suite.addTests(loader.loadTestsFromModule(msaviewer))
suite.addTests(loader.loadTestsFromModule(common_utils))
suite.addTests(loader.loadTestsFromModule(visualiseNetwork))

suite.addTests(loader.loadTestsFromModule(ma_linkage))
suite.addTests(loader.loadTestsFromModule(identify_sequence_set))
suite.addTests(loader.loadTestsFromModule(guidLookup))
suite.addTests(loader.loadTestsFromModule(depictStatus))
suite.addTests(loader.loadTestsFromModule(clusternomenclature))


# initialize a runner, pass it your suite and run it
runner = unittest.TextTestRunner(verbosity=3)
result = runner.run(suite)