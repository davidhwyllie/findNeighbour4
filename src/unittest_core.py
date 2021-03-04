# tests/runner.py
import unittest

# import your test modules
import pycw_client
import hybridComparer
import seqComparer
import preComparer
import mongoStore
import msa
import msaviewer
import NucleicAcid
import visualiseNetwork
import read_config
import ma_linkage
import guidLookup
import depictStatus
import identify_sequence_set
import clusternomenclature

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
suite.addTests(loader.loadTestsFromModule(visualiseNetwork))

suite.addTests(loader.loadTestsFromModule(read_config))
suite.addTests(loader.loadTestsFromModule(ma_linkage))
suite.addTests(loader.loadTestsFromModule(identify_sequence_set))
suite.addTests(loader.loadTestsFromModule(guidLookup))
suite.addTests(loader.loadTestsFromModule(depictStatus))
suite.addTests(loader.loadTestsFromModule(clusternomenclature))

# initialize a runner, pass it your suite and run it
runner = unittest.TextTestRunner(verbosity=3)
result = runner.run(suite)