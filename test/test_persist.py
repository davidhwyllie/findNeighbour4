""" tests Persistence class, which returns either a 
    monogo or rdbms based fn3persistence object.

"""

import unittest
from findn.persistence import Persistence
from findn.mongoStore import fn3persistence
from findn.rdbmsstore import fn3persistence_r


class Test_Persistence(unittest.TestCase):
    """tests persistence class"""

    def setUp(self):
        self.engines = {}
        self.engines["Sqlite"] = "sqlite://"  # in memory sqlite
        self.engines["mongo"] = "mongodb://localhost"  # in memory sqlite

    def runTest(self):
        """yields fn3persistence objects, one for each database server being tested."""
        for engine, config in self.engines.items():
            print(engine, type(self).__name__)

            sf = Persistence()
            pdm = sf.get_storage_object(connString=config, debug=2)

            if engine == "mongo":
                self.assertEqual(pdm.storage_technology, "mongodb")
                self.assertIsInstance(pdm, fn3persistence)
            else:
                self.assertEqual(pdm.storage_technology, "rdbms")
                self.assertIsInstance(pdm, fn3persistence_r)

            # explicitly close connection (required for unittesting)
            pdm.closedown()
