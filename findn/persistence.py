from findn.mongoStore import fn3persistence
from findn.rdbmsstore import fn3persistence_r


class Persistence:
    """returns either a mongodb based, or rdbms based persistence object"""

    def __init__(self):
        """factory object for returning Persistence objects"""
        pass

    def get_storage_object(
        self,
        connString: str,
        dbname: str = "fn3_unittesting",
        debug: int = 0,
        config_settings: dict = {},
        max_neighbours_per_document: int = 100000,
        server_monitoring_min_interval_msec: int = 0,
        verbose=True,
    ):
        if connString.startswith("mongodb"):
            pdm = fn3persistence(
                connString,
                dbname,
                debug,
                config_settings,
                max_neighbours_per_document,
                server_monitoring_min_interval_msec,
            )
        else:
            pdm = fn3persistence_r(
                connString, debug, server_monitoring_min_interval_msec
            )

        if verbose:
            print(
                "Set up data access object {0} called {1} **".format(
                    pdm.storage_technology, dbname
                )
            )

        return pdm
