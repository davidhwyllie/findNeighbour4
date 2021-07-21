
import os
import findn.mongoStore as mongoStore
import findn.rdbmsstore as rdbmsstore

def createPersistence(
    connString: str,
    dbname: str = "fn3_unittesting",
    debug: int = 0,
    config_settings: dict = {},
    max_neighbours_per_document: int = 100000,
    server_monitoring_min_interval_msec: int = 0,
    ):
    if connection_config := os.environ.get("FINDNEIGHBOUR_USE_RDBMS"):
        return rdbmsstore.findNeighbour4(connection_config, debug, server_monitoring_min_interval_msec)
    else:
        return mongoStore.findNeighbour4(connString, dbname, debug, config_settings, max_neighbours_per_document, server_monitoring_min_interval_msec)

