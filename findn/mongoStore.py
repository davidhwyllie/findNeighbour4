#!/usr/bin/env python
""" fnPersistence, a class which provides a storage layer for meta-data and snv distances from the
findneighbour4 system in mongodb

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

 
 """


import bson  # type: ignore
import datetime
import json
import uuid
import pandas as pd  # type: ignore
import logging
import pymongo  # type: ignore
import gridfs  # type: ignore
import pickle
import psutil  # type: ignore
import io
import statistics
import numpy as np
from typing import (
    Any,
    Dict,
    Iterable,
    List,
    NoReturn,
    Optional,
    Set,
    Tuple,
    TypedDict,
    Union,
)

Guid2NeighboursFormat1 = List[Union[str, int]]
Guid2NeighboursFormat3 = Union[str]
Guid2NeighboursFormat4 = Dict[str, Union[str, int]]
Guid2NeighboursFormats = Union[
    Guid2NeighboursFormat1, Guid2NeighboursFormat3, Guid2NeighboursFormat4
]


class RecentDatabaseMonitoringRet(TypedDict, total=False):
    recompression_data: bool
    latest_stats: Dict[str, Union[int, np.float64]]
    trend_stats: List[Dict[str, Any]]


class NPEncoder(json.JSONEncoder):
    """encodes Numpy types as jsonisable equivalents"""

    def default(self, obj: Any) -> Any:
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NPEncoder, self).default(obj)


class fn3persistence:
    """System for persisting results from  large numbers of sequences stored in FindNeighbour 3+.
    Uses Mongodb.

    in the current schema there are the following collections:
    -'config', containing configuration information

    -'refcompressedseq' contains reference compressed sequences. note that this is a gridfs 'filesystem'.  Keys are guids.

    -'clusters' contains a graph of inter-guid links. note that this is a gridfs 'filesystem'.  Keys are names of clustering algorithms.

    -'guid2meta', contains guid -> metadata

    -'guid2neighbour', contains links between guids, including snv
        Here, individuals documents are identified by mongo assigned unique ids.
        Each document contains three keys:
        {'guid':'a1234', 'rstat':'s', 'neighbours':{}}

        Up to max_neighbours_per_document neighbours can be stored per document.

    *max_neighbours_per_document* should be less than 5,000, because there is a max. document size in mongodb.
    In debug mode, it is automatically set to 3.
    if max_neighbours_per_document exist in the document, 'rstat' is set to 'f' (full).
    If there is a single item only, rstat is set to 's' (single); if there are multiple items, it is set to 'm'.

    Indices exist on (i) guid - allowing you to find all the documents contains guid X's neighbours and
                     (ii) guid/rstat combination- allowing one to find guid X's most recent document, useful for addition.

    This class provides methods to access these four entities.

    NOTE:  regarding sharding, the most important collection is guid2neighbour.
    A hashed sharding based on guid should work well when ensuring database scalability.

    """

    # code handling startup and shutdown.
    def __init__(
        self,
        connString: str,
        dbname: str = "fn3_unittesting",
        debug: int = 0,
        config_settings: dict = {},
        max_neighbours_per_document: int = 100000,
        server_monitoring_min_interval_msec: int = 0,
    ) -> None:
        """Creates a connection to a MongoDb database.

        connString : the mongoDb connection string
        dbname: the name of the mongoDb database to use.
        if debug = 0 or 1, the database is opened or created.
        if debug = 2, any existing collections are deleted.
        config_settings: only used on db creation; optional dictionary to note items in the database's config collection.
        """

        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)

        # client calling mongostore should trap for connection errors etc
        self.connString = connString
        self.dbname = dbname
        self.debug = debug
        self._connect()  # will raise ConnectionError if fails

        # can check what exists with connection.database_names()
        self.expected_collections = [
            "server_monitoring",
            "guid2meta",
            "guid2neighbour",
            "config",
            "refcompressedseq.chunks",
            "refcompressedseq.files",
            "clusters.chunks",
            "clusters.files",
            "msa.chunks",
            "msa.files",
            "tree.chunks",
            "tree.files",
            "fnlock",
        ]
        self.expected_clustering_collections = [
            "clusters.chunks",
            "clusters.files",
            "msa.chunks",
            "msa.files",
            "tree.chunks",
            "tree.files",
        ]
        self.storage_technology = "mongodb"
        self.using_sqlite = False
        self.max_neighbours_per_document = max_neighbours_per_document
        self.server_monitoring_min_interval_msec = server_monitoring_min_interval_msec
        self.previous_server_monitoring_data: Dict[str, Any] = {}
        self.previous_server_monitoring_time = None

        # delete any pre-existing data if we are in debug mode.
        if debug == 2:
            self.logger.warning(
                "Debug mode operational [DEBUG={0}]; deleting all data from collections.".format(
                    debug
                )
            )
            self._delete_existing_clustering_data()
            self._delete_existing_data()

            self.max_neighbours_per_document = 3  # used for unittests
        else:
            self.logger.info("Using stored data in mongostore")

        # create indices on guid2neighbours; note will do nothing if index already exists
        ix1 = pymongo.IndexModel(
            [("guid", pymongo.ASCENDING), ("rstat", pymongo.ASCENDING)],
            name="by_guid_full",
        )
        ix2 = pymongo.IndexModel([("rstat", pymongo.ASCENDING)], name="by_rstat")

        self.db["guid2neighbour"].create_indexes([ix1, ix2])

        # create indices on msa and trees; note will do nothing if index already exists
        ix3 = pymongo.IndexModel(
            [("filename", pymongo.ASCENDING), ("uploadDate", pymongo.ASCENDING)],
            name="filename_date",
        )
        self.db["msa.files"].create_indexes([ix3])

        ix3b = pymongo.IndexModel(
            [("filename", pymongo.ASCENDING), ("uploadDate", pymongo.ASCENDING)],
            name="filename_date",
        )
        self.db["tree.files"].create_indexes([ix3b])

        # create indices on guid2meta, allowing recovery of valid and invalid specimens rapidly.
        ix4 = pymongo.IndexModel(
            [('"sequence_meta.DNAQuality.invalid"', pymongo.ASCENDING)],
            name="guid_validity",
        )
        ix5 = pymongo.IndexModel(
            [('"sequence_meta.DNAQuality.propACTG"', pymongo.ASCENDING)],
            name="guid_quality",
        )
        ix5b = pymongo.IndexModel(
            [('"sequence_meta.DNAQuality.examinationDate"', pymongo.ASCENDING)],
            name="examinationDate",
        )

        # note: if additional metadata is added, such as sequence names etc which might be searched for, then we need to add additional indices here.

        self.db["guid2meta"].create_indexes([ix4, ix5, ix5b])

        # create index on server_monitoring insert times
        ix6 = pymongo.IndexModel([("context|time|time_now", pymongo.ASCENDING)])
        self.db["server_monitoring"].create_indexes([ix6])

    def delete_server_monitoring_entries(self, before_seconds: int) -> None:
        """deletes server monitoring entries more than before_seconds ago"""
        now = datetime.datetime.now()
        earliest_allowed = now - datetime.timedelta(seconds=before_seconds)

        earliest_allowed_str = str(earliest_allowed.isoformat())
        self.db["server_monitoring"].delete_many(
            {"context|time|time_now": {"$lt": earliest_allowed_str}}
        )

    def summarise_stored_items(self) -> Dict[str, Any]:
        """counts how many sequences exist of various types"""
        retVal = {}
        collections_present = self.db.list_collection_names()
        for this_collection in self.expected_collections:
            if this_collection in collections_present:
                res = self.db.command("collstats", this_collection)
                for relevant_metric in [
                    "totalIndexSize",
                    "storageSize",
                    "count",
                    "avgObjSize",
                ]:
                    if relevant_metric in res.keys():
                        target_key = "dstats|{0}|{1}".format(
                            this_collection.replace(".", "-"), relevant_metric
                        )
                        retVal[target_key] = res[relevant_metric]
        return retVal

    def connect(self) -> None:
        """test whether the database is connected, and if not, tries to connect.
        if the connection fails, raises pymongo.errors.ConnectionFailure"""
        if not self.is_connected():
            self._connect()

    def _connect(self) -> None:
        """connect to the database"""

        # try to close any existing session, if it exists
        self.closedown()

        # open new client
        self.client = pymongo.MongoClient(self.connString, retryWrites=True)
        self.db = self.client[self.dbname]

        # open gridfs systems
        self.rcs = gridfs.GridFS(self.db, collection="refcompressedseq")
        self.clusters = gridfs.GridFS(self.db, collection="clusters")
        self.monitor = gridfs.GridFS(self.db, collection="monitor")
        self.msa = gridfs.GridFS(self.db, collection="msa")
        self.tree = gridfs.GridFS(self.db, collection="tree")

        # enable sharding at database level
        # self.client.admin.command('enableSharding', self.dbname)

    def is_connected(self) -> bool:
        """Tests whether db is connected cf
        http://api.mongodb.com/python/current/api/pymongo/mongo_client.html"""
        try:
            # The ismaster command is cheap and does not require auth.
            self.client.admin.command("ismaster")
            # success
            return True
        except pymongo.errors.ConnectionFailure:
            return False

    def rotate_log(self) -> None:
        """forces rotation of the mongo log file"""
        self.client.admin.command("logRotate")

    def raise_error(self, token: str) -> NoReturn:
        """raises a ZeroDivisionError, with token as the message.
        useful for unit tests of error logging"""
        raise ZeroDivisionError(token)

    def _delete_existing_data(self) -> None:
        """deletes existing data from the databases"""
        for collection in self.expected_collections:
            self.db[collection].delete_many({})

    def _delete_existing_clustering_data(self) -> None:
        """deletes any clustering data from the databases"""
        for collection in self.expected_clustering_collections:
            self.db[collection].delete_many({})

    def first_run(self) -> bool:
        """if there is no config entry, it is a first-run situation"""
        if self.db.config.find_one({"_id": "config"}) is None:
            return True
        else:
            return False

    def __del__(self) -> None:
        """closes any existing session"""
        self.closedown()

    def closedown(self) -> None:
        """closes any session"""
        # client object has already been destroyed on reaching here
        pass

    # generic routines to handle insertion and read from standard mongodb stores
    def _store(self, collection: str, key: str, object: Dict[str, Any]) -> Any:
        """stores key:object in collection. It is assumed object is a dictionary.  Updates if appropriate."""
        if not isinstance(object, dict):
            raise TypeError(" object{0}  passed must be a dictionary".format(object))
        object["_id"] = key
        res = self.db[collection].replace_one({"_id": key}, object, upsert=True)
        if res.acknowledged is not True:
            raise IOError(
                "Mongo {0} did not acknowledge write of data: {1}".format(
                    self.db, object
                )
            )
        return res

    def _load(self, collection: str, key: str) -> Any:
        """loads object from collection[key]"""
        return self.db[collection].find_one({"_id": key})

    def _load_ids(self, collection: str) -> Set[str]:
        """loads guids from collection"""
        retVal: Set[str] = set()
        for item in self.db[collection].find({}):
            retVal.add(item["_id"])
        return retVal

    def memory_usage(self) -> Dict[str, Union[int, float]]:
        """returns memory usage by current python3 process
        Uses the psutil module, as the resource module is not available in windows.
        """
        memdict = psutil.virtual_memory()._asdict()
        sm = {"server|mstat|" + k: v for k, v in memdict.items()}
        return sm

    # methods for the config collection
    def config_store(self, key: str, object: Dict[str, Any]) -> Any:
        """stores object into config collection
        It is assumed object is a dictionary
        """
        # if "excluded" in object.keys():
        #    del object["excluded"]
        return self._store("config", key, object)

    def config_read(self, key: str) -> Any:
        """loads object from config.
        It is assumed object is a dictionary"""

        return self._load("config", key)

    # methods for the server and database monitoring
    def recent_database_monitoring(
        self, max_reported: int = 100
    ) -> RecentDatabaseMonitoringRet:
        """computes trends in the number of records holding pairs (guid2neighbours) vs. records.
        This ratio is a measure of database health.  Ratios > 100 indicate the database may become very large, and query slowly"""

        db_data = self.recent_server_monitoring(
            selection_field="content|activity|whatprocess",
            selection_string="dbManager",
            max_reported=max_reported,
        )

        res_df = pd.DataFrame.from_dict(db_data)

        retDict: RecentDatabaseMonitoringRet

        if len(res_df.index > 0):
            res_df["storage_ratio"] = res_df["dstats|guid2neighbour|count"] / (
                1 + res_df["dstats|guid2meta|count"]
            )
            res_df["context|time|time_now_dt"] = pd.to_datetime(
                res_df["context|time|time_now"]
            )
            res_df["latest_time"] = max(res_df["context|time|time_now_dt"])
            res_df["interval_seconds"] = (
                res_df["context|time|time_now_dt"] - res_df["latest_time"]
            ).dt.total_seconds()

            desired_cols = set(
                [
                    "_id",
                    "storage_ratio",
                    "dstats|guid2neighbour|count",
                    "dstats|guid2meta|count",
                    "context|time|time_now",
                    "interval_seconds",
                ]
            )
            available_cols = set(res_df.columns.to_list())
            select_cols = desired_cols.intersection(available_cols)
            res_df = res_df[list(select_cols)]  # select what we want
        if len(res_df.index > 0):
            retDict = {
                "recompression_data": True,
                "latest_stats": {"storage_ratio": res_df.at[0, "storage_ratio"]},
                "trend_stats": res_df.to_dict(orient="records"),
            }
        else:
            retDict = {
                "recompression_data": False,
                "latest_stats": {"storage_ratio": 1},
            }  # if there's no data, record as 1 (optimal)

        # store the ratio as 1 if we can't compute it
        if "dstats|guid2meta|count" not in res_df.columns.tolist():
            retDict["latest_stats"]["storage_ratio"] = 1
        elif res_df.at[0, "dstats|guid2meta|count"] == 0:
            retDict["latest_stats"]["storage_ratio"] = 1

        return retDict

    def recent_server_monitoring(
        self,
        max_reported: int = 100,
        selection_field: Optional[str] = None,
        selection_string: Optional[str] = None,
    ) -> List[dict]:
        """returns a list containing recent server monitoring, in reverse order (i.e. tail first).
        The _id field is an integer reflecting the order added.  Lowest numbers are most recent.

        Inputs
        max_reported - return this number of lines, at most.
        selection_field - if not None, will only return lines containing selection_string
                          in the 'selection_field' key of the returned dictionary.
        selection_string -if selection_field is not None, only returns rows if
                          selection_string is present in the 'selection_field' key of the
                          monitoring element. If None, this constraint is ignored.
        """

        if not isinstance(max_reported, int):
            raise TypeError(
                "limit must be an integer, but it is a {0}".format(type(max_reported))
            )
        if not max_reported >= 0:
            raise ValueError("limit must be more than or equal to zero")

        if max_reported == 0:
            return []

        n = 0
        retVal = []

        if selection_field is None:
            formerly_cursor = (
                self.db["server_monitoring"]
                .find({})
                .sort("_id", pymongo.DESCENDING)
                .limit(max_reported)
            )
        else:
            formerly_cursor = (
                self.db["server_monitoring"]
                .find({selection_field: selection_string})
                .sort("_id", pymongo.DESCENDING)
                .limit(max_reported)
            )

        for formerly in formerly_cursor:
            n += 1
            formerly["_id"] = n
            retVal.append(formerly)

        return retVal

    def server_monitoring_store(
        self,
        message: str = "No message provided",
        what: Optional[str] = None,
        guid: Optional[str] = None,
        content: Dict[str, Any] = {},
    ) -> bool:
        """stores content, a dictionary, into the server monitoring log"""
        now = dict(**content)
        if what is not None:
            now["content|activity|whatprocess"] = what
        if guid is not None:
            now["content|activity|guid"] = guid
        now["context|info|message"] = message
        current_time = datetime.datetime.now()
        now["context|time|time_now"] = str(current_time.isoformat())
        now["context|time|time_boot"] = datetime.datetime.fromtimestamp(
            psutil.boot_time()
        ).strftime("%Y-%m-%d %H:%M:%S")

        # should we write this data?  We have the option not to log all messages, to prevent the store getting very full.
        write_content = False
        if self.previous_server_monitoring_time is None:
            write_content = True  # yes if this is the first record written.
        else:
            time_since_last_write = (
                current_time - self.previous_server_monitoring_time
            )  # yes if it's after the server_moni
            t = (
                1000 * float(time_since_last_write.seconds)
                + float(time_since_last_write.microseconds) / 1000
            )
            if t >= self.server_monitoring_min_interval_msec:
                write_content = True

        if write_content:

            self.db["server_monitoring"].insert_one(now)
            self.previous_server_monitoring_time = current_time
            self.previous_server_monitoring_data = now
            return True
        else:
            return False

    # methods for monitor, which store the contents of an html file
    # in a gridFS store.
    def monitor_store(self, monitoring_id: str, html: str) -> str:
        """stores the monitor output string html.  Overwrites any prior object."""
        self.monitor.delete(monitoring_id)
        with io.BytesIO(html.encode("utf-8")) as f:
            id = self.monitor.put(f, _id=monitoring_id, filename=monitoring_id)
            return id

    def monitor_read(self, monitoring_id: str) -> Optional[str]:
        """loads stored string (e.g. html object) from the monitor collection."""
        try:
            res = self.monitor.get(monitoring_id)
        except gridfs.errors.NoFile:
            return None

        if res is None:
            return None
        else:
            return res.read().decode("utf-8")

    # methods for multisequence alignments
    def msa_store(self, msa_token: str, msa: dict) -> Optional[str]:
        """stores the msa object msa under token msa_token."""

        if not isinstance(msa, dict):
            raise TypeError(
                "Can only store dictionary objects, not {0}".format(type(msa))
            )

        res = self.msa.find_one({"_id": msa_token})
        if res is None:
            json_repr = json.dumps(msa).encode("utf-8")
            with io.BytesIO(json_repr) as f:
                self.msa.put(f, _id=msa_token, filename=msa_token)

            return msa_token
        else:
            return None

    def msa_read(self, msa_token: str) -> Optional[dict]:
        """loads object from msa collection.
        It is assumed object is a dictionary"""

        res = self.msa.find_one({"_id": msa_token})
        if res is None:
            return None
        json_repr = json.loads(res.read().decode("utf-8"))
        return json_repr

    def msa_delete(self, msa_token: str) -> None:
        """deletes the msa with token msa_token"""

        self.msa.delete(msa_token)

    def msa_stored_ids(self) -> List[str]:
        """returns a list of msa tokens of all objects stored"""
        return [stored_msa._id for stored_msa in self.msa.find({})]

    def msa_delete_unless_whitelisted(self, whitelist: Iterable[str]) -> None:
        """deletes the msa unless the id is in whitelist"""
        to_delete: Set[str] = set()
        for id in self.msa_stored_ids():
            if id not in whitelist:
                to_delete.add(id)
        for msa_token in to_delete:
            self.msa.delete(msa_token)

    # methods for trees

    def tree_store(self, tree_token: str, tree: dict) -> Optional[str]:
        """stores the tree object tree under token tree_token."""

        if not isinstance(tree, dict):
            raise TypeError(
                "Can only store dictionary objects, not {0}".format(type(tree))
            )

        res = self.tree.find_one({"_id": tree_token})
        if res is None:
            json_repr = json.dumps(tree).encode("utf-8")
            with io.BytesIO(json_repr) as f:
                self.tree.put(f, _id=tree_token, filename=tree_token)

            return tree_token
        else:
            return None

    def tree_read(self, tree_token: str) -> Optional[dict]:
        """loads object from tree collection.
        It is assumed object is a dictionary"""

        res = self.tree.find_one({"_id": tree_token})
        if res is None:
            return None
        json_repr = json.loads(res.read().decode("utf-8"))
        return json_repr

    def tree_delete(self, tree_token: str) -> None:
        """deletes the tree with token tree_token"""

        self.tree.delete(tree_token)

    def tree_stored_ids(self) -> List[str]:
        """returns a list of tree tokens of all objects stored"""
        return [stored_tree._id for stored_tree in self.tree.find({})]

    def tree_delete_unless_whitelisted(self, whitelist: Iterable[str]) -> None:
        """deletes the tree unless the id is in whitelist"""
        to_delete: Set[str] = set()
        for id in self.tree_stored_ids():
            if id not in whitelist:
                to_delete.add(id)
        for tree_token in to_delete:
            self.tree.delete(tree_token)

    # methods for clusters
    def cluster_store(self, clustering_key: str, obj: dict) -> str:
        """stores the clustering object obj.  retains previous version.  To clean these up, call cluster_delete_legacy.

        obj: a dictionary to store
        clustering_key: the name of the clustering, e.g. TBSNP12-graph

        Returns:
        current cluster version

        Note; does not replace previous version, but stores a new one.

        cf. warning in Mongo docs:
        Do not use GridFS if you need to update the content of the entire file atomically.
        As an alternative you can store multiple versions of each file and specify the current version of the file in the metadata.
        You can update the metadata field that indicates “latest” status in an atomic update after uploading the new version of the file,
        and later remove previous versions if needed.
        """

        if not isinstance(obj, dict):
            raise TypeError(
                "Can only store dictionary objects, not {0}".format(type(obj))
            )
        json_repr = json.dumps(obj, cls=NPEncoder).encode("utf-8")

        with io.BytesIO(json_repr) as f:
            id = self.clusters.put(f, filename=clustering_key)
            return id  # this is the current cluster version

    def cluster_read(self, clustering_key: str) -> Optional[dict]:
        """loads object from clusters collection corresponding to the most recent version of
        the clustering, saved with filename = 'clustering_key'.
        """

        cursor = (
            self.clusters.find({"filename": clustering_key})
            .sort("uploadDate", -1)
            .limit(1)
        )
        for res in cursor:
            json_repr = json.loads(res.read().decode("utf-8"))
            return json_repr
        # nothing there
        return None

    def cluster_read_update(
        self, clustering_key: str, current_cluster_version: bson.objectid.ObjectId
    ) -> Optional[dict]:
        """loads object from clusters collection corresponding to the most recent version
        of the clustering, saved with filename = 'clustering_key'.
        it will read only if the current version is different from current_cluster_version; other wise, it returns None
        It is assumed object is a dictionary"""
        latest_version = self.cluster_latest_version(clustering_key)
        if latest_version == current_cluster_version:
            # no update
            return None
        else:
            return self.cluster_read(clustering_key)

    def cluster_latest_version(self, clustering_key: str) -> bson.objectid.ObjectId:
        """returns id of latest version"""
        cursor = (
            self.clusters.find({"filename": clustering_key})
            .sort("uploadDate", -1)
            .limit(1)
        )
        for res in cursor:
            return res._id
        return None

    def cluster_keys(self, clustering_name: Optional[str] = None) -> List[str]:
        """lists  clustering keys beginning with clustering_name.  If clustering_name is none, all clustering keys are returned."""

        cursor = self.clusters.find({})
        filenames = set()
        retVal = []
        for res in cursor:
            filenames.add(res.filename)

        if (
            clustering_name is not None
        ):  # only report keys starting with clustering_name
            retVal = [x for x in sorted(filenames) if x.startswith(clustering_name)]
        else:
            retVal = list(sorted(filenames))
        return retVal

    def cluster_versions(self, clustering_key: str) -> List[bson.objectid.ObjectId]:
        """lists ids and storage dates corresponding to versions of clustering identifed by clustering_key.
        the newest version is first.
        """

        cursor = self.clusters.find({"filename": clustering_key}).sort("uploadDate", -1)
        retVal = []
        for res in cursor:

            retVal.append(res._id)
        return retVal

    def cluster_delete_all(self, clustering_key: str) -> None:
        """delete all clustering objects, including the latest version, stored under clustering_key"""
        ids = self.cluster_versions(clustering_key)
        for this_id in ids:
            self.clusters.delete(this_id)

    def cluster_delete_legacy_by_key(self, clustering_key: str) -> None:
        """delete all clustering objects, except latest version, stored with key clustering_key"""
        ids = self.cluster_versions(clustering_key)
        ids = ids[1:]
        for i, this_id in enumerate(ids):
            logging.info(
                "Removing historical data for {0} {1} / {2}".format(
                    clustering_key, i, len(ids)
                )
            )
            self.clusters.delete(this_id)

    def cluster_delete_legacy(self, clustering_name: str) -> None:
        """delete all clustering objects, except latest version, stored with  clustering_name"""
        clustering_keys = self.cluster_keys(clustering_name=clustering_name)
        for clustering_key in clustering_keys:
            self.cluster_delete_legacy_by_key(clustering_key)

    def refcompressedseq_store(self, guid: str, obj: Any) -> str:
        """stores the sequence object obj with guid guid.
        Issues an error FileExistsError
        if the guid already exists."""
        pickled_obj = pickle.dumps(obj, protocol=2)
        res = self.db.refcompressedseq.files.find_one({"_id": guid}, {"_id": 1})
        if res is not None:  # it exists
            raise FileExistsError("Attempting to overwrite {0}".format(guid))

        id = self.rcs.put(pickled_obj, _id=guid, filename=guid)

        # do a functional test to verify write
        recovered_obj = self.refcompressedsequence_read(guid)
        if not recovered_obj == obj:
            raise IOError(
                "Integrity check failed on reference compressed item write/read for {0}".format(
                    guid
                )
            )
        return id

    def refcompressedsequence_read(self, guid: str) -> Any:
        """loads object from refcompressedseq collection.
        It is assumed object stored is a dictionary"""

        res = self.rcs.find_one({"_id": guid})
        if res is None:
            return None

        return pickle.loads(res.read())

    def refcompressedsequence_guids(self) -> Set[str]:
        """loads guids from refcompressedseq collection."""

        # altered syntax because the .load() syntax previously used loaded > 16MB data and failed with > 600k samples
        res = self.db.refcompressedseq.files.find({}, {"_id": 1})
        guids = list()
        for item in res:
            guids.append(item["_id"])
        return set(guids)

    # methods for guid2meta
    def guid_annotate(self, guid: str, nameSpace: str, annotDict: dict) -> None:
        """adds multiple annotations of guid from a dictionary;
        all annotations go into a namespace.
        creates the record if it does not exist"""

        # check whethere there is an existing metadata object for this

        metadataObj = self.db.guid2meta.find_one({"_id": guid})
        if metadataObj is None:
            # it doesn't exist.  we create a new one.
            metadataObj = {"_id": guid, "sequence_meta": {nameSpace: annotDict}}

        if (
            "sequence_meta" not in metadataObj.keys()
        ):  # this is key is mandatory and is always present
            metadataObj["sequence_meta"] = {}

        # if the namespace does not exist as a subsidiary of sequence_meta, then we create it
        if nameSpace not in metadataObj["sequence_meta"].keys():
            metadataObj["sequence_meta"][nameSpace] = {}

        # we add any annotations to the existing data
        metadataObj["sequence_meta"][nameSpace] = {
            **metadataObj["sequence_meta"][nameSpace],
            **annotDict,
        }

        res = self.db.guid2meta.replace_one({"_id": guid}, metadataObj, upsert=True)
        if res.acknowledged is not True:
            raise IOError(
                "Mongo {0} did not acknowledge write of data: {1}".format(
                    self.db, metadataObj
                )
            )

    def guids(self) -> Set[str]:
        """returns all registered guids"""

        retVal = [x["_id"] for x in self.db.guid2meta.find({}, {"_id": 1})]
        return set(retVal)

    def guids_added_after_sample(self, guid: str) -> Set[str]:
        """ returns all guids added after a sample"""
        print("*** SEARCHING FOR ", guid)
        this_examination_time = self.guid_examination_time(guid)
        if this_examination_time is None:
            return None

        return self.guids_considered_after(addition_datetime = this_examination_time)

    def guids_considered_after(self, addition_datetime: datetime.datetime) -> Set[str]:
        """returns all registered guid added after addition_datetime
        addition_datetime: a date of datetime class."""

        if not isinstance(addition_datetime, datetime.datetime):
            raise TypeError(
                "addition_datetime must be a datetime.datetime value.  It is {0}.  Value = {1}".format(
                    type(addition_datetime), addition_datetime
                )
            )
            
        retVal = [
            x["_id"]
            for x in self.db.guid2meta.find(
                {
                    "sequence_meta.DNAQuality.examinationDate": {
                        "$gt": addition_datetime
                    }
                },
                {"_id": 1},
            )
        ]
        return set(retVal)

    def _guids_selected_by_validity(self, validity: int) -> Set[str]:
        """returns  registered guids, selected on their validity

        0 = guid is valid
        1 = guid is invalid

        """
        if validity not in [0, 1]:
            raise ValueError("Validity must be 0 or 1, not {0}".format(validity))

        retVal = [
            x["_id"]
            for x in self.db.guid2meta.find(
                {"sequence_meta.DNAQuality.invalid": validity}
            )
        ]

        return set(retVal)

    def singletons(
        self, method: str = "approximate", return_top: int = 1000
    ) -> pd.DataFrame:
        """returns guids and the number of singleton records, which
        (if high numbers are present) indicates repacking is needed.
        Inclusion of max_records is important for very large datasets, or the query is slow.

        Parameters:
        method: either 'approximate' or 'exact'.  For ranking samples for repacking, the approximate method is recommended.
        return_top:  the number of results to return.  Set > number of samples to return all records.  If method = 'exact', return_top is ignored and all records are returned.

        The approximate method is much faster than the exact one if large numbers of singletons are present; it randomly downsamples the guid-neighbour pairs in
        the guid2neighbour collection, taking 100k samples only, and the counts the number of singletons in the downsample.  This is therefore a method of ranking
        samples which may benefit from repacking.  This sampling method returns queries in a few milliseconds in testing on large tables.
        This method is not deterministic.

        In the event of failure to successfully generate a sample of sufficient size (set to 100k at present),
        which can happen if there are fewer than singletons, then the exact method will be used as a fallback.  If fallback of this kind occurs, only return_top samples will be returned.

        The exact method computes the exact non-zero number of singletons for each sample.  This requires an index scan which can be
        surprisingly slow, with the query taking > 30 seconds with > table size ~ 3 x10^7.   This method is deterministic.

        Returns:
        a set of sample identifiers ('guids') which contain > min_number_records singleton entries.
        """

        # input validation
        if not isinstance(return_top, int):
            raise TypeError(
                "return_top is {0}; this must be an integer not a {1}".format(
                    return_top, type(return_top)
                )
            )
        if not return_top > 0:
            raise TypeError(
                "return_top is {0}; this must be a non zero positive integer".format(
                    return_top
                )
            )
        if method not in ["exact", "approximate"]:
            raise ValueError(
                "Method must be either 'exact' or 'approximate', not {0}".format(method)
            )

        # use mongodb pipeline
        # shell example: db.guid2neighbour.aggregate([ {$sample: {size: 100000}}, { $match: {rstat:'s'}}, { $sortByCount: "$guid"}, {$limit: 1000}  ] )
        approximate_pipeline = [
            {"$sample": {"size": 100000}},
            {"$match": {"rstat": "s"}},
            {"$sortByCount": "$guid"},
            {"$limit": return_top},
        ]
        exact_pipeline = [{"$match": {"rstat": "s"}}, {"$sortByCount": "$guid"}]
        fallback = False
        if method == "approximate":

            try:
                results = self.db.guid2neighbour.aggregate(
                    approximate_pipeline, allowDiskUse=True
                )
            except pymongo.errors.OperationFailure:
                # occurs if there are very few samples left; a random sample of the required size (in this case 100k) cannot be generated
                method = "exact"
                logging.info(
                    "mongoStore.singletons | unable to generate a random sample of the required size, due to a small numbers of singletons.  Falling back to an exact method."
                )
                fallback = True

        if method == "exact":
            results = self.db.guid2neighbour.aggregate(
                exact_pipeline, allowDiskUse=True
            )

        ret_list = []
        for item in results:
            ret_list.append(item)

        if fallback is True and len(ret_list) > return_top:
            logging.info(
                "Fallback prcesses in place; restricting to top {0}.  Exact search examined {1} samples with singletons".format(
                    return_top, len(ret_list)
                )
            )
            ret_list = ret_list[:return_top]
        ret_df = pd.DataFrame.from_records(ret_list)
        ret_df.rename(columns={"_id": "guid"}, inplace=True)

        if (
            ret_df.empty
        ):  # if empty, the '_id' column is not there, and the rename fails
            ret_df = pd.DataFrame(columns=("guid", "count"))
        ret_df.set_index("guid", drop=True, inplace=True)

        return ret_df

    def guids_valid(self) -> set:
        """return all registered valid guids.

        Validity is determined by the contents of the DNAQuality.invalid field, on which there is an index"""
        return self._guids_selected_by_validity(0)

    def guids_invalid(self) -> set:
        """return all invalid guids

        Validity is determined by the contents of the DNAQuality.invalid field, on which there is an index"""
        return self._guids_selected_by_validity(1)

    def guid_exists(self, guid: str) -> bool:
        """checks the presence of a single guid"""

        res = self.db.guid2meta.find_one({"_id": guid})
        if res is None:
            return False
        else:
            return True

    def guid_valid(self, guid: str) -> int:
        """checks the validity of a single guid

        Parameters:
        guid: the sequence identifier

        Returns
        -1   The guid does not exist
        0    The guid exists and the sequence is valid
        1    The guid exists and the sequence is invalid
        -2    The guid exists, but there is no DNAQuality.valid key"""

        res = self.db.guid2meta.find_one({"_id": guid})
        if res is None:
            return -1
        else:

            try:
                return int(res["sequence_meta"]["DNAQuality"]["invalid"])
            except KeyError:
                return -2

    def guid_examination_time(self, guid: str) -> Optional[datetime.datetime]:
        """returns the examination time for a single guid

        Parameters:
        guid: the sequence identifier

        Returns either
        The examination datetime value for this guid OR
        None if the guid does not exist, or the sequence_meta.DNAQuality.examinationTime key does not exist.
        """

        res = self.db.guid2meta.find_one(
            {"_id": guid}, {"sequence_meta.DNAQuality.examinationDate": 1}
        )
        if res is None:
            return None
        try:
            return res["sequence_meta"]["DNAQuality"]["examinationDate"]
        except KeyError:
            return None

    def guids_considered_after_guid(self, guid: str) -> Set[str]:
        """returns all registered guids added after guid
        guid: a sequence identifier"""

        addition_datetime = self.guid_examination_time(guid)
        if addition_datetime is None:
            raise ValueError("guid is not valid: {0}".format(guid))
        else:
            return self.guids_considered_after(addition_datetime)

    def guid_quality_check(
        self, guid: str, cutoff: Union[float, int]
    ) -> Optional[bool]:
        """Checks whether the quality of one guid exceeds the cutoff.

        If the guid does not exist, returns None.
        If the guid does exist and has quality< cutoff, returns False.
        Otherwise, returns True.
        """

        # test input
        if not type(cutoff) in [float, int]:
            raise TypeError(
                "Cutoff should be either floating point or integer, but it is %s"
                % type(cutoff)
            )
        if not type(guid) == str:
            raise TypeError("The guid passed should be as string, not %s" % str(guid))

        # recover record, compare with quality

        res = self.db.guid2meta.find_one({"_id": guid}, {"sequence_meta": 1})
        if res is None:  # no entry for this guid
            return None
        else:
            try:
                dnaq = res["sequence_meta"]["DNAQuality"]
            except KeyError:
                raise KeyError(
                    "DNA quality is not present in the sequence metadata {0}: {1}".format(
                        guid, res
                    )
                )

            # check the DNA quality metric expected is present
            if "propACTG" not in dnaq.keys():
                raise KeyError(
                    "propACTG is not present in DNAQuality namespace of guid {0}: {1}".format(
                        guid, dnaq
                    )
                )

            # report whether it is larger or smaller than cutoff
            return dnaq["propACTG"] >= cutoff

    def guid2item(
        self, guidList: Optional[List[str]], namespace: str, tag: str
    ) -> Optional[dict]:
        """returns the annotation (such as sequence quality, which is stored as an annotation)
        in namespace:tag for all guids in guidlist.
        If guidList is None, all items are returned.
        An error is raised if namespace and tag is not present in each record.
        """

        retDict = {}
        if guidList is None:
            results = self.db.guid2meta.find({}, {"sequence_meta": 1})
        else:
            results = self.db.guid2meta.find(
                {"_id": {"$in": guidList}}, {"sequence_meta": 1}
            )

        if results is None:  # nothing found
            return None

        for res in results:
            try:
                namespace_content = res["sequence_meta"][namespace]
            except KeyError:
                raise KeyError(
                    "{2} is not present in the sequence metadata {0}: {1}".format(
                        guidList, res, namespace
                    )
                )

            # check the DNA quality metric expected is present
            if tag not in namespace_content.keys():
                raise KeyError(
                    "{2} is not present in {3} namespace of guid {0}: {1}".format(
                        guidList, namespace_content, tag, namespace
                    )
                )

            # return property
            retDict[res["_id"]] = namespace_content[tag]
        return retDict

    def guid2ExaminationDateTime(
        self, guidList: Optional[List[str]] = None
    ) -> Optional[dict]:
        """returns quality scores and examinationDate for all guids in guidlist.  If guidList is None, all results are returned."""

        return self.guid2item(guidList, "DNAQuality", "examinationDate")

    def guid2quality(self, guidList: Optional[List[str]] = None) -> Optional[dict]:
        """returns quality scores for all guids in guidlist (or all samples if guidList is None)
        potentially expensive query if guidList is None."""

        return self.guid2item(guidList, "DNAQuality", "propACTG")

    def guid2propACTG_filtered(self, cutoff: Union[int, float] = 0) -> Dict[str, float]:
        """recover guids which have good quality, > cutoff.
        These are in the majority, so we run a table scan to find these.

        This query is potentially very inefficient- best avoided
        """

        allresults = self.db.guid2meta.find(
            {"sequence_meta.DNAQuality.propACTG": {"$gte": cutoff}},
            {"_id": 1, "sequence_meta.DNAQuality.propACTG": 1},
        )

        retDict = {}
        for item in allresults:
            retDict[item["_id"]] = item["sequence_meta"]["DNAQuality"]["propACTG"]

        return retDict  # note: slightly different from previous api

    def guid2items(
        self, guidList: Optional[List[str]], namespaces: Optional[Set[str]]
    ) -> Optional[Dict[Any, Dict[str, Any]]]:
        """returns all annotations in namespaces, which is a list, as a pandas dataframe.
        If namespaces is None, all namespaces are returned.
        If guidList is None, all items are returned.
        To do this, a table scan is performed - indices are not used.
        """

        retDict = {}
        if guidList is None:
            results = self.db.guid2meta.find({}, {"sequence_meta": 1})
        else:
            results = self.db.guid2meta.find(
                {"_id": {"$in": guidList}}, {"sequence_meta": 1}
            )

        if results is None:  # nothing found
            return None

        for res in results:
            row = {}
            sought_namespaces = set(res["sequence_meta"].keys())
            if namespaces is not None:  # we only want a subset
                sought_namespaces = sought_namespaces.intersection(
                    namespaces
                )  # what we want, intersect what we've got

            for sought_namespace in sought_namespaces:
                for tag in res["sequence_meta"][sought_namespace].keys():
                    col_name = "{0}:{1}".format(sought_namespace, tag)
                    row[col_name] = res["sequence_meta"][sought_namespace][tag]
            retDict[res["_id"]] = row

        return retDict

    def guid_annotations(self) -> Optional[Dict[Any, Dict[str, Any]]]:
        """return all annotations of all guids"""

        return self.guid2items(None, None)  # no restriction by namespace or by guid.

    def guid_annotation(self, guid: str) -> Optional[Dict[Any, Dict[str, Any]]]:
        """return all annotations of one guid"""

        return self.guid2items([guid], None)  # restriction by guid.

    def guid2neighbour_add_links(
        self,
        guid: str,
        targetguids: Dict[str, Dict[str, int]],
        use_update: bool = False,
    ) -> Dict[str, int]:
        """adds links between guid and their neighbours ('targetguids')

        Parameters:
        guid: the 'source' guid for the matches eg 'guid1'
        targetguids: what is guid linked to, eg
        {
                'guid2':{'dist':12},
                'guid3':{'dist':2}
        }
        use_update -  currently ignored, always False.  Setting True yields NotImplementedError

        This stores links in the guid2neighbour collection.
        If use_update = True, will update target documents, adding a new link to the targetguid document.
            {targetguid -> {previousguid: distance1, previousguid2: distance2}} --> {targetguid -> {previousguid: distance1, previousguid2: distance2, guid: distance}
            This approach has many disadvantages
            - multiple database accesses: one to find a document to update, and one to update it - may be hundreds or thousands of database connections for each insert operation
            - approach is (with Mongodb) not inherently atomic.  If failure occurs in the middle of the process, database can be left in an inconsistent state, requiring use of transactions (slower)
            - by contrast, the no_update method (default) requires a single insert_many operation, which is much cleaner and is atomic.
            - the use_update approach is not implemented
        If use_update = False (default), for each new link from targetguid -> guid, a new document will be inserted linking {targetguid -> {guid: distance}}

        The function guid2neighbour_repack() reduces the number of documents
        required to store the same information.

        if annotation is not None, will additionally write an annotation dictionary

        Returns:
        The number of records written
        """

        # find guid2neighbour entry for guid.
        to_insert: List[Dict[str, Any]] = []

        current_m: Optional[
            Dict[str, Any]
        ] = None  # no target record for guid -> multiple targets.  We made a new one.
        for targetguid in targetguids.keys():
            payload = targetguids[targetguid]  # a distance

            ## ADD REVERSE LINK
            # NOTE: for the (new) guid --> targetguids, we can write a single record with many entries
            # however, the other way round, we have to add multiple single samples
            #  targetguid --> guid (reverse link) (singleton)

            if (
                use_update
            ):  # slower: multiple update operations against the database on insert, not atomic see above
                raise NotImplementedError(
                    "Updating method for adding links is not implemented"
                )
            else:
                # insert a new record, which can be processed by update_many
                to_insert.append(
                    {"guid": targetguid, "rstat": "s", "neighbours": {guid: payload}}
                )

            # Now add one record containing many target guids
            # guid --> {many targetguids } (forward link)
            if current_m is not None:
                if (
                    len(current_m["neighbours"].keys())
                    >= self.max_neighbours_per_document
                ):  # need to make another one
                    current_m["rstat"] = "f"  # full
                    to_insert.append(current_m)
                    current_m = None

            # if we don't have a current_m, which is a data structure containing what we're going to write to disc,
            # either because this is the first pass through the loop, or becuase our previous current_m was full (see above)
            # then we make one
            if current_m is None:
                current_m = {"guid": guid, "rstat": "m", "neighbours": {}}  # make one

            # add the element to current_m
            current_m["neighbours"][targetguid] = payload

        # now we've covered all the elements to write
        # do we have a current m with data in it?
        if current_m is not None:
            if len(current_m["neighbours"].keys()) > 0:
                # check if it's full
                if (
                    len(current_m["neighbours"].keys())
                    >= self.max_neighbours_per_document
                ):
                    current_m["rstat"] = "f"  # full
                # check if it's actually a singleton
                if len(current_m["neighbours"].keys()) == 1:
                    current_m["rstat"] = "s"  # just one
                # add to to_insert
                to_insert.append(current_m)

        # when complete, do update
        if len(to_insert) > 0:
            res = self.db.guid2neighbour.insert_many(to_insert, ordered=False)
            if res.acknowledged is not True:
                raise IOError(
                    "Mongo {0} did not acknowledge write of data: {1}".format(
                        self.db, to_insert
                    )
                )
        # check there is a metadata object for the guid
        metadataObj = self.db.guid2meta.find_one({"_id": guid})
        if metadataObj is None:
            # it doesn't exist.  we create a new one.
            metadataObj = {
                "_id": guid,
                "created": {"created_at": datetime.datetime.now().isoformat()},
            }
            res = self.db.guid2meta.insert_one({"_id": guid}, metadataObj)
            if res.acknowledged is not True:
                raise IOError(
                    "Mongo {0} did not acknowledge write of data: {1}".format(
                        self.db, metadataObj
                    )
                )

        return {"records_written": len(to_insert)}

    def _audit_storage(self, guid: str) -> Tuple[set, pd.DataFrame, dict]:
        """returns a pandas data frame containing all neighbours of guid, as stored in mongo, as well as
        a summary of storage statistics and a list of _ids which could be optimised

        Parameters:
        guid : the identifier of a sample to repack

        Returns:
        A tuple containing three elements:
        to_optimise: a set of mongodb record _ids which are either single, not full (m) or contain duplicate guids and so could be optimised
        content_df: a dataframe containing all record _ids, the neighbouring guid, and the snp distance.
        audit_stats: a dictionary containing the following:
            singletons: number of singleton records
            m_records:  number of records which are not completely full
            f_records:  number of full records
            f_records_with_duplicates: number of full records containing the same guid twice
            neighbouring_guids: the number of neighbouring guids
            neighbouring_guids_recordings: number of times a neighbouring guid is recorded.  Maybe larger than neighbouring_guids.

        Designed for internal use
        """

        ## audit the storage of neighbours

        # read all occurrences of each neighbour and its distances into a pandas dataframe
        # note that its is OK for a neighbour to be present more than once
        bytes_per_record: Dict[str, List[int]] = {"s": [], "m": [], "f": []}
        contents = []
        for res in self.db.guid2neighbour.find({"guid": guid}):
            bytes_per_record[res["rstat"]].append(len(str(res)))
            location = res.copy()
            del location["neighbours"]

            for neighbouring_guid in res["neighbours"].keys():
                storage_element = location.copy()
                storage_element["guid"] = neighbouring_guid
                storage_element["dist"] = res["neighbours"][neighbouring_guid]["dist"]
                contents.append(storage_element)
        content_df = pd.DataFrame.from_records(contents)

        # are there any neighbours?
        if len(contents) == 0:
            audit_stats: Dict[str, Union[int, float]] = {
                "singletons": 0,
                "stored_in_m_records": 0,
                "stored_in_f_records": 0,
                "f_records_with_duplicates": 0,
                "neighbouring_guids": 0,
            }
            return (set([]), content_df, audit_stats)

        # there are some neighbours
        # diagnostics - how many neighbours per document
        neighbours_per_document = (
            content_df.groupby(["_id", "rstat"]).size().reset_index(name="counts")
        )

        max_pr = (
            neighbours_per_document.groupby(["rstat"])["counts"]
            .max()
            .reset_index(name="max_neighbours_per_record")
        )
        mean_pr = (
            neighbours_per_document.groupby(["rstat"])["counts"]
            .mean()
            .reset_index(name="mean_neighbours_per_record")
        )
        report_loading = {}
        for ix in max_pr.index:
            report_loading[
                max_pr.at[ix, "rstat"] + "_max_neighbours_per_record"
            ] = max_pr.at[ix, "max_neighbours_per_record"]

        for ix in mean_pr.index:
            report_loading[
                mean_pr.at[ix, "rstat"] + "_mean_neighbours_per_record"
            ] = mean_pr.at[ix, "mean_neighbours_per_record"]

        # compute things to optimise.  These are
        # records containing duplicate guids
        duplicates = content_df.groupby(["guid"]).size().reset_index(name="counts")
        n_guids = len(duplicates.index)
        duplicates = duplicates[duplicates["counts"] > 1]
        duplicate_guids = set(duplicates["guid"])

        # and any singletons
        # and any ms
        # plus any fs which have duplicates in them (these should be very rare, if ever occurring)
        to_optimise = set()
        n_s = 0
        n_m = 0
        n_duplicates = 0
        n_duplicates_f = 0
        n_f = 0
        for ix in content_df.index:
            if (
                content_df.at[ix, "rstat"] in ["s", "m"]
                or content_df.at[ix, "guid"] in duplicate_guids
            ):
                to_optimise.add(content_df.at[ix, "_id"])

            # this section purely to gather audit data
            if content_df.at[ix, "rstat"] == "s":
                n_s += 1
            if content_df.at[ix, "rstat"] == "m":
                n_m += 1
            if content_df.at[ix, "guid"] in duplicate_guids:
                n_duplicates += 1
            if (
                content_df.at[ix, "guid"] in duplicate_guids
                and content_df.at[ix, "rstat"] == "f"
            ):
                n_duplicates_f += 1
            if (
                not content_df.at[ix, "guid"] in duplicate_guids
                and content_df.at[ix, "rstat"] == "f"
            ):
                n_f += 1
        audit_stats = {
            "singletons": n_s,
            "stored_in_m_records": n_m,
            "stored_in_f_records": n_f,
            "f_records_with_duplicates": n_duplicates_f,
            "neighbouring_guids": n_guids,
            "neighbouring_guids_recordings": len(content_df),
        }
        audit_stats.update(report_loading)
        for rstat in ["s", "f", "m"]:
            if len(bytes_per_record[rstat]) > 0:
                audit_stats[rstat + "_mean_bytes_per_record"] = statistics.mean(
                    bytes_per_record[rstat]
                )
                audit_stats[rstat + "_max_bytes_per_record"] = max(
                    bytes_per_record[rstat]
                )
                audit_stats[rstat + "_records"] = len(bytes_per_record[rstat])

        return (to_optimise, content_df, audit_stats)

    class Guid2NeighbourRepackRet(TypedDict):
        guid: str
        finished: str
        pre_optimisation: dict
        actions_taken: Dict[str, int]

    def guid2neighbour_repack(
        self, guid: str, always_optimise: bool = False, min_optimisable_size: int = 1
    ) -> Guid2NeighbourRepackRet:
        """alters the mongodb representation of the links of guid.

        This stores links in the guid2neighbour collection;
        each stored document links one guid to one target.

        This function reduces the number of documents
        required to store the same information.

        Internally, the documents in guid2neighbour are of the form
        {'guid':'guid1', 'rstat':'s', 'neighbours':{'guid2':{'dist':12, ...}}} OR
        {'guid':'guid1', 'rstat':'m', 'neighbours':{'guid2':{'dist':12, ...}, 'guid3':{'dist':5, ...}} OR
        {'guid':'guid1', 'rstat':'f', 'neighbours':{'guid2':{'dist':12, ...}, 'guid3':{'dist':5, ...}}

        The last example occurs when the maximum number of neighbours permitted per record has been reached.

        This operation moves rstat='s' record's neighbours into rstat='m' type records.
        It guarantees that no guid-guid links will be lost, but more than one record of the same link may exist
        during processing.  This does not matter, because guid2neighbours() deduplicates.

        Parameters:
        guid : the identifier of a sample to repack
        always_optimise: consider for optimisation even if there are no 's' (single item) records
        """

        if not always_optimise:
            # determine whether how many one rstat 's' entries for this guid.
            n_s_records = self.db.guid2neighbour.count_documents(
                {"guid": guid, "rstat": "s"}
            )
            if (
                n_s_records < min_optimisable_size
            ):  # either no singletons, or just one: nil to optimise.
                return {
                    "guid": guid,
                    "finished": datetime.datetime.now().isoformat(),
                    "pre_optimisation": {
                        "s_records": n_s_records,
                        "msg": "Below optimisation cutoff",
                    },
                    "actions_taken": {"n_written": 0, "n_deleted": 0},
                }

        # get a summary of how the records are currently stored
        # logging.info("Auditing")
        to_optimise, content_df, audit_stats = self._audit_storage(guid)

        # now we adopt a very simple strategy.
        # We assign all the to_optimise ids ('s' single records, and 'm' partially full records) for destruction,
        # and repartition all of their contents across new rows in the collection.
        # This prevents any updating: only insert and delete operations are used
        current_m: Optional[Dict[str, Any]] = None  # the record to store the records in
        to_write = []
        for ix in content_df.index:
            if (
                content_df.at[ix, "_id"] in to_optimise
            ):  # if this data element needs to be rewritten

                # if current_m, which contains the data we are preparing to write, is full, then we mark it as full, and
                # transfer it to to_write, an array containing things we need to write to disc.
                if current_m is not None:
                    if (
                        len(current_m["neighbours"].keys())
                        >= self.max_neighbours_per_document
                    ):
                        current_m["rstat"] = "f"  # full
                        to_write.append(current_m)
                        current_m = None

                # if we don't have a current_m, which is a data structure containing what we're going to write to disc,
                # either because this is the first pass through the loop, or becuase our previous current_m was full (see above)
                # then we make one
                if current_m is None:
                    current_m = {
                        "guid": guid,
                        "rstat": "m",
                        "neighbours": {},
                    }  # make one

                # add the element to current_m
                current_m["neighbours"][content_df.at[ix, "guid"]] = {
                    "dist": int(content_df.at[ix, "dist"])
                }  # int because stored as numpy.int64, which is not jsonable

        # now we've covered all the elements to rewrite
        # do we have a current m with data in it?
        if current_m is not None:
            if len(current_m["neighbours"].keys()) > 0:
                # check if it's full
                if (
                    len(current_m["neighbours"].keys())
                    >= self.max_neighbours_per_document
                ):
                    current_m["rstat"] = "f"  # full
                # add to to_write
                to_write.append(current_m)

        # write the documents which need writing
        # logging.info("Inserting many")
        self.db.guid2neighbour.insert_many(to_write, ordered=False)

        # delete those processed single records
        # logging.info("Deleting many")
        self.db.guid2neighbour.delete_many({"_id": {"$in": list(to_optimise)}})

        # return a report of what we did
        actions = dict(n_written=len(to_write), n_deleted=len(to_optimise))
        # logging.info("complete")
        return {
            "guid": guid,
            "finished": datetime.datetime.now().isoformat(),
            "pre_optimisation": audit_stats,
            "actions_taken": actions,
        }

    class Guid2NeighboursRet(TypedDict):
        guid: str
        neighbours: List[Guid2NeighboursFormats]

    def guid2neighbours(
        self, guid: str, cutoff: int = 20, returned_format: int = 2
    ) -> Guid2NeighboursRet:
        """returns neighbours of guid with cutoff <=cutoff.
        Returns links either as

        format 1 [[otherGuid, distance],[otherGuid2, distance2],...]
        or as
        format 2 [[otherGuid, distance, N_just1, N_just2, N_either],[],...]
        or as
        format 3 [otherGuid1, otherGuid2, otherGuid3]
        or as
        format 4 [{'guid':otherguid, 'snv':distance}]

            Internally, the documents in guid2neighbour are of the form
            {'guid':'guid1', 'rstat':'s', 'neighbours':{'guid2':{'dist':12, ...}}} OR
            {'guid':'guid1', 'rstat':'m', 'neighbours':{'guid2':{'dist':12, ...}, 'guid3':{'dist':5, ...}} OR
            {'guid':'guid1', 'rstat':'f', 'neighbours':{'guid2':{'dist':12, ...}, 'guid3':{'dist':5, ...}}

            However, irrespective of their internal representation, this function always returns
            exactly one item for each link of 'guid'; duplicates are not possible.
            The last example occurs when the maximum number of neighbours permitted per record has been reached.
        """

        retVal: List[Guid2NeighboursFormats] = []
        formatting: Dict[int, List[str]] = {
            1: ["dist"],
            2: ["dist", "N_just1", "N_just2", "N_either"],
            3: [],
            4: ["dist"],
        }
        desired_fields = formatting[returned_format]
        results = self.db.guid2neighbour.find({"guid": guid})
        reported_already = set()
        for result in results:
            for otherGuid in result["neighbours"].keys():
                if otherGuid not in reported_already:  # exclude duplicates
                    if (
                        result["neighbours"][otherGuid]["dist"] <= cutoff
                    ):  # if distance < cutoff
                        # available_fields = result['neighbours'][otherGuid].keys()
                        reported_fields = {}
                        for desired_field in desired_fields:
                            try:
                                observed = result["neighbours"][otherGuid][
                                    desired_field
                                ]
                            except KeyError:  # doesn't exist
                                observed = None
                            reported_fields[desired_field] = observed

                        returned_data: Guid2NeighboursFormats

                        if returned_format == 1:
                            returned_data = [otherGuid, reported_fields["dist"]]

                        elif returned_format == 2:
                            raise NotImplementedError("format 2 is no longer supported")
                        elif returned_format == 3:
                            returned_data = otherGuid

                        elif returned_format == 4:
                            returned_data = {
                                "guid": otherGuid,
                                "snv": reported_fields["dist"],
                            }

                        else:
                            raise ValueError(
                                "Unable to understand returned_format = {0}".format(
                                    returned_format
                                )
                            )

                        reported_already.add(otherGuid)
                        retVal.append(returned_data)

        # recover the guids
        return {"guid": guid, "neighbours": retVal}

    def _set_lock_status(self, lock_int_id, lock_status):
        """locks or unlocks resources identified by lock_int_id, allowing cross- process sequential processing (e.g. insertion)

        To lock, set lock_status =1 ; to unlock, set lock_status =0
        To return the relevant row, set lock_status to None

        See the acquire_lock() method for more details

        returns:
        True if update succeeded, false if it did not

        Technical notes:
        https://www.mongodb.com/blog/post/how-to-select--for-update-inside-mongodb-transactions

        """
        # make sure there is an entry for this lock
        lock_row = self.db.fnlock.find_one({"_id": lock_int_id})

        # if the row doesn't exist, we add it, with the lock not set.
        if lock_row is None:
            lock_row = dict(
                _id=lock_int_id,
                lock_status=0,
                lock_set_date=datetime.datetime.now(),
                uuid=uuid.uuid4().hex,
            )
            self.db.fnlock.insert_one(lock_row)

        # analyse the record for this row
        with self.client.start_session(causal_consistency=True):
            self.db.fnlock.update_one(
                {"_id": lock_int_id}, {"$set": {"uuid": uuid.uuid4().hex}}
            )  # this places a lock on the result, becuase it has been modified.  Attempts to write by other transactions will yield errors.
            lock_row = self.db.fnlock.find_one({"_id": lock_int_id})
            if lock_status is None:
                retval = lock_row

            if lock_row["lock_status"] == 0 and lock_status == 0:
                # it's already unlocked
                retval = True

            elif lock_row["lock_status"] == 1 and lock_status == 1:
                # it's already locked and we're asked to acquire a lock.  We can't.
                retval = False

            elif lock_row["lock_status"] == 0 and lock_status == 1:
                # it's already unlocked, we can lock
                lock_row["lock_status"] = 1
                lock_row["lock_set_date"] = datetime.datetime.now()
                lock_row["uuid"] = uuid.uuid4().hex
                self.db.fnlock.replace_one({"_id": lock_int_id}, lock_row)
                retval = True

            elif lock_row["lock_status"] == 1 and lock_status == 0:
                lock_row["lock_status"] = 0
                lock_row["lock_set_date"] = datetime.datetime.now()
                lock_row["uuid"] = uuid.uuid4().hex
                self.db.fnlock.replace_one({"_id": lock_int_id}, lock_row)
                retval = True

            return retval

    def lock_status(self, lock_int_id):
        """determine whether a database-based lock is open (0) or closed (1).

        Parameters:
        lock_int_id: an integer identifier to the lock of interest

        Returns:
        0 if the lock is open
        1 if it is locked"""

        return self._set_lock_status(lock_int_id, None)

    def lock(self, lock_int_id):
        """obtains a database-based lock.

        Parameters:
        lock_int_id: an integer identifier to the lock of interest

        Returns:
        True if the lock is acquired
        False if it is not"""

        return self._set_lock_status(lock_int_id, 1)

    def unlock(self, lock_int_id, force=False):
        """obtains a database-based lock.

        Parameters:
        lock_int_id: an integer identifier to the lock of interest
        force: if True, will unlock irrespective of current status, returning True

        Returns:
        True if the lock is acquired
        False if it is not"""

        res = self._set_lock_status(lock_int_id, 0)
        if force:
            return True
        else:
            return res
