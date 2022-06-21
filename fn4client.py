#!/usr/bin/env python3
""" A python3 client for findNeighbour4-server.

Provides a class which allows access to all routes described in:
doc/rest-routes.md for a list.
  
Unit testing:
* launch a test server
python findNeighbour4-server.py

* run fn4Client unit tests
python -m unittest test/test_fn4client

Note: the script run_tests.sh will do this automatically

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@ukhsa.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

 

"""

import requests
import json
import urllib.parse
import logging
import gzip
import io
import pandas as pd

# used for loading fasta files.
from Bio import SeqIO


class fn4Client:
    """python3 API to the findNeighbour3 + -server REST endpoint.

    All endpoints are supported.
    See doc/rest-routes.md for a list.

    """

    def __init__(
        self,
        baseurl="http://127.0.0.1:5020",
    ):

        """ set up logging """
        logging.getLogger()

        # load config
        self.baseurl = baseurl

        # run connection check
        res = self.server_time()
        logging.info("Connection established at server time {0}".format(res["server_time"]))

    def _decode(self, response):
        """checks that a response object has response code in the 200s.
        If it doesn't, raises an error.
        If it does, decodes the object and returns json."""
        response.raise_for_status()
        return json.loads(response.content.decode("utf-8"))

    def _isjson(self, content):
        """ returns true if content parses as json, otherwise false """
        try:
            json.loads(content.decode("utf-8"))
            return True

        except json.decoder.JSONDecodeError:
            return False

    def absurl(self, relpath):
        """ constructs an absolute URL from the relative path requested """
        # removed: self.version+
        absurl = urllib.parse.urljoin(self.baseurl, relpath)
        return absurl

    def getpost(self, relpath, method="GET", payload=None, timeout=None):
        """issues GET or POST against url.  returns a response object
        will raise errors if generated, but does not raise errors if a valid
        response object is returned, even if it has a status_code indicating the call failed
        """

        with requests.Session() as session:
            session.trust_env = False

            url = self.absurl(relpath)
            if method == "GET":
                response = session.get(url=url, timeout=timeout)

            elif method == "POST":
                response = session.post(url=url, data=payload, timeout=timeout)

            else:
                raise NotImplementedError("either GET or POST is required as a method.  was passed {0}".format(method))

        if response.status_code >= 500:
            response.raise_for_status()  # raise error if there's an error.  Sub 500 (404 etc) are let through and handled by the client routines
        return response

    def get(self, relpath, timeout=None):
        """issues GET against url.  returns a response object.
        raises errors if fails"""
        response = self.getpost(relpath=relpath, timeout=timeout, method="GET")
        return response

    def post(self, relpath, payload, timeout=None):
        """ issues  POST against url.  returns a response object. """
        response = self.getpost(relpath=relpath, payload=payload, timeout=timeout, method="POST")
        return response

    def mirror(self, payload, timeout=0.2):
        """returns payload via server round trip.
        payload much be a dictionary, with key- value pairs where values are strings.
        Useful for testing server"""
        r = self.getpost(relpath="/api/v2/mirror", timeout=timeout, payload=payload, method="POST")
        rd = self._decode(r)
        return rd

    def server_config(self, timeout=None):
        """ returns server config as a dictionary """
        return self._decode(self.getpost("/api/v2/server_config", timeout=timeout, method="GET"))

    def server_time(self, timeout=None):
        """ returns server time as an isoformat string"""
        return self._decode(self.getpost("/api/v2/server_time", timeout=timeout, method="GET"))

    def nucleotides_excluded(self, timeout=0.2):
        """ returns the nucleotides excluded (i.e. the mask) """
        return self._decode(self.getpost("/api/v2/nucleotides_excluded", timeout=timeout, method="GET"))

    def guids(self, timeout=None):
        """ returns all guids in the server """
        return self._decode(self.getpost("/api/v2/guids", method="GET", timeout=timeout))

    def annotations(self, timeout=None):
        """ returns all guids and their annotations as a pandas dataframe"""
        retVal = self._decode(self.getpost("/api/v2/annotations", timeout=timeout, method="GET"))
        return pd.DataFrame.from_dict(retVal, orient="index")

    def clustering(self, timeout=None):
        """ return clustering pipelines available """
        return self._decode(self.getpost("/api/v2/clustering", timeout=timeout, method="GET"))

    def guids_and_examination_times(self, timeout=1):
        return self._decode(self.getpost("/api/v2/guids_and_examination_times", timeout=timeout, method="GET"))

    def guid_exists(self, guid, timeout=None):
        """ returns True or False depending on whether the guid exists """
        if not isinstance(guid, str):
            raise TypeError("guid {0} passed must be a string, not a {1}".format(guid, type(guid)))
        return self._decode(self.getpost("/api/v2/{0}/exists".format(guid), timeout=timeout, method="GET"))

    def sequence(self, guid, timeout=None):
        """ returns masked sequence of an existing guid """
        if not isinstance(guid, str):
            raise TypeError("guid {0} passed must be a string, not a {1}".format(guid, type(guid)))
        res = self.getpost("/api/v2/{0}/sequence".format(guid), timeout=timeout, method="GET")
        if res.status_code == 404:
            return None
        else:
            return self._decode(res)

    def change_id(self, clustering_algorithm, timeout=None):
        """ returns the current change_id associated with clustering_algorithm """
        if not isinstance(clustering_algorithm, str):
            raise TypeError("clustering_algorithm must be str not {0}".format(type(clustering_algorithm)))
        return self._decode(
            self.getpost("/api/v2/clustering/{0}/change_id".format(clustering_algorithm), timeout=timeout, method="GET")
        )

    def guids2clusters(self, clustering_algorithm, after_change_id=None, timeout=None):
        """ returns a guid2cluster lookup """
        if not isinstance(clustering_algorithm, str):
            raise TypeError("clustering_algorithm must be str not {0}".format(type(clustering_algorithm)))
        if after_change_id is None:
            res = self._decode(
                self.getpost(
                    "/api/v2/clustering/{0}/guids2clusters".format(clustering_algorithm), timeout=timeout, method="GET"
                )
            )
            return pd.DataFrame.from_records(res)

        elif isinstance(after_change_id, int):
            res = self._decode(
                self.getpost(
                    "/api/v2/clustering/{0}/guids2clusters/after_change_id/{1}".format(clustering_algorithm, after_change_id),
                    timeout=timeout,
                    method="GET",
                )
            )
            return pd.DataFrame.from_records(res)
        else:
            raise TypeError("after must be None or an integer, not {0}".format(type(after_change_id)))

    def clusters(self, clustering_algorithm, timeout=None):
        """ returns a clusters for a given clustering_algorithm """
        if not isinstance(clustering_algorithm, str):
            raise TypeError("clustering_algorithm must be str not {0}".format(type(clustering_algorithm)))

        res = self._decode(
            self.getpost("/api/v2/clustering/{0}/clusters".format(clustering_algorithm), timeout=timeout, method="GET")
        )
        return pd.DataFrame.from_records(res["members"])

    def cluster_members(self, clustering_algorithm, timeout=None):
        """ synonym for clusters """
        return self.clusters(clustering_algorithm, timeout=timeout)

    def cluster_summary(self, clustering_algorithm, timeout=None):
        """ returns a clusters and counts in each cluster for a given clustering_algorithm """
        if not isinstance(clustering_algorithm, str):
            raise TypeError("clustering_algorithm must be str not {0}".format(type(clustering_algorithm)))

        res = self._decode(
            self.getpost("/api/v2/clustering/{0}/summary".format(clustering_algorithm), timeout=timeout, method="GET")
        )
        return pd.DataFrame.from_records(res["summary"])

    def cluster_ids(self, clustering_algorithm, timeout=None):
        """ returns a cluster_ids for a given clustering_algorithm """
        if not isinstance(clustering_algorithm, str):
            raise TypeError("clustering_algorithm must be str not {0}".format(type(clustering_algorithm)))

        res = self._decode(
            self.getpost("/api/v2/clustering/{0}/cluster_ids".format(clustering_algorithm), timeout=timeout, method="GET")
        )
        return res

    def network(self, clustering_algorithm, cluster_id, timeout=None):
        """ returns a cytoscape compatible network """
        if not isinstance(clustering_algorithm, str):
            raise TypeError("clustering_algorithm must be str not {0}".format(type(clustering_algorithm)))
        if not isinstance(cluster_id, int):
            raise TypeError("cluster_id must be int not {0}".format(type(cluster_id)))

        res = self._decode(
            self.getpost(
                "/api/v2/clustering/{0}/{1}/network".format(clustering_algorithm, cluster_id), timeout=timeout, method="GET"
            )
        )
        return res

    def server_memory_usage(self, nrows=100, timeout=None):
        if not isinstance(nrows, int):
            raise TypeError("nrows must be integer not {0}".format(type(nrows)))
        retVal = self._decode(self.getpost("/api/v2/server_memory_usage/{0}".format(nrows), timeout=timeout, method="GET"))
        return pd.DataFrame.from_records(retVal)

    def server_database_usage(self, nrows=100, timeout=None):
        if not isinstance(nrows, int):
            raise TypeError("nrows must be integer not {0}".format(type(nrows)))
        targeturl = "/api/v2/server_database_usage/{0}".format(nrows)
        retVal = self._decode(self.getpost(targeturl, timeout=timeout, method="GET"))
        if "trend_stats" in retVal.keys():
            retVal["trend_stats"] = pd.DataFrame.from_dict(retVal["trend_stats"], orient="columns")
        return retVal

    def guids_with_quality_over(self, cutoff=0, timeout=None):
        if not type(cutoff) in [int, float]:
            raise TypeError("cutoff must be float or int not {0}".format(type(cutoff)))
        return self._decode(self.getpost("/api/v2/guids_with_quality_over/{0}".format(cutoff), timeout=timeout, method="GET"))

    def guid2neighbours(self, guid, threshold, quality_cutoff=None, timeout=None):
        """ returns a guid2cluster lookup """
        if not isinstance(guid, str):
            raise TypeError("guid must be str not {0}".format(type(guid)))
        if not isinstance(threshold, int):
            raise TypeError("threshold must be int not {0}".format(type(threshold)))
        if quality_cutoff is None:
            quality_cutoff = 0  # no cutoff
        if isinstance(quality_cutoff, float) or isinstance(quality_cutoff, int):
            return self._decode(
                self.getpost(
                    "/api/v2/{0}/neighbours_within/{1}/with_quality_cutoff/{2}".format(guid, threshold, quality_cutoff),
                    timeout=timeout,
                    method="GET",
                )
            )
        else:
            raise TypeError("unhandled: quality_cutoff must be None or float, not {0}".format(type(quality_cutoff)))

    def msa_cluster(self, clustering_algorithm, cluster_id, output_format="json", timeout=None):
        """ does MSA on a cluster """

        res = self.getpost(
            "/api/v2/multiple_alignment_cluster/{0}/{1}/{2}".format(clustering_algorithm, cluster_id, output_format),
            timeout=timeout,
            method="GET",
        )
        if output_format == "json":
            retDict = self._decode(res)
            return pd.DataFrame.from_dict(retDict, orient="index")
        else:
            return res.content

    def msa(self, guids, output_format="json-records", what="N", timeout=None):
        """performs msa

        valid values for 'what', which determines how the p-values are computed, are
        M
        N
        N_or_M
        """
        guidstring = ";".join(guids)
        payload = {"guids": guidstring, "output_format": output_format, "what": what}
        res = self.post("/api/v2/multiple_alignment/guids", payload=payload, timeout=timeout)
        if output_format in ["json"]:
            return res.json()
        elif output_format in ["json-records"]:
            return pd.DataFrame.from_records(res.json())
        else:
            return res.content

    def reset(self, timeout=30):
        """ resets the server to a state with no data """
        return self.post("/api/v2/reset", payload={}, timeout=timeout)

    def insert(self, guid, seq, timeout=None):
        """inserts a sequence seq with guid"""

        # check input
        if not isinstance(guid, str):
            raise TypeError("guid {0} passed must be a string, not a {1}".format(guid, type(guid)))
        if not isinstance(seq, str):
            raise TypeError("sequence passed must be a string, not a {0}".format(type(seq)))
        return self.post("/api/v2/insert", payload={"guid": guid, "seq": seq}, timeout=timeout)

    def read_fasta_file(self, fastafile):
        """reads the content of a fasta file into memory.
        returns a dictionary {seqid:(first part of defline), seq:(nucleic acid), content:(entire file content)}.
        Supports both .gz and uncompressed files transparently.
        Does not support multi-fasta files.  Will raise an error if such are detected.
        """
        # first determine whether it is a .gz file or not; read into RAM.
        if fastafile.endswith(".gz"):
            # we decompress it on the fly.
            with gzip.open(fastafile, "r") as f:
                content = f.read().decode("utf-8")
        else:
            with open(fastafile, "rt") as f:
                content = f.read()

        # use BioPython3 SeqIO library to read the file.
        nFiles = 0
        with io.StringIO(content) as f:
            for record in SeqIO.parse(f, "fasta"):
                nFiles += 1
                if nFiles > 1:  # that's a multifasta, and we don't support that
                    raise ValueError(
                        "Multifasta file is present in {0}.  Multifasta files are not supported".format(fastafile)
                    )
                else:
                    res = {"seq": str(record.seq), "seqid": str(record.id), "content": content}
                    return res
        raise IOError("no content parsed from result of length {0}".format(len(content)))
