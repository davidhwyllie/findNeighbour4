#!/usr/bin/env python
"""  findNeighbour4 is a  server providing relatedness information for bacterial genomes via a Restful API.
See documentation for full details of its functionality.

There are unit tests for the server component.  To run them:

# starting a test RESTFUL server
nohup pipenv run python3 findNeighbour4_server.py &

# And then  launching unit tests with
pipenv run python3 -m unittest test/test_server.py
"""

# import libraries
import os
import requests
import json
import warnings
import datetime
import pandas as pd
import markdown
import codecs

# only used for unit testing
from Bio import SeqIO
import unittest
from urllib.parse import urljoin as urljoiner

# default parameters for unit testing only.
RESTBASEURL = "http://127.0.0.1:5020"
ISDEBUG = True
LISTEN_TO = "127.0.0.1"  # only local addresses


def isjson(content):
    """returns true if content parses as json, otherwise false. used by unit testing."""
    try:
        json.loads(content.decode("utf-8"))
        return True

    except json.decoder.JSONDecodeError:
        return False


def tojson(content):
    """json dumps, formatting dates as isoformat"""

    def converter(o):
        if isinstance(o, datetime.datetime):
            return o.isoformat()
        else:
            return json.JSONEncoder.default(o)

    return json.dumps(content, default=converter)


def do_GET(relpath):
    """makes a GET request  to relpath.
    Used for unit testing."""

    url = urljoiner(RESTBASEURL, relpath)
    # print("GETing from: {0}".format(url))

    session = requests.Session()
    session.trust_env = False

    # print out diagnostics
    # print("About to GET from url {0}".format(url))
    response = session.get(url=url, timeout=None)

    # print("Result:")
    # print("code: {0}".format(response.status_code))
    # print("reason: {0}".format(response.reason))
    try:
        "text: {0}".format(response.text[:100])

    except UnicodeEncodeError:
        # which is what happens if you try to display a gz file as text, which it isn't
        warnings.warn(
            "Response cannot be coerced to unicode ? a gz file.  The response content had {0} bytes.".format(
                len(response.text)
            )
        )
        warnings.warn("headers: {0}".format(response.headers))

    session.close()
    return response


def do_POST(relpath, payload):
    """makes a POST request  to relpath.
    Used for unit testing.
    payload should be a dictionary"""

    url = urljoiner(RESTBASEURL, relpath)

    # print out diagnostics
    # print("POSTING to url {0}".format(url))
    if not isinstance(payload, dict):
        raise TypeError("not a dict {0}".format(payload))
    response = requests.post(url=url, data=payload)

    # print("Result:")
    # print("code: {0}".format(response.status_code))
    # print("reason: {0}".format(response.reason))
    # print("content: {0}".format(response.content))

    return response


def render_markdown(md_file):
    """render markdown as html"""
    with codecs.open(md_file, mode="r", encoding="utf-8") as f:
        text = f.read()
        html = markdown.markdown(text, extensions=["tables"])
    return html


class test_reset(unittest.TestCase):
    """tests route /api/v2/reset and /guids"""

    def runTest(self):
        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))  # get all the guids

        guid_to_insert = "guid_{0}".format(n_pre + 1)

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq)

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})

        self.assertTrue(isjson(content=res.content))
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(info, "Guid {0} inserted.".format(guid_to_insert))

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_post = len(json.loads(str(res.text)))  # get all the guids

        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_post_reset = len(json.loads(str(res.text)))  # get all the guids

        self.assertTrue(n_post > 0)
        self.assertTrue(n_post_reset == 0)


class test_guids(unittest.TestCase):
    """tests routes  /guids, /valid_guids and /invalid_guids"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        self.assertEqual(0, len(json.loads(str(res.text))))  # get all the guids

        guid_to_insert = "valid"

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq)

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})

        seq2 = "".join("N" * len(seq))
        guid_to_insert = "invalid"
        res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq2})

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        self.assertEqual(
            set(["valid", "invalid"]), set(json.loads(str(res.text)))
        )  # get all the guids

        relpath = "/api/v2/valid_guids"
        res = do_GET(relpath)
        self.assertEqual(
            set(["valid"]), set(json.loads(str(res.text)))
        )  # get all the guids

        relpath = "/api/v2/invalid_guids"
        res = do_GET(relpath)
        self.assertEqual(
            set(["invalid"]), set(json.loads(str(res.text)))
        )  # get all the guids


class test_guid_validity(unittest.TestCase):
    """tests routes  /validity_guids"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        self.assertEqual(0, len(json.loads(str(res.text))))  # get all the guids

        guid_to_insert = "valid_guid"

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq)

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})

        seq2 = "".join("N" * len(seq))
        guid_to_insert = "invalid_guid"
        res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq2})

        relpath = "/api/v2/invalid_guid/valid"
        res = do_GET(relpath)
        valid_code = json.loads(str(res.text))
        self.assertEqual(valid_code, 1)

        relpath = "/api/v2/valid_guid/valid"
        res = do_GET(relpath)
        valid_code = json.loads(str(res.text))
        self.assertEqual(valid_code, 0)

        relpath = "/api/v2/missing_guid/valid"
        res = do_GET(relpath)
        valid_code = json.loads(str(res.text))
        self.assertEqual(valid_code, -1)


class test_cl2network(unittest.TestCase):
    """tests return of a change_id number from clustering engine"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        # add four samples, two mixed
        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                originalseq = list(str(record.seq))
        guids_inserted = list()
        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))  # get all the guids

        for i in range(1, 4):

            seq = originalseq
            if i % 2 == 0:
                is_mixed = True
                guid_to_insert = "mixed_{0}".format(n_pre + i)
            else:
                is_mixed = False
                guid_to_insert = "nomix_{0}".format(n_pre + i)
            # make i mutations at position 500,000

            offset = 500000
            for j in range(i):
                mutbase = offset + j
                ref = seq[mutbase]
                if is_mixed is False:
                    if not ref == "T":
                        seq[mutbase] = "T"
                    if not ref == "A":
                        seq[mutbase] = "A"
                if is_mixed is True:
                    seq[mutbase] = "N"
            seq = "".join(seq)
            guids_inserted.append(guid_to_insert)

            relpath = "/api/v2/insert"
            res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})
            self.assertEqual(res.status_code, 200)

        # run the clustering engine.

        os.system("pipenv run python3 findNeighbour4_clustering.py")

        # do tests
        relpath = "/api/v2/clustering/SNV12_ignore/cluster_ids"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        retVal = json.loads(res.text)
        # plot the cluster with the highest clusterid
        relpath = "/api/v2/clustering/SNV12_ignore/{0}/network".format(max(retVal))
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        jsonresp = json.loads(str(res.text))
        self.assertTrue(isinstance(jsonresp, dict))
        self.assertTrue("elements" in jsonresp.keys())

        # plot the cluster with the highest clusterid
        res = None
        relpath = "/api/v2/clustering/SNV12_ignore/{0}/minimum_spanning_tree".format(
            max(retVal)
        )
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        jsonresp = json.loads(str(res.text))
        self.assertTrue(isinstance(jsonresp, dict))
        self.assertTrue("elements" in jsonresp.keys())


class test_msa_2(unittest.TestCase):
    """tests route /api/v2/multiple_alignment/guids, with additional samples."""

    def runTest(self):
        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))  # get all the guids

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                originalseq = list(str(record.seq))
        inserted_guids = ["guid_ref"]
        seq = "".join(originalseq)
        res = do_POST("/api/v2/insert", payload={"guid": "guid_ref", "seq": seq})

        for k in range(0, 1):
            # form one clusters
            for i in range(0, 3):
                guid_to_insert = "msa2_{1}_guid_{0}".format(n_pre + k * 100 + i, k)
                inserted_guids.append(guid_to_insert)
                muts = 0
                seq = originalseq
                # make i mutations at position 500,000
                offset = 500000
                if k == 1:
                    for j in range(1000000, 1000100):  # make 100 mutants at position 1m
                        mutbase = offset + j
                        ref = seq[mutbase]
                        if not ref == "T":
                            seq[mutbase] = "T"
                        if not ref == "A":
                            seq[mutbase] = "A"
                        muts += 1

                for j in range(i):
                    mutbase = offset + j
                    ref = seq[mutbase]
                    if not ref == "T":
                        seq[mutbase] = "T"
                    if not ref == "A":
                        seq[mutbase] = "A"
                    muts += 1
                seq = "".join(seq)

                # print("Adding TB sequence {2} of {0} bytes with {1} mutations relative to ref.".format(len(seq), muts, guid_to_insert))
                self.assertEqual(len(seq), 4411532)  # check it's the right sequence

                relpath = "/api/v2/insert"
                res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})
                self.assertTrue(isjson(content=res.content))
                info = json.loads(res.content.decode("utf-8"))
                self.assertEqual(info, "Guid {0} inserted.".format(guid_to_insert))

        relpath = "/api/v2/multiple_alignment/guids"
        payload = {"guids": ";".join(inserted_guids), "output_format": "html"}
        res = do_POST(relpath, payload=payload)
        self.assertFalse(isjson(res.content))
        self.assertEqual(res.status_code, 200)
        self.assertTrue(b"</table>" in res.content)

        payload = {"guids": ";".join(inserted_guids), "output_format": "json"}
        res = do_POST(relpath, payload=payload)
        self.assertTrue(isjson(res.content))
        self.assertEqual(res.status_code, 200)
        self.assertFalse(b"</table>" in res.content)
        d = json.loads(res.content.decode("utf-8"))
        expected_keys = set(
            [
                "variant_positions",
                "invalid_guids",
                "valid_guids",
                "expected_p1",
                "sample_size",
                "df_dict",
                "what_tested",
                "outgroup",
                "creation_time",
                "fconst",
            ]
        )
        self.assertEqual(set(d.keys()), set(expected_keys))

        payload = {"guids": ";".join(inserted_guids), "output_format": "json-records"}
        res = do_POST(relpath, payload=payload)
        self.assertTrue(isjson(res.content))
        self.assertEqual(res.status_code, 200)
        self.assertFalse(b"</table>" in res.content)
        d = json.loads(res.content.decode("utf-8"))

        payload = {"guids": ";".join(inserted_guids), "output_format": "fasta"}
        res = do_POST(relpath, payload=payload)
        self.assertFalse(isjson(res.content))
        self.assertEqual(res.status_code, 200)

        payload = {"guids": ";".join(inserted_guids), "output_format": "json-fasta"}
        res = do_POST(relpath, payload=payload)
        self.assertTrue(isjson(res.content))
        self.assertEqual(res.status_code, 200)
        retVal = json.loads(res.content.decode("utf_8"))
        self.assertTrue(isinstance(retVal, dict))
        self.assertEqual(set(retVal.keys()), set(["fasta"]))

        relpath = "/api/v2/multiple_alignment/guids"
        payload = {
            "guids": ";".join(inserted_guids),
            "output_format": "html",
            "what": "N",
        }
        res = do_POST(relpath, payload=payload)
        self.assertFalse(isjson(res.content))
        self.assertEqual(res.status_code, 200)
        self.assertTrue(b"</table>" in res.content)

        relpath = "/api/v2/multiple_alignment/guids"
        payload = {
            "guids": ";".join(inserted_guids),
            "output_format": "html",
            "what": "M",
        }
        res = do_POST(relpath, payload=payload)
        self.assertFalse(isjson(res.content))
        self.assertEqual(res.status_code, 200)
        self.assertTrue(b"</table>" in res.content)

        relpath = "/api/v2/multiple_alignment/guids"
        payload = {
            "guids": ";".join(inserted_guids),
            "output_format": "html",
            "what": "N_or_M",
        }
        res = do_POST(relpath, payload=payload)
        self.assertFalse(isjson(res.content))
        self.assertEqual(res.status_code, 200)
        self.assertTrue(b"</table>" in res.content)

        relpath = "/api/v2/multiple_alignment/guids"
        payload = {
            "guids": ";".join(inserted_guids),
            "output_format": "interactive",
            "what": "N_or_M",
        }
        res = do_POST(relpath, payload=payload)
        self.assertFalse(isjson(res.content))
        self.assertEqual(res.status_code, 200)
        self.assertTrue(b"</html>" in res.content)

        relpath = "/api/v2/multiple_alignment/guids"
        payload = {
            "guids": ";".join(inserted_guids),
            "output_format": "json-records",
            "what": "N",
        }
        res = do_POST(relpath, payload=payload)
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(res.content))
        d = json.loads(res.content.decode("utf-8"))
        df = pd.DataFrame.from_records(d)
        self.assertEqual(df.loc[df.index[0], "what_tested"], "N")

        relpath = "/api/v2/multiple_alignment/guids"
        payload = {
            "guids": ";".join(inserted_guids),
            "output_format": "json-records",
            "what": "M",
        }
        res = do_POST(relpath, payload=payload)
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(res.content))
        d = json.loads(res.content.decode("utf-8"))
        df = pd.DataFrame.from_records(d)

        self.assertEqual(df.loc[df.index[0], "what_tested"], "M")

        relpath = "/api/v2/multiple_alignment/guids"
        payload = {
            "guids": ";".join(inserted_guids),
            "output_format": "json-records",
            "what": "N_or_M",
        }
        res = do_POST(relpath, payload=payload)
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(res.content))
        d = json.loads(res.content.decode("utf-8"))
        df = pd.DataFrame.from_records(d)
        self.assertEqual(df.loc[df.index[0], "what_tested"], "N_or_M")

        relpath = "/api/v2/multiple_alignment/guids"
        payload = {
            "guids": ";".join(inserted_guids),
            "output_format": "json-records",
            "what": "N",
        }
        res = do_POST(relpath, payload=payload)
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(res.content))
        d = json.loads(res.content.decode("utf-8"))
        df = pd.DataFrame.from_records(d)
        self.assertEqual(df.loc[df.index[0], "what_tested"], "N")

        relpath = "/api/v2/multiple_alignment/guids"
        payload = {
            "guids": ";".join(inserted_guids),
            "output_format": "html",
            "what": "X",
        }
        res = do_POST(relpath, payload=payload)
        self.assertEqual(res.status_code, 404)

        # Do clustering
        # "Doing clustering")
        os.system("pipenv run python3 findNeighbour4_clustering.py")

        relpath = "/api/v2/clustering/SNV12_ignore/guids2clusters"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        retVal = json.loads(str(res.text))
        self.assertTrue(isinstance(retVal, list))
        res = json.loads(res.content.decode("utf-8"))
        cluster_id = None

        for item in res:
            if item["guid"] in inserted_guids:
                cluster_id = item["cluster_id"]
        # print("Am examining cluster_id",cluster_id)
        self.assertTrue(cluster_id is not None)
        relpath = (
            "/api/v2/multiple_alignment_cluster/SNV12_ignore/{0}/json-records".format(
                cluster_id
            )
        )
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(res.content))
        d = json.loads(res.content.decode("utf-8"))
        df = pd.DataFrame.from_records(d)
        self.assertEqual(df.loc[df.index[0], "what_tested"], "M")

        relpath = (
            "/api/v2/multiple_alignment_cluster/SNV12_include/{0}/json-records".format(
                cluster_id
            )
        )
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(res.content))
        d = json.loads(res.content.decode("utf-8"))
        df = pd.DataFrame.from_records(d)
        self.assertEqual(df.loc[df.index[0], "what_tested"], "M")

        relpath = "/api/v2/multiple_alignment_cluster/SNV12_exclude/{0}/fasta".format(
            cluster_id
        )
        res = do_GET(relpath)
        self.assertFalse(isjson(res.content))
        self.assertEqual(res.status_code, 200)


class test_msa_1(unittest.TestCase):
    """tests route /api/v2/multiple_alignment/guids, with additional samples."""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})
        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))  # get all the guids
        # print("There are {0} existing samples".format(n_pre))
        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                originalseq = list(str(record.seq))
        inserted_guids = []
        for i in range(0, 3):
            guid_to_insert = "msa1_guid_{0}".format(n_pre + i)
            inserted_guids.append(guid_to_insert)

            seq = originalseq
            # make i mutations at position 500,000
            offset = 500000
            for j in range(i):
                mutbase = offset + j
                ref = seq[mutbase]
                if not ref == "T":
                    seq[mutbase] = "T"
                if not ref == "A":
                    seq[mutbase] = "A"
            seq = "".join(seq)

            # print("Adding TB sequence {2} of {0} bytes with {1} mutations relative to ref.".format(len(seq), i, guid_to_insert))
            self.assertEqual(len(seq), 4411532)  # check it's the right sequence

            relpath = "/api/v2/insert"
            res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})
            self.assertEqual(res.status_code, 200)

            self.assertTrue(isjson(content=res.content))
            info = json.loads(res.content.decode("utf-8"))
            self.assertEqual(info, "Guid {0} inserted.".format(guid_to_insert))

        relpath = "/api/v2/multiple_alignment/guids"
        payload = {"guids": ";".join(inserted_guids), "output_format": "html"}
        res = do_POST(relpath, payload=payload)
        self.assertFalse(isjson(res.content))
        self.assertEqual(res.status_code, 200)
        self.assertTrue(b"</table>" in res.content)

        payload = {"guids": ";".join(inserted_guids), "output_format": "json"}
        res = do_POST(relpath, payload=payload)
        self.assertTrue(isjson(res.content))
        self.assertEqual(res.status_code, 200)
        self.assertFalse(b"</table>" in res.content)
        d = json.loads(res.content.decode("utf-8"))
        expected_keys = set(
            [
                "variant_positions",
                "invalid_guids",
                "valid_guids",
                "expected_p1",
                "sample_size",
                "df_dict",
                "what_tested",
                "outgroup",
                "creation_time",
                "fconst",
            ]
        )
        self.assertEqual(set(d.keys()), set(expected_keys))

        payload = {"guids": ";".join(inserted_guids), "output_format": "fasta"}
        res = do_POST(relpath, payload=payload)
        self.assertFalse(isjson(res.content))
        self.assertEqual(res.status_code, 200)


class test_server_config(unittest.TestCase):
    """tests route v2/server_config"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/server_config"
        res = do_GET(relpath)

        self.assertTrue(isjson(content=res.content))

        config_dict = json.loads(res.content.decode("utf-8"))
        self.assertTrue("PRECOMPARER_PARAMETERS" in config_dict.keys())
        self.assertEqual(res.status_code, 200)


class test_server_memory_usage(unittest.TestCase):
    """tests route /api/v2/server_memory_usage"""

    def runTest(self):

        relpath = "/api/v2/server_memory_usage"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content=res.content))

        res = json.loads(res.content.decode("utf-8"))
        self.assertTrue(isinstance(res, list))


class test_server_database_usage(unittest.TestCase):
    """tests route /api/v2/server_database_usage"""

    def runTest(self):

        relpath = "/api/v2/server_database_usage"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content=res.content))

        res = json.loads(res.content.decode("utf-8"))
        self.assertTrue(isinstance(res, dict))
        self.assertEqual(res["latest_stats"]["storage_ratio"], 1)


class test_snpceiling(unittest.TestCase):
    """tests route /api/v2/snpceiling"""

    def runTest(self):
        res = "/api/v2/reset"
        relpath = "/api/v2/snpceiling"
        res = do_POST(relpath, payload={})

        res = do_GET(relpath)
        self.assertTrue(isjson(content=res.content))
        config_dict = json.loads(res.content.decode("utf-8"))
        self.assertTrue("snpceiling" in config_dict.keys())
        self.assertEqual(res.status_code, 200)


class test_server_time(unittest.TestCase):
    """tests route /api/v2/server_time"""

    def runTest(self):
        relpath = "/api/v2/server_time"
        res = do_GET(relpath)
        # print(res)
        self.assertTrue(isjson(content=res.content))
        config_dict = json.loads(res.content.decode("utf-8"))
        self.assertTrue("server_time" in config_dict.keys())
        self.assertEqual(res.status_code, 200)


class test_server_name(unittest.TestCase):
    """tests route /api/v2/server_name"""

    def runTest(self):
        relpath = "/api/v2/server_name"
        res = do_GET(relpath)

        self.assertTrue(isjson(content=res.content))
        config_dict = json.loads(res.content.decode("utf-8"))
        self.assertTrue("server_name" in config_dict.keys())
        self.assertEqual(res.status_code, 200)


class test_get_all_guids_1(unittest.TestCase):
    """tests route /api/v2/guids"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        self.assertTrue(isjson(content=res.content))
        guidlist = json.loads(str(res.content.decode("utf-8")))
        self.assertTrue(isinstance(guidlist, list))
        self.assertEqual(res.status_code, 200)
        ## TODO: insert guids, check it doesn't fail.


class test_guids_with_quality_over_1(unittest.TestCase):
    """tests route /api/v2/guids_with_quality_over"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/guids_with_quality_over/0.7"
        res = do_GET(relpath)
        self.assertTrue(isjson(content=res.content))
        guidlist = json.loads(res.content.decode("utf-8"))
        self.assertTrue(isinstance(guidlist, list))
        self.assertEqual(res.status_code, 200)


class test_get_all_guids_examination_time_1(unittest.TestCase):
    """tests route /api/v2/guids_and_examination_times"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/guids_and_examination_times"
        res = do_GET(relpath)
        self.assertTrue(isjson(content=res.content))
        guidlist = json.loads(res.content.decode("utf-8"))

        self.assertTrue(isinstance(guidlist, dict))
        self.assertEqual(res.status_code, 200)

        #  test that it actually works
        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))  # get all the guids

        guid_to_insert = "guid_{0}".format(n_pre + 1)

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq)

        # print("Adding TB reference sequence of {0} bytes".format(len(seq)))
        self.assertEqual(len(seq), 4411532)  # check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})
        self.assertEqual(res.status_code, 200)

        self.assertTrue(isjson(content=res.content))
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(info, "Guid {0} inserted.".format(guid_to_insert))

        relpath = "/api/v2/guids_and_examination_times"
        res = do_GET(relpath)
        et = len(json.loads(res.content.decode("utf-8")))
        self.assertTrue(et > 0)


class test_get_matching_guids_1(unittest.TestCase):
    """tests route /api/v2/guids_beginning_with"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        #  get existing guids
        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))  # get all the guids

        guid_to_insert = "guid_{0}".format(n_pre + 1)

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq)

        # print("Adding TB reference sequence of {0} bytes".format(len(seq)))
        self.assertEqual(len(seq), 4411532)  # check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})
        self.assertEqual(res.status_code, 200)

        self.assertTrue(isjson(content=res.content))
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(info, "Guid {0} inserted.".format(guid_to_insert))

        relpath = "/api/v2/guids_beginning_with/{0}".format(guid_to_insert)
        res = do_GET(relpath)
        self.assertEqual(json.loads(res.content.decode("utf-8")), [guid_to_insert])


class test_annotations_1(unittest.TestCase):
    """tests route /api/v2/annotations"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/annotations"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content=res.content))
        inputDict = json.loads(res.content.decode("utf-8"))
        self.assertTrue(isinstance(inputDict, dict))
        guiddf = pd.DataFrame.from_dict(inputDict, orient="index")  # , orient='index'
        self.assertTrue(isinstance(guiddf, pd.DataFrame))


class test_exist_sample(unittest.TestCase):
    """tests route /api/v2/guid/exists"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/non_existent_guid/exists"
        res = do_GET(relpath)

        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content=res.content))
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(type(info), bool)
        self.assertEqual(info, False)


class test_clusters_sample(unittest.TestCase):
    """tests route /api/v2/guid/clusters"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        # what happens if there is nothing there
        relpath = "/api/v2/non_existent_guid/clusters"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content=res.content))

        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(type(info), list)
        self.assertEqual(info, [])
        # add one
        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))  # get all the guids

        guid_to_insert = "guid_{0}".format(n_pre + 1)

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq)

        # print("Adding TB reference sequence of {0} bytes".format(len(seq)))
        self.assertEqual(len(seq), 4411532)  # check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content=res.content))
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(info, "Guid {0} inserted.".format(guid_to_insert))

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        guids_loaded = json.loads(res.content.decode("utf-8"))
        n_post = len(guids_loaded)
        # print("*** GUIDS LOADED ", guids_loaded)
        self.assertEqual(n_pre + 1, n_post)

        relpath = "/api/v2/{0}/clusters".format(guid_to_insert)
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)  # clustering hasn't happened yet
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(info, [])

        # Do clustering
        os.system("pipenv run python3 findNeighbour4_clustering.py")

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        guids_loaded = json.loads(res.content.decode("utf-8"))
        n_post = len(guids_loaded)
        # print("*** GUIDS RELOADED ", guids_loaded)
        self.assertEqual(n_pre + 1, n_post)

        relpath = "/api/v2/{0}/clusters".format(guid_to_insert)
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)  # clustering has happened

        self.assertTrue(isjson(content=res.content))
        cluster_list = res.content.decode("utf-8")
        info = json.loads(cluster_list)
        # print("CLUSTERINFO",info)
        self.assertEqual(len(info), 4)
        # print(info)


class test_clusters_what(unittest.TestCase):
    """tests implementation of 'what' value, stored in clustering results object"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        # add one
        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))  # get all the guids

        guid_to_insert = "guid_{0}".format(n_pre + 1)

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq)

        # print("Adding TB reference sequence of {0} bytes".format(len(seq)))
        self.assertEqual(len(seq), 4411532)  # check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content=res.content))
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(info, "Guid {0} inserted.".format(guid_to_insert))

        # Do clustering
        os.system("pipenv run python3 findNeighbour4_clustering.py")

        # what happens if there is nothing there
        relpath = "/api/v2/non_existent_guid/clusters"
        res = do_GET(relpath)

        self.assertTrue(isjson(content=res.content))
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(info, [])


class test_annotation_sample(unittest.TestCase):
    """tests route /api/v2/guid/annotation"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/non_existent_guid/annotation"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 404)
        self.assertTrue(isjson(content=res.content))
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(type(info), dict)


class test_algorithms(unittest.TestCase):
    """tests return of a change_id number"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/clustering"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        retDict = json.loads(str(res.text))
        self.assertEqual(
            set(retDict["algorithms"]),
            set(["SNV12_ignore", "SNV12_include", "SNV12_exclude", "SNV12_include_n"]),
        )


class test_what_tested(unittest.TestCase):
    """tests return of what is tested"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        # add one
        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))  # get all the guids

        guid_to_insert = "guid_{0}".format(n_pre + 1)

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq)

        # print("Adding TB reference sequence of {0} bytes".format(len(seq)))
        self.assertEqual(len(seq), 4411532)  # check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content=res.content))
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(info, "Guid {0} inserted.".format(guid_to_insert))

        # Do clustering
        os.system("pipenv run python3 findNeighbour4_clustering.py")

        relpath = "/api/v2/clustering/SNV12_exclude/what_tested"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        retDict = json.loads(str(res.text))
        self.assertEqual(
            retDict, {"clustering_algorithm": "SNV12_exclude", "what_tested": "M"}
        )

        relpath = "/api/v2/clustering/SNV12_include/what_tested"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        retDict = json.loads(str(res.text))
        self.assertEqual(
            retDict, {"clustering_algorithm": "SNV12_include", "what_tested": "M"}
        )

        relpath = "/api/v2/clustering/SNV12_ignore/what_tested"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        retDict = json.loads(str(res.text))
        self.assertEqual(
            retDict, {"clustering_algorithm": "SNV12_ignore", "what_tested": "M"}
        )


class test_g2c(unittest.TestCase):
    """tests return of guid2clusters data structure"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        # add two
        relpath = "/api/v2/guids"
        for i in range(3):

            n_pre = len(json.loads(str(res.text)))  # get all the guids
            res = do_GET(relpath)
            guid_to_insert = "guid_{0}".format(n_pre + 1)

            inputfile = "COMPASS_reference/R39/R00000039.fasta"
            with open(inputfile, "rt") as f:
                for record in SeqIO.parse(f, "fasta"):
                    seq = str(record.seq)

            self.assertEqual(len(seq), 4411532)  # check it's the right sequence

            relpath = "/api/v2/insert"
            res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})
            self.assertEqual(res.status_code, 200)
            self.assertTrue(isjson(content=res.content))
            info = json.loads(res.content.decode("utf-8"))
            self.assertEqual(info, "Guid {0} inserted.".format(guid_to_insert))

        # Do clustering
        os.system("pipenv run python3 findNeighbour4_clustering.py")

        relpath = "/api/v2/clustering/SNV12_ignore/guids2clusters"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        retVal = json.loads(str(res.text))
        self.assertTrue(isinstance(retVal, list))

        relpath = "/api/v2/clustering/SNV12_ignore/guids2clusters/after_change_id/0"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        retVal = json.loads(str(res.text))
        self.assertEqual(len(retVal), 3)

        self.assertEqual(retVal[0]["clusterSize"], 3)

        self.assertTrue(isinstance(retVal, list))


class test_clusters2cnt(unittest.TestCase):
    """tests return of guid2clusters data structure"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        # add one
        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))  # get all the guids

        guid_to_insert = "guid_{0}".format(n_pre + 1)

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq)

        # print("Adding TB reference sequence of {0} bytes".format(len(seq)))
        self.assertEqual(len(seq), 4411532)  # check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content=res.content))
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(info, "Guid {0} inserted.".format(guid_to_insert))

        # Do clustering
        os.system("pipenv run python3 findNeighbour4_clustering.py")

        relpath = "/api/v2/clustering/SNV12_ignore/clusters"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        retVal = json.loads(str(res.text))
        self.assertTrue(isinstance(retVal, dict))
        self.assertEqual(set(retVal.keys()), set(["summary", "members"]))

        relpath = "/api/v2/clustering/SNV12_ignore/members"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        retVal = json.loads(str(res.text))
        self.assertTrue(isinstance(retVal, dict))
        self.assertEqual(set(retVal.keys()), set(["members"]))
        # print(retVal)
        relpath = "/api/v2/clustering/SNV12_ignore/summary"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        retVal = json.loads(str(res.text))
        self.assertTrue(isinstance(retVal, dict))
        self.assertEqual(set(retVal.keys()), set(["summary"]))


class test_cluster2cnt1(unittest.TestCase):
    """tests return of guid2clusters data structure"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/clustering/SNV12_ignore/0"  # doesn't exist
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 404)

        # get existing clusterids
        relpath = "/api/v2/clustering/SNV12_ignore/clusters"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        retVal = json.loads(str(res.text))
        self.assertTrue(isinstance(retVal, dict))
        valid_cluster_ids = set()
        for item in retVal["members"]:
            valid_cluster_ids.add(item["cluster_id"])

        for this_cluster_id in valid_cluster_ids:
            relpath = "/api/v2/clustering/SNV12_ignore/{0}".format(
                this_cluster_id
            )  # may exist
            res = do_GET(relpath)
            self.assertEqual(res.status_code, 200)

            retVal = json.loads(str(res.text))
            self.assertTrue(isinstance(retVal, dict))
            self.assertTrue(len(retVal["summary"]) == 1)
            self.assertTrue(len(retVal["members"]) > 0)
            self.assertEqual(set(retVal.keys()), set(["summary", "members"]))
            break


class test_g2cl(unittest.TestCase):
    """tests return of a change_id number"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/clustering/SNV12_ignore/cluster_ids"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        retVal = json.loads(str(res.text))
        self.assertTrue(isinstance(retVal, list))


class test_g2ca(unittest.TestCase):
    """tests return of a change_id number"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/clustering/SNV12_ignore/guids2clusters"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)

        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))  # get all the guids

        guid_to_insert = "guid_{0}".format(n_pre + 1)

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq)

        # print("Adding TB reference sequence of {0} bytes".format(len(seq)))
        self.assertEqual(len(seq), 4411532)  # check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content=res.content))
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(info, "Guid {0} inserted.".format(guid_to_insert))

        # Do clustering
        os.system("pipenv run python3 findNeighbour4_clustering.py")

        relpath = "/api/v2/clustering/SNV12_ignore/guids2clusters"

        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        retVal = json.loads(str(res.text))
        self.assertTrue(isinstance(retVal, list))


class test_change_id(unittest.TestCase):
    """tests return of a change_id number"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/clustering/SNV12_ignore/change_id"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        retDict = json.loads(str(res.text))
        self.assertEqual(
            set(retDict.keys()), set(["change_id", "clustering_algorithm"])
        )
        self.assertEqual(retDict["clustering_algorithm"], "SNV12_ignore")

        relpath = "/api/v2/clustering/SNV12_ignore/change_id"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        retDict = json.loads(str(res.text))
        self.assertEqual(
            set(retDict.keys()), set(["change_id", "clustering_algorithm"])
        )
        self.assertEqual(retDict["clustering_algorithm"], "SNV12_ignore")

        relpath = "/api/v2/clustering/not_exists/change_id"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 404)


class test_compare_two(unittest.TestCase):
    """tests route /api/v2/{guid1}/{guid2}/exact_distance"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq)

        # make variant
        vseq = list(seq)
        for i in range(100):
            if not vseq[100000 + i] == "A":
                vseq[100000 + i] = "A"
            else:
                vseq[100000 + i] = "T"
        vseq = "".join(vseq)
        # print("Adding TB reference sequence of {0} bytes".format(len(seq)))
        self.assertEqual(len(seq), 4411532)  # check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload={"guid": "one", "seq": seq})
        self.assertEqual(res.status_code, 200)

        res = do_POST(relpath, payload={"guid": "two", "seq": seq})
        self.assertEqual(res.status_code, 200)

        res = do_POST(relpath, payload={"guid": "three", "seq": vseq})
        self.assertEqual(res.status_code, 200)

        relpath = "/api/v2/one/two/exact_distance"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)

        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual({"guid1": "one", "guid2": "two", "dist": 0}, info)
        relpath = "/api/v2/one/three/exact_distance"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)

        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual({"guid1": "one", "guid2": "three", "dist": 100}, info)

        relpath = "/api/v2/one/four/exact_distance"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 404)


class test_insert_1(unittest.TestCase):
    """tests route /api/v2/insert"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))  # get all the guids

        guid_to_insert = "guid_{0}".format(n_pre + 1)

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq)

        # print("Adding TB reference sequence of {0} bytes".format(len(seq)))
        self.assertEqual(len(seq), 4411532)  # check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(content=res.content))
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(info, "Guid {0} inserted.".format(guid_to_insert))

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_post = len(json.loads(res.content.decode("utf-8")))
        self.assertEqual(n_pre + 1, n_post)

        # check if it exists
        relpath = "/api/v2/{0}/exists".format(guid_to_insert)
        res = do_GET(relpath)
        self.assertTrue(isjson(content=res.content))
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(type(info), bool)
        self.assertEqual(res.status_code, 200)
        self.assertEqual(info, True)


class test_insert_10(unittest.TestCase):
    """tests route /api/v2/insert, with additional samples.
    Also provides a set of very similar samples, testing recompression code."""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))  # get all the guids

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                originalseq = list(str(record.seq))

        for i in range(1, 10):
            guid_to_insert = "guid_{0}".format(n_pre + i)

            seq = originalseq
            # make i mutations at position 500,000
            offset = 500000
            for j in range(i):
                mutbase = offset + j
                ref = seq[mutbase]
                if not ref == "T":
                    seq[mutbase] = "T"
                if not ref == "A":
                    seq[mutbase] = "A"
            seq = "".join(seq)

            # print("Adding TB sequence {2} of {0} bytes with {1} mutations relative to ref.".format(len(seq), i, guid_to_insert))
            self.assertEqual(len(seq), 4411532)  # check it's the right sequence

            relpath = "/api/v2/insert"
            res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})
            self.assertTrue(isjson(content=res.content))
            info = json.loads(res.content.decode("utf-8"))
            self.assertEqual(info, "Guid {0} inserted.".format(guid_to_insert))

            relpath = "/api/v2/guids"
            res = do_GET(relpath)
            n_post = len(json.loads(res.content.decode("utf-8")))
            self.assertEqual(n_pre + i, n_post)

            # check if it exists
            relpath = "/api/v2/{0}/exists".format(guid_to_insert)
            res = do_GET(relpath)
            self.assertTrue(isjson(content=res.content))
            info = json.loads(res.content.decode("utf-8"))
            self.assertEqual(type(info), bool)
            self.assertEqual(res.status_code, 200)
            self.assertEqual(info, True)


class test_insert_10a(unittest.TestCase):
    """tests route /api/v2/insert, with additional samples.
    Also provides a set of very similar samples, testing mixture addition."""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))  # get all the guids

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                originalseq = list(str(record.seq))

        for i in range(1, 10):
            guid_to_insert = "guid_{0}".format(n_pre + i)

            seq = originalseq
            # make i mutations at position 500,000
            offset = 500000
            for j in range(i):
                mutbase = offset + j
                ref = seq[mutbase]
                if not ref == "T":
                    seq[mutbase] = "T"
                if not ref == "A":
                    seq[mutbase] = "M"
            seq = "".join(seq)

            # print("Adding TB sequence {2} of {0} bytes with {1} mutations relative to ref.".format(len(seq), i, guid_to_insert))
            self.assertEqual(len(seq), 4411532)  # check it's the right sequence

            relpath = "/api/v2/insert"
            res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})
            self.assertTrue(isjson(content=res.content))
            info = json.loads(res.content.decode("utf-8"))
            self.assertEqual(info, "Guid {0} inserted.".format(guid_to_insert))

            relpath = "/api/v2/guids"
            res = do_GET(relpath)
            n_post = len(json.loads(res.content.decode("utf-8")))
            self.assertEqual(n_pre + i, n_post)

            # check if it exists
            relpath = "/api/v2/{0}/exists".format(guid_to_insert)
            res = do_GET(relpath)
            self.assertTrue(isjson(content=res.content))
            info = json.loads(res.content.decode("utf-8"))
            self.assertEqual(type(info), bool)
            self.assertEqual(res.status_code, 200)
            self.assertEqual(info, True)


# @unittest.skip("skipped; to investigate if this causes a timeout")
class test_insert_60(unittest.TestCase):
    """tests route /api/v2/insert, with additional samples.
    Also provides a set of very similar samples, testing recompression code."""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))  # get all the guids

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                originalseq = list(str(record.seq))
        guids_inserted = list()
        for i in range(1, 40):

            seq = originalseq
            if i % 5 == 0:
                is_mixed = True
                guid_to_insert = "mixed_{0}".format(n_pre + i)
            else:
                is_mixed = False
                guid_to_insert = "nomix_{0}".format(n_pre + i)
            # make i mutations at position 500,000

            offset = 500000
            for j in range(i):
                mutbase = offset + j
                ref = seq[mutbase]
                if is_mixed is False:
                    if not ref == "T":
                        seq[mutbase] = "T"
                    if not ref == "A":
                        seq[mutbase] = "A"
                if is_mixed is True:
                    seq[mutbase] = "N"
            seq = "".join(seq)
            guids_inserted.append(guid_to_insert)
            if is_mixed:
                pass
                # print("Adding TB sequence {2} of {0} bytes with {1} mutations relative to ref.".format(len(seq), i, guid_to_insert))
            else:
                # print("Adding mixed TB sequence {2} of {0} bytes with {1} Ns relative to ref.".format(len(seq), i, guid_to_insert))
                pass

            self.assertEqual(len(seq), 4411532)  # check it's the right sequence

            relpath = "/api/v2/insert"
            res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})
            self.assertEqual(res.status_code, 200)

            self.assertTrue(isjson(content=res.content))
            info = json.loads(res.content.decode("utf-8"))
            self.assertEqual(info, "Guid {0} inserted.".format(guid_to_insert))

            relpath = "/api/v2/guids"
            res = do_GET(relpath)
            n_post = len(json.loads(res.content.decode("utf-8")))
            self.assertEqual(n_pre + i, n_post)

            # check if it exists
            relpath = "/api/v2/{0}/exists".format(guid_to_insert)
            res = do_GET(relpath)
            self.assertTrue(isjson(content=res.content))
            info = json.loads(res.content.decode("utf-8"))
            self.assertEqual(type(info), bool)
            self.assertEqual(res.status_code, 200)
            self.assertEqual(info, True)

        # check: is everything there?
        for guid in guids_inserted:
            relpath = "/api/v2/{0}/exists".format(guid)
            res = do_GET(relpath)
            self.assertTrue(isjson(content=res.content))
            info = json.loads(res.content.decode("utf-8"))
            self.assertEqual(type(info), bool)
            self.assertEqual(res.status_code, 200)
            self.assertEqual(info, True)

        # is everything clustered?
        relpath = "/api/v2/clustering/SNV12_ignore/guids2clusters"
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)
        retVal = json.loads(str(res.text))
        self.assertTrue(isinstance(retVal, list))

        # generate MSA
        relpath = "/api/v2/multiple_alignment/guids"
        payload = {"guids": ";".join(guids_inserted), "output_format": "json-records"}
        res = do_POST(relpath, payload=payload)
        self.assertEqual(res.status_code, 200)
        self.assertTrue(isjson(res.content))
        d = json.loads(res.content.decode("utf-8"))
        pd.DataFrame.from_records(d)

        # print("running mixed checks:")
        for item in retVal:
            if "mixed_" in item["guid"]:
                # print(item['guid'], item['is_mixed'])
                # self.assertTrue(item['is_mixed'])
                pass


class test_mirror(unittest.TestCase):
    """tests route /api/v2/mirror"""

    def runTest(self):

        relpath = "/api/v2/mirror"
        payload = {"guid": "1", "seq": "ACTG"}
        res = do_POST(relpath, payload=payload)
        res_dict = json.loads(res.content.decode("utf-8"))
        self.assertEqual(payload, res_dict)


class test_neighbours_within_1(unittest.TestCase):
    """tests route /api/v2/guid/neighbours_within/"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/non_existent_guid/neighbours_within/12"
        res = do_GET(relpath)
        self.assertTrue(isjson(content=res.content))
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(type(info), dict)
        self.assertEqual(res.status_code, 404)


class test_neighbours_within_2(unittest.TestCase):
    """tests route /api/v2/guid/neighbours_within/"""

    def runTest(self):
        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = (
            "/api/v2/non_existent_guid/neighbours_within/12/with_quality_cutoff/0.5"
        )
        res = do_GET(relpath)
        self.assertTrue(isjson(content=res.content))
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(type(info), dict)
        self.assertEqual(res.status_code, 404)


class test_neighbours_within_3(unittest.TestCase):
    """tests route /api/v2/guid/neighbours_within/"""

    def runTest(self):

        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/non_existent_guid/neighbours_within/12/with_quality_cutoff/0.5/in_format/1"
        res = do_GET(relpath)

        self.assertTrue(isjson(content=res.content))
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(type(info), dict)
        self.assertEqual(res.status_code, 404)


class test_neighbours_within_4(unittest.TestCase):
    """tests route /api/v2/guid/neighbours_within/"""

    def runTest(self):

        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/non_existent_guid/neighbours_within/12/with_quality_cutoff/0.5/in_format/1"
        res = do_GET(relpath)

        self.assertTrue(isjson(content=res.content))
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(type(info), dict)
        self.assertEqual(res.status_code, 404)


class test_neighbours_within_5(unittest.TestCase):
    """tests route /api/v2/guid/neighbours_within/"""

    def runTest(self):

        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/non_existent_guid/neighbours_within/12/in_format/1"
        res = do_GET(relpath)
        # print(res)
        self.assertTrue(isjson(content=res.content))
        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(type(info), dict)
        self.assertEqual(res.status_code, 404)


class test_neighbours_within_6(unittest.TestCase):
    """tests all the /api/v2/guid/neighbours_within methods using test data"""

    def runTest(self):

        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq)

        # generate variants
        variants = {}
        for i in range(4):
            guid_to_insert = "guid_insert_{0}".format(n_pre + i + 1)
            vseq = list(seq)
            vseq[100 * i] = "A"
            vseq = "".join(vseq)
            variants[guid_to_insert] = vseq

        for guid_to_insert in variants.keys():

            # print("Adding mutated TB reference sequence called {0}".format(guid_to_insert))
            relpath = "/api/v2/insert"

            res = do_POST(
                relpath,
                payload={"guid": guid_to_insert, "seq": variants[guid_to_insert]},
            )
            self.assertTrue(isjson(content=res.content))
            info = json.loads(res.content.decode("utf-8"))
            self.assertTrue("inserted" in info)

            # check if it exists
            relpath = "/api/v2/{0}/exists".format(guid_to_insert)
            res = do_GET(relpath)
            self.assertTrue(isjson(content=res.content))
            info = json.loads(res.content.decode("utf-8"))
            self.assertEqual(type(info), bool)
            self.assertEqual(res.status_code, 200)
            self.assertEqual(info, True)

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_post = len(json.loads(res.content.decode("utf-8")))
        self.assertEqual(n_pre + 4, n_post)

        test_guid = min(variants.keys())
        # print("Searching for ",test_guid)

        search_paths = [
            "/api/v2/{0}/neighbours_within/1",
            "/api/v2/{0}/neighbours_within/1/with_quality_cutoff/0.5",
            "/api/v2/{0}/neighbours_within/1/with_quality_cutoff/0.5/in_format/1",
            "/api/v2/{0}/neighbours_within/1/with_quality_cutoff/0.5/in_format/3",
            "/api/v2/{0}/neighbours_within/1/with_quality_cutoff/0.5/in_format/4",
            "/api/v2/{0}/neighbours_within/1/in_format/1",
            "/api/v2/{0}/neighbours_within/1/in_format/3",
            "/api/v2/{0}/neighbours_within/1/in_format/4",
        ]

        for search_path in search_paths:

            url = search_path.format(test_guid)
            res = do_GET(url)
            self.assertTrue(isjson(content=res.content))
            info = json.loads(res.content.decode("utf-8"))
            self.assertEqual(type(info), list)
            guids_found = set()
            for item in info:
                if isinstance(item, list):
                    guids_found.add(item[0])
                elif isinstance(item, dict):
                    guids_found.add(item["guid"])
                elif isinstance(item, str):
                    guids_found.add(item)
                else:
                    self.fail("Unknown class returned {0}".format(type(item)))
            recovered = guids_found.intersection(variants.keys())
            self.assertEqual(len(recovered), 3)
            self.assertEqual(res.status_code, 200)


class test_sequence_1(unittest.TestCase):
    """tests route /api/v2/*guid*/sequence"""

    def runTest(self):

        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        n_pre = len(json.loads(str(res.text)))  # get all the guids

        guid_to_insert = "guid_{0}".format(n_pre + 1)

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq)

        # print("Adding TB reference sequence of {0} bytes".format(len(seq)))
        self.assertEqual(len(seq), 4411532)  # check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})
        self.assertEqual(res.status_code, 200)

        relpath = "/api/v2/{0}/sequence".format(guid_to_insert)
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)

        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(info["guid"], guid_to_insert)
        self.assertEqual(info["invalid"], 0)
        self.assertEqual(info["masked_dna"].count("N"), 557291)


class test_sequence_2(unittest.TestCase):
    """tests route /api/v2/*guid*/sequence"""

    def runTest(self):

        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/{0}/sequence".format("no_guid_exists")
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 404)


class test_sequence_3(unittest.TestCase):
    """tests route /api/v2/*guid*/sequence"""

    def runTest(self):

        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        # print(res)
        n_pre = len(json.loads(res.content.decode("utf-8")))  # get all the guids

        guid_to_insert = "guid_{0}".format(n_pre + 1)

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq)
        seq = "N" * 4411532
        # print("Adding TB reference sequence of {0} bytes with {1} Ns".format(len(seq), seq.count('N')))
        self.assertEqual(len(seq), 4411532)  # check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload={"guid": guid_to_insert, "seq": seq})
        self.assertEqual(res.status_code, 200)

        relpath = "/api/v2/{0}/sequence".format(guid_to_insert)
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)

        info = json.loads(res.content.decode("utf-8"))
        # print(info)
        self.assertEqual(info["guid"], guid_to_insert)
        self.assertEqual(info["invalid"], 1)


class test_sequence_4(unittest.TestCase):
    """tests route /api/v2/*guid*/sequence"""

    def runTest(self):

        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        # print(res)
        n_pre = len(json.loads(res.content.decode("utf-8")))  # get all the guids

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq2 = str(record.seq)

        guid_to_insert1 = "guid_{0}".format(n_pre + 1)
        guid_to_insert2 = "guid_{0}".format(n_pre + 2)

        seq1 = "N" * 4411532
        # print("Adding TB reference sequence of {0} bytes with {1} Ns".format(len(seq1), seq1.count('N')))
        self.assertEqual(len(seq1), 4411532)  # check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload={"guid": guid_to_insert1, "seq": seq1})
        self.assertEqual(res.status_code, 200)

        # print("Adding TB reference sequence of {0} bytes with {1} Ns".format(len(seq2), seq2.count('N')))
        self.assertEqual(len(seq2), 4411532)  # check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload={"guid": guid_to_insert2, "seq": seq2})
        self.assertEqual(res.status_code, 200)

        relpath = "/api/v2/{0}/sequence".format(guid_to_insert1)
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)

        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(info["guid"], guid_to_insert1)
        self.assertEqual(info["invalid"], 1)
        relpath = "/api/v2/{0}/sequence".format(guid_to_insert2)
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)

        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(info["guid"], guid_to_insert2)
        self.assertEqual(info["invalid"], 0)


class test_sequence_5(unittest.TestCase):
    """tests route /api/v2/*guid*/sequence"""

    def runTest(self):

        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "/api/v2/guids"
        res = do_GET(relpath)
        # print(res)
        n_pre = len(json.loads(res.content.decode("utf-8")))  # get all the guids

        inputfile = "COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq2 = str(record.seq)

        guid_to_insert1 = "guid_{0}".format(n_pre + 1)
        guid_to_insert2 = "guid_{0}".format(n_pre + 2)

        seq1 = "R" * 4411532
        # print("Adding TB reference sequence of {0} bytes with {1} Rs".format(len(seq1), seq1.count('R')))
        self.assertEqual(len(seq1), 4411532)  # check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload={"guid": guid_to_insert1, "seq": seq1})
        self.assertEqual(res.status_code, 200)

        # print("Adding TB reference sequence of {0} bytes with {1} Ns".format(len(seq2), seq2.count('N')))
        self.assertEqual(len(seq2), 4411532)  # check it's the right sequence

        relpath = "/api/v2/insert"
        res = do_POST(relpath, payload={"guid": guid_to_insert2, "seq": seq2})
        self.assertEqual(res.status_code, 200)

        relpath = "/api/v2/{0}/sequence".format(guid_to_insert1)
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)

        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(info["guid"], guid_to_insert1)
        self.assertEqual(info["invalid"], 1)
        relpath = "/api/v2/{0}/sequence".format(guid_to_insert2)
        res = do_GET(relpath)
        self.assertEqual(res.status_code, 200)

        info = json.loads(res.content.decode("utf-8"))
        self.assertEqual(info["guid"], guid_to_insert2)
        self.assertEqual(info["invalid"], 0)


class test_nucleotides_excluded(unittest.TestCase):
    """tests route /api/v2/nucleotides_excluded"""

    def runTest(self):

        relpath = "/api/v2/reset"
        res = do_POST(relpath, payload={})

        relpath = "api/v2/nucleotides_excluded"
        res = do_GET(relpath)
        resDict = json.loads(res.text)
        self.assertTrue(isinstance(resDict, dict))
        self.assertEqual(set(resDict.keys()), set(["exclusion_id", "excluded_nt"]))
        self.assertEqual(res.status_code, 200)
