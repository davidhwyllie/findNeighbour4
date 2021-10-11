#!/usr/bin/env python3
""" compares bacterial sequences using catwalk (an external distance provider)

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

"""

import datetime
import hashlib
import json
import copy
import random
import numpy as np
from scipy.stats import binom_test
import pandas as pd
import logging
import progressbar
import time
from collections import Counter
from findn.persistence import Persistence
from findn.mongoStore import fn3persistence
from findn.rdbmsstore import fn3persistence_r
from findn.msa import MSAResult
from catwalk.pycw_client import CatWalk


# connection to mongodb on localhost; used for unittesting
UNITTEST_MONGOCONN = "mongodb://localhost"


class NeighbourStorageFailureError(Exception):
    """what we stored isn't what we got back"""

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


class NoCWParametersProvidedError(Exception):
    """no catwalk parameters provided"""

    def __init__(self):
        pass


class cw_seqComparer:
    def __init__(
        self,
        reference,
        maxNs,
        snpCeiling,
        excludePositions=set(),
        preComparer_parameters={},
        PERSIST=UNITTEST_MONGOCONN,
        unittesting=False,
        disable_insertion=False,
    ):

        """a module for determining relatedness between two sequences.

         * it stores data in catwalk, a separate relatedness engine
         * it does not persist the entire data set with RAM, but accesses this using an PERSISTENCE (database access) object

        It exposes a similar interface to py_seqComparer.

        Note that unlike py_seqComparer and earlier xComparer versions, this class is stateless (state is held in catwalk and in a database) and data is not held in RAM.
        As s

         Parameters:
         ===========
         reference is a string consisting of the reference sequence.  This is required because, as a data compression technique,  only differences from the reference are stored.

         maxNs: If the number of Ns are more than maxNs, no data from the sequence is stored.

         snpCeiling:  Results > snpCeiling are not returned or stored.

         excludePositions contains a zero indexed set of bases which should not be considered at all in the sequence comparisons.  Any bases which are always N should be added to this set.  Not doing so will substantially degrade the algorithm's performance.

         preComparer_parameters:
            Must include the following:
            selection_cutoff: (int)
            snp distances more than this are not of interest epidemiologically and are not reported

            uncertain_base (str): one of 'M', 'N', or 'N_or_M'
            which bases to regard as uncertain in MixPore (assessment of mixed bases) computations.  Default to M.

            over_selection_cutoff_ignore_factor:
            [currently ignored]

            catWalk_parameters:  parameters for catWalk function.  If empty {}, an error is raised.
                To run catWalk, the following keys are needed:
                    cw_binary_filepath : path to the catwalk binary.  If present, the CW_BINARY_FILEPATH environment variable will be used instead.  To set this, edit or create an .env file next to the Pipfile (if using a virtual environment)
                    reference_filepath : path to the relevant fasta reference sequence
                    reference_name: a human readable version of the reference's name
                    mask_filepath : path to a mask file, which consists of zero indexed positions to exclude
                    bind_host: the host the catwalk is running on
                    bind_port: the port the catwalk is running on
                    An example would look like:
                                    catWalk_parameters ={'cw_binary_filepath':None,
                                    'reference_name':"h37rv",
                                    'reference_filepath':"reference/TB-ref.fasta",
                                    'mask_filepath':"reference/TB-exclude-adaptive.txt",
                                    'bind_host':"127.0.0.1",
                                    'bind_port':5999}

                Catwalk also requires additional parameters,  These are supplied:
                    max_distance : maximum distance to report.  Catwalk does not store distances, so high distances can be requested, albeit at the cost of slightly slower computation.  max_distance is set to selection_cutoff.

         PERSIST: either a mongo/rdbms connection string, or an instance of a Persistence object

         unittesting: if True, will remove any stored data from the associated database

         disable_insertion:  catwalk is not used; can only access data already present in the database.

         - to run unit tests, do
         python3 -m unittest cw_seqComparer
        """

        # reference based compression relative to reference 'compressed_sequence' have keys as follows:
        self.distances_are_exact = (
            True  # this class does not deliver estimated distances
        )
        self.compressed_sequence_keys = set(["invalid", "A", "C", "G", "T", "N", "M"])

        # snp distances more than this will not be reported
        self.snpCeiling = snpCeiling

        # sequences with more than maxNs Ns will be considered invalid and their details (apart from their invalidity) will not be stored.
        self.maxNs = maxNs

        # check composition of the reference.
        self.reference = str(reference)  # if passed a Bio.Seq object, coerce to string.
        self.reference_list = list(self.reference)  # a list, one element per nt
        letters = Counter(self.reference)
        if len(set(letters.keys()) - set(["A", "C", "G", "T"])) > 0:
            raise TypeError(
                "Reference sequence supplied contains characters other than ACTG: {0}".format(
                    letters
                )
            )

        # load the excluded bases
        self.excluded = excludePositions

        # define what is included
        self.included = set(range(len(self.reference))) - self.excluded

        # open a persist object
        if isinstance(PERSIST, str):
            pm = Persistence()
            self.PERSIST = pm.get_storage_object(PERSIST, verbose=True)

        elif isinstance(PERSIST, fn3persistence) | isinstance(
            PERSIST, fn3persistence_r
        ):
            self.PERSIST = PERSIST  # passed a persistence object
        else:
            raise TypeError(
                "Expected a sqlalchemy connection string or an Persistence object but got {0}; can't proceed".format(
                    type(PERSIST)
                )
            )

        # attempt to load preComparer_parameters (which are used to parameterise catwalk) from the mongo / rdbms
        stored_preComparer_parameters = self.PERSIST.config_read("preComparer")
        if stored_preComparer_parameters is not None:
            try:
                del stored_preComparer_parameters["_id"]
            except KeyError:
                pass  # not present in rdbms
            preComparer_parameters = stored_preComparer_parameters

        else:
            # store them, if any are supplied
            if len(preComparer_parameters.keys()) > 0:
                self.PERSIST.config_store("preComparer", preComparer_parameters)

        # store whether any insertion is necessary
        self.disable_insertion = disable_insertion
        if self.disable_insertion:
            logging.info("cw_seqComparer is running in read-only mode")

        # start catwalk
        if "catWalk_parameters" not in preComparer_parameters.keys():
            raise NoCWParametersProvidedError()

        self.catWalk_parameters = preComparer_parameters["catWalk_parameters"]
        self.uncertain_base = preComparer_parameters["uncertain_base"]
        self.catWalk_parameters["max_distance"] = preComparer_parameters[
            "selection_cutoff"
        ]
        if len(self.catWalk_parameters) == 0:
            raise NoCWParametersProvidedError()

        self.catWalk = CatWalk(**self.catWalk_parameters, unittesting=unittesting)

        self.repopulate_all()  # reload samples
        logging.info(
            "CatWalk server started, operating with uncertain bases representing {0}".format(
                self.uncertain_base
            )
        )

        # estimates for msa
        self.estimated_p1 = None
        self.estimated_p1_sample_size = 0

        if unittesting:
            if (
                not disable_insertion
            ):  # if insertion (and deletion) is disabled, we can't restart
                self.restart_empty()

    def guids(self):
        """returns samples in the catwalk object"""
        return self.catWalk.sample_names()

    def repopulate_all(self):
        """repopulates catwalk with all samples from database
        Note: will not add a sample twice"""

        # eheck what needs to be done
        guids = self.PERSIST.guids_valid()
        already_added = self.catWalk.sample_names()
        to_add = set(guids) - set(already_added)

        bar = progressbar.ProgressBar(max_value=len(to_add))
        logging.info(
            "Repopulating catwalk.  Already present: {0}; remaining {1}".format(
                len(already_added), len(to_add)
            )
        )
        for i, guid in enumerate(to_add):
            bar.update(i)
            self.repopulate(guid)

        bar.finish()
        logging.info(
            "Repopulation of catwalk complete; added in total {0}".format(len(to_add))
        )

    def restart_empty(self):
        """restarts the catwalk server, with no data in it;
        also deletes all data from the database store."""

        self.PERSIST._delete_existing_data()
        self.catWalk.stop()  # terminate the catwalk instance
        time.sleep(1)
        self.catWalk.start()  # restart empty
        time.sleep(1)
        self.repopulate_all()

    def repopulate(self, guid):
        """saves a reference compressed object into catwalk,
        having loaded from the db using its reference guid.
        Note: will not add a sample twice.

        parameters:
        guid: the identity (name) of the sequence

        """
        # load
        obj = self.PERSIST.refcompressedsequence_read(guid)

        # store object in the precomparer
        self.catWalk.add_sample_from_refcomp(guid, obj)
        return

    def persist(self, obj, guid, annotations={}, unittesting_omit_link=False):
        """saves obj, a reference compressed object generated by .compress(),
        in the PERSISTENCE store.

        also saves object in the preComparer, and ensures the preComparer's settings
        are in sync with those on in the db. this enables external software, such as that calibrating preComparer's performance, to update preComparer settings while  the server is operating.

        note that this call should be made serially - not simulataneously from different threads.

        parameters:
        obj: a reference compressed object generated by .compress()
            {'A':set([1,2,3,4]} which would represent non-reference A at positions 1-4.
            The following keys are possible:
            A,C,G,T,N,M, invalid.
            A,C,G,T,N keys map to sets representing the positions where the reference differs by the respective keys
            M is a dictionary e.g. {1:'y',2:'R'} where the keys of the dictionary represent the positions where the base is mixed and the value is the IUPAC code for the mixture.
            invalid is 0 of the sequence is valid (i.e. comparisons should be made) and 1 if it is not (no comparison are made with the sequence).
            This format is that used and stored by findNeighbour4.

        guid: the identity (name) of the sequence
        annotations: a dictionary containing optional annotations, in the format
            {namespace:{key:value,...}
            e.g.
            {'Sequencing':{'Instrument':'MiSeq','Library':'Nextera'}}
            The persistence will store a validity key 'invalid' in the annotations DNAQuality namespace, overwriting if it exists.

        unittesting_omit_link: do not add links.  Do not use this; it is only used for generating conditions where the links are missing for unittesting the verify_insertion method

        """
        # check that insertion is allowed
        if self.disable_insertion:
            raise NotImplementedError(
                "Persistence is not permitted via this cw_seqComparer"
            )

        if (
            "invalid" not in obj.keys()
        ):  # assign it a valid status if we're not told it's invalid
            obj["invalid"] = 0

        # check it is not invalid
        isinvalid = obj["invalid"] == 1

        # construct an object to be stored in catWalk
        # determine the object's composition, i.e. numbers of bases differing from reference.
        # this is not required for SNV computation, but it is required for composition-based mixture testing (mixPORE etc).
        # computing this here is very lightweight computationally
        obj_composition = {"A": 0, "C": 0, "G": 0, "T": 0, "M": 0, "N": 0, "invalid": 0}
        if isinvalid:
            obj_composition["invalid"] = 1

        # construct an object to store
        for key in set(["A", "C", "G", "T"]) - set(obj.keys()):  # what is missing
            obj[key] = set()  # add empty set if no key exists.

        # create a smaller object to store in Catwalk; M and N are combined as Us
        smaller_obj = {}
        for item in ["A", "C", "G", "T"]:
            smaller_obj[item] = copy.deepcopy(obj[item])

        # store uncertain bases - either Us, Ns, or both
        # if there are no M or N keys, add them
        if "N" not in obj.keys():
            obj["N"] = set()
        if "M" not in obj.keys():
            obj["Ms"] = set()
        else:
            obj["Ms"] = set(obj["M"].keys())  # make a set of the positions of us

        if self.uncertain_base == "M":
            obj["U"] = obj["Ms"]
        elif self.uncertain_base == "N":
            obj["U"] = obj["N"]
        elif self.uncertain_base == "N_or_M":
            obj["U"] = obj["N"].union(obj["Ms"])
        else:
            raise KeyError(
                "Invalid uncertain_base: got {0}".format(self.uncertain_base)
            )

        del obj["Ms"]  # no longer needed

        # add the uncertain bases to the smaller object for storage in catwalk
        smaller_obj["U"] = copy.deepcopy(obj["U"])
        key_mapping = {
            "A": "A",
            "C": "C",
            "G": "G",
            "T": "T",
            "U": "N",
        }  # catwalk uses N, not U, for unknown. map specifies this.

        # add to database
        try:
            self.PERSIST.refcompressedseq_store(guid, obj)
        except FileExistsError:
            return ["Sample {0} already exists".format(guid)]  # already present

        # add annotations
        if "DNAQuality" in annotations.keys():

            annotations["DNAQuality"]["invalid"] = int(isinvalid)
        else:
            annotations["DNAQuality"] = {"invalid": int(isinvalid)}

        # add all annotations
        if annotations is not None:
            for nameSpace in annotations.keys():

                self.PERSIST.guid_annotate(
                    guid=guid, nameSpace=nameSpace, annotDict=annotations[nameSpace]
                )

        if not isinvalid:  # we only store valid sequences in catWalk
            to_catwalk = {}
            for key in smaller_obj.keys():
                to_catwalk[key_mapping[key]] = list(
                    smaller_obj[key]
                )  # make a dictionary for catwalk

            retVal = self.catWalk.add_sample_from_refcomp(guid, to_catwalk)  # add it
        else:
            return 1

        if (
            retVal == 200
        ):  # the sample was already there; no need to add neighbours; 201 is returned if insert is successful
            return obj["invalid"]

        # add links if the sample is valid and newly inserted
        links = {}
        loginfo = ["inserted {0} into catwalk; isinvalid = {1}".format(guid, isinvalid)]

        mcompare_result = self.mcompare(guid)  # compare guid against all

        for i, item in enumerate(mcompare_result["neighbours"]):  # all against all
            (guid1, guid2, dist) = item
            if not guid1 == guid2:
                link = {"dist": dist}
            if dist is not None:
                if link["dist"] <= self.snpCeiling:
                    links[guid2] = link

        ## now persist links in database.
        # insert any links found, and the datetime we added the statement.

        # insert
        if unittesting_omit_link is False:
            self.PERSIST.guid2neighbour_add_links(guid=guid, targetguids=links)

        return loginfo

    def load(self, guid):
        """returns a reference compressed object into a variable
        Note: this function loads a reference compressed json, stored relative to the reference.
        """
        return self.PERSIST.refcompressedsequence_read(guid)

    def verify_insertion(self, guid):
        """checks that the sample identified by guid has been correctly inserted

        this is called if a web worker doing insertion crashes,
        leaving the database in a potentially inconsistent state.

        such a crash is readily detectable because the insertion lock is left locked, so the sample identifier involved in the
        failure can be discerned.

        parameters:
        guid:       the sample identifier

        returns:
        1        if the record is in place
        0       if the record is not in place.

        Logic
        =====
        1. Check whether guid is present in (a) Catwalk (b) the database

        Actions taken:
        STEP 1
        Present in catwalk & database --> proceed to next stage
        Present in catwalk but not database --> delete from catwalk
        Present in database but not catwalk --> add to catwalk
        Present in neither --> return 0

        STEP 2:
        recover all edges from the database associated with guid
        compare with edges in catwalk
        add any missing edges

        returns a dictionary
        {'status':x, 'n_updated': y}

        x reflects the situation found
        y reflects the number of links updated

        x =

        'Neither' Present in neither
        'Catwalk only' Present in catwalk but not database
        'Database not catwalk' Present in database but not catwalk
        'Both' Present in catwalk & database

        """

        # STEP 1
        exists_in_catwalk = guid in self.catWalk.sample_names()
        exists_in_database = self.PERSIST.guid_exists(guid)

        if exists_in_catwalk is False and exists_in_database is False:
            logging.info(
                "Verify insertion of {0}.  Present in neither catwalk nor database".format(
                    guid
                )
            )
            return {"status": "Neither", "n_updated": 0}
        elif exists_in_catwalk is True and exists_in_database is False:
            # delete from catwalk
            self.catWalk.remove_sample(guid)
            logging.info(
                "Verify insertion of {0}.  Present in catwalk, not present in database.  Deleted from catwalk".format(
                    guid
                )
            )
            return {"status": "Catwalk only", "n_updated": 0}
        elif exists_in_catwalk is False and exists_in_database is True:
            # load from database, add to catwalk
            logging.info(
                "Verify insertion of {0}.  Not present in catwalk,  present in database.  Repopulating catwalk".format(
                    guid
                )
            )

            self.repopulate(guid)
            status = "Database, not catwalk"
        else:
            # otherwise exists_in_catwalk and exists_in_database: we need take no action except checking links
            logging.info(
                "Verify insertion of {0}.  Present in catwalk,  present in database. ".format(
                    guid
                )
            )
            status = "Both"

        # step 2: get all existing links from database
        registered_links = self.PERSIST.guid2neighbours(
            guid=guid, cutoff=1e12, returned_format=3
        )
        registered_links = set(registered_links["neighbours"])

        # get all links from catwalk
        links = {}
        mcompare_result = self.mcompare(guid)  # compare guid against all

        for i, item in enumerate(mcompare_result["neighbours"]):  # all against all
            (guid1, guid2, dist) = item
            if not guid1 == guid2:
                link = {"dist": dist}
            if dist is not None:

                if guid2 not in registered_links:
                    links[guid2] = link

        # if there are extra links to add
        if len(links.keys()) > 0:
            logging.info(
                "Missing links found.  Added {0} links".format(len(links.keys()))
            )
            self.PERSIST.guid2neighbour_add_links(guid=guid, targetguids=links)
        else:
            logging.info("Links verified.")
        return {"status": status, "n_updated": len(links.keys())}

    def mcompare(self, guid):
        """performs comparison of one guid with
        all valid guids, which are also stored samples.

        input:
            guid: the guid to compare
        output:
            a dictionary, including as keys
            'invalid': 0 if valid, 1 if invalid.
            'neighbours' : a list of neighbours and their distances.  If invalid, an empty list
            'timings': a dictionary including numbers analysed and timings for comparisons.
        """

        t1 = datetime.datetime.now()

        # check validity
        #  scores are
        # validity -1   The guid does not exist
        #   0    The guid exists and the sequence is valid
        #   1    The guid exists and the sequence is invalid

        validity_score = self.PERSIST.guid_valid(guid)

        if validity_score == -1:
            raise KeyError("mcompare: Sample {0} does not exist".format(guid))

        elif validity_score == 0:
            # sample is mcompare using catwalk
            cw_neighbours = self.catWalk.neighbours(guid, self.snpCeiling)
            neighbours = []
            for match, dist in cw_neighbours:
                neighbours.append([guid, match, dist])
            invalid = 0

        elif validity_score == 1:
            neighbours = []
            invalid = 1
        else:
            raise ValueError(
                "PERSIST.guid_valid({0}) returned an invalid score: 0,1,-1 expected got {1}".format(
                    guid, validity_score
                )
            )

        t2 = datetime.datetime.now()
        i1 = t2 - t1
        n1 = len(neighbours)

        if n1 > 0:
            rate1 = 1e3 * i1.total_seconds() / n1
        else:
            rate1 = 0

        timings = {
            "py_preComparer_msec_per_comparison": rate1,
            "candidates": n1,
            "matches": len(neighbours),
        }

        return {"invalid": invalid, "neighbours": neighbours, "timings": timings}

    def summarise_stored_items(self):
        """counts how many sequences exist in catwalk."""
        retVal = {}

        # call the catWalk server, and return the dictionary including status information in key-value format.
        # the keys must be pipe delimited and should be for the format
        # 'server|catWalk|{key}'.  The values must be scalars, not lists or sets.
        # the dictionary thus manipulated should be added to retVal.
        cw_status = self.catWalk.info()
        for item in cw_status:
            key = "server|catwalk|{0}".format(item)
            retVal[key] = cw_status[item]

        return retVal

    def excluded_hash(self):
        """returns a string containing the number of nt excluded, and a hash of their positions.
        This is useful for version tracking & storing patterns of masking."""
        list_excluded = sorted(list(self.excluded))
        len_l = len(list_excluded)
        h = hashlib.md5()
        h.update(json.dumps(list_excluded).encode("utf-8"))
        md5_l = h.hexdigest()
        return "Excl {0} nt [{1}]".format(len_l, md5_l)

    def uncompress(self, compressed_sequence):
        """returns a sequence from a compressed_sequence"""
        if "invalid" in compressed_sequence.keys():
            if compressed_sequence["invalid"] == 1:
                raise ValueError(
                    "Cannot uncompress an invalid sequence, because the sequence it is not stored {0}".format(
                        compressed_sequence.keys()
                    )
                )

        seq = list(self.reference)

        # mark all positions excluded as N
        for x in self.excluded:
            seq[x] = "N"

        # mark any bases definitively called as whatever they are called as
        for item in ["A", "C", "T", "G", "N"]:
            for x in compressed_sequence[item]:  # these are different from reference;
                seq[x] = item

        # add any mixed bases
        for x in compressed_sequence[
            "M"
        ].keys():  # the 'M' option includes iupac coding
            seq[x] = compressed_sequence["M"][x]
        return "".join(seq)

    def uncompress_guid(self, guid):
        """uncompresses a saved sequence"""
        seq = self.load(guid)
        return self.uncompress(seq)

    def compress(self, sequence):
        """reads a string sequence and extracts position - genome information from it.
        returns a dictionary consisting of zero-indexed positions of non-reference bases.

        """
        if not len(sequence) == len(self.reference):
            raise TypeError(
                "sequence must of the same length as reference; seq is {0} and ref is {1}".format(
                    len(sequence), len(self.reference)
                )
            )
        if len(self.reference) == 0:
            raise TypeError("reference cannot be of zero length")

        # we consider - characters to be the same as N
        sequence = sequence.replace("-", "N")

        # we only record differences relative to to refSeq.
        # anything the same as the refSeq is not recorded.
        # a dictionary, M, records the mixed base calls.
        diffDict = {
            "A": set([]),
            "C": set([]),
            "T": set([]),
            "G": set([]),
            "N": set([]),
            "M": {},
        }

        for i in self.included:  # for the bases we need to compress
            if not sequence[i] == self.reference[i]:  # if it's not reference
                if sequence[i] in ["A", "C", "T", "G", "N"]:
                    diffDict[sequence[i]].add(i)  # if it's a definitively called base
                else:
                    # we regard it as a code representing a mixed base.  we store the results in a dictionary
                    diffDict["M"][i] = sequence[i]

        # check how many Ns
        if len(diffDict["N"]) + len(diffDict["M"].keys()) > self.maxNs:
            # we store it, but not with sequence details if is invalid
            diffDict = {"invalid": 1}
        else:
            diffDict["invalid"] = 0
        return diffDict

    def compressed_sequence_hash(self, compressed_sequence):
        """returns a string containing a hash of a compressed object.
        Used for identifying compressed objects, including consensus sequences.
        """
        keys = sorted(compressed_sequence.keys())
        serialised_compressed_sequence = ""
        for key in keys:
            if isinstance(compressed_sequence[key], set):
                list_compressed = sorted(list(compressed_sequence[key]))
            else:
                list_compressed = compressed_sequence[key]
            serialised_compressed_sequence = (
                serialised_compressed_sequence + key + ":" + str(list_compressed) + ";"
            )
        h = hashlib.md5()
        h.update(serialised_compressed_sequence.encode("utf-8"))
        md5 = h.hexdigest()
        return md5

    def estimate_expected_proportion(self, seqs):
        """computes the median Ns for seqs, a list of referencecompressed sequence objects.
        Returns None if
        * the length of the sequences in seqs are 0, or
        * there are <= 3 seqs
        """
        if len(seqs) < 3:
            return None
        if len(seqs[0]) == 0:
            return None
        Ns = []
        for seq in seqs:
            Ns.append(seq.count("N"))
        return np.median(Ns) / len(seqs[0])

    def estimate_expected_unk(
        self,
        sample_size=30,
        exclude_guids=set(),
        unk_type="N",
        what="median",
        sites=None,
    ):
        """computes the median numbers of unk_type (N,M,or N_or_N) for sample_size guids, randomly selected from all guids except for exclude_guids.

        Parameters
        sample_size: how many items to sample when determining expected_unk.  If None, uses all samples
        exclude_guids: an iterable of guids not to analyse
        unk_type: one of 'N' 'M' 'N_or_M'

        If sites are None, whole sequence is analysed.  If sites is a list/set, only these positions are analysed

        Return
        either the median unk_type for the sample, or None if fewer than sample_size items are present.

        """

        if unk_type not in ["N", "M", "N_or_M"]:
            raise KeyError(
                "unk_type can be one of 'N' 'M' 'N_or_M', not {0}".format(unk_type)
            )

        # step 0: find all valid guids
        valid_guids = set(self.PERSIST.guids_valid())
        guids = valid_guids - set(exclude_guids)

        guids = list(guids)
        if sample_size is not None:
            if sample_size > len(guids):
                return None

        guids = random.sample(guids, sample_size)

        if isinstance(sites, set):
            sites = list(sites)

        composition_dict = {}
        for guid in guids:
            rcs = self.load(guid)

            # should not happen
            if rcs is None:
                raise KeyError("Tries to load {0} but got None ".format(guid))

            seq = self.uncompress(rcs)
            if sites is not None:
                if len(sites) > 0:
                    start_list = list(seq)
                    seq = list(start_list[i] for i in sites)
                    seq = "".join(seq)

            composition_onesample = Counter(seq)
            for unk_nt in ["N", "M"]:
                if unk_nt not in composition_onesample.keys():
                    composition_onesample[unk_nt] = 0
            composition_onesample["N_or_M"] = (
                composition_onesample["N"] + composition_onesample["M"]
            )
            composition_dict[guid] = composition_onesample

        composition = pd.DataFrame.from_dict(composition_dict, orient="index")
        if len(composition.index) == 0:
            return None

        report_value = True
        if sample_size is not None:
            if len(composition.index) < sample_size:
                report_value = False

        if report_value:
            if what == "median":
                return np.median(composition[unk_type])
            elif what == "mean":
                return np.mean(composition[unk_type])
            elif what == "sum":
                return np.sum(composition[unk_type])
            elif what == "sum_and_median":
                return np.sum(composition[unk_type]), np.median(composition[unk_type])
            else:
                raise ValueError("Invalid what passed : {0}".format(what))
        else:
            return None

    def estimate_expected_unk_sites(
        self, unk_type, sites, what="median", sample_size=30
    ):
        """estimates the median unk_type for all guids,  at the positions in sites().

        Parameters
        unk_type: one of 'N' 'M' 'N_or_M'
        sites: the sites to consider, zero indexed
        sample_size: the minimum number of samples from which to report an estimate

        Returns
        either the estimated median unk_type for the sample, or None if fewer than sample_size items are present.

        """
        if len(sites) > 0:
            return self.estimate_expected_unk(
                unk_type=unk_type, what=what, sites=sites, sample_size=sample_size
            )
        else:
            return self.estimate_expected_unk(
                unk_type=unk_type, what=what, sites=None, sample_size=sample_size
            )

    def update_p1_estimate(self, sample_size, uncertain_base_type):
        """estimates p1 (see MSA for description) from existing data using a sample of size sample_size,
        if it has not already been recorded
        returns: Nothing
        side effects: sets estimated_p1 and estimated_p1_sample_size"""

        # Estimate expected N or M as median(observed N or Ms),
        # which is a valid thing to do if the proportion of mixed samples is low.

        if (
            self.estimated_p1 is not None
            and self.estimated_p1_sample_size >= sample_size
        ):
            # already computed, no need to repeat.
            return

        expected_N1 = self.estimate_expected_unk(
            sample_size=sample_size,
            unk_type=uncertain_base_type,
        )
        if expected_N1 is None:
            self.estimated_p1 = None
        else:
            self.estimated_p1 = expected_N1 / (len(self.reference) - len(self.excluded))

        self.estimated_p1_sample_size = sample_size

    def raise_error(self, token):
        """raises a ZeroDivisionError, with token as the message.
        useful for unit tests of error logging"""
        raise ZeroDivisionError(token)

    def multi_sequence_alignment(
        self,
        guids,
        sample_size=30,
        expected_p1=None,
        uncertain_base_type="N",
        outgroup=None,
    ):
        """computes a multiple sequence alignment containing only sites which vary between guids.

            Parameters:
            guids: an iterable of guids to analyse

            sample_size is the number of samples to randomly sample to estimate the expected number of N or Ms in
            the population of sequences currently in the server.  From this, the routine computes expected_p1,
            which is expected_expected_N/M the length of sequence.
            if expected_p1 is supplied, then such sampling does not occur.

            uncertain_base_type: the kind of base which is to be analysed, either N,M, or N_or_M
            outgroup: the outgroup sample, if any.  Not used in computations, but stored in the output. ** NOT IMPLEMENTED AT PRESENT **

            Output:
            A MSAResult object.

            Statistical approach:
            A key use of the function is to compute statistics comparing the frequency of mixed (N,N,N_or_M) bases between variant bases in the alignment, and whose which are not.
            This includes statistical testing of mixed base frequencies.

            The p values reported are derived from exact, one-sided binomial tests as implemented in python's scipy.stats.binom_test().

            TEST 1:
            This tests the hypothesis that the number of Ns in the *alignment*
            is GREATER than those expected from the expected_N in the population of whole sequences.

            Does so by comparing the observed number of Ns in the alignment (alignN),
            given the alignment length (4 in the above case) and an expectation of the proportion of bases which will be N.
            The expected number of Ns is estimated by
            i) randomly sampling sample_size guids from those stored in the server and
            observing the number of Ns per base across the genome.  The estimate_expected_unk() function performs this.
            ii) randomly sampling sample_size guids from those stored in the server and
            observing the number of Ns per base across the relevant  genome.  The estimate_expected_unk() function performs this.

            This approach determines the median number of Ns in valid sequences, which (if bad samples with large Ns are rare)
            is a relatively unbiased estimate of the median number of Ns in the good quality samples.

            If there  are not enough samples in the server to obtain an estimate, p_value is not computed, being
            reported as None.

            TEST 2:
            This tests the hypothesis that the number of Ns in the *alignment*
            is GREATER than those expected from the expected_N in the population of whole sequences
            *at the bases examined in the alignment*.
            This might be relevant if these particular bases are generally hard to call.

            Does so by comparing the observed number of Ns in the alignment (alignN),
            given the alignment length (4 in the above case) and an expectation of the proportion of bases which will be N.
            The expected number of Ns is estimated by randomly sampling sample_size guids from those stored in the server and
            observing the number of Ns per base at the relevant sites.  The estimate_expected_unk_sites() function performs this.

            This approach determines the median number of Ns in valid sequences, which (if bad samples with large Ns are rare)
            is a relatively unbiased estimate of the median number of Ns in the good quality samples.

            If there  are not enough samples in the server to obtain an estimate, p_value is not computed, being
            reported as None.

            ## TEST 3:
        tests whether the proportion of Ns in the alignment is greater
            than in the bases not in the alignment, for this sequence.
        This is the test published.

            ## TEST 4:
            tests whether the proportion of Ns in the alignment  for this sequence
            is greater than the median proportion of Ns in the alignment for all other sequences.
            This is a sensible option if the test is performed per-cluster, and the proportion of Ns in the
            aligned sequences differs markedly by cluster.

            This test is computed if there are four or more samples in the cluster, and if the alignment length is non-zero.

        """

        # -1 validate input

        if sample_size is None:
            sample_size = 100

        if expected_p1 is not None:
            if expected_p1 < 0 or expected_p1 > 1:
                raise ValueError("Expected_p1 must lie between 0 and 1")
        else:
            # estimate expected_p1
            if (
                self.estimated_p1_sample_size > sample_size
            ):  # we have an estimate derived from a sufficient sample
                expected_p1 = self.estimated_p1
            else:
                # we need to construct an estimate
                self.update_p1_estimate(
                    sample_size=sample_size, uncertain_base_type=uncertain_base_type
                )

        valid_guids = set(guids).intersection(set(self.catWalk.sample_names()))
        return self._msa(
            valid_guids=valid_guids,
            invalid_guids=[],
            expected_p1=self.estimated_p1,
            sample_size=self.estimated_p1_sample_size,
            uncertain_base_type=uncertain_base_type,
            outgroup=outgroup,
        )

    def _msa(
        self,
        valid_guids,
        invalid_guids,
        expected_p1,
        sample_size,
        uncertain_base_type="N",
        outgroup=None,
    ):
        """perform multisequence alignment and significance tests.
        It assumes that the relevant sequences (valid_guids) are in-ram in this cw_seqComparer object as part of the seqProfile dictionary.

        Parameters:
        valid_guids:  the guids in valid_guids are those on which the msa is computed.
        invalid_guids: passed to the function for information only.

        uncertain_base_type:
        the software performs four significance tests
        on the number of uncertain_bases in an alignment containing only the variant bases between variant_guids.

        If uncertain_base_type == 'N' analyses Ns.
        If uncertain_base_type == 'M' analyses mixed bases.
        If uncertain_base_type == 'N_or_M' analyses both combined

        expected_p1:  The expected proportion of Ns or Ms (as specified by uncertain_base_type) is expected_p1.


        Returns:
        MSAResult object.

        Statistical note:
        The p values reported are derived from exact, one-sided binomial tests as implemented in pythons scipy.stats.binom_test().

        TEST 1:
        This tests the hypothesis that the number of Ns in the *alignment*
        is GREATER than those expected from the expected_N in the population of whole sequences.

        Does so by comparing the observed number of Ns in the alignment (alignN),
        given the alignment length (4 in the above case) and an expectation of the proportion of bases which will be N.
        The expected number of Ns is estimated by
        i) randomly sampling sample_size guids from those stored in the server and
        observing the number of Ns per base across the genome.  The estimate_expected_unk() function performs this.
        ii) randomly sampling sample_size guids from those stored in the server and
        observing the number of Ns per base across the relevant  genome.  The estimate_expected_unk() function performs this.

        This approach determines the median number of Ns in valid sequences, which (if bad samples with large Ns are rare)
        is a relatively unbiased estimate of the median number of Ns in the good quality samples.

        If there  are not enough samples in the server to obtain an estimate, p_value is not computed, being
        reported as None.

        TEST 2:
        This tests the hypothesis that the number of Ns in the *alignment*
        is GREATER than those expected from the expected_N in the population of whole sequences
        *at the bases examined in the alignment*.
        This might be relevant if these particular bases are generally hard to call.

        Does so by comparing the observed number of Ns in the alignment (alignN),
        given the alignment length (4 in the above case) and an expectation of the proportion of bases which will be N.
        The expected number of Ns is estimated by randomly sampling sample_size guids from those stored in the server and
        observing the number of Ns per base at the relevant sites.  The estimate_expected_unk_sites() function performs this.

        This approach determines the median number of Ns in valid sequences, which (if bad samples with large Ns are rare)
        is a relatively unbiased estimate of the median number of Ns in the good quality samples.

        If there  are not enough samples in the server to obtain an estimate, p_value is not computed, being
        reported as None.

        TEST 3: tests whether the proportion of Ns in the alignment is greater
        than in the bases not in the alignment, for this sequence.

        TEST 4: tests whether the proportion of Ns in the alignment  for this sequence
        is greater than the median proportion of Ns in the alignment for all other sequences.
        This is a sensible option if the test is performed per-cluster, and the proportion of Ns in the
        aligned sequences differs markedly by cluster.

        This test is computed if there are two or more samples in the cluster.

        """

        # step 1: define all positions in these sequences for which there is a non-reference position
        # we ignore Ms and Ns in this analysis
        # we also compute the total number of Ms or Ns in each guid

        nrps = {}
        guid2all = {"N": {}, "M": {}, "N_or_M": {}}
        guid2align = {"N": {}, "M": {}, "N_or_M": {}}
        seq_profile = {}

        for guid in valid_guids:
            seq_profile[guid] = self.load(guid)

        for guid in valid_guids:
            for unk_type in ["N", "M"]:

                unks = len(seq_profile[guid][unk_type])
                guid2all[unk_type][guid] = unks

            guid2all["N_or_M"][guid] = guid2all["N"][guid] + guid2all["M"][guid]

            for base in ["A", "C", "T", "G"]:
                positions = seq_profile[guid][base]
                for position in positions:
                    if (
                        position not in nrps.keys()
                    ):  # if it's non-reference, and we've got no record of this position
                        nrps[
                            position
                        ] = set()  # then we generate a set of bases at this position
                    nrps[position].add(
                        base
                    )  # either way add the current non-reference base there

        # step 2: for the non-reference called positions,
        # record whether there's a reference base there.
        for guid in valid_guids:
            for position in nrps.keys():
                psn_accounted_for = 0
                for base in ["A", "C", "T", "G"]:
                    if position in seq_profile[guid][base]:
                        psn_accounted_for = 1
                if psn_accounted_for == 0:  # no record of non-ref seq here
                    if not (
                        position in seq_profile[guid]["N"]
                        or position in seq_profile[guid]["M"].keys()
                    ):  # not M or N so non-reference
                        # it is reference; this guid has no record of a variant base at this position, so it must be reference.
                        nrps[position].add(
                            self.reference[position]
                        )  # add reference psn

        # step 3: find those which have multiple bases at a position
        variant_positions = set()
        for position in nrps.keys():
            if len(nrps[position]) > 1:
                variant_positions.add(position)

        # step 3b: define constant base frequencies - i.e. frequencies outside the region of variation.  This is used by iqtree, see http://www.iqtree.org/doc/Command-Reference
        constant_positions = list(self.included - variant_positions)
        constant_nt = [self.reference_list[x] for x in constant_positions]

        fconst_dict = Counter(constant_nt)
        fconst = [0, 0, 0, 0]
        for i, base in enumerate(["A", "C", "G", "T"]):
            if base in fconst_dict.keys():
                fconst[i] = fconst_dict[base]

        # step 4: determine the sequences of all bases.
        ordered_variant_positions = sorted(list(variant_positions))
        guid2seq = {}
        guid2mseq = {}
        guid2msa_seq = {}
        guid2msa_mseq = {}
        for guid in valid_guids:
            guid2seq[guid] = []
            guid2mseq[guid] = []
            for position in ordered_variant_positions:  # positions of variation
                this_base = self.reference[position]
                this_base_m = self.reference[position]
                for base in ["A", "C", "T", "G", "N", "M"]:
                    if not base == "M":
                        positions = seq_profile[guid][base]
                    else:
                        positions = seq_profile[guid][base].keys()

                    if position in positions:
                        this_base = base
                        this_base_m = base
                        if base == "M":
                            this_base_m = seq_profile[guid]["M"][position]
                guid2seq[guid].append(this_base)  # include M for mix if present
                guid2mseq[guid].append(this_base_m)  # include iupac code for mix
            guid2msa_seq[guid] = "".join(guid2seq[guid])
            guid2msa_mseq[guid] = "".join(guid2mseq[guid])

        # step 5: determine the expected_p2 at the ordered_variant_positions:
        expected_N2 = self.estimate_expected_unk_sites(
            sites=set(ordered_variant_positions), unk_type=uncertain_base_type
        )
        if expected_N2 is None:
            expected_p2 = None
        elif len(ordered_variant_positions) == 0:
            expected_p2 = None
        else:
            expected_p2 = expected_N2 / len(ordered_variant_positions)

        # step 6: perform Binomial tests on all samples
        if len(valid_guids) > 0:
            guid2pvalue1 = {}
            guid2pvalue2 = {}
            guid2pvalue3 = {}
            guid2pvalue4 = {}
            guid2observed_p = {}
            guid2expected_p1 = {}
            guid2expected_p2 = {}
            guid2expected_p3 = {}
            guid2expected_p4 = {}
            for guid in valid_guids:

                # compute p value 1.  This tests the hypothesis that the number of Ns in the *alignment*
                # is GREATER than those expected from the expected_N in the population of whole sequences.
                for unk_type in ["N", "M"]:
                    guid2align[unk_type][guid] = guid2msa_seq[guid].count(unk_type)
                guid2align["N_or_M"][guid] = (
                    guid2align["N"][guid] + guid2align["M"][guid]
                )

                if (
                    expected_p1 is None
                ):  # we don't have an expectation, so we can't assess the first binomial test;
                    p_value1 = None
                    observed_p = None
                elif (
                    len(guid2msa_seq[guid]) == 0
                ):  # we don't have any information to work with
                    p_value1 = None
                    observed_p = None
                else:
                    observed_p = guid2align[uncertain_base_type][guid] / len(
                        guid2msa_seq[guid]
                    )
                    p_value1 = binom_test(
                        guid2align[uncertain_base_type][guid],
                        len(guid2msa_seq[guid]),
                        expected_p1,
                        alternative="greater",
                    )

                guid2pvalue1[guid] = p_value1
                guid2observed_p[guid] = observed_p
                guid2expected_p1[guid] = expected_p1

                # compute p value 2.  This tests the hypothesis that the number of Ns or Ms in the *alignment*
                # is GREATER than those expected from the expected_N in the population of whole sequences
                # at these sites.
                if (
                    expected_p2 is None
                ):  # we don't have an expectation, so we can't assess the binomial test;
                    p_value2 = None
                elif (
                    len(guid2msa_seq[guid]) == 0
                ):  # we don't have any information to work with
                    p_value2 = None
                else:

                    p_value2 = binom_test(
                        guid2align[uncertain_base_type][guid],
                        len(guid2msa_seq[guid]),
                        expected_p2,
                        alternative="greater",
                    )
                guid2pvalue2[guid] = p_value2
                guid2expected_p2[guid] = expected_p2

                # compute p value 3.  This tests the hypothesis that the number of Ns or Ms in the alignment of THIS SEQUENCE
                # is GREATER than the number of Ns or Ms not in the alignment  IN THIS SEQUENCE
                # based on sequences not in the alignment

                expected_p3 = (
                    guid2all[uncertain_base_type][guid]
                    - guid2align[uncertain_base_type][guid]
                ) / (len(self.reference) - len(guid2msa_seq[guid]))
                p_value = binom_test(
                    guid2align[uncertain_base_type][guid],
                    len(guid2msa_seq[guid]),
                    expected_p3,
                    alternative="greater",
                )
                guid2pvalue3[guid] = p_value
                guid2expected_p3[guid] = expected_p3

                # compute p value 4.
                # tests whether the proportion of Ns in the alignment  for this sequence
                # is greater than the median proportion of Ns in the alignment for all other sequences.
                # This is a sensible option if the test is performed per-cluster, and the proportion of Ns in the
                # aligned sequences differs markedly by cluster.
                # This test is computed if there are two or more samples in the cluster.

                # get the sequences which are not guid;
                other_seqs = []
                for this_guid in guid2msa_seq.keys():
                    if not guid == this_guid:
                        other_seqs.append(guid2msa_seq[this_guid])
                expected_p4 = self.estimate_expected_proportion(other_seqs)
                if expected_p4 is None:
                    p_value = None
                else:
                    p_value = binom_test(
                        guid2align[uncertain_base_type][guid],
                        len(guid2msa_seq[guid]),
                        expected_p4,
                        alternative="greater",
                    )
                guid2pvalue4[guid] = p_value
                guid2expected_p4[guid] = expected_p4

            # assemble dataframe
            df1 = pd.DataFrame.from_dict(guid2msa_seq, orient="index")
            df1.columns = ["aligned_seq"]
            df1["aligned_seq_len"] = len(ordered_variant_positions)

            df1_iupac = pd.DataFrame.from_dict(guid2msa_mseq, orient="index")
            df1_iupac.columns = ["aligned_mseq"]

            df2n = pd.DataFrame.from_dict(guid2all["N"], orient="index")
            df2n.columns = ["allN"]
            df3n = pd.DataFrame.from_dict(guid2align["N"], orient="index")
            df3n.columns = ["alignN"]

            df2m = pd.DataFrame.from_dict(guid2all["M"], orient="index")
            df2m.columns = ["allM"]
            df3m = pd.DataFrame.from_dict(guid2align["M"], orient="index")
            df3m.columns = ["alignM"]

            df2mn = pd.DataFrame.from_dict(guid2all["N_or_M"], orient="index")
            df2mn.columns = ["allN_or_M"]
            df3mn = pd.DataFrame.from_dict(guid2align["N_or_M"], orient="index")
            df3mn.columns = ["alignN_or_M"]

            df4 = pd.DataFrame.from_dict(guid2pvalue1, orient="index")
            df4.columns = ["p_value1"]
            df5 = pd.DataFrame.from_dict(guid2pvalue2, orient="index")
            df5.columns = ["p_value2"]
            df6 = pd.DataFrame.from_dict(guid2pvalue3, orient="index")
            df6.columns = ["p_value3"]
            df7 = pd.DataFrame.from_dict(guid2pvalue4, orient="index")
            df7.columns = ["p_value4"]
            df8 = pd.DataFrame.from_dict(guid2observed_p, orient="index")
            df8.columns = ["observed_proportion"]
            df9 = pd.DataFrame.from_dict(guid2expected_p1, orient="index")
            df9.columns = ["expected_proportion1"]
            df10 = pd.DataFrame.from_dict(guid2expected_p3, orient="index")
            df10.columns = ["expected_proportion2"]
            df11 = pd.DataFrame.from_dict(guid2expected_p3, orient="index")
            df11.columns = ["expected_proportion3"]
            df12 = pd.DataFrame.from_dict(guid2expected_p4, orient="index")
            df12.columns = ["expected_proportion4"]
            df12["what_tested"] = uncertain_base_type

            df = df1.merge(df1_iupac, left_index=True, right_index=True)
            df = df.merge(df2n, left_index=True, right_index=True)
            df = df.merge(df3n, left_index=True, right_index=True)
            df = df.merge(df2m, left_index=True, right_index=True)
            df = df.merge(df3m, left_index=True, right_index=True)
            df = df.merge(df2mn, left_index=True, right_index=True)
            df = df.merge(df3mn, left_index=True, right_index=True)
            df = df.merge(df4, left_index=True, right_index=True)
            df = df.merge(df5, left_index=True, right_index=True)
            df = df.merge(df6, left_index=True, right_index=True)
            df = df.merge(df7, left_index=True, right_index=True)
            df = df.merge(df8, left_index=True, right_index=True)
            df = df.merge(df9, left_index=True, right_index=True)
            df = df.merge(df10, left_index=True, right_index=True)
            df = df.merge(df11, left_index=True, right_index=True)
            df = df.merge(df12, left_index=True, right_index=True)

            retDict = {
                "variant_positions": ordered_variant_positions,
                "invalid_guids": invalid_guids,
                "valid_guids": valid_guids,
                "expected_p1": expected_p1,
                "sample_size": sample_size,
                "df_dict": df.to_dict(orient="index"),
                "what_tested": uncertain_base_type,
                "outgroup": outgroup,
                "creation_time": datetime.datetime.now().isoformat(),
                "fconst": fconst_dict,
            }

            return MSAResult(**retDict)

        else:
            return None
