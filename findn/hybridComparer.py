#!/usr/bin/env python3
""" compares bacterial sequences, using either catwalk (an external distance provider), python code achieving the same functions

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
import datetime
import hashlib
import json
import copy
import random
import numpy as np
from scipy.stats import binom_test
import pandas as pd
import logging

from collections import Counter
from findn.mongoStore import fn3persistence
from findn.preComparer import preComparer  # catwalk enabled
from findn.msa import MSAResult

# connection to mongodb on localhost; used for unittesting
UNITTEST_MONGOCONN = "mongodb://localhost"


class NeighbourStorageFailureError(Exception):
    """ what we stored isn't what we got back"""

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


class hybridComparer:
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
         It is a hybrid in that
         * it stores a subset of data in RAM and uses this to pre-screen sequences.  This is delivered by a preComparer class.
         * it does not persist the entire data set with RAM, but accesses this using an fn3persistence object, connecting to the database found at mongo_connstring.

         Otherwise, it tries to seek the same interface as seqComparer.

         Parameters:
         ===========
         reference is a string consisting of the reference sequence.  This is required because, as a data compression technique,  only differences from the reference are stored.

         maxNs: If the number of Ns are more than maxNs, no data from the sequence is stored.

         snpCeiling:  Results > snpCeiling are not returned or stored.

         excludePositions contains a zero indexed set of bases which should not be considered at all in the sequence comparisons.  Any bases which are always N should be added to this set.  Not doing so will substantially degrade the algorithm's performance.

         preComparer_parameters: parameters passed to the preComparer.  Used to judge which sequences do not need further analysis.
         Must include the following:
         selection_cutoff: (int)
         snp distances more than this are not of interest epidemiologically and are not reported

         uncertain_base (str): one of 'M', 'N', or 'N_or_M'
         which bases to regard as uncertain in the computation.  Default to M.
         if 'N_or_M', will store all us and Ns and will give the same result as the seqComparer module.

         over_selection_cutoff_ignore_factor:
         SNP distances more than over_selection_cutoff_ignore_factor * selection_cutoff do not need to be further analysed.  For example, if a SNP cutoff was 20, and over_selection_cutoff_ignore_factor is 5, we can safely consider with SNV distances > 100 (=20*5) as being unrelated.

         catWalk_parameters:  parameters for catWalk function.  If empty {}, catWalk is not run.
         PERSIST: either a mongo connection string, or an instance of fn3persistence

         unittesting: if True, will remove any stored data from the associated database [Careful!]
        David Wyllie, September 2019

         - to run unit tests, do
         python3 -m unittest hybridComparer
        """

        # reference based compression relative to reference 'compressed_sequence' have keys as follows:

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
            raise TypeError("Reference sequence supplied contains characters other than ACTG: {0}".format(letters))

        # load the excluded bases
        self.excluded = excludePositions

        # define what is included
        self.included = set(range(len(self.reference))) - self.excluded

        # initialise pairwise sequences for comparison.
        self._refresh()

        if isinstance(PERSIST, str):
            self.PERSIST = fn3persistence(PERSIST)
        else:
            self.PERSIST = PERSIST  # passed a persistence object

        if unittesting:
            self.PERSIST._delete_existing_data()

        # attempt to load preComparer_parameters from the mongoStore
        stored_preComparer_parameters = self.PERSIST.config_read("preComparer")
        if stored_preComparer_parameters is not None:
            del stored_preComparer_parameters["_id"]
            preComparer_parameters = stored_preComparer_parameters
        else:
            # store them, if any are supplied
            if len(preComparer_parameters.keys()) > 0:
                self.PERSIST.config_store("preComparer", preComparer_parameters)

        # store whether any insertion is necessary
        self.disable_insertion = disable_insertion
        if self.disable_insertion:
            preComparer_parameters["catWalk_parameters"] = {}  # no catwalk
            logging.info("HybridComparer is running in read-only mode")

        self.pc = preComparer(**preComparer_parameters)

    def repopulate_sample(self, n=100):
        """saves a sample of reference compressed objects in the precomparer,
        for the purpose of assessing sequence composition etc.

        parameters:
        n the number of samples to load

        """

        guids = self.PERSIST.guids() - self.pc.guids()
        if n < len(guids):
            guids = random.sample(guids, n)
        for guid in guids:
            obj = self.PERSIST.refcompressedsequence_read(guid)

            # store object in the precomparer
            self.pc.persist(obj, guid)  # store in the preComparer
        return

    def repopulate(self, guid):
        """saves a reference compressed object in the precomparer,
        as loaded from the mongo db using its reference guid.

        parameters:
        guid: the identity (name) of the sequence


        """
        # load
        obj = self.PERSIST.refcompressedsequence_read(guid)

        # store object in the precomparer
        self.pc.persist(obj, guid)  # store in the preComparer
        return

    def update_precomparer_parameters(self):
        """ ensures precomparer settings are up to date.  Precomparer settings may be optimised by external processes, such as preComparer_calibrator. """

        stored_preComparer_parameters = self.PERSIST.config_read("preComparer")
        if stored_preComparer_parameters is not None:  # if settings exist on disc
            del stored_preComparer_parameters["_id"]  # remove the mongoDb id
            if len(stored_preComparer_parameters.keys()) > 0:
                if self.pc.check_operating_parameters(**stored_preComparer_parameters) is False:
                    self.pc.set_operating_parameters(**stored_preComparer_parameters)
        return ()

    def persist(self, obj, guid, annotations={}):
        """saves obj, a reference compressed object generated by .compress(),
        in the mongo db.

        also saves object in the preComparer, and ensures the preComparer's settings
        are in sync with those on in the mongoStore. this enables external software, such as that calibrating preComparer's performance, to update preComparer settings while  the server is operating.

        note that this call should be made serially - not synchronously from different threads.

        parameters:
        obj: a reference compressed object generated by .compress()
        guid: the identity (name) of the sequence
        annotations: a dictionary containing optional annotations, in the format
            {namespace:{key:value,...}
            e.g.
            {'Sequencing':{'Instrument':'MiSeq','Library':'Nextera'}}
            The persistence will store a validity key 'invalid' in the annotations DNAQuality namespace, overwriting if it exists.

        """
        # check that insertion is allowed
        if self.disable_insertion:
            raise NotImplementedError("Persistence is not permitted via this hybridComparer")

        # check preComparer settings are up to date
        # self.update_precomparer_parameters()           # no longer necessary

        # add sequence and its links to disc
        try:
            self.PERSIST.refcompressedseq_store(guid, obj)  # store in the mongostore
        except FileExistsError:
            return {"guid": guid, "result": "already present"}  # exists - we don't need to re-add it.

        pc_obj = copy.deepcopy(obj)
        # store object in the precomparer
        is_invalid = self.pc.persist(pc_obj, guid)  # store in the preComparer.
        # in the current implementation, this either loads CatWalk or an in memory python relatedness data structure doing the same thing.
        # If the object exists already, nothing will happen.  This does not write persistently to a db.

        # add links
        links = {}
        loginfo = ["inserted {0} into precomparer; is_invalid = {1}".format(guid, is_invalid)]

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
        self.PERSIST.guid2neighbour_add_links(guid=guid, targetguids=links)

        # make sure the invalid annotation is set correctly
        if "DNAQuality" in annotations.keys():
            annotations["DNAQuality"]["invalid"] = is_invalid
        else:
            annotations["DNAQuality"] = {"invalid": is_invalid}

        # add all annotations
        if annotations is not None:
            for nameSpace in annotations.keys():
                self.PERSIST.guid_annotate(guid=guid, nameSpace=nameSpace, annotDict=annotations[nameSpace])

        # if we found neighbours, report mcompare timings to log
        msg = "mcompare|"
        for key in mcompare_result["timings"].keys():
            msg = msg + " {0} {1}; ".format(key, mcompare_result["timings"][key])
        loginfo.append(msg)

        return loginfo

    def load(self, guid):
        """returns a reference compressed object into RAM.
        Note: this function loads stored on disc/db relative to the reference.
        """
        return self.PERSIST.refcompressedsequence_read(guid)

    def _refresh(self):
        """ empties any sequence data from ram in any precomparer object existing. """
        if hasattr(self, "pc"):
            self.pc.seqProfile = {}
        self.seqProfile = {}

    def mcompare(self, guid, guids=None):
        """performs comparison of one guid with
        all guids, which are also stored samples.

        input:
            guid: the guid to compare
            guids: guids to compare with; if blank, compares with all.
        output:
            a dictionary, including as keys
            'neighbours' : a list of neighbours and their distances
            'timings': a dictionary including numbers analysed and timings for comparisons.
        """

        t1 = datetime.datetime.now()

        # if guids are not specified, we do all vs all
        if guids is None:
            guids = set(self.pc.seqProfile.keys())

        if guid not in self.pc.seqProfile.keys():
            raise KeyError(
                "Asked to compare {0}  but guid requested has not been stored in the preComparer.  call .persist() on the sample to be added before using mcompare.".format(
                    guid
                )
            )

        # mcompare using preComparer
        candidate_guids = set()

        neighbours = []
        exact_comparison = False
        for match in self.pc.mcompare(guid, guids):

            if self.pc.distances_are_exact:  # then the precomparer is computing an exact distance
                exact_comparison = True
                if match["dist"] <= self.snpCeiling:

                    neighbours.append([match["guid1"], match["guid2"], match["dist"]])
            else:  # the precomparer is computing an estimated distance, and more detailed comparison is needed.
                if not match["no_retest"]:
                    candidate_guids.add(match["guid2"])

        t2 = datetime.datetime.now()

        if not exact_comparison:
            # need to do second phase computation to determine neighbours
            guids = list(set(guids))

            # load guid
            load_time = 0
            seq1 = self.load(guid)
            for key2 in candidate_guids:
                if not guid == key2:
                    l1 = datetime.datetime.now()
                    seq2 = self.load(key2)
                    l2 = datetime.datetime.now()
                    i4 = l2 - l1
                    load_time = load_time + i4.total_seconds()
                    comparison = self.countDifferences(guid, key2, seq1, seq2, cutoff=self.snpCeiling)
                    if comparison is not None:  # is is none if one of the samples is invalid

                        (guid1, guid2, dist) = comparison
                        if dist <= self.snpCeiling:
                            neighbours.append([guid1, guid2, dist])

        t3 = datetime.datetime.now()
        i1 = t2 - t1
        i2 = t3 - t2
        i3 = t3 - t1
        n1 = len(self.pc.seqProfile.keys())
        n2 = len(candidate_guids)
        if n1 > 0:
            rate1 = 1e3 * i1.total_seconds() / n1
        else:
            rate1 = 0
        if n2 > 0:
            rate2 = 1e3 * i2.total_seconds() / n2
            rate3 = 1e3 * load_time / n2
        else:
            rate2 = 0
            rate3 = 0

        timings = {
            "preComparer_msec_per_comparison": rate1,
            "seqComparer_msec_per_comparison": rate2,
            "preCompared": n1,
            "candidates": n2,
            "matches": len(neighbours),
            "total_sec": i3.total_seconds(),
            "seqComparer_msec_per_sequence_loaded": rate3,
            "catWalk_enabled": self.pc.catWalk_enabled,
            "preComparer_distances_are_exact": self.pc.distances_are_exact,
        }

        return {"neighbours": neighbours, "timings": timings}

    def summarise_stored_items(self):
        """ counts how many sequences exist of various types in the preComparer and seqComparer objects """

        return self.pc.summarise_stored_items()

    def iscachedinram(self, guid):
        """ returns true or false depending whether we have a  copy of the limited refCompressed representation of a sequence (name=guid) in RAM"""
        return self.pc.iscachedinram(guid)

    def guidscachedinram(self):
        """ returns all guids with full sequence profiles in RAM """
        return self.pc.guidscachedinram()

    def _delta(self, x):
        """ returns the difference between two numbers in a tuple x """
        return x[1] - x[0]

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
        """ returns a sequence from a compressed_sequence """
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
        for x in compressed_sequence["M"].keys():  # the 'M' option includes iupac coding
            seq[x] = compressed_sequence["M"][x]
        return "".join(seq)

    def uncompress_guid(self, guid):
        """ uncompresses a saved sequence """
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
        diffDict = {"A": set([]), "C": set([]), "T": set([]), "G": set([]), "N": set([]), "M": {}}

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

    def _setStats(self, i1, i2):
        """compares either:
        * two sets (if i1 or i2 is a set)
        OR
        * the keys of the dictionaries i1 or i2 is a dictionary

        returns
        * the number of elements in i1
        * the number of elements in i2
        * the number of elements in the union of i1 and i2
        * i1
        * i2
        * the union of set1 and set2

        """
        if isinstance(i1, set):
            retVal1 = i1
        elif isinstance(i1, dict):
            retVal1 = set(i1.keys())
        if isinstance(i2, set):
            retVal2 = i2
        elif isinstance(i2, dict):
            retVal2 = set(i2.keys())

        retVal = retVal2 | retVal1

        return (len(retVal1), len(retVal2), len(retVal), retVal1, retVal2, retVal)

    def countDifferences(self, key1, key2, seq1, seq2, cutoff=None):
        """compares seq1 with seq2.

        Ms (uncertain bases) are ignored in snp computations.


        """
        #  if cutoff is not specified, we use snpCeiling
        if cutoff is None:
            cutoff = self.snpCeiling

        nDiff = 0
        if seq1 is None or seq2 is None:
            return None

        if seq1["invalid"] == 1 or seq2["invalid"] == 1:
            return None

        # compute positions which differ;
        differing_positions = set()
        for nucleotide in ["C", "G", "A", "T"]:

            # we do not consider differences relative to the reference if the other nucleotide is an N or M
            nonN_seq1 = seq1[nucleotide] - (seq2["N"] | set(seq2["M"].keys()))
            nonN_seq2 = seq2[nucleotide] - (seq1["N"] | set(seq1["M"].keys()))
            differing_positions = differing_positions | (nonN_seq1 ^ nonN_seq2)

        nDiff = len(differing_positions)

        if nDiff <= cutoff:
            return (key1, key2, nDiff)
        else:
            return (key1, key2, nDiff)

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
            serialised_compressed_sequence = serialised_compressed_sequence + key + ":" + str(list_compressed) + ";"
        h = hashlib.md5()
        h.update(serialised_compressed_sequence.encode("utf-8"))
        md5 = h.hexdigest()
        return md5

    def remove(self, guid):
        """ removes guid  from preComparer """
        self.pc.remove(guid)

    def estimate_expected_proportion(self, seqs):
        """computes the median Ns for seqs, a list.
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

    def estimate_expected_unk(self, sample_size=30, exclude_guids=set(), unk_type="N", what="median"):
        """computes the median numbers of unk_type (N,M,or N_or_N) for sample_size guids, randomly selected from all guids except for exclude_guids.

                Parameters
                sample_size: how many items to sample when determining expected_unk.  If None, uses all samples
                exclude_guids: an iterable of guids not to analyse
                unk_type: one of 'N' 'M' 'N_or_M'
        e
                Return
                either the median unk_type for the sample, or None if fewer than sample_size items are present.

                Note: this is a very fast method, if there are sequences in the preComparer.
        """

        if unk_type not in ["N", "M", "N_or_M"]:
            raise KeyError("unk_type can be one of 'N' 'M' 'N_or_M'")
        current_composition = copy.copy(self.pc.composition)  # can be changed by flask, so duplicate it
        composition = pd.DataFrame.from_dict(current_composition, orient="index")  # preComparer maintains a composition list
        composition.drop(exclude_guids, inplace=True)  # remove the ones we want to exclude

        if len(composition) == 0:

            return None
        if sample_size is not None:
            guids = composition.index.tolist()
            np.random.shuffle(guids)
            not_these_guids = guids[sample_size:]  # if there are more guids than sample_size, drop these
            composition.drop(not_these_guids, inplace=True)

        composition["N_or_M"] = composition["N"] + composition["M"]
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

    def estimate_expected_unk_sites(self, unk_type, sites=set(), sample_size=30):
        """estimates the median unk_type for all guids,  at the positions in sites().

        Parameters
        unk_type: one of 'N' 'M' 'N_or_M'
        sites: the sites to consider, zero indexed
        sample_size: the minimum number of samples from which to report an estimate

        Returns
        either the estimated median unk_type for the sample, or None if fewer than sample_size items are present.

        Note:
        this is an exact method.  Computation is slower as samples are loaded from disc.  Typical load is about 10msec/seq, or 300msec per 30 sample.
        """

        # get samples
        if sample_size <= len(self.pc.seqProfile.keys()):
            to_test = random.sample(self.pc.seqProfile.keys(), sample_size)
        else:
            return None

        observed = []

        for guid in to_test:

            obj = self.PERSIST.refcompressedsequence_read(guid)

            if obj is None:
                # that's an error, should never happen.
                raise KeyError(
                    "A guid {0} found in the preComparer, selected in to_test {1}, was not found in the PERSIST object (call returned None).  Internal software error.".format(
                        guid, to_test
                    )
                )

            try:
                N_sites = set(obj["N"]).intersection(sites)
            except KeyError:
                N_sites = set()
            try:
                M_sites = set(obj["M"].keys()).intersection(sites)
            except KeyError:
                M_sites = set()

            relevant = 0

            if unk_type in ["M", "N_or_M"]:
                relevant = relevant + len(M_sites)
            if unk_type in ["N", "N_or_M"]:
                relevant = relevant + len(N_sites)

            observed.append(relevant)

        return np.median(observed)

    def multi_sequence_alignment(self, guids, sample_size=30, expected_p1=None, uncertain_base_type="N", outgroup=None):
        """computes a multiple sequence alignment containing only sites which vary between guids.

            Parameters:
            guids: an iterable of guids to analyse

            sample_size is the number of samples to randomly sample to estimate the expected number of N or Ms in
            the population of sequences currently in the server.  From this, the routine computes expected_p1,
            which is expected_expected_N/M the length of sequence.
            if expected_p1 is supplied, then such sampling does not occur.

            uncertain_base_type: the kind of base which is to be analysed, either N,M, or N_or_M
            outgroup: the outgroup sample, if any.  Not used in computations, but stored in the output

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
        if expected_p1 is not None:
            if expected_p1 < 0 or expected_p1 > 1:
                raise ValueError("Expected_p1 must lie between 0 and 1")
        if sample_size is None:
            sample_size = 30

        # step 0: find all valid guids
        valid_guids = []
        invalid_guids = []

        self.remove_all_temporary_seqs()  # remove all full seqProfiles from this object
        for guid in guids:
            seq = self.load(guid)
            if seq is None:
                raise KeyError("Sequence missing {0}".format(guid))

            if seq["invalid"] == 0:
                self.seqProfile[guid] = seq
                valid_guids.append(guid)

            else:

                invalid_guids.append(guid)

        # Estimate expected N or M as median(observed N or Ms),
        # which is a valid thing to do if the proportion of mixed samples is low.
        # this is a very fast method which won't materially slow down the msa
        if expected_p1 is None:
            expected_N1 = self.estimate_expected_unk(
                sample_size=sample_size, exclude_guids=invalid_guids, unk_type=uncertain_base_type
            )
            if expected_N1 is None:
                expected_p1 = None
            else:
                expected_p1 = expected_N1 / len(self.reference)
        else:
            expected_N1 = np.floor(expected_p1 * len(self.reference))

        return self._msa(valid_guids, invalid_guids, expected_p1, sample_size, uncertain_base_type, outgroup=outgroup)

    def _msa(self, valid_guids, invalid_guids, expected_p1, sample_size, uncertain_base_type="N", outgroup=None):
        """perform multisequence alignment and significance tests.
        It assumes that the relevant sequences (valid_guids) are in-ram in this hybridComparer object as part of the seqProfile dictionary.

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

        for guid in valid_guids:

            if self.seqProfile[guid]["invalid"] == 1:
                raise TypeError("Invalid sequence {0} passed in valid_guids".format(guid))

            for unk_type in ["N", "M"]:

                unks = len(self.seqProfile[guid][unk_type])
                guid2all[unk_type][guid] = unks

            guid2all["N_or_M"][guid] = guid2all["N"][guid] + guid2all["M"][guid]

            for base in ["A", "C", "T", "G"]:
                positions = self.seqProfile[guid][base]
                for position in positions:
                    if position not in nrps.keys():  # if it's non-reference, and we've got no record of this position
                        nrps[position] = set()  # then we generate a set of bases at this position
                    nrps[position].add(base)  # either way add the current non-reference base there

        # step 2: for the non-reference called positions,
        # record whether there's a reference base there.
        for guid in valid_guids:
            for position in nrps.keys():
                psn_accounted_for = 0
                for base in ["A", "C", "T", "G"]:
                    if position in self.seqProfile[guid][base]:
                        psn_accounted_for = 1
                if psn_accounted_for == 0:  # no record of non-ref seq here
                    if not (
                        position in self.seqProfile[guid]["N"] or position in self.seqProfile[guid]["M"].keys()
                    ):  # not M or N so non-reference
                        # it is reference; this guid has no record of a variant base at this position, so it must be reference.
                        nrps[position].add(self.reference[position])  # add reference psn

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
                        positions = self.seqProfile[guid][base]
                    else:
                        positions = self.seqProfile[guid][base].keys()

                    if position in positions:
                        this_base = base
                        this_base_m = base
                        if base == "M":
                            this_base_m = self.seqProfile[guid]["M"][position]
                guid2seq[guid].append(this_base)  # include M for mix if present
                guid2mseq[guid].append(this_base_m)  # include iupac code for mix
            guid2msa_seq[guid] = "".join(guid2seq[guid])
            guid2msa_mseq[guid] = "".join(guid2mseq[guid])

        # step 5: determine the expected_p2 at the ordered_variant_positions:
        expected_N2 = self.estimate_expected_unk_sites(sites=set(ordered_variant_positions), unk_type=uncertain_base_type)
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
                guid2align["N_or_M"][guid] = guid2align["N"][guid] + guid2align["M"][guid]

                if expected_p1 is None:  # we don't have an expectation, so we can't assess the first binomial test;
                    p_value1 = None
                    observed_p = None
                elif len(guid2msa_seq[guid]) == 0:  # we don't have any information to work with
                    p_value1 = None
                    observed_p = None
                else:
                    observed_p = guid2align[uncertain_base_type][guid] / len(guid2msa_seq[guid])
                    p_value1 = binom_test(
                        guid2align[uncertain_base_type][guid], len(guid2msa_seq[guid]), expected_p1, alternative="greater"
                    )

                guid2pvalue1[guid] = p_value1
                guid2observed_p[guid] = observed_p
                guid2expected_p1[guid] = expected_p1

                # compute p value 2.  This tests the hypothesis that the number of Ns or Ms in the *alignment*
                # is GREATER than those expected from the expected_N in the population of whole sequences
                # at these sites.
                if expected_p2 is None:  # we don't have an expectation, so we can't assess the binomial test;
                    p_value2 = None
                elif len(guid2msa_seq[guid]) == 0:  # we don't have any information to work with
                    p_value2 = None
                else:

                    p_value2 = binom_test(
                        guid2align[uncertain_base_type][guid], len(guid2msa_seq[guid]), expected_p2, alternative="greater"
                    )
                guid2pvalue2[guid] = p_value2
                guid2expected_p2[guid] = expected_p2

                # compute p value 3.  This tests the hypothesis that the number of Ns or Ms in the alignment of THIS SEQUENCE
                # is GREATER than the number of Ns or Ms not in the alignment  IN THIS SEQUENCE
                # based on sequences not in the alignment

                expected_p3 = (guid2all[uncertain_base_type][guid] - guid2align[uncertain_base_type][guid]) / (
                    len(self.reference) - len(guid2msa_seq[guid])
                )
                p_value = binom_test(
                    guid2align[uncertain_base_type][guid], len(guid2msa_seq[guid]), expected_p3, alternative="greater"
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
                        guid2align[uncertain_base_type][guid], len(guid2msa_seq[guid]), expected_p4, alternative="greater"
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

    def remove_all_temporary_seqs(self):
        """ empties any sequence data from ram in the .seqProfile dictionary; these are used for MSAs.  Does not affect the preComparer. """
        self.seqProfile = {}

    def raise_error(self, token):
        """raises a ZeroDivisionError, with token as the message.
        useful for unit tests of error logging"""
        raise ZeroDivisionError(token)
