#!/usr/bin/env python3
""" manages multisequence alignments 

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


import pandas as pd
import datetime
from findn.identify_sequence_set import IdentifySequenceSet
from findn.msaviewer import DepictMSA


class MSAStore:
    """stores multisequence alignments, both in ram and in a database"""

    def __init__(self, PERSIST, in_ram_persistence_time=60):
        """sets up an MSA cache, which stores MSAs on disc and in ram

        Parameters:
        PERSIST is an fn3persistence object
        in_ram_persistence_time the time MSAs are cached in ram, in seconds, by the findNeighbour4 server
        """
        self.cache = {}
        self.PERSIST = PERSIST
        self.iss = IdentifySequenceSet()
        self.in_ram_persistence_time = in_ram_persistence_time
        self.tree = None  # placeholder, not used
        self.comment = None  # placeholder, not used

    def get_token(self, what, has_outgroup, guids):
        """computes a token representing a set of guids with or without an outgroup"""
        return self.iss.make_identifier("msa", what, has_outgroup, guids)

    def existing_tokens(self):
        """returns a set of all tokens on disc"""
        return set(self.PERSIST.msa_stored_ids())

    def is_in_ram(self, token):
        """whether the msa identifed by token is cached in ram"""
        return token in self.cache.keys()

        return token in self.cache.keys()

    def persist(self, token, msa_result):
        """stores a MSAResult object in mongo identified by token"""
        self.PERSIST.msa_store(token, msa_result.serialise())
        self.cache_in_ram(token, msa_result)

    def load(self, token):
        """recovers an MSAResult identified by token.  Caches it in ram"""
        persisted_dict = self.PERSIST.msa_read(token)
        if persisted_dict is not None:
            msa_result = MSAResult(**persisted_dict)
            self.cache_in_ram(token, msa_result)
            return msa_result
        else:
            return None

    def cache_in_ram(self, token, msa_result):
        """stores a MSAResult object in ram identified by token"""
        self.purge_ram()
        now = datetime.datetime.now()
        expiry_time = now + datetime.timedelta(seconds=self.in_ram_persistence_time)
        if token not in self.cache.keys():
            self.cache[token] = {"expiry_time": expiry_time, "MSA": msa_result}
        else:
            self.cache[token]["expiry_time"] = expiry_time  # update expiry time

    def purge_ram(self):
        """removes any in-ram objects which are time expired"""
        now = datetime.datetime.now()
        to_delete = set()
        for token in self.cache.keys():
            if now > self.cache[token]["expiry_time"]:
                to_delete.add(token)
        for token in to_delete:
            del self.cache[token]

    def unpersist(self, whitelist):
        """removes any on-disc msa results the ids of which are not in the whitelist"""
        self.PERSIST.msa_delete_unless_whitelisted(whitelist)


class MSAResult:
    """a representation of a multisequence alignment"""

    def __init__(
        self,
        variant_positions,
        invalid_guids,
        valid_guids,
        expected_p1,
        sample_size,
        df_dict,
        what_tested,
        outgroup,
        creation_time,
        fconst,
    ):

        """a representation of a multisequence alignment
        The multisequence alignment is generates by SeqComparer.multi_sequence_alignment(), not by this class.

        Parameters:
                    'variant_positions': positions of variation in msa
                    'invalid_guids': invalid guids associated with, but not included in, msa
                    'valid_guids': guids in msa
                    'expected_p1': expected_p1 (see ._msa in hybridComparer for discussion)
                    'sample_size': sample size used to estimate p1
                    'df_dict': a dictionary serialising the pandas dataframe containing the msa
                    'what_tested': what (M/N/N_or_M) was tested
                    'outgroup': either none, or a guid used as the outgroup.
                    'creation_time': timestamp the msa was created
                    'fconst': a dictionary containing the number of constant sites (A,C,G,T) in the reference genome but OUTSIDE the alignment.
                              Can be an empty dictionary {} if this concept is not relevant.  May be required by maximal likelihood tree gen software, e.g. iqTree.
        """

        self.variant_positions = variant_positions
        self.alignment_length = len(variant_positions)
        self.invalid_guids = invalid_guids
        self.valid_guids = valid_guids
        self.expected_p1 = expected_p1
        self.sample_size = sample_size
        self.df = pd.DataFrame.from_dict(df_dict, orient="index")
        self.fconst = fconst
        self.what_tested = what_tested
        self.outgroup = outgroup
        self.creation_time = creation_time
        self.annotation = None

        iss = IdentifySequenceSet()
        self.token = iss.make_identifier(
            "msa", self.what_tested, outgroup is not None, valid_guids
        )

    def msa_df(self):
        """return the msa as a dataframe"""
        return self.df

    def msa_html(self):
        """returns the msa as an html table"""
        return self.df.to_html()

    def msa_fasta(self):
        """returns the msa as fasta"""
        fasta = ""
        for guid in self.df.index:
            fasta = fasta + ">{0}\n{1}\n".format(guid, self.df.loc[guid, "aligned_seq"])
        return fasta

    def msa_dict(self):
        """return the msa as a dictionary"""
        return self.df.to_dict(orient="index")

    def msa_interactive_depiction(self):
        """returns an interactive depiction of the msa object.
        This is an html file string, with embedded data and javascript.
        rendering it will produce an interactive visualisation"""
        dep_msa = DepictMSA(
            self.df, positions_analysed=sorted(list(self.variant_positions))
        )
        return dep_msa.render_msa()

    def serialise(self):
        """return the entire object in a serialisable format (dictionary)"""
        return {
            "variant_positions": self.variant_positions,
            "invalid_guids": self.invalid_guids,
            "valid_guids": self.valid_guids,
            "expected_p1": self.expected_p1,
            "sample_size": self.sample_size,
            "df_dict": self.df.to_dict(orient="index"),
            "what_tested": self.what_tested,
            "outgroup": self.outgroup,
            "creation_time": self.creation_time,
            "fconst": self.fconst,
        }
