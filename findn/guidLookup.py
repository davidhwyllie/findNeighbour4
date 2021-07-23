#!/usr/bin/env python
""" 
modules to search for guids (sampleIds) rapidly

used for autocompletion in web front end
fast bisection based search is used

Unittests:
pipenv run python3 -m unittest test/test_guid_lookup.py

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

import bisect
import datetime


class guidSearcher:
    """stores a list of sample numbers (such as guids) and provides
    methods to add more, and to search them rapidly"""

    def __init__(self):
        """initialise a guidSearcher.

        Example usage:

        gs = guidSearcher()
        gs.add('b1')
        gs.add('b3')
        gs.add('b3')
        gs.add('b2')
        gs.add('a1')
        gs.add('c1')
        print(gs.guids)


        """
        self.guids = []

    def add(self, guid):
        """add a string ('guid') into an ordered list"""
        self._add(guid)

    def _add(self, guid):
        """add a string ('guid') into an ordered list"""
        insertion_point = bisect.bisect_left(self.guids, guid)

        if not isinstance(guid, str):
            raise TypeError("added item must be a string, not a {0}".format(type(guid)))

        # test whether the item exists
        already_exists = False
        try:
            if guid == self.guids[insertion_point]:
                already_exists = True
        except IndexError:
            # doesn't exist
            pass
        if not already_exists:
            self.guids.insert(insertion_point, guid)

    def search(self, search_string, max_returned=30, return_subset=True):
        """
        Searches guids stored in the guidSearcher object for search_string
        search_string  the substring in self.guids sought at the beginning of the string
        max_returned   the maximum number of matches returned
        return_subset  whether to return max_returned records if more than max_returned records are found.
        """
        return self._search(search_string, max_returned, return_subset)

    def _search(self, search_string, max_returned, return_subset):
        """
        Searches guids stored in the guidSearcher object for search_string
        search_string  the substring in self.guids sought at the beginning of the string
        max_returned   the maximum number of matches returned
        return_subset  whether to return max_returned records if more than max_returned records are found.
        """
        search_start = bisect.bisect_left(self.guids, search_string)
        retVal = []

        for i in range(max_returned + 1):
            search_point = i + search_start
            if search_point > len(self.guids) - 1:  # got to the end
                break
            else:
                if self.guids[search_point].startswith(search_string):
                    retVal.append(self.guids[search_point])
            if len(retVal) == max_returned and return_subset is True:
                return retVal
            elif len(retVal) == max_returned and return_subset is False:
                # don't return partial lists of matches
                return []

        return retVal


class guidDbSearcher(guidSearcher):
    """maintains a list of sample numbers (such as guids), keeping the list synchronised with an underlying
    PERSISTENCE object, which is a database"""

    def __init__(self, PERSIST, recheck_interval_seconds=60):
        """create a guidDbSearcher object, which recovers its guids from the fn3/4 persistence object PERSIST

        parameters:
        PERSIST: an fn3/4 persistence object
        recheck_interval_seconds:  integer, how many seconds to wait before rechecking PERSIST for new samples."""
        self.PERSIST = PERSIST
        self.recheck_interval_seconds = recheck_interval_seconds
        self.last_check_time = datetime.datetime.now()
        self.next_recheck_time = self.last_check_time + datetime.timedelta(
            seconds=recheck_interval_seconds
        )

        self.guids = sorted(self.PERSIST.guids())  # return a list

        if not (recheck_interval_seconds, int):
            raise TypeError(
                "recheck_interval_seconds must be an integer, not a {0}.  Value={1}".format(
                    type(recheck_interval_seconds), recheck_interval_seconds
                )
            )

    def add(self, guid):
        """not implemented: to add a sample to a guidDbSearcher, add it to the underlying PERSIST object (database)"""
        raise NotImplementedError(
            "guidDbSearcher objects do not allow direct addition of samples.  Just ensure the database"
        )

    def search(self, search_string, max_returned=30, return_subset=True):
        """
        Searches guids stored in the guidSearcher object for search_string
        search_string  the substring in self.guids sought at the beginning of the string
        max_returned   the maximum number of matches returned
        return_subset  whether to return max_returned records if more than max_returned records are found.
        """

        # test whether we need to check for more samples
        if datetime.datetime.now() > self.next_recheck_time:
            new_guids = self.PERSIST.guids_considered_after(self.last_check_time)
            self.last_check_time = datetime.datetime.now()
            self.next_recheck_time = self.last_check_time + datetime.timedelta(
                seconds=self.recheck_interval_seconds
            )
            for this_guid in new_guids:
                self._add(this_guid)

        return self._search(search_string, max_returned, return_subset)
