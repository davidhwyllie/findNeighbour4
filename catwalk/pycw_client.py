"""
A catwalk client for python

Tested on Linux with python 3.9
Uses linux-specific commands to run the server, so is not expected to work on Windows.

unittests:
pipenv run python -m unittest test/test_pycw_client.py

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

"""

import shlex
import time
import json
import requests
import logging
import os
import psutil
import uuid
import warnings

from findn.rdbmsstore import fn3persistence_r
from findn.mongoStore import fn3persistence

class CatWalkNoneSequenceError(Exception):
    """sequence is None; this cannot be loaded"""

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


class CatWalkServerInsertError(Exception):
    """insert failed"""

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


class CatWalkServerDeleteError(Exception):
    """delete failed"""

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


class CatWalkServerDidNotStartError(Exception):
    """the catwalk server did not start"""

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


class CatWalkStartupLockNotAcquiredError(Exception):
    """the acquisition of a lock to prevent race conditions did not succeed"""

    def __init__(self):
        pass


class CatWalkBinaryNotAvailableError(Exception):
    """no catwalk binary"""

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


class CatWalkMultipleServersRunningError(Exception):
    """multiple catwalk servers with identical specification are running"""

    def __init__(self, message):
        self.message = message


class CatWalk:
    """start, stop, and communicate with a CatWalk server"""

    def __init__(
        self,
        cw_binary_filepath,
        reference_name,
        reference_filepath,
        mask_filepath,
        max_n_positions,
        bind_host,
        bind_port,
        identity_token=None,
        unittesting=False,
        lockmanager=None,
    ):
        """
        Start the catwalk process in the background, if it not running.
        Parameters:
        cw_binary_filepath: where the catwalk binary is
        reference_name:     parameter passed to catwalk binary;  the name of the catwalk instance
        reference_filepath: parameter passed to catwalk binary: where there reference fasta file is
        mask_filepath:      parameter passed to catwalk binary: where the file containing any mask is
        max_n_positions:    parameter passed to catwalk binary: maximum number of Ns to permit ; sequences with higher numbers are not used in comparisons
        bind_host:          parameter passed to catwalk binary: the ip on which catwalk is running.  normally 'localhost'
        bind_port:          parameter passed to catwalk binary: the port on which catwalk is running

        identity_token: a string identifying the process.  If not provided, a random guid is generated
        unittesting: if True, will shut down and restart (empty) any catwalk on bind_port on creation

        lockmanager: if None (default) no lockmanager is used during catwalk startup/shutdown operations.  This is appropriate for situations in which only one software is interacting with catwalk; however,
                     if multiple processes are interacting with catwalk it is possible that race conditions can result such that multiple catwalks can be started, which only one is allowed.
                     to address this, a lockmanager has to be provided.
                     at present, the only lockmanager allowed is a findneighbour4 PERSISTENCE object, as instantiated by Persistence.get_storage_object() method.
                     Please see findn/persistence.py for details about how this should be instantiated.
        """

        # a message used to inform if a catwalk server application could not be found
        no_catwalk_exe_message = """
The catWalk client could not find a CatWalk server application.  This is because either:
i) the cw_binary_filepath parameter was None
ii) the above, and there is no CW_BINARY_FILEPATH environment variable.  To set this, you insert a line like
CW_BINARY_FILEPATH=/path/to/cw/executable
in either
.bashrc - if you're not using a virtual environment
.env    - a file in the same directory as the PipFile, if you are using a virtual environment.
          This file is not committed into the repository, so you'll have to create it once in your installation.
"""
        # if cw_binary_filepath is "", we regard it as not specified (None)
        if cw_binary_filepath == "":
            cw_binary_filepath = None

        # if cw_binary_filepath is None, we check the env. variable CW_BINARY_FILEPATH and use that if present.
        if cw_binary_filepath is None:
            if "CW_BINARY_FILEPATH" in os.environ:
                cw_binary_filepath = os.environ["CW_BINARY_FILEPATH"]
            else:
                raise CatWalkBinaryNotAvailableError(
                    expression=None, message=no_catwalk_exe_message
                )
        if not os.path.exists(cw_binary_filepath):
            raise FileNotFoundError(
                "Was provided a cw_binary_filepath, but there is no file there {0}".format(
                    cw_binary_filepath
                )
            )

        # store parameters
        self.bind_host = bind_host
        self.bind_port = int(bind_port)  # has to be an integer
        self.cw_url = "http://{0}:{1}".format(bind_host, self.bind_port)
        self.cw_binary_filepath = cw_binary_filepath
        self.reference_filepath = reference_filepath
        self.mask_filepath = mask_filepath
        self.max_n_positions = int(max_n_positions)
        self.reference_name = reference_name
        self.instance_stem = "CatWalk-PORT-{0}".format(self.bind_port)
        if identity_token is None:
            identity_token = str(uuid.uuid1())
        self.instance_name = "{0}-MAXN-{1}-{2}".format(
            self.instance_stem, self.max_n_positions, identity_token
        )

        # check the lockmanager passed
        self.PERSIST = lockmanager
        if self.PERSIST is not None:
            if not isinstance(self.PERSIST, (fn3persistence, fn3persistence_r)):
                raise TypeError(
                    "lockmanager must be either fn3persistence or fn3persistence_r object"
                )

        # start up if not running

        # if we are unittesting and a server is running, we shut it down
        if unittesting and self.server_is_running():
            self.stop()  # removes any data from server  and any other running cws

        # try to start up.  Will have no effect if the server is already running.
        self.start()

    def _running_servers(self):
        """returns details of running servers matching the details of this server

        There should be either 0 or 1 of these only"""

        servers = []
        for proc in psutil.process_iter():
            if "cw_server" in proc.name():
                cmdline_parts = proc.cmdline()
                for i, cmdline_part in enumerate(proc.cmdline()):
                    if cmdline_part == "--instance_name":
                        if cmdline_parts[i + 1].startswith(self.instance_stem):
                            servers.append(cmdline_parts)
        return servers

    def server_is_running(self):
        """returns true if the relevant process is running, otherwise false.

        The alternative strategy, returning true if a response is received by the server,
        can result in reporting false if the server is busy"""

        servers = self._running_servers()
        if len(servers) == 0:
            return False
        elif len(servers) == 1:
            return True
        else:
            # raise an error. See also issue #117
            raise CatWalkMultipleServersRunningError(
                message="{0} servers with specification {1} detected.  This is not permitted and should not occur.  Servers running are:{2}".format(
                    len(servers), self.instance_stem, servers
                )
            )  # there cannot be multiple servers running

    def start(self):
        """starts a catwalk process in the background"""

        # don't startup if it's running
        if self.server_is_running():
            return

        cw_binary_filepath = shlex.quote(self.cw_binary_filepath)
        instance_name = shlex.quote(self.instance_name)
        reference_filepath = shlex.quote(self.reference_filepath)
        mask_filepath = shlex.quote(self.mask_filepath)

        # prevent potential race conditions where different processes try to start catwalk
        n_tries = 0
        lock_acquired = False
        if self.PERSIST is not None:
            while n_tries < 6 and lock_acquired is False:
                n_tries = n_tries + 1
                lock_acquired = self.PERSIST.lock(2, "startup_catwalk")
                if (
                    lock_acquired is False
                ):  # true if an catwalk startup lock was acquired
                    logging.info(
                        "Catwalk startup lock could not be acquired, try = {0}/6.  Waiting 2 seconds".format(
                            n_tries
                        )
                    )
                    time.sleep(2)

            if lock_acquired is False:
                # lock acquisition failed, indicative of a problem as yet unrecognised
                info_msg = """A lock to prevent race conditions on starting catwalk could not be acquired.
                    Tried 6 times at 2 second intervals.
                    This may arise because
                    (i) some failure of catwalk to start has occurred
                    (ii) the lock is held inappropriately after a crash. 
                    The findNeighbour4_lockmonitor should unlock it automatically in 90 seconds.  
                    If needed, you can reset the lock as follows:  fn4_configure <path to config file> --drop_locks"""
                logging.warning(info_msg)
                raise CatWalkStartupLockNotAcquiredError()

        if self.server_is_running() is False:
            cmd = f"nohup {cw_binary_filepath} --instance_name {instance_name}  --bind_host {self.bind_host} --bind_port {self.bind_port} --reference_filepath {reference_filepath}  --mask_filepath {mask_filepath} --max_n_positions {self.max_n_positions} > cw_server_nohup.out &"
            logging.info("Attempting startup of CatWalk server : {0}".format(cmd))

            os.system(
                cmd
            )  # synchronous: will return when it has started; runs via nohup

            time.sleep(1)  # short break to ensure it has started

        # check there is exactly one running.
        # will raise an error otherwise
        if self.server_is_running() is False:
            raise CatWalkServerDidNotStartError()

        info = self.info()  # functional test: test responsiveness
        if info is None:
            raise CatWalkServerDidNotStartError()
        else:
            logging.info("Catwalk server running: {0}".format(info))

        if self.PERSIST is not None:
            # release lock
            self.PERSIST.unlock(2)  # release the lock if it was acquired

    def stop(self):
        """stops a catwalk server launched by with this specification, if running.
        The server is not killed by pid, so if a different process (but with the same identifier)
        is running started by a different program, this will be shut down.

        If more than one catwalk server with the same id is running (this is an error condition)
        then all will be shut down.

        The server process is identified by the instance_stem, which is of the for Catwalk-PORT-XXXX-hash
        where the hash is the hash on an 'identity token' passed to the constructor, see above."""

        for proc in psutil.process_iter():
            if "cw_server" in proc.name():
                cmdline_parts = proc.cmdline()
                for i, cmdline_part in enumerate(proc.cmdline()):
                    if cmdline_part == "--instance_name":
                        if cmdline_parts[i + 1].startswith(self.instance_stem):
                            proc.kill()
        if self.server_is_running():
            warnings.warn(
                "Attempt to shutdown a catwalk process with name cw_server and --instance_name beginning with {0} failed.  It may be that another process (with a different instance name) is running on port {1}.  Review running processes (ps x) and kill residual processes on this port  manually if appropriate".format(
                    self.instance_stem, self.bind_port
                )
            )
        if self.PERSIST is not None:
            # release lock
            self.PERSIST.unlock(2)  # release the lock if it was acquired

    def stop_all(self):
        """stops all catwalk servers"""
        nKilled = 0
        for proc in psutil.process_iter():
            if "cw_server" in proc.name():
                proc.kill()
                nKilled += 1
        if nKilled > 0:
            warnings.warn(
                "Catwalk client.stop_all() executed. Kill instruction issues on {0} processes.  Beware, this will kill all cw_server processes on the server, not any specific one".format(
                    nKilled
                )
            )
        if self.PERSIST is not None:
            # release lock
            self.PERSIST.unlock(2)  # release the lock if it was acquired

    def info(self):
        """
        Get status information from catwalk
        """
        target_url = "{0}/info".format(self.cw_url)

        r = requests.get(target_url)
        r.raise_for_status()  # report errors

        return r.json()

    def _filter_refcomp(self, refcomp):
        """examines the keys in a dictionary, refcomp, and only lets through keys with a list
        This will remove Ms (linked to a dictionary) and invalid keys which are linked to an integer.
        These keys are not required by catwalk
        """
        refcompressed = {}
        for key in refcomp.keys():
            if isinstance(refcomp[key], set):
                refcompressed[key] = list(refcomp[key])
            elif isinstance(refcomp[key], list):
                refcompressed[key] = refcomp[key]
            else:
                pass  # drop everything else, such as invalid
        return refcompressed

    def add_sample_from_refcomp(self, name, refcomp):
        """
        Add a reference compressed (dict with ACGTN keys and list of positions as values) sample.
        Note, if the sample already exists, it will not be added twice.
        The json dict must have all keys: ACGTN, even if they're empty

        Returns:
        status code
        201 = added successfully
        200 = was already present
        """

        # note particular way of creating json, but catwalk accepts this (refcomp has to be in '')
        # cannot json serialise sets; use lists instead
        if refcomp is None:
            # raise error
            raise CatWalkNoneSequenceError(
                "Asked to reload catwalk with {0} but the refcomp was None".format(name)
            )

        refcompressed = self._filter_refcomp(refcomp)
        payload = {"name": name, "refcomp": json.dumps(refcompressed), "keep": True}

        r = requests.post(
            "{0}/add_sample_from_refcomp".format(self.cw_url), json=payload
        )
        r.raise_for_status()
        if r.status_code not in [200, 201]:
            raise CatWalkServerInsertError(
                message="Failed to insert {0}; return code was {1}".format(name, r.text)
            )
        return r.status_code

    def remove_sample(self, name):
        """deletes a sample called name"""

        r = requests.get("{0}/remove_sample/{1}".format(self.cw_url, name))
        r.raise_for_status()
        if r.status_code not in [200]:
            raise CatWalkServerDeleteError(
                message="Failed to delete {0}; return code was {1}".format(name, r.text)
            )
        return r.status_code

    def neighbours(self, name, distance=None):
        """get neighbours.  neighbours are recomputed on demand.

        Parameters:
        name:  the name of the sample to search for
        distance: the maximum distance reported.  if distance is not supplied, 99 is used.
        """
        if not distance:
            logging.warning("no distance supplied. Using 99")
            distance = 99

        distance = int(distance)  # if a float, url contstruction may fail

        r = requests.get("{0}/neighbours/{1}/{2}".format(self.cw_url, name, distance))
        r.raise_for_status()
        j = r.json()
        return [(sample_name, int(distance_str)) for (sample_name, distance_str) in j]

    def sample_names(self):
        """get a list of samples in catwalk"""
        r = requests.get("{0}/list_samples".format(self.cw_url))
        r.raise_for_status()
        return r.json()
