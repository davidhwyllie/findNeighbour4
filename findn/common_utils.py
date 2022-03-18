#!/usr/bin/env python
""" methods for reading the findNeighbour config file 

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

 
"""
import os
import glob
import json
import datetime
from Bio import SeqIO
import logging
from findn.persistence import Persistence
from pathlib import Path
from typing import Union

PathLike = Union[str, Path]
ConfigLike = Union[str, dict]


class EnvWriter:
    """a utility for reading .env files, and for adding and removing elements from them."""

    def __init__(self, env_file=".env"):
        """reads a .env file  and prepares to manipulate it"""
        self.env_file = env_file
        if not os.path.exists(self.env_file):
            raise FileNotFoundError(".env file {0} was not found".format(self.env_file))

        self.env_vars = {}
        with open(env_file, "rt") as f:
            for i, x in enumerate(f.readlines()):
                x = x.strip()
                if len(x) > 0:  # skip blank lines
                    if "=" not in x:
                        raise ValueError(
                            "Line #{0}  ({1}, length={2}) does not contain '='".format(
                                i, x, len(x)
                            )
                        )

                    key, value = x.split("=")
                    self.env_vars[key.strip()] = value.strip()

    def number_of_env_vars(self):
        """returns the number of environment variables read from file"""
        return len(self.env_vars.keys())

    def set_env_var(self, key, value):
        """set environment variable key to value value.  Note: if the value needs quoting when exported to file, do this e.g. ev.set_env_var("new_key", "'/path/to/file'")"""
        logging.info("Reset environment variable {0} = {1}".format(key, value))
        self.env_vars[key] = value

    def del_env_var(self, key):
        """deletes environment variable key if present"""
        if key in self.env_vars.keys():
            del self.env_vars[key]

    def save_changes(self):
        with open(self.env_file, "wt") as f:
            for key in self.env_vars.keys():
                value = self.env_vars[key]
                f.write("{0}={1}\n".format(key, value))


class ConfigManager:
    """reads, and where approporiate modifies from environmental variables containing secret
    configuration parameters etc, a findNeighbour config dictionary

    Can read from file, or (where configuration is present) from database.

    Background
    ==========
    CONFIG files are json files with format similar to the below:


    An example CONFIG is below:

    {
    "DESCRIPTION":"A test server operating in ../unittest_tmp, only suitable for testing",
    "IP":"127.0.0.1",
    "INPUTREF":"reference/TB-ref.fasta",
    "EXCLUDEFILE":"reference/TB-exclude.txt",
    "DEBUGMODE":0,
    "SERVERNAME":"TBSNP",
    "FNPERSISTENCE_CONNSTRING":"mongodb://127.0.0.1",
    "MAXN_STORAGE":100000,
    "MAXN_PROP_DEFAULT":0.70,
    "PRECOMPARER_PARAMETERS":{},
    "LOGFILE":"../unittest_tmp/logfile.log",
    "LOGLEVEL":"INFO",
    "SNPCEILING": 20,
    "SERVER_MONITORING_MIN_INTERVAL_MSEC":0,
    "SENTRY_URL":"https://c******************@sentry.io/1******",
    "CLUSTERING":{'SNV12_ignore' :{'snv_threshold':12, 'mixed_sample_management':'ignore', 'mixture_criterion':'pvalue_1', 'cutoff':0.001},
                    'SNV12_include':{'snv_threshold':12, 'mixed_sample_management':'include', 'mixture_criterion':'pvalue_1', 'cutoff':0.001}
                    },
    "LISTEN_TO":"127.0.0.1"
    }

    where the database connection binds to
    FNPERSISTENCE_CONNSTRING
    Note: if a FNPERSISTENCE_CONNSTRING environment variable is present, then the value of this will take precedence over any values in the config file.
    This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuration file.

    related to what monitoring the server uses
    SERVER_MONITORING_MIN_INTERVAL_MSEC (optional)

    related to error handling
    SENTRY_URL (optional)
    Note: if a FN_SENTRY URL environment variable is present, then the value of this will take precedence over any values in the config file.
    This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
    PERSIST is a storage object needs to be supplied.
    """

    def __init__(self, config_fpath):
        """creates a config manager, which reads data from either
        - the fn3persistence object PERSIST (a class accessing a mongodb, or None) OR
        - on disc config files, except in first-run situations"""

        self.config_fpath = config_fpath
        self.CONFIG = None

    def _enforce_key_presence(self, key_dict: dict, required_keys: dict) -> None:
        """check that the required keys are in config file"""
        missing_keys = set(required_keys.keys()) - set(key_dict.keys())
        if len(missing_keys) > 0:
            raise KeyError(f"Keys: {missing_keys} are required but were not found.")

    def _enforce_config_object_type(self, config_like: ConfigLike) -> dict:
        """require that the config object is valid json, or a dictionary"""
        CONFIG = config_like
        if isinstance(CONFIG, str):
            CONFIG = json.loads(CONFIG)  # assume JSON string; convert.

        if not isinstance(CONFIG, dict):
            raise KeyError(
                f"Configuration object must be either a dictionary or a JSON string encoding a dictionary, but is {type(CONFIG)}"
            )
        return CONFIG

    def read_config(self, not_debug_mode=False):
        """reads a configuration dictionary from file, with persistence to disc in first-run situations

        returns: configuration dictionary

        not_debug_mode is required for unit testing of clustering and lockmanagement, but should not otherwise be needed

        """

        # read the results from disc
        required_keys = {
            "IP": True,
            "REST_PORT": True,
            "DEBUGMODE": True,
            "LOGFILE": True,
            "MAXN_PROP_DEFAULT": True,
        }
        self.CONFIG = self._read_config_from_file(
            self.config_fpath, required_keys=required_keys
        )

        # check precomparer parameters
        if len(self.CONFIG["PRECOMPARER_PARAMETERS"]) > 0:
            # these are supplied
            observed_keys = set(list(self.CONFIG["PRECOMPARER_PARAMETERS"].keys()))
            expected_keys = set(
                [
                    "selection_cutoff",
                    "over_selection_cutoff_ignore_factor",
                    "uncertain_base",
                ]
            )

            missing = expected_keys - observed_keys
            if len(missing) > 0:
                raise KeyError(
                    "Precomparer parameters were supplied, but the required keys were not found . Missing are {0}".format(
                        missing
                    )
                )

        debug_status = self.CONFIG["DEBUGMODE"]
        if not_debug_mode:
            debug_status = 0

        pm = Persistence()
        self.PERSIST = pm.get_storage_object(
            dbname=self.CONFIG["SERVERNAME"],
            connString=self.CONFIG["FNPERSISTENCE_CONNSTRING"],
            debug=debug_status,
            verbose=False,
        )

        do_not_persist_keys = set(
            [
                "IP",
                "SERVERNAME",
                "FNPERSISTENCE_CONNSTRING",
                "LOGFILE",
                "LOGLEVEL",
                "REST_PORT",
                "SENTRY_URL",
                "SERVER_MONITORING_MIN_INTERVAL_MSEC",
            ]
        )

        if self.PERSIST.first_run():
            # we don't persist some things; some are secret, others might change
            self._first_run(do_not_persist_keys)

        # load the result from database
        stored_config = self.PERSIST.config_read("config")

        stored_config["excludePositions"] = set(stored_config["excludePositions"])
        for (
            key
        ) in (
            do_not_persist_keys
        ):  # update stored config with any of the do_not_persist_keys
            if key in self.CONFIG.keys():
                stored_config[key] = self.CONFIG[key]
        self.CONFIG = stored_config

        # set min logging interval if not supplied
        if "SERVER_MONITORING_MIN_INTERVAL_MSEC" not in self.CONFIG.keys():
            self.CONFIG["SERVER_MONITORING_MIN_INTERVAL_MSEC"] = 0

        # ensure the log directory exists
        # check the log directory exists.  raises an error if it is not there
        self.logdir = os.path.dirname(self.CONFIG['LOGFILE'])
        if not os.path.exists(self.logdir):
            raise FileNotFoundError("Logdir does not exist at {0}.  You must create this directory".format(self.logdir))
        
        # create a storage directory; check it is writeable.
        # two subdirectories exist for multisequence fasta files (mfa) and 
        # reference compressed directories (rcs)
        self.localcache = os.path.join(self.logdir, 'localcache', self.CONFIG['SERVERNAME'])
        self.mfacache = os.path.join(self.localcache, 'mfa')
        self.rcscache = os.path.join(self.localcache, 'rcs')
        
        os.makedirs(self.localcache, exist_ok=True)
        os.makedirs(self.mfacache, exist_ok=True)
        os.makedirs(self.rcscache, exist_ok=True)

        return self.CONFIG

    def delete_existing_data(self):
        """deletes all existing data from the database
        Note: this routine is in place for unit testing purposes; it should not normally be called in other settings"""
        pm = Persistence()
        self.PERSIST = pm.get_storage_object(
            dbname=self.CONFIG["SERVERNAME"],
            connString=self.CONFIG["FNPERSISTENCE_CONNSTRING"],
            debug=2,
            verbose=False,
        )

        # delete the contents of any catwalk backupdir
        for cw_backup_file in glob.glob(os.path.join(self.localcache, "rc_*.tar")):
            os.unlink(cw_backup_file)

    def _first_run(self, do_not_persist_keys):
        """first run actions.  Stores CONFIG to database, minus any keys in do_not_persist_keys"""

        # create a config dictionary
        config_settings = {}

        # store start time
        config_settings["createTime"] = datetime.datetime.now()

        # store description
        config_settings["description"] = self.CONFIG["DESCRIPTION"]

        # store clustering settings
        self.clustering_settings = self.CONFIG["CLUSTERING"]
        config_settings["clustering_settings"] = self.clustering_settings

        # store precomparer settings
        self.PERSIST.config_store("preComparer", self.CONFIG["PRECOMPARER_PARAMETERS"])

        # load the excluded bases
        excluded = set()
        if self.CONFIG["EXCLUDEFILE"] is not None:
            with open(self.CONFIG["EXCLUDEFILE"], "rt") as f:
                rows = f.readlines()
            for row in rows:
                excluded.add(int(row))

        logging.info("Noted {0} positions to exclude.".format(len(excluded)))

        # load reference
        with open(self.CONFIG["INPUTREF"], "rt") as f:
            for r in SeqIO.parse(f, "fasta"):
                config_settings["reference"] = str(r.seq)

        # persist other config settings.
        for item in self.CONFIG.keys():
            if item not in do_not_persist_keys:
                config_settings[item] = self.CONFIG[item]
        config_settings["excludePositions"] = list(sorted(excluded))

        self.PERSIST.config_store("config", config_settings)

    def _read_config_from_file(
        self, config_fpath: PathLike, required_keys: dict = dict()
    ) -> dict:
        """read a config file

        input:      config_fpath, a path to a configuration file
        returns:    the config file as a dictionary"""

        used_fpath = Path(config_fpath)
        if not used_fpath.exists():
            raise FileNotFoundError(f"Config file {used_fpath} not found")

        with used_fpath.open("r") as f:
            CONFIG = f.read()

        CONFIG = self._enforce_config_object_type(CONFIG)
        self._enforce_key_presence(CONFIG, required_keys)
        CONFIG = self._environment_variables_override_config_defaults(CONFIG)

        return CONFIG

    def _environment_variables_override_config_defaults(
        self, config_like: ConfigLike
    ) -> dict:
        """if environment variables are present containing connection strings, use these in precedence to results in CONFIG"""
        CONFIG = config_like
        # determine whether a FNPERSISTENCE_CONNSTRING environment variable is present,
        # if so, the value of this will take precedence over any values in the config file.
        # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
        if os.environ.get("FNPERSISTENCE_CONNSTRING") is not None:
            CONFIG["FNPERSISTENCE_CONNSTRING"] = os.environ.get(
                "FNPERSISTENCE_CONNSTRING"
            )
            logging.info("Set mongodb connection string  from environment variable")
        else:
            logging.info("Using mongodb connection string from configuration file.")

        # determine whether a FN_SENTRY_URLenvironment variable is present,
        # if so, the value of this will take precedence over any values in the config file.
        # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
        if os.environ.get("FN_SENTRY_URL") is not None:
            CONFIG["SENTRY_URL"] = os.environ.get("FN_SENTRY_URL")
            logging.info("Set Sentry connection string from environment variable")
        else:
            logging.info("Using Sentry connection string from configuration file.")

        # if the CW_BINARY_FILEPATH is present, replace whatever is in the config file
        if os.environ.get("CW_BINARY_FILEPATH") is not None:

            try:
                CONFIG["PRECOMPARER_PARAMETERS"]["catWalk_parameters"][
                    "cw_binary_filepath"
                ] = os.environ.get("CW_BINARY_FILEPATH")
            except KeyError:
                pass  # no key

        return CONFIG

    def validate_server_config(
        self, config_like: ConfigLike, required_keys: dict = dict()
    ) -> None:
        CONFIG = self._enforce_config_object_type(config_like)
        self._enforce_key_presence(CONFIG, required_keys)
