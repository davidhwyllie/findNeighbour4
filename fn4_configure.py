#!/usr/bin/env python
""" 
Configures the environment in which findneighbour4_server.py operates, which is 
- required before running with gunicorn
- not required if running using the flask internal development server werkzeug.

It can perform the following actions:
- configure the .env file such that findneighbour4_server, run with gunicorn, will load the correct configuration file
- preload catwalk with data.  This is a good idea because otherwise each gunicorn thread will try to do so independently.
- delete data from the findNeighbour4 database.  This is only a good idea if you are testing/benchmarking software. 
- unconfigure the .env file, removing references to any configuration file.  
(Following unconfiguring, running findNeighbour4_server with gunicorn will only be possible in debug mode, with the unittesting config, and using only one worker.)

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

"""

# import libraries
import os
import logging
import argparse
import multiprocessing
from findn.persistence import Persistence
from findNeighbour4_server import findNeighbour4
from findn.common_utils import ConfigManager, EnvWriter

if __name__ == "__main__":
    # command line usage.  Pass the location of a config file as a single argument.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""Configures the environment in which findneighbour4_server.py operates, which is 
- required before running with gunicorn
- not required if running using the flask internal development server werkzeug.

It can perform the following actions:
- configure the .env file such that findneighbour4_server, run with gunicorn, will load the correct configuration file
- preload catwalk with data.  This is a good idea because otherwise each gunicorn thread will try to do so independently.
- delete data from the findNeighbour4 database.  This is only a good idea if you are testing/benchmarking software. 
- unconfigure the .env file, removing references to any configuration file. 
- drop the insert semaphore (a flag which prevents insertion of more than one sample at a time, but which could stay raised in some kinds of server crash) 
(Following unconfiguring, running findNeighbour4_server with gunicorn will only be possible in debug mode, with the unittesting config, and using only one worker.)

Returns a suitable command to run gunicorn 

Example usage: 
============== 
# show command line options 
pipenv run python fn4_configure.py --help  

# normal usage; provide the path to config file in .env file
pipenv run python fn4_configure.py /path/to/config_file.json --set  

""",
    )
    parser.add_argument(
        "path_to_config_file",
        type=str,
        action="store",
        nargs="?",
        help="the path to the configuration file",
        default=None,
    )
    parser.add_argument(
        "--set",
        action="store_true",
        help="set the FN4_SERVER_CONFIG_FILE in .env to path_to_config_file",
    )
    parser.add_argument(
        "--clear",
        action="store_true",
        help="clear (remove) the FN4_SERVER_CONFIG_FILE in .env",
    )
    parser.add_argument(
        "--remove_data_as_part_of_debugging",
        action="store_true",
        help="delete any data in the underlying databases",
    )
    parser.add_argument(
        "--drop_locks",
        action="store_true",
        help="drop the insert and catwalk startup locks (these prevents insertion of more than one sample at a time, or startup of > 1 catwalk at a time, but which could stay raised in some rare kinds of server crash)",
    )
    parser.add_argument(
        "--prepopulate_catwalk",
        action="store_true",
        help="prepopulate catwalk with data from database",
    )
    parser.add_argument(
        "--startup",
        action="store_true",
        help="perform all necessary preparations to launching a gunicorn based server; write startup command to STDOUT",
    )
    parser.add_argument(
        "--shutdown",
        action="store_true",
        help="write shutdown command to STDOUT",
    )
    parser.add_argument(
        "--n_workers",
        action="store",
        type=int,
        default=multiprocessing.cpu_count(),
        help="the number of workers to be used by gunicorn",
    )
    args = parser.parse_args()

    ############################ LOAD CONFIG ######################################
    config_file = args.path_to_config_file
    cfm = ConfigManager(config_file)
    CONFIG = (
        cfm.read_config()
    )  # tests that it can be read.  Will raise a suitable error if it cannot

    pm = Persistence()
    PERSIST = pm.get_storage_object(
        dbname=CONFIG["SERVERNAME"],
        connString=CONFIG["FNPERSISTENCE_CONNSTRING"],
        debug=0,
        verbose=False,
    )

    # read env file present at .env
    ev = EnvWriter()

    logging.info(
        "Read .env file with the following keys: {0}".format(ev.env_vars.keys())
    )
    if args.set or args.startup:
        logging.info("Setting FN4_SERVER_CONFIG_FILE to {0}".format(config_file))
        ev.set_env_var(
            "FN4_SERVER_CONFIG_FILE", "'{0}'".format(args.path_to_config_file)
        )
        ev.save_changes()

    if args.drop_locks or args.startup:
        logging.info("Dropping insert and catwalk startup semaphores")
        PERSIST.unlock(1, force=True)
        PERSIST.unlock(2, force=True)

    if args.prepopulate_catwalk or args.startup:

        fn4 = findNeighbour4(CONFIG, PERSIST)
        fn4.prepopulate_catwalk()

    if args.remove_data_as_part_of_debugging:
        logging.warning("Deleting any data in underlying server databases")
        cfm.delete_existing_data()

    if args.clear:
        logging.info("Clearing FN4_SERVER_CONFIG_FILE env var")
        ev.del_env_var("FN4_SERVER_CONFIG_FILE")
        ev.save_changes()

    # construct command to run server with gunicorn.
    # construct the required global variables
    LISTEN_TO = "127.0.0.1"  # localhost by default
    if "LISTEN_TO" in CONFIG.keys():
        LISTEN_TO = CONFIG["LISTEN_TO"]

    if args.startup:
        error_output_file = os.path.join(
            cfm.logdir,
            "gunicorn_error_logging_{0}_{1}.log".format(
                CONFIG["SERVERNAME"], CONFIG["REST_PORT"]
            ),
        )
        access_output_file = os.path.join(
            cfm.logdir,
            "gunicorn_access_logging_{0}_{1}.log".format(
                CONFIG["SERVERNAME"], CONFIG["REST_PORT"]
            ),
        )
        nohup_output_file = os.path.join(
            cfm.logdir,
            "gunicorn_nohup_logging_{0}_{1}.log".format(
                CONFIG["SERVERNAME"], CONFIG["REST_PORT"]
            ),
        )

        startup_cmd = "nohup pipenv run gunicorn wsgi:app --bind {1}:{2} --log-level info  --workers {0} --error-logfile {3} --access-logfile {4} --timeout 90  > {5} &".format(
            args.n_workers,
            LISTEN_TO,
            CONFIG["REST_PORT"],
            error_output_file,
            access_output_file,
            nohup_output_file,
        )
        logging.info(
            "Configure finished.  Startup command returned to STDOUT {0}".format(
                startup_cmd
            )
        )
        print(startup_cmd)

    if args.shutdown:
        shutdown_cmd = (
            'pkill -f "gunicorn wsgi:app --bind {0}:{1} --log-level info"'.format(
                LISTEN_TO, CONFIG["REST_PORT"]
            )
        )
        logging.info(
            "Configure finished.  Shutdown command returned to STDOUT {0}".format(
                shutdown_cmd
            )
        )
        print(shutdown_cmd)
