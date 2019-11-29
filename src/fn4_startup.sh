#!/bin/bash

# script to startup a findNeighbour server and accessory software (monitor, clustering, dbmanager).
# expects a configuration file to be provided as a single parameter.


if [ $# -gt 0 ]; then
    echo "findNeighbour4 startup script; using config file $1"
else
    echo "Your command line contains no arguments.  A config file must be passed"
    exit 1
fi

if [ -f $1 ]; then
    echo "startup script exists."
else
    echo "The startup script was not found."
    exit 1
fi

# set sentry url for bug tracking
FN_SENTRY_URL="https://49ebb508a10f428aaf82f9e1b6b11def@sentry.io/110110"

# startup catwalk server - need to tet whther it is running
#nohup /home/phe.gov.uk/david.wyllie/catwalk/src/cw_server --instance_name test --reference_name ../reference/TB-ref.fasta --reference_filepath ../reference/TB-ref.fasta --mask_name ../reference/TB-exclude-adaptive.txt --mask_filepath ../reference/TB-exclude-adaptive.txt --max_distance 20 > nohup_catwalk.out &


echo "Starting server"
nohup pipenv run python3 findNeighbour4-server.py $1 > nohup_server.out &
echo "Starting dbmanager"
nohup pipenv run python3 findNeighbour4-dbmanager.py $1 > nohup_dbmanager.out &
echo "Starting monitor"
nohup pipenv run python3 findNeighbour4-monitor.py $1 > nohup_monitor.out &
echo "Starting clustering"
nohup pipenv run python3 findNeighbour4-clustering.py $1 > nohup_clustering.out &


echo "Startup complete"

