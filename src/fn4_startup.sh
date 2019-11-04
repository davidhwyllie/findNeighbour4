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

echo "Starting server"
nohup pipenv run python3 findNeighbour4-server.py $1 > nohup_server.out &
echo "Starting dbmanager"
nohup pipenv run python3 findNeighbour4-dbmanager.py $1 > nohup_dbmanager.out &
echo "Starting monitor"
nohup pipenv run python3 findNeighbour4-monitor.py $1 > nohup_monitor.out &
echo "Starting clustering"
nohup pipenv run python3 findNeighbour4-clustering.py $1 > nohup_clustering.out &


echo "Startup complete"

