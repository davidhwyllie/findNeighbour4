#!/bin/bash

# script to startup a findNeighbour server.
# expects a configuration file to be provided as a single parameter.


if [ $# -gt 0 ]; then
    echo "findNeighbour3 startup script; using config file $1"
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
nohup pipenv run python3 findNeighbour3-server.py $1 &
echo "Starting dbmanager"
nohup pipenv run python3 findNeighbour3-dbmanager.py $1 &
echo "Starting monitor"
nohup pipenv run python3 findNeighbour3-monitor.py $1 &

echo "Startup complete"

