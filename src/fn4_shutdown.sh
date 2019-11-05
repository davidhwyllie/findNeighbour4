#!/bin/bash

# script to shutdown a findNeighbour server and accessory software (monitor, clustering, dbmanager).
# expects a configuration file to be provided as a single parameter.


if [ $# -gt 0 ]; then
    echo "findNeighbour4 shutdown script; using config file $1"
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

echo "Shutdown planned. processes running are:"
ps -x | grep ".py $1"
echo "--------------------"
echo "Stopping server"
pkill -f "findNeighbour4-server.py $1" 
echo  "Stopping dbmanager"
pkill -f "python3 findNeighbour4-dbmanager.py $1" 
echo "Stopping monitor"
pkill -f "python3 findNeighbour4-monitor.py $1" 
echo "Stopping clustering"
pkill -f "python3 findNeighbour4-clustering.py $1" 

echo "Shutdown attempted.  see below: there should be no processes running "
ps -x | grep ".py $1"
echo "--------------------"
