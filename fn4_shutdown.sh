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
    echo "config file exists."
else
    echo "The config file was not found."
    exit 1
fi


echo "Shutdown planned. Targeted processes currently running are:"
echo "--------------------"	
ps -x | grep ".[p]y $1"
echo "--------------------"

echo  "Stopping dbmanagers"
pkill -f "findNeighbour4_dbmanager.py --recompress_subset 01 $1" 
pkill -f "findNeighbour4_dbmanager.py --recompress_subset 23 $1" 
pkill -f "findNeighbour4_dbmanager.py --recompress_subset 45 $1" 
pkill -f "findNeighbour4_dbmanager.py --recompress_subset 67 $1" 
pkill -f "findNeighbour4_dbmanager.py --recompress_subset 89 $1" 
pkill -f "findNeighbour4_dbmanager.py --recompress_subset ab $1" 
pkill -f "findNeighbour4_dbmanager.py --recompress_subset cd $1" 
pkill -f "findNeighbour4_dbmanager.py --recompress_subset ef $1" 

echo "Stopping server"
pkill -f "findNeighbour4_server.py $1" 

echo "Stopping monitor"
pkill -f "findNeighbour4_monitor.py $1" 
echo "Stopping clustering"
pkill -f "findNeighbour4_clustering.py $1" 

echo "Shutdown attempted.  see below: there should be no findNeighbour processes running "
echo "Targeted processes still running are as follows:"
echo "--------------------"
ps -x | grep ".[p]y $1"
echo "--------------------"

echo "Targeting any CatWalk server operating on the relevant port"
BIND_PORT=`cat $1 | grep -o -P '(?<=bind_port":).*(?=, )'`
BIND_PORT_LENGTH=${#BIND_PORT}
if [ "$BIND_PORT_LENGTH" -eq 4 ]
then
   echo "Examination of the config file indicates a CatWalk server may be running on port $BIND_PORT"
   echo "--------------------"
   ps -x | grep "bind_port $BIND_PORT"
   echo "--------------------"
   echo "Stopping any CatWalk process .."
   pkill -f "bind_port $BIND_PORT"
   echo "Shutdown attempted.  see below; there should be no catWalk process running on $BIND_PORT"
   echo "--------------------"
   ps -x | grep "bind_port $BIND_PORT"
   echo "--------------------"
fi
echo "Shutdown completed"

