#!/bin/bash

# script to shutdown a findNeighbour server and accessory software (monitor, clustering, dbmanager).
# expects a configuration file to be provided as a single parameter.

HELPMESSAGE="Your command line contains no arguments. A config file must be passed as the first parameter.  If you do not want to shutdown catwalk (perhaps because restarting it slow) pass --leave_catwalk_running as a second parameter."
if [ $# -gt 0 ]; then
    echo "findNeighbour4 shutdown script; using config file $1"
else
    echo $HELPMESSAGE
    exit 1
fi

if [ -f $1 ]; then
    echo "config file exists."
else
    echo "A config file was not found at $1."
    echo $HELPMESSAGE
    exit 1
fi

SHUTDOWN_CATWALK="YES"      
if [ $# -eq 2 ]; then
    echo "findNeighbour4 shutdown script has 2 parameters; got option $2"
   if [ "$2" = "--leave_catwalk_running" ]; then
        echo "Will leave catwalk running."
        SHUTDOWN_CATWALK="NO"
   else
        echo "Your command line contains two arguments; the only valid second argument is --leave_catwalk_running.  Instead, got $2"
        exit 1
    fi
fi

echo "Shutdown planned. Targeted processes currently running are:"
echo "--------------------"	
ps -x | grep ".[p]y $1"
ps -x | grep "path_to_config_file $1"
echo "--------------------"

echo  "Stopping database managers"
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

echo "Stopping lockmanager"
pkill -f "findNeighbour4_lockmanager.py --path_to_config_file $1" 

echo "Stopping local store manager"
pkill -f "findNeighbour4_lsmanager.py --path_to_config_file $1" 

echo "Shutdown attempted.  see below: there should be no findNeighbour processes running "
echo "Targeted processes still running are as follows:"
echo "--------------------"
ps -x | grep ".[p]y $1"
ps -x | grep "path_to_config_file $1"
echo "--------------------"

echo "Targeting any gunicorn based servers"

# checksum the config file
MD5CHECKSUM=`md5sum $1 | cut -d' ' -f1`

# launch script file
STOPSCRIPT="fn4server_stop_${MD5CHECKSUM}.sh"

echo "Stopping server  worker processes"
pipenv run python3 fn4_configure.py $1 --shutdown --n_workers 8 > $STOPSCRIPT
chmod +x $STOPSCRIPT

echo "running $STOPSCRIPT"
./$STOPSCRIPT $1 


if [ $SHUTDOWN_CATWALK = "YES" ]; then

    BIND_PORT=`cat $1 | grep -o -P '(?<=bind_port":).{4}.*'`
    BIND_PORT=${BIND_PORT:0:4}
    BIND_PORT_LENGTH=${#BIND_PORT}
    echo "Targeting any CatWalk server operating on the relevant port ${BIND_PORT}; initiating shutdown"

    if [ "$BIND_PORT_LENGTH" -eq 4 ]; then
        echo "Examination of the config file indicates a CatWalk server may be running on port $BIND_PORT"
        echo "--------------------"
        ps -x | grep "$BIND_PORT"
        echo "--------------------"
        echo "Stopping any CatWalk process .."
        pkill -f "$BIND_PORT"
        echo "Shutdown attempted.  see below; there should be no catWalk process running on $BIND_PORT"
        echo "--------------------"
        ps -x | grep "$BIND_PORT"
        echo "--------------------"
    else
        echo "Failed to parse the config file.  Identified a BIND_PORT of $BIND_PORT, but this is not the right length"
    fi  
else
    echo "Left catwalk running, as instructed by --leave_catwalk_running"
fi


echo "Shutdown completed"

