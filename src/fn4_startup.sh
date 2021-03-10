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

NOHUP_LOGGING=1
if [ $# -gt 1 ]; then
    if [ $2=="NO_NOHUP_LOGGING" ]; then
	NOHUP_LOGGING=0
	echo "Logging of STDERR and STDOUT from server processes disabled (recommended in production)"
    else
	echo "Second parameter $2 passed, but this is not understood.  To disable logging of STDERR and STDOUT from server processes  (as recommended in production) use NO_NOHUP_LOGGING as the second parameter after the config filename."
	exit 1
    fi
fi

echo "NOHUP_LOGGING is $NOHUP_LOGGING (1= enabled)"

# checksum the config file
MD5CHECKSUM=`md5sum $1 | cut -d' ' -f1`

# get the output directory from the config file
LOGDIR=`python3 get_log_dir_from_config_file.py $1`


# set sentry url for bug tracking: FN_SENTRY_URL
# not set here. Should be set in a .env file in the pipenv

MONLOG="${LOGDIR}nohup_fn4_monitor_${MD5CHECKSUM}.out"
SRVLOG="${LOGDIR}nohup_fn4_server_${MD5CHECKSUM}.out"
MANLOG="${LOGDIR}nohup_fn4_dbmanager_${MD5CHECKSUM}.out"
CLUSTLOG="${LOGDIR}nohup_fn4_clustering_${MD5CHECKSUM}.out"

if [ $NOHUP_LOGGING -eq 0 ]; then
	MANLOG="/dev/null"
	SRVLOG="/dev/null"
	MONLOG="/dev/null"
	CLUSTLOG="/dev/null"
fi

echo "Starting server"
#nohup pipenv run python3 findNeighbour4_server.py $1 > $SRVLOG &
sleep 5
echo "Starting dbmanager instance 1"
nohup pipenv run python3 findNeighbour4_dbmanager.py --recompress_subset 0123 $1 > $MANLOG &
sleep 5
echo "Starting dbmanager instance 2"
nohup pipenv run python3 findNeighbour4_dbmanager.py --recompress_subset 4567 $1 > $MANLOG &
sleep 5
echo "Starting dbmanager instance 3"
nohup pipenv run python3 findNeighbour4_dbmanager.py --recompress_subset 89ab $1 > $MANLOG &
sleep 5
echo "Starting dbmanager instance 4"
nohup pipenv run python3 findNeighbour4_dbmanager.py --recompress_subset cdef $1 > $MANLOG &
sleep 5

echo "Starting monitor"
#nohup pipenv run python3 findNeighbour4_monitor.py $1 > $MONLOG &
sleep 5
echo "Starting clustering"
#nohup pipenv run python3 findNeighbour4_clustering.py $1 > $CLUSTLOG &
sleep 5

echo "Startup complete.  Output files containing STDOUT and STDERR output are being written to $LOGDIR."
echo "Server processes making their own log in the database. Nohup output is going to the following locations:"
echo $MONLOG
echo $MANLOG
echo $SRVLOG
echo $CLUSTLOG

if [ $NOHUP_LOGGING -eq 0 ]; then
	echo "Either use the parameter NO_NOHUP_LOGGING  (recommended) or (2) arrange log rotation of the nohup output using the linux logrotate command, see: https://support.rackspace.com/how-to/understanding-logrotate-utility/"
fi
	 
exit 0
