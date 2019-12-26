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

echo "NOHUP_LOGGING $NOHUP_LOGGING"

# checksum the config file
MD5CHECKSUM=`md5sum $1 | cut -d' ' -f1`

# get the output directory from the config file
LOGDIR=`python3 get_log_dir_from_config_file.py $1`


# set sentry url for bug tracking
FN_SENTRY_URL="https://49ebb508a10f428aaf82f9e1b6b11def@sentry.io/110110"

if [ $NOHUP_LOGGING -eq 1 ]; then
	echo "Starting server"
	nohup pipenv run python3 findNeighbour4-server.py $1 > ${LOGDIR}nohup_fn4_server_${MD5CHECKSUM}.out &
	echo "Starting dbmanager"
	nohup pipenv run python3 findNeighbour4-dbmanager.py $1 > ${LOGDIR}nohup_fn4_dbmanager_${MD5CHECKSUM}.out &
	echo "Starting monitor"
	nohup pipenv run python3 findNeighbour4-monitor.py $1 > ${LOGDIR}nohup_fn4_monitor_${MD5CHECKSUM}.out &
	echo "Starting clustering"
	nohup pipenv run python3 findNeighbour4-clustering.py $1 > ${LOGDIR}nohup_fn4_clustering_${MD5CHECKSUM}.out &
        sleep 1
	echo "Startup complete.  Output files containing STDOUT and STDERR output are being written to $LOGDIR."
	echo "Server processes make their own logs; retaining the STDOUT and STDERR of the server processes should not be necessary, although we do so currently for debugging purposes.  This will yield very large log files in production.  In a production system, you should either (i) use .fn4_startup.sh {LOGFILE} NO_NOHUP_LOGGING  (recommended) or (2) arrange log rotation of the nohup output using the linux logrotate command, see: https://support.rackspace.com/how-to/understanding-logrotate-utility/"

else

	echo "Starting server"
	nohup pipenv run python3 findNeighbour4-server.py $1 > /dev/null &
	echo "Starting dbmanager"
	nohup pipenv run python3 findNeighbour4-dbmanager.py $1 > /dev/null &
	echo "Starting monitor"
	nohup pipenv run python3 findNeighbour4-monitor.py $1 > /dev/null &
	echo "Starting clustering"
	nohup pipenv run python3 findNeighbour4-clustering.py $1 > /dev/null &
        sleep 1
	echo "Startup complete.  Autorotating server log files are being written to $LOGDIR. Nohup output is not retained.  This is the recommended production arrangement, as large log files do not have to be managed external to findNeighbour4."

fi
	 
exit 0
