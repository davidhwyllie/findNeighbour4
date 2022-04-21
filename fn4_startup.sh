#!/bin/bash

# script to startup a findNeighbour server and accessory software (monitor, clustering, dbmanager).
# expects a configuration file to be provided as a single parameter.

#A component of the findNeighbour4 system for bacterial relatedness monitoring
#Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
#repo: https://github.com/davidhwyllie/findNeighbour4

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as published
#by the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.  See see <https://www.gnu.org/licenses/>.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

set -e      # terminate if errors

# set version
echo "Setting software version from git repo."
VERSION=`python3 setup.py --version`
rm version.py -f
touch version.py
echo "version = '$VERSION'" > version.py

# start server
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

# launch script file
LAUNCHSCRIPT="fn4server_launch_${MD5CHECKSUM}.sh"

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

echo "Server processes making their own log in the database. Nohup output is going to the following locations:"
echo $MONLOG
echo $MANLOG
echo $SRVLOG
echo $CLUSTLOG

echo "Starting server with prespecified worker processes (remove --n_workers X from fn4_startup.sh to autopick number of workers)"
pipenv run python3 fn4_configure.py $1 --startup --n_workers 8 > $LAUNCHSCRIPT
chmod +x $LAUNCHSCRIPT

echo "running $LAUNCHSCRIPT to start gunicorn based server"
echo " -------------------------------------- running launch script to start gunicorn -----------------------------------"
echo $LAUNCHSCRIPT
cat $LAUNCHSCRIPT
echo " -------------------------------------- --------------------------------------- -----------------------------------"

## this process launches multiple catwalks -it's a bug
./$LAUNCHSCRIPT $1 > $SRVLOG 


sleep 0.5
echo "Starting lockmanager"
nohup pipenv run python3 findNeighbour4_lockmanager.py --path_to_config_file $1 --max_run_time 90 > $MANLOG &

sleep 0.5
echo "Starting dbmanager instance 1"
nohup pipenv run python3 findNeighbour4_dbmanager.py --recompress_subset 01 $1 > $MANLOG &
sleep 0.5
echo "Starting dbmanager instance 2"
nohup pipenv run python3 findNeighbour4_dbmanager.py --recompress_subset 23 $1 > $MANLOG &
sleep 0.5
echo "Starting dbmanager instance 3"
nohup pipenv run python3 findNeighbour4_dbmanager.py --recompress_subset 45 $1 > $MANLOG &
sleep 0.5
echo "Starting dbmanager instance 4"
nohup pipenv run python3 findNeighbour4_dbmanager.py --recompress_subset 67 $1 > $MANLOG &
sleep 0.5
echo "Starting dbmanager instance 5"
nohup pipenv run python3 findNeighbour4_dbmanager.py --recompress_subset 89 $1 > $MANLOG &
sleep 0.5
echo "Starting dbmanager instance 6"
nohup pipenv run python3 findNeighbour4_dbmanager.py --recompress_subset ab $1 > $MANLOG &
sleep 0.5
echo "Starting dbmanager instance 7"
nohup pipenv run python3 findNeighbour4_dbmanager.py --recompress_subset cd $1 > $MANLOG &
sleep 0.5
echo "Starting dbmanager instance 8"
nohup pipenv run python3 findNeighbour4_dbmanager.py --recompress_subset ef $1 > $MANLOG &
sleep 0.5

echo "Starting monitor [disabled until issue #71 is resolved]"
#nohup pipenv run python3 findNeighbour4_monitor.py $1 > $MONLOG &
#sleep 0.5

echo "Starting clustering"
nohup pipenv run python3 findNeighbour4_clustering.py $1 > $CLUSTLOG &
sleep 0.5

echo "Starting localstore manager"
nohup pipenv run python3 findNeighbour4_lsmanager.py --path_to_config_file $1  > $MANLOG &

echo "Startup complete.  Output files containing STDOUT and STDERR output are being written to $LOGDIR."
echo "Server processes making their own log in the database. Nohup output is going to the following locations:"
echo $MONLOG
echo $MANLOG
echo $SRVLOG
echo $CLUSTLOG

if [ $NOHUP_LOGGING -eq 0 ]; then
	echo "Either use the parameter NO_NOHUP_LOGGING  (recommended) or (2) arrange log rotation of the nohup output using the linux logrotate command, see: https://support.rackspace.com/how-to/understanding-logrotate-utility/"
fi

rm $LAUNCHSCRIPT    # temporary file	 
exit 0
