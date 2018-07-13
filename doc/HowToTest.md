
Description of how to test and use the system
=============================================

Python version
--------------
This application does not work with python 2.7.  
It has been tested with python 3.5 and 3.7, and with Mongodb v 3.6.1 and 4.0

Dependencies
------------

* Python 3.5+.   
You may need to install pip3 with: 
```
sudo apt-get install python3-pip
```

Then install the following packages:
 requests  
 logging  
 hashlib  
 queue  
 threading  
 pymongo  
 pandas as pd  
 numpy as np  
 flask  
 psutil    
 BioPython  
 
Database backend
----------------
The server requires a mongodb database to work.
This server has been tested both with a local mongodb database and with a free cloud instance of mongodb, Mongo Atlas.


Protocol for configuring a clean Linux ubuntu 14.04 instance
============================================================
Note that this protocol does not use a virtual environment.
Note: this has not been tested with findneighbour3, only findNeighbour2.
```
sudo apt-get update  
sudo apt-get upgrade  
sudo apt-get install git  
sudo apt-get install python3  
sudo apt-get install -y python3-pip --force-yes  
sudo apt-get install build-essential libssl-dev libffi-dev python3-dev  
sudo apt-get install python3-numpy  

# optionally inform git of the proxy's location, depending whether there is one
git config --global http.proxy http://[ip of proxy]

# clone repository
git clone https://github.com/davidhwyllie/findNeighbour3.git


After this, please follow the below steps.


Start the server
-----------------

To start the server, go to the findNeighbour3 *src* folder and run the command:


```
python3 findNeighbour3-server.py
```
Note: This application doesn't work with python2, so be sure to use python 3.
This will try to start the webserver with a default configuration, in debug mode.
** Debug mode means, amongst other things, that all existing data will be wiped on server restart.  This is good for unittesting, but not in many other settings.  You need to edit the config file (see below) for your needs.**.

If it fails to start, it's probably due to missing dependencies (which it will report).  Install them, then try again.  When it works, terminate the server, and kill any remaining process.

The more general form for starting the server is:
```
nohup python3 findNeighbour3-server.py {configfile.json} &  
```

* If {configfile.json} is not provided, then it will use a default config file, config/default_config.json  
This is suitable for testing. It expects a mongodb running on localhost on the default port.
**It is unsuitable for production.  A warning is emitted if the server is running with this default configuration.**  


Unit tests
----------

At the moment, some kinds of unit testing assume a server is running.  Unit tests don't start the server.
You will need to do.  After this, you can run unit tests:

```

# starting a test RESTFUL server
nohup python3 findNeighbour3-server.py &

# And then (e.g. in a different terminal) launching unit tests with
python3 -m unittest findNeighbour3-server
# all should pass

# you can also test the internal classes used by findNeighbour2; all should pass
python3 -m unittest  seqComparer  
python3 -m unittest  clustering  
python3 -m unittest  mongoStore  

```
All should pass.
Now kill the webserver

```
ps -x | grep findNeighbour3-server
# kill  servers:
kill -9 <pid>
```
  
Using the web server
--------------------
You need to start the web server with a sensible configuration, e.g. something like
```
# start the restful server
nohup python3 findNeighbour3-server.py config/tbproduction.json &
```

									
		An example CONFIG is below:
		
		{			
		"DESCRIPTION":"A test server operating in ../unittest_tmp, only suitable for testing",
		"IP":"127.0.0.1",
		"INPUTREF":"../reference/TB-ref.fasta",
		"EXCLUDEFILE":"../reference/TB-exclude.txt",
		"DEBUGMODE":0,
		"SERVERNAME":"TBSNP",
		"FNPERSISTENCE_CONNSTRING":"mongodb://127.0.0.1",
		"MAXN_STORAGE":100000,
		"SNPCOMPRESSIONCEILING":250,
		"MAXN_PROP_DEFAULT":0.70,
		"LOGFILE":"../unittest_tmp/logfile.log",
		"LOGLEVEL":"INFO",
		"SNPCEILING": 20,
		"GC_ON_RECOMPRESS":1,
		"RECOMPRESS_FREQUENCY":5,
		"CLUSTERING":{'SNV12_ignore' :{'snv_threshold':12, 'mixed_sample_management':'ignore'},
		              'SNV12_include':{'snv_threshold':12, 'mixed_sample_management':'include'}
					  }
		}
    
		CONFIG contains Configuration parameters relevant to the reference based compression system which lies
		at the core of the server.  More details on these are below.
		
		  INPUTREF:       the path to fasta format reference file.
		  EXCLUDEFILE:    a file containing the zero-indexed positions in the supplied sequences which should be ignored in all cases.
			           			Typically, this is because the software generating the mapped fasta file has elected not to call these regions,
                            in any samples, e.g. because of difficulty mapping to these regions.
                            Such regions can occupy up 5- 20% of the genome and it is important for efficient working of this software
                            that these regions are supplied for exclusion on sequence loading.  Not doing so will slow loading, and markedly increase
                            memory requirements, but will not alter the results produced.
      DEBUGMODE:      False by default.  If true, will delete any samples in the backend data store on each run.
      SERVERNAME:     the name of the server (used for display purposes only)
			FNPERSISTENCE_CONNSTRING: a valid mongodb connection string. if shard keys are set, the 'guid' field is suitable key.
      MAXN_STORAGE:   The maximum number of Ns in the sequence <excluding those defined in > EXCLUDEFILE which should be indexed.
                            Other files, e.g. those with all Ns, will be tagged as 'invalid'.  Although a record of their presence in the database
                            is kept, they are not compared with other sequences.
			MAXN_PROP_DEFAULT: if the proportion not N in the sequence exceeds this, the sample is analysed, otherwise considered invalid.
			LOGFILE:        the log file used
			LOGLEVEL:		default logging level used by the server.  Valid values are DEBUG INFO WARNING ERROR CRITICAL
			SNPCEILING: 	links between guids > this are not stored in the database
			GC_ON_RECOMPRESS: if 'recompressing' sequences to a local reference, something the server does automatically, perform
			                a full mark-and-sweep gc at this point.  This setting alters memory use and compute time, but not the results obtained.
			RECOMPRESS_FREQ: if recompressable records are detected, recompress every RECOMPRESS_FREQ th detection (e.g. 5).
							Trades off compute time with mem usage.  This setting alters memory use and compute time, but not the results obtained.
			CLUSTERING:		a dictionary of parameters used for clustering.  In the below example, there are two different
							clustering settings defined, one named 'SNV12_ignore' and the other 'SNV12_include.
							{'SNV12_ignore' :{'snv_threshold':12, 'mixed_sample_management':'ignore'},
							'SNV12_include':{'snv_threshold':12, 'mixed_sample_management':'include'}
							}
							Each setting is defined by two parameters:
							snv_threshold: clusters are formed if samples are <= snv_threshold from each other
							mixed_sample_management: this defines what happens if mixed samples are detected.
								Suppose there are three samples, A,B and M.  M is a mixture of A and B.
								A and B are > snv_threshold apart, but their distance to M is zero.
								If mixed_sample_management is
								'ignore', one cluster {A,B,M} is returned
								'include', two clusters {A,M} and {B,M}
								'exclude', three clusters are returns {A},{B},{C}
		

		Some of these settings are read when the server is first-run, stored in a database, and the server will not
		change the settings on re-start even if the config file is changed.  Examples are:
		SNPCEILING
		MAXN_PROP_DEFAULT
		EXCLUDEFILE
		INPUTREF
		CLUSTERING
		These settings cannot be changed because they alter the way that the data is stored; if you want to change
		the settings, the data will have to be re-loaded. 
		
		However, most other settings can be changed and will take effect on server restart.  These include:
		server location
		IP
		SERVERNAME
		REST_PORT
		
		internal logging	
		LOGFILE
		LOGLEVEL
		
		where the database connection binds to
		FNPERSISTENCE_CONNSTRING
		
		related to internal server memory management:
		GC_ON_RECOMPRESS
		RECOMPRESS_FREQUENCY
		SNPCOMPRESSIONCEILING

Edit the config file as appropriate.

Database backend
----------------
Mongodb is required.
Provided the findNeighbour server connection has sufficient priviledges, no configuration or pre-creation of databases is needed.


Benchmarking
============
To follow.  The machine used to do the benchmarking was as follows:


Services available
==================
These are listed [here](rest-routes.md).


Ensuring the services restart with the server
=============================================
(thanks to Hemanth M for this)

* create a file “start_fn3” in /etc/init.d which runs a shell script as below, 
```
#!/bin/sh -e
su -l <entity_under_which_server_runs> -c "sh <directory>/startup_fn3.sh"
```

The shell script, located in <directory> in turn starts the required web services:

```
#!/bin/bash
cd <directory into which project cloned>/src
nohup python3 findNeighbour3-server.py <options> &

```

Run the following to change the file permissions of the script and init.d file 
sudo chmod +x /etc/init.d/startup_fn3.sh

Update the run levels to start our script at boot, by running the following command.
sudo update-rc.d /etc/init.d/startup_fn3.sh defaults

Test the changes by restarting the server.
