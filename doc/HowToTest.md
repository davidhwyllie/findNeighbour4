
Description of how to test and use the system
=============================================

Python version
--------------
This application does not work with python 2.7.  
It has been tested with python 3.5 and 3.7, and with Mongodb v 3.6.1 and 4.0

Operating system
----------------
The code has been tested on Linux (Ubuntu 16.04, 18.04), Windows 7 and Windows 10.  The below instruction include linux-specific commands.

Dependencies
------------
* Python >= 3.5.  On Windows there are MSI installers available.  See below for Linux.  
You may need to install pip3 with: ```sudo apt-get install python3-pip```


Install pipenv, as follows:  
```pip3 install pipenv```  

If you do not have root/admin access, you can install this locally:  
```pip3 install pipenv --user```

Then install dependencies.  The [recommended method](dependencies.md) uses a virtual environment.  

Database backend
----------------
The server requires a mongodb database to work.
[These instructions](mongoinstall.md) describe installation of mongodb.
This server has been tested both with a local mongodb database and with a free cloud instance of mongodb, Mongo Atlas.


Protocol for configuring a clean Linux ubuntu 14.04 instance
-----------------------

Note: **at present this has not been tested with findneighbour3, only findNeighbour2.**
```
sudo apt-get update  
sudo apt-get upgrade  
sudo apt-get install git  
sudo apt-get install python3  
sudo apt-get install -y python3-pip --force-yes  
sudo apt-get install build-essential libssl-dev libffi-dev python3-dev  
sudo apt-get install python3-numpy  
```

Optionally, set a proxy: inform git of the proxy's location, depending whether there is one:
```git config --global http.proxy http://[ip of proxy]```

Then clone repository:
```git clone https://github.com/davidhwyllie/findNeighbour3.git```

After this, please follow the below steps.

Virtual environments
--------------------
It is recommended, but not essential, to use a virtual environment.
This [section](dependencies.md) describes how to set this up.
To run with a virtual environment, preface command with ```pipenv run ..```
e.g.
```pipenv run python3 findNeighbour3-server.py```.

The below commands will run without a virtual environment.   

Start the server
-----------------
To start the server, go to the findNeighbour3 *src* folder and run the command:  
```python3 findNeighbour3-server.py```  

Note: This application doesn't work with python2, so be sure to use python 3.
This will try to start the webserver with a default configuration, in debug mode.  
**Important**: *Debug mode means, amongst other things, that all existing data will be wiped on server restart.  This is good for testing, but not in most other settings.  You need to edit the config file (see below) for your needs.*

If the server fails to start, it's probably due to one of the following:
* mongodb not being operational (a ```pymongo.errors.ServerSelectionTimeOutError``` will indicate this; in Windows, check in *Services* that the service is running; in linux, a command like ```sudo systemctl start mongod``` will be needed), or
* missing dependencies (which it will report).  Install them, then try again.  When it works, terminate the server, and kill any remaining process.

The more general form for starting the server is:
```
nohup python3 findNeighbour3-server.py {configfile.json} &  
```

If {configfile.json} is omitted, then it will use a default config file, config/default_test_config.json.  This is suitable for unit testing, and other kinds of one-off tests. It expects a mongodb running on localhost on the default port.
It is **unsuitable for production**, because:  
1 it runs the flask webserver in debug mode, which is insecure   
2 all data is wiped on server restart.   
A warning is emitted if the server is running with this default configuration.  

Unit tests
----------

```
# you can test the internal classes used by findNeighbour3; all should pass if a mongodb server is operational on local host
python3 -m unittest  NucleicAcid
python3 -m unittest  seqComparer  
python3 -m unittest  clustering  
python3 -m unittest  mongoStore   # requires mongodb server on localhost

More complex testing requires a findNeighbour3 server running.
Unit tests don't start the server. You will need to do.  After this, you can run unit tests.  

# starting a test RESTFUL server
nohup python3 findNeighbour3-server.py &

# And then (e.g. in a different terminal, in windows) launching unit tests as below.
#
# Note: unittesting is changes the data in the server.
# Do not run unittests against a production server.
# In the below configuration, the unittests will run against a
# separate instance of the server used for debugging, called 'fn3_unittesting'

python3 -m unittest findNeighbour3-server 

```
All should pass.
Now kill the webserver

```
ps -x | grep findNeighbour3-server
# kill  servers:
kill -9 <pid>
```

Integration tests
-----------------
see [here](integration.md)

Demonstrations
--------------
see [here](demos.md)

Using the web server
--------------------
You need to start the web server with a sensible configuration, e.g. something like
```nohup python3 findNeighbour3-server.py config/tbproduction.json & ```

The parameter is a json file containing a number of important parameters:
```
INPUTREF:       the path to fasta format reference file.
EXCLUDEFILE:    a file containing the zero-indexed positions in the supplied sequences which should be ignored in all cases.
                Typically, this is because the software generating the mapped fasta file has elected not to call these regions,
                in any samples, e.g. because of difficulty mapping to these regions.
                Such regions can occupy up 5- 20% of the genome and it is important for efficient working of this software
                that these regions are supplied for exclusion on sequence loading.  Not doing so will slow loading, and markedly increase
                memory requirements, but will not alter the results produced.
DEBUGMODE:      Controls operation of the server:

                DEBUGMODE =                                       0       1        2
                Run server                                        Y       N        N
                Run server in debug mode (errors reported)        N       Y        Y
                Create Database if don't exist                    Y       Y        Y
                Delete all data on startup                        N       N        Y

If true, will delete any samples in the backend data store on each run.
SERVERNAME:     the name of the server.  Used as the name of mongodb database which is bound to the server.
FNPERSISTENCE_CONNSTRING: a valid mongodb connection string. if shard keys are set, the 'guid' field is suitable key.  if you don't want to put this is
                in a configuration file, leave the config. file value as "" and  make an environment variable with this name;
                and the server will use the value in that instead.

MAXN_STORAGE:   The maximum number of Ns in the sequence <excluding those defined in > EXCLUDEFILE which should be indexed.
                Other files, e.g. those with all Ns, will be tagged as 'invalid'.  Although a record of their presence in the database
                is kept, they are not compared with other sequences.
MAXN_PROP_DEFAULT: if the proportion not N in the sequence exceeds this, the sample is analysed, otherwise considered invalid.
LOGFILE:        the log file used
LOGLEVEL:		default logging level used by the server.  Valid values are DEBUG INFO WARNING ERROR CRITICAL
SNPCEILING: 	links between guids > this are not stored in the database
GC_ON_RECOMPRESS: if 'recompressing' sequences to a local reference, something the server does automatically, perform
                a full mark-and-sweep gc at this point.  This setting alters memory use and compute time, but not the results obtained.
RECOMPRESS_FREQUENCY: if recompressable records are detected, recompress every RECOMPRESS_FREQ th detection (e.g. 5).
                Trades off compute time with mem usage.  This setting alters memory use and compute time, but not the results obtained.
REPACK_FREQUENCY: concerns how the matrix is stored in mongodb.
                if REPACK_FREQ=0, there will be one document for every non-empty matrix cell.
			             if REPACK_FREQ>0, then if a guid has REPACK_FREQ-1 neighbours, then a 'repack' operation
							         occurs.  This transfers multiple matrix cells into one mongodb document: essentially, part or all of a row
							         will be packed into a single document.  This reduces query times, but the repack operation slows inserts.
							         Repacking doesn't alter the results at all, and could be performed independently of inserts.
CLUSTERING:		a dictionary of parameters used for clustering.  In the below example, there are two different
                clustering settings defined, one named 'SNV12_ignore' and the other 'SNV12_include.
                {'SNV12_ignore' :{'snv_threshold':12, 'mixed_sample_management':'ignore', 'mixture_criterion':'pvalue_1', 'cutoff':0.001},
           'SNV12_include':{'snv_threshold':12, 'mixed_sample_management':'include', 'mixture_criterion':'pvalue_1', 'cutoff':0.001}
 }
                Each setting is defined by four parameters:
                snv_threshold: clusters are formed if samples are <= snv_threshold from each other
                mixed_sample_management: this defines what happens if mixed samples are detected.
                    Suppose there are three samples, A,B and M.  M is a mixture of A and B.
                    A and B are > snv_threshold apart, but their distance to M is zero.
                    If mixed_sample_management is
                    'ignore', one cluster {A,B,M} is returned
                    'include', two clusters {A,M} and {B,M}
                    'exclude', three clusters are returns {A},{B},{C}
                mixture_criterion: sensible values include 'p_value1','p_value2','p_value3' but other output from  seqComparer._msa() is also possible.
                     these p-values arise from three different tests for mixtures.  Please see seqComparer._msa() for details.
                cutoff: samples are regarded as mixed if the mixture_criterion is less than or equal to this value.
```
	
An example CONFIG is below:

```
	{
"DESCRIPTION":"A test server operating in on localhost for unit testing using mapped TB data",
"IP":"127.0.0.1",
"INPUTREF":"../reference/TB-ref.fasta",
"EXCLUDEFILE":"../reference/TB-exclude-adaptive.txt",
"DEBUGMODE":2,
"SERVERNAME":"fn3_unittesting",      
"FNPERSISTENCE_CONNSTRING":"mongodb://localhost",
"MAXN_STORAGE":130000,
"RECOMPRESS_FREQUENCY":5,
"REPACK_FREQUENCY":1,
"GC_ON_RECOMPRESS":1,
"SNPCOMPRESSIONCEILING":250,
"MAXN_PROP_DEFAULT":0.85,
"LOGFILE":"../unittest_tmp/logfile_unittesting.log",
"LOGLEVEL":"DEBUG",	
"SNPCEILING": 20,
"REST_PORT":5000,
"CLUSTERING":{"SNV12_ignore":{"snv_threshold":12,"mixed_sample_management":"ignore","mixture_criterion":"p_value1","cutoff":0.001},
              "SNV12_include":{"snv_threshold":12,"mixed_sample_management":"include","mixture_criterion":"p_value1","cutoff":0.001}
             }
}

```
Notes on altering the CONFIG file

```
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
REPACK_FREQUENCY

```

Database backend
----------------
Mongodb is required.
Provided the findNeighbour server connection has sufficient priviledges, no configuration or pre-creation of databases is needed.


Benchmarking
============
To follow.  The machine used to do the benchmarking was as follows:


Multiple instances of findNeighbour3
----------------------------------------
You can run multiple instances of findNeighbour3 (e.g. multiple different organisms) on the same physical server.
However, you cannot run multiple instances on the same port.
The API is not parameterised by 'instance' or 'organism' etc.
* One port, one server, one config file.
* Use different 'SERVERNAME' settings for each server.  This name becomes the name of the backend mongodb database.
* Be aware that unittesting is destructive.  The database named in the CONFIG file used for unittesting (currently fn3_unittesting) will be recreated.


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
