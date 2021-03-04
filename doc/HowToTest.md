
Description of how to test and use the system
=============================================

Python version
--------------
This application does not work with python 2.7.  
It has been tested with python 3.5 and 3.7, and with Mongodb v 3.6.1 and 4.0

Operating system
----------------
The code has been tested on Linux (Ubuntu 16.04, 18.04).  
It is not expected to run on Windows, as it uses unix-specific concepts and libraries.
The below instruction include linux-specific commands.

Dependencies
------------
* Python >= 3.5.  
You may need to install pip3 with: ```sudo apt-get install python3-pip```

Note, Ubuntu 18 users, there is a [known issue](https://github.com/pypa/pipenv/issues/2122) and you may need to unisntall pip first, e.g.  
```sudo python3 -m pip uninstall pip && sudo apt install python3-pip --reinstall```


Database backend
----------------
The server requires a mongodb database to work.  Tested with v. 4.02 and 4.2.
[These instructions](mongoinstall.md) describe installation of mongodb.
This server has been tested both with a local mongodb database and with a free cloud instance of mongodb, Mongo Atlas.

Clone git repository
-----------------------
Optionally, set a proxy: inform git of the proxy's location, depending whether there is one:
```git config --global http.proxy http://[ip of proxy]```

Clone repository:
```git clone https://github.com/davidhwyllie/findNeighbour4.git```

After this, please follow the below steps.

Virtual environments and dependencies
-------------------------------------
It is strongly recommended, but not essential, to use a virtual environment.
A pipenv Pipfile is provided which specifies dependencies.  See also [section](dependencies.md).

In the root of your directory, you will need to create a .env file.  This sets environment variables required for the code to run.
Here is an example:
```
FN_SENTRY_URL="https://********************.ingest.sentry.io/*****"
CW_BINARY_FILEPATH="/home/phe.gov.uk/david.wyllie/catwalk/cw_server"
IQTREE="/home/phe.gov.uk/david.wyllie/software/iqtree-2.1.2-Linux/bin/iqtree2"
```
All are optional, but the CW_BINARY_FILEPATH is required if (as is strongly recommended) you are using the catwalk comparison engine as part of findNeighbour4.
FN_SENTRY_URL is an optional url for the sentry.io (error logging) service.  If present and Sentry.io (small free or paid service) configured, error logging will be collated there.  
This is very useful for collating  & debugging server side errors.   If considering this, be aware that if identifiable data is in the server, errors trapped may be sent to the Sentry.io service.  
IQTREE is an optional path to the IQTREE executable, a tree drawing software.


To run with a virtual environment, preface command with ```pipenv run ..```
e.g.
```pipenv run python3 findNeighbour4_server.py```.

The below commands will run without a virtual environment if relevant packages are installed at a machine level, 
and relevant environment variables present..   

To start the server
-------------------
The easiest way to start the server is go to the findNeighbour4 *src* folder and run the command:  
```./fn4_startup.sh {configfile}```  where {configfile} is the path to a configuration file. 
This will launch the server, database monitor, server monitor, and clustering systems.  The latter three are recommended, but not essential.

However, in a first run setting you can check it starts manually
```pipenv run python3 findNeighbour4_server.py```

Note: This application doesn't work with python2, so be sure to use python 3.  It also won't work on Windows, because the clustering system uses linux-specific packages.  These packages (networkit) are supposed to be compatible with WSL, but we have not tested this.
This will try to start the webserver with a default configuration which is useful for unit testing, in debug mode.  **Important**: *Debug mode means, amongst other things, that all existing data in the fn3_unittesting database will be wiped on server restart.   This is good for testing, but not in most other settings.  You need to edit the config file (see below) for your needs.*


If the server fails to start, it's probably due to one of the following:
* mongodb not being operational (a ```pymongo.errors.ServerSelectionTimeOutError``` will indicate this; in Windows, check in *Services* that the service is running; in linux, a command like ```sudo systemctl start mongod``` will be needed), or
* missing dependencies (which it will report).  Install them, then try again.  When it works, terminate the server, and kill any remaining process.

The more general form for starting the server is:
```
nohup python3 findNeighbour4_server.py {configfile.json} &  
```
or, if using a virtual environment

```
nohup pipenv run python3 findNeighbour4_server.py {configfile.json} &  
```  

If {configfile.json} is omitted, then it will use a default config file, config/default_test_config.json.  This is suitable for unit testing, and other kinds of one-off tests. It expects a mongodb running on localhost on the default port.
It is **unsuitable for production**, because:  
1 it runs the flask webserver in debug mode, which is insecure   
2 all data is wiped on server restart because it is running in debug mode
3 it enables the /restart and /status endpoints, which respectively wipe the database and disclose confidential information about the server config, because it is running in debug mode.   
A warning is emitted if the server is running with this default configuration.  

Unit tests
----------

Do not run unit tests on a production server with findNeighbour server instances running.
Side effects are not expected and have not been observed in the most recent versions, but it is better not to take the chance.

Complete testing can be achieved with
```
./unittest_all.sh
```
however a stepwise approach is recommended after initial installation.
Automated testing of everything except the server itself can be performed with
```
pipenv run python3 unittest_core.py
```
Failures initially are likely due to missing dependencies.

If you wish to run unittests individually, you can do 
```
# you can test the internal classes used by findNeighbour4; all should pass if a mongodb server is operational on local host
python3 -m unittest  NucleicAcid
python3 -m unittest  mongoStore   # requires mongodb server on localhost
```

etc or with a virtual environment, do 
```
pipenv run python3 -m unittest  NucleicAcid seqComparer clustering mongoStore
```


More complex testing requires a findNeighbour4 server running.  Note that unit tests don't start the server. You will need to do.  After this, you can run unit tests.  

```
# starting a test RESTFUL server
nohup python3 findNeighbour4_server.py &          # if run in a shell without nohup, this will also work, but you will need to run the tests in a different shell
```
or if using a virtual environment 
```
# starting a test RESTFUL server
nohup pipenv run python3 findNeighbour4_server.py &
```  

And then (e.g. in a different termina) launching unit tests as below.
Note: unittesting is changes the data in the server.
Do not run unittests against a production server, although side effects from testing are not expected; unittest use a separate database and separate port.
In the below configuration, the unittests will run against a
separate instance of the server used for debugging, called 'fn3_unittesting'
 
```
# just tests the python client
pipenv run python3 -m unittest fn4client  
# test the server
pipenv run python3 -m unittest findNeighbour4_server 
```

All tests should pass.
Now kill the webserver

```
ps -x | grep findNeighbour4_server
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
```nohup python3 findNeighbour4_server.py config/tbproduction.json & ```

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
PRECOMPARER_PARAMETERS:
               describe how the Catwalk server is to be operated, see example below.
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
"DESCRIPTION":"PHE fn4 covid test server",
"IP":"127.0.0.1",
"INPUTREF":"../reference/nc_045512.fasta",
"EXCLUDEFILE":"../reference/covid-exclude.txt",
"DEBUGMODE":0,
"SERVERNAME":"PHE_covid_5",
"FNPERSISTENCE_CONNSTRING":"mongodb://localhost",
"MAXN_STORAGE":20000,
"MAXN_PROP_DEFAULT":0.70,
"LISTEN_TO":"0.0.0.0",
"LOGFILE":"/srv/data/mixfiles/log/phe_covid.log",
"LOGLEVEL":"INFO",
"SNPCEILING": 5,
"REST_PORT":5023,
"PRECOMPARER_PARAMETERS":{"selection_cutoff":5,"uncertain_base":"N_or_M",
"over_selection_cutoff_ignore_factor":1,
"catWalk_parameters":{"cw_binary_filepath":"","reference_name":"covid5snp",
"reference_filepath":"../reference/nc_045512.fasta","mask_filepath":"../reference/covid-exclude.txt", "bind_port":5024, "bind_host":"localhost"}},

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

```

Database backend
----------------
Mongodb is required.
Provided the findNeighbour server connection has sufficient priviledges, no configuration or pre-creation of databases is needed.


Benchmarking
============
To follow.  


Multiple instances of findNeighbour4
----------------------------------------
You can run multiple instances of findNeighbour4 (e.g. multiple different organisms) on the same physical server.
However, you cannot run multiple instances on the same port.
The API is not parameterised by 'instance' or 'organism' etc.
* One port, one server, one config file.
* Use different 'SERVERNAME' settings for each server.  This name becomes the name of the backend mongodb database.
* Be aware that unittesting is potentially destructive if done with the wrong settings.  The database named in the CONFIG file used for unittesting (currently fn3_unittesting) will be recreated.  Use the default CONFIG file for unittesting.


Services available
==================
These are listed [here](rest-routes.md).


Ensuring the services restart with the server
=============================================
(thanks to Hemanth M for this)

* create a file “fn4_start.sh” in /etc/init.d which runs a shell script as below, 
```
#!/bin/sh -e
su -l <entity_under_which_server_runs> -c "sh <directory>/fn4_startup.sh {configfile.json} {options}"
```

The shell script, located in <directory> in turn starts the required web services:

Note: set file permissions of the script and init.d file, e.g.
sudo chmod +x /etc/init.d/fn4_start.sh

Update the run levels to start our script at boot, by running the following command.
sudo update-rc.d /etc/init.d/fn4_start.sh defaults

Test the changes by restarting the server.
