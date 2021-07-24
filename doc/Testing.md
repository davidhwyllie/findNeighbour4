Testing
=============================================
instructions on how to install dependencies are [here](HowToTest.md)

To run with a virtual environment, which is strongly recommended, preface command with   
```pipenv run ..```  
e.g.  
```pipenv run python3 findNeighbour4_server.py```.

## To start the server

The easiest way to start the server is go to the findNeighbour4 root folder and run the command:  
```./fn4_startup.sh {configfile}```  
where {configfile} is the path to a configuration file. 
This will launch the server, database monitor, server monitor, and clustering systems.  The latter three are recommended, but not essential.  Not all of these are required if you are using an RDBMS backend, but the unnecessary software will sense this and shut down.

However, in a first run setting you can check it starts manually

```pipenv run python3 findNeighbour4_server.py config/default_test_config.json```

*Note: This application doesn't work with python2, so be sure to use python 3; also  this software won't work on Windows.*

This will try to start the webserver with a default configuration which is useful for unit testing, in debug mode, using mongodb.    

**Important**: *Debug mode means, amongst other things, that all existing data in the fn3_unittesting database will be wiped on server restart.   This is good for testing, but not in most other situations.  You need to edit the config file (see below) for your needs.*

If you don't want to use mongodb, you need to [configure an RDBMS](database_credentials.md) first.  TO date, we have only tested Oracle ADW.  You will need to edit the config/default_test_config_rdbms.json file so the FNPERSISTENCE_CONNSTRING key in the json file matches your configuration, see above.

*Note: You cannot use SQLite as a backend for the server, because the server is multithreaded, and SQLite doesn't support connections from more than 1 thread.  Althought the storage engine (findn/rdbmsstore.py) reads / writes to sqlite fine, it isn't suitable for backing a web application.  If you are interested in testing a database backend other than Oracle, please [see here](not_oracle.md).*


If the server fails to start, it's probably due to one of the following:
* mongodb, or the database, not being operational/accessible (a ```pymongo.errors.ServerSelectionTimeOutError``` or similar will indicate this.

 ```sudo systemctl start mongod``` 
 
 may help.

* missing dependencies (which it will report).  Install them, then try again.  When it works, terminate the server, and kill any remaining process.  If you find dependencies we haven't documented, or encounter issues, please let us know.

The more general form for starting the server is:
```
nohup python3 findNeighbour4_server.py configfile.json &  
```
or, if using a virtual environment

```
nohup pipenv run python3 findNeighbour4_server.py configfile.json &  
```  

**Important**  
* If {configfile.json} is omitted, then it will use a default config file, config/default_test_config.json.  This is suitable for unit testing, and other kinds of one-off tests. It expects a mongodb running on localhost on the default port.
It is **unsuitable for production**, because:  
1 it runs the flask webserver in debug mode, which is insecure   
2 all data is wiped on server restart because it is running in debug mode
3 it enables the /restart and /status endpoints, which respectively wipe the database and disclose confidential information about the server config, because it is running in debug mode.   
A warning is emitted if the server is running with this default configuration.  


## Unit tests

We do not recommend running unit tests on a production server with production findNeighbour server instances running.  Unit testing shuts down & starts up servers.  Side effects (shutting down other servers than the ones unit testing) is not expected, but it is better not to take the chance.

Complete testing can be achieved with
```
./run_test.sh
```
*Note: this script starts test servers, which are required for some unittesting, and then runs pytest to run all available tests.  You may have to run chmod +x run_test.sh first.*

If you wish to run unittests individually, you can do 
```
# you can test the internal classes used by findNeighbour4; all should pass if a mongodb server is operational on local host.  Here is an example:

```
pipenv run python3 -m unittest  test/test_seq2dict.py
```

## Demonstrations

see [here](demos.md) and [here](integration.md)
TODO: check all these run

Setting up your own server
==========================
You need to start the web server with a sensible configuration, e.g. something like
```nohup python3 findNeighbour4_server.py config/tbproduction.json & ```

The parameter is a json file containing a number of important parameters.

**Important **  
 You can copy and edit the example config files provided for your own server applications.  The config file tells the findNeighbour and catwalk servers what ports to operation on.  Do not run production servers on ports 5020 or 5021, nor catwalk servers on ports 599*.  Servers listening on these ports are assumed by the unit testing infrastructure to be test servers, and will be killed, and alternative servers started up, if unit testing is carried out.   

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
