
Description of how to test and use the system
=============================================

Python version
--------------
This application does not work with python 2.7.  
It has been tested with python 3.5+ with SQLite and postgresql backends.  
If you have two versions of python (e.g. python for 2.7, and python3 for python3) ensure you use the correct version.


Dependencies
------------

* Python 3.4+.   
You may need to install pip3 with: 
```
sudo apt-get install python3-pip
```

Then install the following packages:
* BioPython  
* sqlalchemy  
* psutil  
* queue
* numpy
* pandas
* banyan
* psycopg2 (if using postgres as database)
* pymysql (if using mysql as database)
* flask (if using restful endpoints)

Database backend
----------------
The server will run using SQLite without any additional configuration.
Please use this for initial testing; if you wish to use a more scalable database,
please see the section on configuration files below.


Protocol for configuring a clean Linux ubuntu 14.04 instance
============================================================
Note that this protocol does not use a virtual environment.

```
sudo apt-get update  
sudo apt-get upgrade  
sudo apt-get install git  
sudo apt-get install python3  
sudo apt-get install -y python3-pip --force-yes  
sudo apt-get install build-essential libssl-dev libffi-dev python3-dev  
sudo apt-get install python3-numpy  
sudo apt-get install postgres postgres-dev-all  
sudo apt-get install postgres-client postgres-contrib   

# use of --proxy depends on whether the server is using one
sudo pip3 install sqlalchemy --proxy http://[ip of proxy] 
sudo pip3 install psutil --proxy http://[ip of proxy] 
sudo pip3 install numpy --proxy http://[ip of proxy] 
sudo pip3 install pandas --proxy http://[ip of proxy]
sudo pip3 install banyan --proxy http://[ip of proxy] 
sudo pip3 install psycopg2 --proxy http://[ip of proxy] 
sudo pip3 install biopython --proxy http://[ip of proxy] 
sudo pip3 install flask --proxy http://[ip of proxy]

# optionally inform git of the proxy's location, depending whether there is one
git config --global http.proxy http://[ip of proxy]

# clone repository
git clone https://github.com/davidhwyllie/findNeighbour2.git

# configure postgres; assumes you're logged in as 'user'
sudo -u postgres createuser user
sudo -u postgres createuser ew2  
sudo -u postgres createdb user   # may be unnecessary if you installed postgres, cf. https://stackoverflow.com/questions/17633422/psql-fatal-database-user-does-not-exist

sudo -u postgres psql
```

The following commands are issued via the psql terminal.

```
CREATE ROLE ew2_daemon;

CREATE DATABASE "TB_ew2_edges"
  WITH OWNER = ew2_daemon
       ENCODING = 'UTF8'
       TABLESPACE = pg_default
       LC_COLLATE = 'en_GB.UTF-8'
       LC_CTYPE = 'en_GB.UTF-8'
       CONNECTION LIMIT = -1;

COMMENT ON DATABASE "TB_ew2_edges"
  IS 'Edge list maintained by the EW2 server';
  
CREATE DATABASE "TB_ew2_seqProps"
  WITH OWNER = ew2_daemon
       ENCODING = 'UTF8'
       TABLESPACE = pg_default
       LC_COLLATE = 'en_GB.UTF-8'
       LC_CTYPE = 'en_GB.UTF-8'
       CONNECTION LIMIT = -1;

COMMENT ON DATABASE "TB_ew2_seqProps"
  IS 'Meta-data (apart from edges) related to EW2 server sequences';
  
-- assign permissions to the EW2 daemon.  This assigns permission on all the tables;
-- more restricted permissions are possible
ALTER DEFAULT PRIVILEGES 
    GRANT INSERT, SELECT, UPDATE, DELETE, TRUNCATE, REFERENCES, TRIGGER ON TABLES
    TO ew2_daemon;
    
ALTER USER "ew2_daemon" WITH PASSWORD "ew2";
ALTER ROLE "ew2_daemon" WITH LOGIN;
```


After this, please follow the below steps.


Start the server
-----------------

To start the server, go to the findNeighbour2 *src* folder and run the command:


```
python3 webservice-server.py
```
Note: This application doesn't work with python2, so be sure to use python 3.

This will try to start the webserver with a default configuration using SQLite.  If it fails to start, it's probably due to missing dependencies (which it will report).  Install them, then try again.  When it works, terminate the server, and kill any remaining process.

The more general form for starting the server is:
```
nohup python3 webservice-server.py {configfile.json} &  
```

* If {configfile.json} is not provided, then it will use a default config file, config/default_config.json  
This is suitable for testing. It runs a test instance on localhost using an SQLite backend in ../unittest_Tmp/.
**It is unsuitable for production.  A warning is emitted if the server is running with this default configuration.**  


Unit tests
----------

At the moment, some kinds of unit testing assume a server is running.  Unit tests don't start the server.
You will need to do.  After this, you can run unit tests:

```
# start the xmlrpc server
nohup python3 webservice-server.py &

# starting a test RESTFUL server
nohup python3 webservice-server-rest.py &

# And then (e.g. in a different terminal) launching unit tests with
python3 -m unittest webservice-server-rest
# all should pass

# you can also test the internal classes used by findNeighbour2; all should pass
python3 -m unittest  seqComparer  
python3 -m unittest  ewsnpstore  
python3 -m unittest  FN  

```
All should pass.
Now kill the webserver

```
ps -x | grep webservice-server
# kill both xmlrpc and restful servers:
kill -9 <pid>
```
  
Using the technology without a web server
-----------------------------------------
If you are interested in doing this (which is unlikely)
please see the example in testdrive.py.
This uses all the components of EW2 without the web server interface.

Using the web server
--------------------
You need to start the web server with a sensible configuration, e.g. something like
```
# start the xmlrpc server
nohup python3 webservice-server.py config/tbproduction.json
# start the restful server
nohup python3 webservice-server-rest.py config/tbproduction.json
```

The json file should look something like this:
```
{
"DESCRIPTION":"The production server used for the PHE relatedness test",
"PORT":8185,
"REST_PORT":8186,
"IP":"IP of the machine that run this service",
"INPUTREF":"../reference/reference.fasta",
"PERSISTENCEDIR":"/home/dwyllie/data/relatednesstest/TB_SERVER/persist",
"EXCLUDEFILE":"../reference/TB.txt",
"SNPDIR":"/home/dwyllie/data/relatednesstest/TB_SERVER/snp",
"DEBUGMODE":0,
"SERVERNAME":"TBSNP",      
"EDGEDB_CONNSTRING":"postgresql+psycopg2://ew2:*******@localhost:5432/ew2_edges",
"FNPERSISTENCE_CONNSTRING":"postgresql+psycopg2://ew2:******@localhost:5432/ew2_seqProps",
"MAXN_STORAGE":100000,
"NCOMPRESSIONCUTOFF":100000,
"MAXN_PROP_DEFAULT":0.85,
"LOGFILE":"/home/dwyllie/data/relatednesstest/TB_SERVER/log/logfile.log",
"LOGLEVEL":"INFO",
"SNPCEILING":20
}
```

Notes on the configuration file are as follows:

Tag | Meaning | Example
--- | --- | ---
DESCRIPTION| A user readable description of the server |The production server used for the PHE relatedness test
PORT| The port the xml rpc server is communicating on  |8185
RESTPORT| The port the RESTful server is communicating on  |8186
IP| The IP of the server |127.0.0.1
INPUTREF| The reference sequence to which the sequence has been mapped, in fasta format |../reference/reference.fasta
PERSISTENCEDIR| The directory to which the reference compressed sequences are stored.  This directory is read by the server on startup, to load sequences into RAM. | home/dwyllie/data/relatednesstest/TB_SERVER/persist
EXCLUDEFILE| A text file containing a single column containing the zero-indexed positions of bases in the reference sequence which cannot be reliably called (in the COMPASS Pipeline, these are marked as 'always-N'). However, the base does not have to be 'N' in the provided sequence; if (for example) regions are subsequently found to required 'extra masking' for relatedness studies, these positions just need to be present in the exclusion file.  Although the server will operate without this information, it is important to provide it as these bases are not considered by the algorithm, which makes for *much* smaller 'deltas' relative to the reference, and so reduced memory footprint and faster operation.  The contents of the exclusion file should not be changed if data is present in the server; doing so will result in incorrect comparisons between samples.|../reference/TB.txt
SNPDIR| The directory in which the SNP database is placed, if SQLite is being used.  |/home/dwyllie/data/relatednesstest/TB_SERVER/snp
DEBUGMODE| Whether to operate in production mode (0) or debug mode (1), in which only 500 samples will be stored by the server.    |0
SERVERNAME| A short code name for the server.  Not used programmatically.   |TBSNP     
EDGEDB_CONNSTRING| An sqlAlchemy connection string for connection to the database holding edges.  For sqlite, this all that is needed; the database will be created if it does not exist. For postgres, the database must exist, and the ew2 login must be able to write to it | postgresql+psycopg2://ew2:*******@localhost:5432/ew2_edges [for postgres];  OR   sqlite:///<<DEFAULT>>/findNeighbour.db  [for sqlite]
FNPERSISTENCE_CONNSTRING| n sqlAlchemy connection string for connection to the database holding meta-data on sequences, e.g. GC content etc.  These statistics are calculated automatically by the server, at at present there are no methods for adding arbitrary meta-data, although this would be easy to do.    For sqlite, this all that is needed; the database will be created if it does not exist. For postgres, the database must exist, and the ew2 login must be able to write to it. | postgresql+psycopg2://ew2:******@localhost:5432/ew2_seqProps [for postgres]; sqlite:///<<DEFAULT>>/{0}.db [for sqlite]
MAXN_STORAGE| Sequences with more than this number of 'N's in the sequence will be recorded in EW2, but will be tagged as 'invalid' and reference-compressed representations of them will not be stored in memory. These Ns do not include any 'excluded' positions, as present in EXCLUDEFILE |100000
NCOMPRESSIONCUTOFF| It is recommended that this is set to equal to MAXN_STORAGE.  Sequences with less than NCOMPRESSIONCUTOFF Ns have Ns stored in a different way.  Instead of storing the positions of Ns as single integers e.g. {1,2,3,4,5} is positions 1-5 are Ns, it stores ranges {(1,5)}.  This markedly reduces memory usage in some settings (up to 3-5 fold) but it slows down the server's addition times up to 10x, despite use of in-memory C++ search trees (Banyan).  Therefore, this option is not recommended. |100000
MAXN_PROP_DEFAULT| By default (unless a different value is given to the API), don't report edges if the sequences are < MAXN_PROP_DEFAULT Ns. |0.85
LOGFILE| The location the server logs to |/home/dwyllie/data/relatednesstest/TB_SERVER/log/logfile.log
LOGLEVEL| The logging level used. |INFO
SNPCEILING | The maximum SNP distance to report | 20

Edit this as appropriate.

Database backend
----------------
If you wish to use a database, you will need to create two databases for each server.
* Their names must match those in the configuration string in the config.json
* Their permissions must be set correctly so that the user connecting to the database (in this case, the 'ew2' user, which is part of the 'ew2_' group.)
has relevant permissions.
* SQL configuring two postgres databases is shown below.

```
-- Database: "TB_ew2_edges"
-- DROP DATABASE "TB_ew2_edges";

CREATE DATABASE "TB_ew2_edges"
  WITH OWNER = ew2_daemon
       ENCODING = 'UTF8'
       TABLESPACE = pg_default
       LC_COLLATE = 'en_GB.UTF-8'
       LC_CTYPE = 'en_GB.UTF-8'
       CONNECTION LIMIT = -1;

COMMENT ON DATABASE "TB_ew2_edges"
  IS 'Edge list maintained by the EW2 server';
  
  
-- Database: "TB_ew2_seqProps"
-- DROP DATABASE "TB_ew2_seqProps";

CREATE DATABASE "TB_ew2_seqProps"
  WITH OWNER = ew2_daemon
       ENCODING = 'UTF8'
       TABLESPACE = pg_default
       LC_COLLATE = 'en_GB.UTF-8'
       LC_CTYPE = 'en_GB.UTF-8'
       CONNECTION LIMIT = -1;

COMMENT ON DATABASE "TB_ew2_seqProps"
  IS 'Meta-data (apart from edges) related to EW2 server sequences';
  

-- alter 'me' to your user name if you want to admin the database 
ALTER DEFAULT PRIVILEGES 
    GRANT INSERT, SELECT, UPDATE, DELETE, TRUNCATE, REFERENCES, TRIGGER ON TABLES
    TO me;

-- assign permissions to the EW2 daemon.  This assigns permission on all the tables;
-- more restricted permissions are possible
ALTER DEFAULT PRIVILEGES 
    GRANT INSERT, SELECT, UPDATE, DELETE, TRUNCATE, REFERENCES, TRIGGER ON TABLES
    TO ew2_daemon;

```

Typical performance
-------------------
Benchmarking results are below.[#2]


Server | Processes Used | Dataset | Masking | MaxNStorage | Min Prop N | SNPCeiling | Memory [#3] | Add 1 sample [#1] | Restart server [#2] | On disc storage [#3] | Find all guids [#4] | Detailed comparison (rpt all variant positions, incl. N) [#4] | Find pairwise edges [#5] | Download all edges < 12 snp [#4] | Check server time [#5] | Get all guids & examination time [#5] | Get all guids and annotations [#5]    
--- | --- | --- | --- | --- | --- | --- | --- | --- | ---  | --- | --- | --- | --- | --- | --- | --- | --- | ---   
1 | 1 | TB, n= 15985 | Ambiguous + 329714 'always N' sites | 330k | 0.85 | 20 | 23.5G | 2.23s | 250 s | 2.4G | 0.19s | 0.09s /comparions | < 2msec | 60 s [198,584 edges found ] | 2 ms | 0.3s | 80s   
2 | 1 | NGON, n= 2455 | Ambiguous. No 'always N' sites. | 1000k | 0.7 | 500 | 19.3G | 2.95s | 105s | 2.9G | 0.03s | 0.03s / comparison |  < 2msec | 63s [126048 edges found] | 2msec | 0.05s | 15s   
3 | 1 | SALM, n= 5380 | Ambiguous + 51897 'always N' sites. | 1000k | 0.66 | 20 | 7.4G | 1.77s | 30s  | 0.84s | 0.03s | 0.06s /comparison | < 2msec |  30s [264209 edges found] | 2msec |  0.08s | 22s   


The machine used to do the benchmarking was as follows:


Property | Result
--- | ---
Server | Baremetal server running Ubuntu 16.04 LTS    
RAM | 128G RAM DDR3
Cpu | 2 x Intel Xeon 2.4GHz
Cores used | 1 for server, 1 for database
Database |  Postgres on localhost.
Database disc | Same as storage disc
Storage disc | Storage of fasta and reference compressed data on 12 TB RAID5 device (4 x 4TB discs, SATA, 7200 rpm)
Note #1 | python3 push_samples.py {configfile.json}
Note #2 | Checked by manual server restart
Note #3 | Checked manually
Note #4 | python pull_edges.py {configfile.json}
Note #5 | python pull_data.py {configfile.json}


Services available
==================
These are listed [here](endpoints.md).


Ensuring the services restart with the server
=============================================
(thanks to Hemanth M for this)

* create a file “start_fn2” in /etc/init.d which runs a shell script as below, 
```
#!/bin/sh -e
su -l <entity_under_which_server_runs> -c "sh <directory>/startup_fn2.sh"
```

The shell script, located in <directory> in turn starts the required web services:

```
#!/bin/bash
cd <directory into which project cloned>/src
nohup python3 webservice-server.py <options> &
nohup python3 webservice-server-rest.py <options> &
```

Run the following to change the file permissions of the script and init.d file 
sudo chmod +x /etc/init.d/startup_fn2.sh

Update the run levels to start our script at boot, by running the following command.
sudo update-rc.d /etc/init.d/startup_fn2.sh defaults

Test the changes by restarting the server.
