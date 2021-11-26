Installing dependencies
=============================================
if you've done this, you can go straght to [testing](Testing.md)

Python version
--------------
Requires python 3.8 + on Linux

Operating system
----------------
The code has been tested on Linux (Ubuntu 20.04).  The server uses a small number of linux specific packages.  It is not expected to work on Windows. 
The below instructions include linux-specific commands.

Dependencies
------------
* Python >= 3.8.  
You may need to install pip3 & pipenv with: 
```
sudo apt-get install python3-pip
pip3 install pipenv
```

Database backend
----------------
You can either use 
* Mongo 4.4
* A relational database; only Oracle ADW has been tested to date.

### Mongo
* Mongo 4.4 (note: tested at very large scale only on this version)
[These instructions](mongoinstall.md) describe installation of mongodb.
This server has been tested both with a local mongodb database and with a free cloud instance of mongodb, Mongo Atlas.

### Oracle ADW
* For details of how to set up an Oracle Autonomous datawarehouse, please see Oracle's cloud documentation.
In addition
* You need to install [dependencies]((https://cx-oracle.readthedocs.io/en/latest/user_guide/initialization.html)) of the python Oracle_cx (database connection) module. 
* You need to set LD_LIBRARY_PATH in the python virtual environment, see below.
See [details](configuring_oracle_connections.md)  

Later
* You need to configure a user for the database correctly
* You need to provide credentials in a configuration file used by FindNeighbour4.
Read on, and also [here](configuring_oracle_connections.md) how to do these steps.  

Running multiple findneighbour4 servers on the same machine
-----------------------------------------------------------
You can run multiple findneighbour4 servers on the same physical or virtual machine.  To do so, you need to arrange the following:
* Each server should be started from its own virtual environment.  Do not use try to use a common virtual environment for all your servers; each findneighbour4 instance writes instance-specific information to a .env file, and there can only be one such in each virtual environment.
* Ensure that the ports (both for the server and catwalk) differ between instances.  You cannot run two findneighbour servers on the same port.  This is configured in a config/ file (see below)

Clone git repository
---------------------
Optionally, set a proxy: inform git of the proxy's location, depending whether there is one:
```
git config --global http.proxy http://[ip of proxy]
```

Clone repository:
```
git clone https://github.com/davidhwyllie/findNeighbour4.git
```

After this, please follow the below steps.
Set version from git repository
------------
```
VERSION=`python3 setup.py --version`
rm version.py -f
touch version.py
echo "version = '$VERSION'" > version.py
```


Virtual environments and dependencies
-------------------------------------
It is strongly recommended, but not essential, to use a virtual environment.
A pipenv Pipfile is provided which specifies dependencies.  See also [here](dependencies.md).

```
cd /mydir/fn4       # or whatever you've installed into
pipenv install -e . --skip-lock   # put fn4 packages in virtualenv; this is essential for effective imports

```

Catwalk
--------
Catwalk is an external component which is required by findneighbour4.  
Installation Example:

```
mkdir external_software

cd external_software
wget https://nim-lang.org/download/nim-1.4.8-linux_x64.tar.xz
tar -xf nim-1.4.8-linux_x64.tar.xz

NIMDIR="`pwd`/nim-1.4.8/bin"

# Appending to PATH
export PATH="$NIMDIR:$PATH"

# check nimble is accessible on path
echo "Nimble version:"
nimble --version

# clone catwalk relatedness engine
git clone https://github.com/dvolk/catwalk.git

cd catwalk
nimble -y build -d:release -d:danger -d:no_serialisation
# add path to executable to .env file 
echo CW_BINARY_FILEPATH=\"`pwd`/cw_server\" > ../../.env

```

.env file
---------

In the root of your directory, you will need to create a .env file.  This sets environment variables required for the code to run.   
Because you are using a virtual environment, it has its own environment variables.

Here is an example if using an Oracle database:
```
CW_BINARY_FILEPATH="/home/phe.gov.uk/david.wyllie/catwalk/cw_server"
IQTREE="/home/phe.gov.uk/david.wyllie/software/iqtree-2.1.2-Linux/bin/iqtree2"
FASTTREE_DIR="/data/software/fasttree"
LD_LIBRARY_PATH="/data/software/instantclient_21_1"
PCA_CONNECTION_CONFIG_FILE='/data/credentials/dbconfig.json'
DB_CONNECTION_CONFIG_FILE='/data/credentials/dbconfig.json'
```

Here is an example if using a Mongo database:
```
CW_BINARY_FILEPATH="/home/phe.gov.uk/david.wyllie/catwalk/cw_server"
IQTREE="/home/phe.gov.uk/david.wyllie/software/iqtree-2.1.2-Linux/bin/iqtree2"
FASTTREE_DIR="/data/software/fasttree"
```

### Sentry.io connection details
```
FN_SENTRY_URL="https://********************.ingest.sentry.io/*****"
```

FN_SENTRY_URL is an optional url for the [sentry.io](sentry.io) error logging service.  If present and Sentry.io (small free or paid service) configured, error logging will be collated there.    This is very useful for collating  & debugging server side errors.   If considering this, be aware that if identifiable data is in the server, errors trapped may include such and may be sent to the Sentry.io service. 

### Catwalk executable location
The CW_BINARY_FILEPATH is required
```
CW_BINARY_FILEPATH="/home/phe.gov.uk/david.wyllie/catwalk/cw_server"
```

### Variables pointing to Database credentials
If using Oracle databases, two keys are required.  
```
PCA_CONNECTION_CONFIG_FILE='/data/credentials/dbconfig.json'
DB_CONNECTION_CONFIG_FILE='/data/credentials/dbconfig.json'
```
These are as above, and point to database credentials used for Oracle or other secure database connections.  The credentials file referred must exist, or the server will not start up.
More detail on this arrangement is [here](database_credentials.md).


### Oracle toolset location
```
LD_LIBRARY_PATH="/data/software/instantclient_21_1"
```
If you are using an Oracle database, you need to set LD_LIBRARY_PATH to point to where you've installed your Oracle connection software.  This has to be set in the .env file, not globally.  This is [required](https://cx-oracle.readthedocs.io/en/latest/user_guide/initialization.html) for Oracle_cx drivers.


### optional phylogenetics tools
```
IQTREE="/home/phe.gov.uk/david.wyllie/software/iqtree-2.1.2-Linux/bin/iqtree2"
FASTTREE_DIR="/data/software/fasttree"
```
These are the locations of the iqtree and fasttree executable, which are used for tree reconstruction by some utility/experimental scripts in the server.  They are not required for server functioning.

Now you can [try the server out](Testing.md).
