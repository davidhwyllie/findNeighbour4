# Configuration an Oracle ADW connection

These notes are about how to configure a findNeighbour4 connection to an [Oracle Autonomous Database](https://blogs.oracle.com/oraclemagazine/getting-started-with-autonomous).
There are several steps.

### Install dependencies
See [instructions](https://cx-oracle.readthedocs.io/en/latest/user_guide/installation.html)
You will need to do something like
```
        wget https://download.oracle.com/otn_software/linux/instantclient/211000/instantclient-basic-linux.x64-21.1.0.0.0.zip
```
and unzip into a directory.

### Set LD_LIBRARY_PATH
You have to set the [LD_LIBRARY_PATH variable](https://www.oracle.com/database/technologies/instant-client/linux-x86-64-downloads.html).
You have add this to the .env file in the project root, e.g. 
```
LD_LIBRARY_PATH="/software/instantclient_21_1"
```

### Provide configuration details
For Oracle connections, you have to provide these in a configuration file, the format of which is findneighbour specific.
You have tell findneighbour where the database connections are.  Add these lines to your .env file; the file pointed to should not be under source control, and should be readable by the findNeighbour4 software.
```
DB_CONNECTION_CONFIGFILE="/secret/config.json"
PCADB_CONNECTION_CONFIGFILE="/secret/config.json"
```
It's fine for both variables to point to the same file.

### Format of the configuration file
The /secret/config.json file should look something like this.
```
        {
            'prod':{'DBTYPE':'oracle',
                    'ENGINE_NAME':''oracle+cx_oracle://PROD:97bxxxxxxxxX@(description: .........)))',
                    'TNS_ADMIN':'/secrets/oracle/pca_prod'
                    }
        }
```
*NOTE: as per normal json conventions, escape quotes (i.e. \" not " around the certificate name, otherwise SSL connections will fail)*
*Note: it's fine to have multiple entries in this file, e.g. for production and testing, e.g.*
```
 {
            'prod':{'DBTYPE':'oracle',
                    'ENGINE_NAME':''oracle+cx_oracle://PROD:97bxxxxxxxxX@(description: .........)))',
                    'TNS_ADMIN':'/secrets/oracle/pca_prod'
                    },

            'dev':{'DBTYPE':'oracle',
                    'ENGINE_NAME':''oracle+cx_oracle://PROD:97bxxxxxxxxX@(description: .........)))',
                    'TNS_ADMIN':'/secrets/oracle/pca_prod'
                    }
        }
```
* Note: the DBTYPE and ENGINE_NAME keys are essential, and can be used to store connection details to non-Oracle database too.  If you want to do this, you can omit the TNS_ADMIN key*  

The content of the config.json file is Oracle provided credentials.  Follow these steps:

1. Download your OCI wallet, & unzip it somewhere.  The OCI wallet can be downloaded (if you have sufficient permissions) from the OCI cloud web portal.
2. Set the TNS_ADMIN entry in the configuration file to point to this directory.  You have to do this - it's not optional; if you've set it globally, when run in a virtual environment, findneighbour4 won't see it.
3. Edit the WALLET_LOCATION in the sqlnet.ora file to point to the relevant directory pointed to by TNS_ADMIN; e.g. if you've got your wallet in /data/credentials/oci_test, set WALLET_LOCATION = (SOURCE = (METHOD = file) (METHOD_DATA = (DIRECTORY="/data/credentials/oci_test")))
4. Create a user with relevant privileges via the ADW front end or otherwise (see below)
5. Set the ENGINE_NAME env var.  An example of this is as below (parts are redacted)
oracle+cx_oracle://scott:tigerX22@(description= (retry_count=20)(retry_delay=3)(address=(protocol=tcps)(port=1522)(host=host.oraclecloud.com))(connect_data=(service_name=redacted))(security=(ssl_server_cert_dn="redacted")))

This service credentials required can be found in the tnsnames.ora file.  They are used to construct a connection URL compatible with the [SQLAlchemy ](https://www.sqlalchemy.org/) module used to connect by FindNeighbour4.  For more details see: [constructing the url](https://stackoverflow.com/questions/14140902/using-oracle-service-names-with-sqlalchemy/35215324); [dialects](https://docs.sqlalchemy.org/en/14/dialects/oracle.html#dialect-oracle-cx_oracle-connect); [alternatives](https://stackoverflow.com/questions/37471892/using-sqlalchemy-dburi-with-oracle-using-external-password-store); 
     

The user name (scott, in the above) and password (tigerX22 in the above) you set in the next step.

## Configuring interactions with external OCI databases

Your application will need to run as a user (we'll call it PCADB) will need some priviledges granted.
The exact privileges required involving creating, dropping tables & indexes, as well as inserting and deleting data.
For more details, see [here](https://blogs.oracle.com/sql/how-to-create-users-grant-them-privileges-and-remove-them-in-oracle-database) and [here](https://docs.oracle.com/en/database/oracle/oracle-database/19/sqlrf/GRANT.html#GUID-20B4E2C0-A7F8-4BC8-A5E8-BE61BDC41AC3).

```
CREATE USER PCADB IDENTIFIED BY 'MyPassword1234!';
GRANT CONNECT TO PCADB;
GRANT CREATE SESSION TO PCADB;
GRANT CREATE SEQUENCE TO PCADB;
GRANT CREATE TABLE TO PCADB;
GRANT CREATE SYNONYM TO PCADB;
ALTER USER PCADB DEFAULT TABLESPACE DATA quota unlimited on DATA;
