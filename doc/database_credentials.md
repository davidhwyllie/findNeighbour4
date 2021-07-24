Database credentials
=======================
FindNeighbour4 requires somewhere to store credentials for database connectivity.
FindNeighbour4 uses a configuration file, set at startup, to store various details of the server's functioning.  This is always required.  In addition, sometimes a separate credentials file is also needed.

## FNPERSISTENCE_CONNSTRING in configuration json file 
Here is part of the findNeighbour4 configuration file:
```
{
"DESCRIPTION":"A test server operating in on localhost for unit testing using mapped TB data and mongo",
"IP":"127.0.0.1",
"INPUTREF":"reference/TB-ref.fasta",
"EXCLUDEFILE":"reference/TB-exclude-adaptive.txt",
"DEBUGMODE":2,
"SERVERNAME":"fn3_unittesting",      
"FNPERSISTENCE_CONNSTRING":"mongodb://127.0.0.1",
......
}
```
The FNPERSISTENCE_CONNSTRING can contain connection details.  In the above example, it's connecting to mongodb on localhost with no authentication.  It will create a database called fn3_unittesting (as set in SERVERNAME).

It could also look something like
```
"mysql+pymysql://user:password@localhost:3306/test_db"
```
In this case, test_db is the name of the databas to connect to, on a mysql database on localhost.
This is fine, provided all the connection information requires is in the connection string.
If it isn't (as with Oracle databases), then a different process is used.  


## Separate credentials JSON

The alternative is a three phase process
 - Set an environment variable pointing at a credentials file
 - Store a json dictionary, the key to which is a particular configuration (e.g. 'unittest', 'production', 'development'), and which contains all the necessary credentials, in a defined location.
 - Store the key to the dictionary entry (e.g. 'production') in he FN3PERSISTENCE_CONNSTRING.

If you're connecting to Oracle, you have to use this mehtod.  For details, [see here](configuring_oracle_connections.md)