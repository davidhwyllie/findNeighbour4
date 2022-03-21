Local storage of sequence data
===============================

findneighbour4 by default keeps a copy of the reference compressed sequence data on the server on which findneighbour4 runs.

Why does it do this?
--------------------
When large numbers (> 100k samples) are present in the server, if the server has to load the reference compressed sequence data from the database then database access can be slow, delay operationsa and consume large amounts of network bandwidth.

It only does this under two circumstances:
* When it is restarting: if catwalk (a relatedness system used by findneighbour4) is not running, then it will populate catwalk from sequences in the database
* When performing PCA on sequence data

Sequence data, once added to findneighbour4, is immutable.  Therefore, findneighbour4 can keep a local copy of the data in the database.  By default it does, keeps it up to date, and uses it when performing restarts and PCA.

The storage location is determined from the specified LOGFILE in the config file supplied on startup.
```
./fn4_startup.sh config.json
```

findneighbour4 uses the directory in which LOGFILE is specified to create a directory called 
```localcache/{SERVERNAME}```.  So if the config.json file looks like this
```
{
"DESCRIPTION":"A test server operating in on localhost for unit testing using mapped TB data",
"IP":"127.0.0.1",
"INPUTREF":"reference/TB-ref.fasta",
"EXCLUDEFILE":"reference/TB-exclude-adaptive.txt",
"DEBUGMODE":2,
"SERVERNAME":"fn4server1",      
..
"LOGFILE":"/data/fn4storage/fn4/fn4_server1.log",
...
}
```
then it will create a directory ```/data/fn4storage/fn4/localcache/fn4server1```

Inside this will be two subdirectories
```\rcs```
```\mfa```

This directory can contain large amounts of data, so you should ensure there is sufficient space.  It does not need to be backed up, because it will be automatically regenerated.

