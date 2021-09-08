Storage of distance matrix when using findNeighbour4
--------------------------------------------

If we wish to store a sparse matrix of pairwise differences in a database, there are at least two ways we can do it.  
Suppose our samples are identified by guid1, guid2, guid3, .. guid_n.  We wish to store the distances between these samples.  

### Method 1:  We store the cells in the matrix.    
Thus, if guid1 and guid2 are 3 SNV apart, we store  
{'guid1' : {'guid2':3} }
This type of insert can occur independently of the content of the database, and is very fast.  This format can be considered write-optimised.  
findNeighbour4 always stores the samples in this format on insertion.  In the RDBMS data store, data stays in this format.  In the Mongodb backend, data is reformatted, see below

### Method 2: We store the rows of the matrix.    
for example, if guid1 is related to guid2 by 3 SNV and to guid3 by 5 SNV, we store  
{'guid1':{'guid2':3, 'guid3':5}}  
This format can have many fewer documents (entries) in the database.  If guid1 had 1,000 close neighbours, we'd be storing one document, not 1,000.
Since mongodb works best if indices are kept in RAM, the in-ram index size will be correspondingly smaller.  This format can be considered read-optimised.  Converting the write-optimised format to the read-optimised one doesn't alter the results produced by findNeighbour3, but it does   
* reduce database size, sometimes markedly  
* take time  
 
This conversion can occur:  
* never.  For small numbers of samples, this is OK.  
* immediately after insertion.  This will slow up insertion, but keep the mongodb as small as possible.  
* at some point after insertion when the server is idle.  This is the recommended option and runs by default with a separate processes repacking in the background.  The software doing this is ```findNeighbour4-dbmanager.py```.  and will run on startup automatically if you run the ./fn4_startup.sh script.  See [here](../doc/HowToTest.md) for more details.
* This conversion does not occur with the RDBMS backend.  ```findNeighbour4-dbmanager``` instances will automatically terminate if you are using an RDBMS backend.

