Storage of distance matrix in findNeighbour3
--------------------------------------------

If we wish to store a sparse matrix of pairwise differences in a database, there are at least two ways we can do it.  
Suppose our samples are identified by guid1, guid2, guid3, .. guid_n.  We wish to store the distances between these samples.  

### Method 1:  We store the cells in the matrix.    
Thus, if guid1 and guid2 are 3 SNV apart, we store  
{'guid1' : {'guid2':3} }
This type of insert can occur independently of the content of the database, and is very fast.  This format can be considered write-optimised.  
findNeighbour3 always stores the samples in this format on insertion.  

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
* at some point after insertion when the server is idle.  This is the recommended option - repacking after insert is disabled, and a separate process does repacking in the background.  The software doing this is ```findNeighbour3-dbmanager.py```.  See [here](../doc/HowToTest.md) for more details.

### REPACK_FREQ
This behaviour is controlled by the REPACK_FREQ setting in the server config file.

If REPACK_FREQ=0, there will be one document for every non-empty matrix cell.  This is the recommended setting 
* where mongodb storage size is not limiting 
* for real-world applications where the samples come in batches, and new sample insertion into the server has to be as fast as possible
* when the findNeighbour3_dbmanager application is running in the background

if REPACK_FREQ>0, then if a guid has REPACK_FREQ-1 neighbours, then a 'repack' operation occurs.
For example, if REPACK_FREQ=1, a repack occurs after every insert.
This minimises mongodb size but increases the time needed for insertion.

Repacking doesn't alter the results obtained at all.
