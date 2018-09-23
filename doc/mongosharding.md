Sharding mongodb
================

By default, findNeighbour3 makes suitable indices for sharding the mongodb.
However, sharding has to be set up when the mongodb is configured.  
This has not been tested so far.  

There are two steps:  
1. Enable sharding on the database
2. Specify sharding for each collection to be sharded.  
There are only two collections which are worth sharding in the schema used:   
* guid2neighbour  
This contains the SNV matrix.  The _id key, which is the guid, is a suitable sharding key.  
* refcompressedseq.chunks  
This contains the reference compressed sequences.  The files_id_1_n_1 key is a suitable sharding key.  


