# Demonstrations of the server using real and simulated data

# large_matrix1
simulates TB samples for which a complete distance matrix will be stored

* Scenario tested  
A large, complete matrix is stored.  No repacking is used.

* Background  
An open question is whether the storage engine (now mongodb) is capable of storing a very large matrix with reasonable efficiency.  This is one of the big unknowns.
In the current configuration, the server stores a sparse matrix of pairwise distances.  In the worst case, it has to store a complete matrix.  If the requirement is to be able to store 10^6 samples, that’s at worst about 10^12 cells.
There are databases which happily store billions of records (like mongo) and the current server uses such.

When it inserts, it writes a document for every cell in the distance matrix.  The content is essentially  
**{guid1}:{guid2:{'snv':3}}** 

There is also an option to ‘repack’ these records.
This can be linked to insertion, switched on or off.  This converts the records, initially representing matrix cells, into a ‘row’, or part of a row if the document would be too big.  
**{guid1}:{guid2:{“snv”:3}, guid3:{“snv”:6} … }**   
If we store (say) 10,000 (10^4) cells from each row/column in a single mongodb document (this is possible, within the limits of mongo document sizes), that leaves us with about 10^8 documents.  This is well within the capabilities of mongodb.  However, the repacking operations may be expensive.
It could also be run as a separate process, were it to be slow; it can occur simultaneously to insertion.

In this test,  
repacking is switched off.  
each sequence differs by only 1nt from the reference, so memory usage is very low  
there is no clustering enabled.  

* Outcome measures    
Insertion time  
Time to read all the neighbours of one sample  

* How outcome measures are recorded  
They are written to files by the script  

* How to run the test
``` python findNeighbour3-server.py ../demos/large_matrix_1/config/config.json ```  
-- run the software adding samples to the server
-- edit line 29 to reflect the size of the test required.
-- if there are already (say 1000) samples in the server,
running the test with a higher max_sequences will increase the number of samples up to max_sequences.
``` python demo_large_matrix_1.py ```

* How the output is analysed
Nothing written yet

# AC587  
a collection of 43 mapped samples containing TB data, as well as 38 control TB samples.  The latter are added before the 43 related samples, as they are used by the server to estimate expected N frequencies in real data.
To run the demo:
- make sure mongodb is running
- from the src directory  
-- start the server  
``` python findNeighbour3-server.py ../demos/AC587/config/config.json ```  
-- run the software adding samples to the server  
``` python demo_ac587.py ```

TODO: depict combinations of clusters over time with both algorithms
check that mixed samples go into both.  ** IMPORTANT ** haven't observed this happening, except in unit tests
a gif with the clusters on the L and the MSA on the R might work