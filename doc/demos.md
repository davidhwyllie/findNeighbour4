# Demonstrations of the server

##Using simulated data

### large_matrix1
simulates samples for which a complete distance matrix will be stored

#### Scenario tested  
A large matrix has to be stored, perhaps because all the samples are very closely related or because the clinically relevant SNV threshold is very high.

#### Background    
An important question is whether the storage engine (mongodb) is capable of storing a very large matrix with reasonable efficiency.  
In the current configuration, the server stores a sparse matrix of pairwise distances.  In the worst case, it has to store a complete matrix.  If the requirement is to be able to store 10^6 samples, that’s at worst about 10^12 cells.
This is a large number of cells to write.

When it inserts, the server writes a document for every cell in the distance matrix.  The content is essentially  
**{guid1}:{guid2:{'snv':3}}** 

There is also an option to ‘repack’ these records.
This can be linked to insertion, switched on or off.  This converts the records, initially representing matrix cells, into a ‘row’, or part of a row if the document would be too big.  
**{guid1}:{guid2:{“snv”:3}, guid3:{“snv”:6} … }**   
If we store (say) 10,000 (10^4) cells from each row/column in a single mongodb document (this is possible, within the limits of mongo document sizes), that leaves us with about 10^8 documents.  This is well within the capabilities of mongodb.  However, the repacking operations may be expensive.
It could also be run as a separate process, were it to be slow; it can occur simultaneously to insertion.

In this test,  
* repacking immediately after insertion is switched off.
* insertion occurs in batches, with a pause after each batch of samples.  This simulates the real life scenario when
* sample addition occurs in batches
* each sequence differs by only 1nt from the reference, so memory usage is very low  
there is no clustering enabled.  

#### Outcome measures    
Insertion time  
Time to read all the neighbours of one sample  
domes
#### How outcome measures are recorded  
They are written to files by the script  

#### How to run the test

To run the test, start up a server, e.g.  
```python findNeighbour3-server.py ../demos/large_matrix_1/config/config.json```

Optionally launch one or more findNeighbour3-dbmanager processes  
```python findNeighbour3-dbmanager.py ../demos/large_matrix_1/config/config.json```  
Then run the test.

**Note**: at present, the server which the test runs against isn't configurable.
It runs against the server running on localhost:5020.  Minor changes to the test script will change this; 
the client url needs to be be passed to the call instantiating the fn3Client().

The below inserts until 500 samples are present in the server.
The below inserts 100 samples, then pauses for 1 hour.  
Set the third parameter to 0 to avoid pausing.  
```python demo_large_matrix.py 500 ../demos/large_matrix_1/output  100 3600```

If we now do this, then 250 more will be added
```python demo_large_matrix.py 750 ../demos/large_matrix_1/output  100 3600```

#### How the output is analysed  

This will analyse all output from the above:
```Rscript demo_depict_timings.R ../demos/large_matrix_1/output```


### lss_tb
adds simulated TB samples, groups of which have evolved from common ancestors

#### Scenario tested
A large number of samples are derived from sequencing of TB over a number of years.
When sequencing started on a large scale (ca. 2015), a number of TB strains ('ancestor sequences') were
circulating.  Over time, these transmitted and their descendents were isolated.
Sequences corresponding to this scenario are simulated, and added to the server in a random order.

No repacking is used following insertion, but findNeighbour3-dbmanager can be run in the background to achieve repacking.
No clustering is enabled.

#### Background  
This test is designed to test the speed of insertion and storage requirements of a set of samples simulating those encountered with a highly clonal pathogen, such as TB.

#### Outcome measures    
Insertion time  
Time to read all the neighbours of one sample  

#### How outcome measures are recorded  
They are written to files by the script  

#### How to run the test
Simulate the samples to be added:
The above will create 50,000 samples; adjust parameters as needed to create more or fewer samples.
Parameter 1 is the number of ancestor sequences.
Parameter 4 is the number of children of each ancestor.  Increasing this to 1000 will produce 1M samples.
``` python make_large_sequence_set.py 1000 500 3 50 1e-8 ../output/lss_tb/seqs```   

To run the test, start up a server, e.g.
```python findNeighbour3-server.py ../demos/lss_tb/config/config.json```

Optionally launch one or more findNeighbour3-dbmanager processes  
```python findNeighbour3-dbmanager.py ../demos/lss_tb/config/config.json```
Then run the test.


**Note**: at present, the server which the test runs against isn't configurable.
It runs against the server running on localhost:5020.  Minor changes to the config file will cange this
The client url needs to be be passed to the call instantiating the fn3Client().

The below inserts until 500 samples are present in the server.
The below inserts 100 samples, then pauses for 1 hour.  
Set the third parameter to 0 to avoid pausing.
```python demo_large_matrix.py 500 ../demos/large_matrix_1/output  100 3600```

If we now do this, then 250 more will be added
```python demo_large_matrix.py 750 ../demos/large_matrix_1/output  100 3600```

* How the output is analysed  

This will analyse all output from the above:
```Rscript demo_depict_timings.R ../demos/large_matrix_1/output```

## Using real data

### AC587  
a collection of 43 mapped samples containing TB data, as well as 38 control TB samples.  The latter are added before the 43 related samples, as they are used by the server to estimate expected N frequencies in real data.
To run the demo:
- make sure mongodb is running
- from the src directory  
-- start the server  
``` python findNeighbour3-server.py ../demos/AC587/config/config.json ```  
-- run the software adding samples to the server  
``` python demo_ac587.py ```


