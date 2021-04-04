#!/usr/bin/env python
""" performs a large scale test of TB sequence storage.

adds simulated TB samples, groups of which have evolved from common ancestors
This program does not generate the simulated sequences; this is done by make_large_sequence_set.py

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
The below will create 100,000 samples; adjust parameters as needed to create more or fewer samples.
Parameter 1 is the number of ancestor sequences.
Parameter 4 is the number of children of each ancestor.  Increasing this to 1000 will produce 1M samples.
Note that these samples are stored on disc.
``` python make_large_sequence_set.py 1000 500 3 100 1e-8 ../demos/lss_tb/seqs```   

To run the test, start up a server, e.g.
```python findNeighbour3-server.py ../demos/lss_tb/config/config.json```

Optionally launch one or more findNeighbour3-dbmanager processes  
```python findNeighbour3-dbmanager.py ../demos/lss_tb/config/config.json```
Then run the test.


**Note**: at present, the server which the test runs against isn't configurable.
It runs against the server running on localhost:5020.  Minor changes to the config file will cange this
The client url needs to be be passed to the call instantiating the fn4Client().

The below inserts until 500 samples are present in the server.
The below inserts 100 samples, then pauses for 1 hour.  
Set the third parameter to 0 to avoid pausing.
```python demo_lss_tb.py 500 ../demos/lss_tb/seqs ../demos/lss_tb/output  100 3600```

If we now do this, then 250 more will be added
```python demo_large_matrix.py 750 ../demos/large_matrix_1/output  100 3600```

* How the output is analysed  

This will analyse all output from the above:
```Rscript demo_depict_timings.R ../demos/large_matrix_1/output```


"""
import copy
import os
import glob
import random
import datetime
import argparse
import pathlib
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from fn4client import fn4Client
    
if __name__ == '__main__':
    
    # maximum number to add
    parser = argparse.ArgumentParser(description='Generate and add to the server groups of sequences which have evolved from each other')
    parser.add_argument('max_sequences', type=int, nargs=1,
                        help='sequences will be added until max_sequences exist in the server.')
    parser.add_argument('inputdir', type=str, nargs=1,
                        help='input fasta files will be read from the inputdir')
    parser.add_argument('outputdir', type=str, nargs=1,
                        help='output will be written to the outputdir')
    parser.add_argument('pause_after', type=int, nargs=1,
                        help='insertion will pause after adding pause_after sequences (set to zero to never pause)')
    parser.add_argument('pause_duration', type=int, nargs=1,
                        help='pause_duration in seconds')

    args = parser.parse_args()
    max_sequences = args.max_sequences[0]
    pause_after = args.pause_after[0]
    pause_duration = args.pause_duration[0]
    outputdir = os.path.abspath(args.outputdir[0])
    inputdir = os.path.abspath(args.inputdir[0])

    # make the directories if they don't exist
    p = pathlib.Path(outputdir)
    p.mkdir(parents=True, exist_ok=True)
    p = pathlib.Path(inputdir)
    p.mkdir(parents=True, exist_ok=True)
        
    # determine input files
    inputfiles = glob.glob(os.path.join(inputdir,'*.fasta'))
    random.shuffle(inputfiles)      # read them in order
    if len(inputfiles)<max_sequences:
        raise ValueError("Asked to add {0} sequences, but only {1} are available in the input directory {2}".format(max_sequences, len(inputfiles), inputdir))
    else:
        inputfiles = inputfiles[0:max_sequences]

    print("opening connection to fn3 server")
    fn4c = fn4Client(baseurl = "http://127.0.0.1:5020")

    # determine all masked positions
    excluded_positions = fn4c.nucleotides_excluded()

    # determine how many samples there are currently in the server.
    nSamples = len(fn4c.guids())
    print("There are {0} existing samples.  Adding more ..".format(nSamples))

    # create output file with header line
    outputfile = os.path.join(outputdir, 'timings_{0}.tsv'.format(nSamples))
    nAdded_this_batch = 0
    with open(outputfile, 'w+t') as f:
        output_line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format('nSamples', 's_insert', 'e_insert', 'd_insert', 's_read', 'e_read', 'd_read')
        f.write(output_line) 
        while nAdded_this_batch < max_sequences:   
        
            if pause_after >0:
                if nAdded_this_batch % pause_after==0:
                     print("Insert paused; will resume after {0} seconds.".format(pause_duration))
                     time.sleep(pause_duration)
        
            # read samples
            inputfile = inputfiles[nAdded_this_batch]
            nSamples +=1  
            nAdded_this_batch +=1      
            with open(inputfile, 'rt') as f_in:
                for record in SeqIO.parse(f_in,'fasta'):               
                    seq = str(record.seq)
                    guid = str(record.id)
       
                    # add
                    print("Inserting", guid, "(samples in this batch = ",nAdded_this_batch,"); will pause every ",pause_after)
                    stime1 = datetime.datetime.now()
                    resp = fn4c.insert(guid=guid, seq=seq)  
                   
                    etime1 = datetime.datetime.now()
                    delta1= etime1-stime1
                    print("Insert yielded status code {0} after {1}".format(resp.status_code, delta1))
    
                    # recover neighbours of guid
                    stime2 = datetime.datetime.now()
                    # check it exists
                    if not fn4c.guid_exists(guid):
                        print("Guid {0} was not inserted".format(guid))   
                    else:
                        neighbours = fn4c.guid2neighbours(guid, threshold=10000000)
                        etime2 = datetime.datetime.now()
                        delta2 = etime2 - stime2
                        print("Recovered {1} neighbours of {0}".format(guid, len(neighbours)))            
                        output_line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(nSamples, stime1, etime1, delta1, stime2, etime2, delta2)
                        f.write(output_line)
                        f.flush()
            # delete the source file - keep disc space down
            try:
                os.unlink(inputfile)
                print("Fasta input file deleted")
            except PermissionError:
                print("Fasta input file not deleted, as locked")
        print("Have added {0} sequences, stopping.".format(nSamples))
        exit(0)              
