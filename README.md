# Abstract
findNeighbour3 is a server application for investigating bacterial relatedness using reference-mapped data.
Accessible via RESTful webservices, findNeighbour3 maintains a sparse distance matrix in a database
for a set of sequences.

It has the following features:
* Allows incremental addition of new sequences to a collection via RESTful web services.  
* Automatically masks a pre-defined set of nucleotides in the reference mapped data.  Such sites are ignored in pairwise sequence computations.
* Maintains a sparse distance matrix of bacterial sequences using reference mapped sequence data.  Very large matrices can be efficiently and transparently stored.
* Uses a highly compressed sequence representation, relying on compression to local reference, having first applied compression to the reference sequence to which mapping occurred.  This *double delta* technique aids storage of large numbers of sequences in RAM.
* Tracks memory usage, logging to database, during routine operation.
* Returns pairwise distance matrices.
* Returns multiple sequence alignments.
* Automatically performs clustering to a range of SNV thresholds.
* Can detect, and appropriately cluster sequences in the presence of, inter-sample mixtures.
* Allows rapid additions (< 2 sec per sample adding to a 30,000 sequence collection)
* Allows queries identifying similar sequences, cluster members, and multisequence alignments with  millisecond response times.
* Would readily allow attachment of arbitrary metadata to each sequence, but the front end for this is not implemented.

It was produced as part of the [Modernising Medical Microbiology](http://modmedmicro.nsms.ox.ac.uk/) initiative, together with [Public Health England](https://www.gov.uk/government/organisations/public-health-england).

# Implementation and Requirements
findNeighbour3 is written entirely in python3.  
It operates on Windows and Linux environments.    
It uses mongodb as a storage layer.

# Access
The server can be accessed via RESTful web services from any language.
A python client (fn3client), which calls the REST endpoints and converts output into python objects, is also provided.

# Comparison with findNeighbour2
findNeighbour3 is a development of [findNeighbour2](https://github.com/davidhwyllie/findNeighbour2).
findNeighbour3's RESTful API is backwards compatible with that of findNeighbour2, but offers increased functionality.  
There are the following other differences:
* It uses additional compression (*double delta*), resulting in it needing about 30-50% of the memory required by findNeighbour2.
* It uses mongodb, not relational databases, for persistent storage.
* Queries are much faster for large numbers of samples
* It performs clustering.
* It is 'mixture-aware' and implements an approach for detecting mixed samples.
* Dependencies on linux-specific packages have been removed.
* It does not use any storage in a filesystem, except for logging.
* Internally, it has been refactored into four components, managing the web server, in-memory storage, on-disc storage, and clustering.
* It is only accessible via a RESTful endpoint.  The xmlrpc API included with findNeighbour2 has been removed.

# More information
[How to test it](doc/HowToTest.md)  
[Endpoints](doc/rest-routes.md)
[Demonstrations using real and simulated data](doc/demos.md)

# Publication
A publication describing findNeighbour3 is planned.

A publication describing findNeighbour2 is in BMC Bioinformatics:
*BugMat and FindNeighbour: command line and server applications for investigating bacterial relatedness*
DOI : 10.1186/s12859-017-1907-2 (https://dx.doi.org/10.1186/s12859-017-1907-2)

# Large test data sets
Test data sets of *N. meningitidis*, *M. tuberculosis* and *S. enterica* data are available to download [here](https://ora.ox.ac.uk/objects/uuid:82ce6500-fa71-496a-8ba5-ba822b6cbb50).  These are .tar.gz files, to a total of 80GB.

