# Abstract
findNeighbour4 is a server application for investigating bacterial relatedness using reference-mapped data.
Accessible via RESTful webservices, findNeighbour4 maintains a sparse distance matrix in a database
for a set of sequences.


It has the following features:
* Allows incremental addition of new sequences to a collection via RESTful web services.  
* Automatically masks a pre-defined set of nucleotides in the reference mapped data.  Such sites are ignored in pairwise sequence computations.
* Maintains a sparse distance matrix of bacterial sequences using reference mapped sequence data.  Very large matrices can be efficiently and transparently stored.
* Returns pairwise distance matrices.
* Returns multiple sequence alignments.
* [Detects mixtures of different sequences](https://www.biorxiv.org/content/10.1101/681502v1).
* Automatically performs clustering to a range of SNV thresholds.
* Can detect, and appropriately cluster sequences in the presence of, inter-sample mixtures. [experimental feature - do not use in production]
* Allows queries identifying similar sequences, cluster members, and multisequence alignments with  millisecond response times.
* Tracks memory usage, logging to database, during routine operation.
* Allow attachment of arbitrary metadata to each sequence, but the front end for this is not implemented.

Compared with findNeighbours 2 and 3, it works much faster (>100x faster inserts) and uses about 5% of the RAM.  
This is achieved using using a compiled component, CatWalk, to store and compare sequences.  The component was developed by Denis Volk (Modernising Medical Microbiology, Oxford) - publication planned.
It was produced as part of the [Modernising Medical Microbiology](http://modmedmicro.nsms.ox.ac.uk/) initiative, together with [Public Health England](https://www.gov.uk/government/organisations/public-health-england).

# Front end
There is a front end, *findNeighbour4 monitor*.  (this also works with findNeighbour4; the API of the two servers is identical).    Although not required to run or use findNeighbour4 effectively, it helps to visualise server status and supports ad hoc queries.  In particular, it allows selecting and browsing of samples and clusters of samples in the server, including multisequence alignment, mixture detection, and depiction of their relationships.  

Note that the findNeighbour4 startup script (fn4_startup.sh) does not startup the web front end.  The web front end startup script requires root priviledges (loads and runs a docker image) but the findNeighbour4 server does not.

![findNeighbour4 monitor example page](https://davidhwyllie.github.io/FNMFINDNEIGHBOUR3/img/startup.PNG)  
The *findNeighbour4 monitor* is easy to use and to install.  See [details](doc/frontend.md).  
findNeighbour4 itself is accessed by [web services](doc/rest-routes.md). In general, these return json objects.

# Implementation and Requirements
findNeighbour4 is written entirely in python3.  
It operates on Windows and Linux environments.    
It uses mongodb as a storage layer.

# Access
The server can be accessed via RESTful web services from any language.
A python client (fnclient), which calls the REST endpoints and converts output into python objects, is also provided.

# Memory and disc usage
This depends on the kind of sequences stored.  The server has been tested with SARS-CoV-2 genomes & with  *M. tuberculosis*:

# Comparison with findNeighbour2
findNeighbour4 is a development of [findNeighbour2](https://github.com/davidhwyllie/findNeighbour2).
findNeighbour4's RESTful API is backwards compatible with that of findNeighbour2, but offers increased functionality.  
There are the following other differences:
* It uses much less RAM and is much faster, due to use of a specialised CatWalk component, and other changes.  
* It uses mongodb, not relational databases, for persistent storage.
* Queries are much faster for large numbers of samples
* It performs clustering.
* It is 'mixture-aware' and implements an approach for detecting mixed samples.
* Clustering requires the linux specific *networkit* library.
* It does not use any storage in a filesystem, except for logging.
* Internally, it has been refactored into four components, managing the web server, in-memory storage, on-disc storage, and clustering.
* It is only accessible via a RESTful endpoint.  The xmlrpc API included with findNeighbour2 has been removed.

# Comparison with findNeighbour3
The server has been much more heavily tested, and is much faster.  Production use of findNeighbour3 is not recommended.

# More information
[Set up and unit testing](doc/HowToTest.md)  
[Endpoints](doc/rest-routes.md)  
[Demonstrations using real and simulated data](doc/demos.md)  
[Integration tests](doc/integration.md)

# Publications
A publication describing findNeighbour4 implementation & performance is planned.  
A publication describing findNeighbour2 is in BMC Bioinformatics:  
*BugMat and FindNeighbour: command line and server applications for investigating bacterial relatedness*
DOI : 10.1186/s12859-017-1907-2 (https://dx.doi.org/10.1186/s12859-017-1907-2)  

The nature of the mixPORE (mixture detection algorithm) provided by the server, and its application to *M. tuberculosis* mixture detection is described [here](https://www.biorxiv.org/content/10.1101/681502v1).

# Large test data sets
For the detection of mixtures, please see the additional test data sets [here](doc/demos_real.md).
