[![Test findNeighbour4](https://github.com/davidhwyllie/findNeighbour4/actions/workflows/test_fn4_with_coverage.yml/badge.svg?branch=master)](https://github.com/davidhwyllie/findNeighbour4/actions/workflows/test_fn4_with_coverage.yml)

# Abstract
findNeighbour4 is a server application for investigating bacterial relatedness using reference-mapped data.
Accessible via RESTful webservices, findNeighbour4 maintains a sparse distance matrix in a database
for a set of sequences.  The objective is to be able to rapidly identify likely transmission events.

It has the following features:
* findNeighbour4 itself is accessed by [web services](doc/rest-routes.md). In general, these return json objects.
* Allows incremental addition of new sequences to a collection via RESTful web services.  
* Automatically masks a pre-defined set of nucleotides in the reference mapped data.  Such sites are ignored in pairwise sequence computations.
* Maintains a sparse distance matrix of bacterial sequences using reference mapped sequence data.  Very large matrices can be efficiently and transparently stored.
* Returns multiple sequence alignments.
* [Detects mixtures of different sequences](https://www.biorxiv.org/content/10.1101/681502v1).
* Optionally, automatically performs clustering to a range of SNV thresholds.
* Can detect, and appropriately cluster sequences in the presence of, inter-sample mixtures. [experimental feature - manuscript in progress]
* Allows queries identifying similar sequences, cluster members, and multisequence alignments with  millisecond response times.
* Tracks memory usage, logging to database, during routine operation.
*
It is much faster and uses much less RAM than [previous findNeighbours versions](cf_previous_versions.md). This is achieved using using a compiled component, [CatWalk](https://github.com/dvolk/catwalk.git), to store and compare sequences, developed by Denis Volk.

It was produced as part of the [Modernising Medical Microbiology](http://modmedmicro.nsms.ox.ac.uk/) initiative, together with [Public Health England](https://www.gov.uk/government/organisations/public-health-england).

# Front end
There is a front end, *findNeighbour4 monitor*.  Although not required to run or use findNeighbour4 effectively, it helps to visualise server status and supports ad hoc queries.  In particular, it allows selecting and browsing of samples and clusters of samples in the server, including multisequence alignment, mixture detection, and depiction of their relationships.  This was developed by [Trien Do](https://github.com/TrienDo).

Note that the findNeighbour4 startup script (fn4_startup.sh) does not startup the web front end.  The web front end startup script (fn4_frontend_startup.sh) requires root priviledges to load and run a docker image, whereas the findNeighbour4 server does not.

![findNeighbour4 monitor example page](https://davidhwyllie.github.io/FNMFINDNEIGHBOUR3/img/startup.PNG)  

# Implementation and Requirements
findNeighbour4 is written entirely in python3.  
It operates in Linux only and has only been tested on Ubuntu 20.04.   
To store the relationship between samples, it can use either
- Mongodb (tested with v. 4.4.4 +) 
- A relational database system.  To date we have only tested the Oracle Autonomous Datawarehouse (ADW). [See also Installation details](doc/HowToTest.md) 

# Access
The server can be accessed via RESTful web services from any language.
A python client (fn4client), which calls the REST endpoints and converts output into python objects, is also provided.

# Throughput
The server runs multiple workers, using gunicorn.  The maximum query rate supported needs quantification for each individual stack; it is likely to dependent on database performance.

# Memory usage
This depends on the kind of sequences stored.  
RAM requirements with 810,000 SARS-CoV-2 genomes were about 2GB.  
With ~ 50,000 *M. tuberculosis* genomes, about 10G RAM was used.

# Disc usage
This also depends on the kind of sequences stored.  
A mongodb database holding links between 630,000 SARS-CoV-2 genomes was about 100GB.  
A similar database holding links between ~ 50,000 *M. tuberculosis* genomes, was about 50G.

# More information
[Set up and unit testing](doc/HowToTest.md)  
[Front end installation](doc/frontend.md).  
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
