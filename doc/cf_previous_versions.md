
Comparisons with previous versions 
==================================================

## Comparison with findNeighbour2
findNeighbour4 is a development of [findNeighbour2](https://github.com/davidhwyllie/findNeighbour2).
findNeighbour4's RESTful API is backwards compatible with that of findNeighbour2, but offers increased functionality.  

There are the following other differences:
* It uses much less RAM and is much faster, due to use of a specialised CatWalk component, and many other changes.  
* It supports parallel server processes allowing concurrent queries, using the gunicorn library.
* It uses mongodb or relational databases, for persistent storage.
* Queries are much faster for large numbers of samples
* It performs clustering using single nucleotide distances.  This is an experimental feature, only tested using M. tuberculosis genomes.
* It is 'mixture-aware' and implements an approach for detecting mixed samples.
* Clustering requires the linux specific *networkit* library, which is linux specific.
* It does not use any storage in a filesystem, except for logging.
* It is only accessible via a RESTful endpoint.  The xmlrpc API included with findNeighbour2 has been removed.

## Comparison with findNeighbour3
The server 
- is much faster (>100x faster inserts) and uses about 5% of the RAM.  This is achieved using using a compiled component, [CatWalk](https://gitea.mmmoxford.uk/dvolk/catwalk.git), to store and compare sequences, developed by Denis Volk.
- can use RDBMS backends
- can be run as multiple WSGI processes


## Contributors
It was produced as part of the [Modernising Medical Microbiology](http://modmedmicro.nsms.ox.ac.uk/) initiative, together with [Public Health England](https://www.gov.uk/government/organisations/public-health-england).

