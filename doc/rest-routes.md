Recommended REST endpoints
==========================

Conventions
-----------
* *italics* describe parameters which can be supplied.
* *guid* refers to the identifier of a sequence
* All endpoints support GET unless otherwise specified.

List routes available
---------------------
/

Testing the server
------------------
* returns the dictionary posted to the server. Can be used for testing network connectivity.   
/api/v2/mirror  (requires POST)

Describing server configuration
---------------------------------
* Describe server config.  Disabled if not in debug mode.  
/api/v2/server_config
* Return the [last *nrows* of the] server's log of internal memory usage  
/api/v2/server_memory_usage  
/api/v2/server_memory_usage/*nrows*
* Return server time  
/api/v2/server_time
* List nucleotides masked (ignored) by the server in distance computations  
/api/v2/nucleotides_excluded  

Insert into server   
-------------------
/api/v2/insert requires POST; see docs for details  

Search/describe all sequences in the server, each identified by a guid
-----------------------------------------------------------------------
* list all guids (sequence identifiers) in the server  
/api/v2/guids
* list all guids with quality (proportion of Ns in the sequence) over *cutoff*  
/api/v2/guids_with_quality_over/*cutoff*
* list all guids and their examination (i.e. insertion) time  
/api/v2/guids_and_examination_times
* describe annotations (e.g. quality) for all sequences  
/api/v2/annotations

Describe properties/neighbours of a single sequence, identified by a guid
-------------------------------------------------------------------------
* test whether it exists  
/api/v2/*guid*/exists

* specifies threshold, uses default quality cutoff and output format  
/api/v2/*guid*/neighbours_within/*threshold*

* specifying quality cutoff  
uses default output format, as specified in MAXN_PROP_DEFAULT in config file  
/api/v2/*guid*/neighbours_within/*threshold*/with_quality_cutoff/*cutoff*

* specify quality cutoff and output format  
/api/v2/*guid*/neighbours_within/*threshold*/in_format/*returned_format*

Recover masked sequences
------------------------
* recover masked sequences for *guid*  
/api/v2/*guid*/sequence

Multiple sequence alignments
----------------------------
* return multiple sequence alignment for an arbitrary set of sequences, either in json or html format.
/api/v2/multiple_alignment/guids   requires POST; see docs for details

* return multiple sequence alignments of members of cluster  
/api/v2/multiple_alignment/*clustering_algorithm*/*cluster_id*/*output_format*

Clustering
----------
* List the clustering settings operating  
/api/v2/clustering  
* Return the change_id, an incrementing integer which rises are changes are made to clusters  
/api/v2/clustering/*clustering_algorithm*/change_id  
* Return a guid -> cluster lookup
/api/v2/clustering/*clustering_algorithm*/guids2clusters
/api/v2/clustering/*clustering_algorithm*/guids2clusters/after_change_id/*change_id*
