Frontend
--------
[/ui/info](/ui/info) server status page   
**/ui/server_status/*absdelta*/*stats_type*/*nrows* **  visualise server disc and application memory usage 

Search/describe sequences in the server
-----------------------------------------------------------------------
[/api/v2/guids](/api/v2/guids)  list all guids (sequence identifiers) in the server  
**/api/v2/guids_beginning_with/*startstr* **  list all guids starting with *startstr*.  Very fast algorithm, suitable for on-keypress prediction of matching guids.  Only up to 30 results are returned.  If more than 30 records match, an empty list is returned.  
**/api/v2/guids_with_quality_over/*cutoff* ** list all guids with quality (proportion of Ns in the sequence) over *cutoff*    
[/api/v2/guids_and_examination_times](/api/v2/guids_and_examination_times) list all guids and their examination (i.e. insertion) time   
[/api/v2/annotations](/api/v2/annotations) describe annotations (e.g. quality) for all sequences   

Describe properties/neighbours of a single sequence, identified by a guid
-------------------------------------------------------------------------
**/api/v2/*guid*/exists**  test whether it exists  
**/api/v2/*guid*/annotation**  return metadata for the guid  
**/api/v2/*guid*/neighbours_within/*threshold* ** specifies threshold, uses default quality cutoff and output format   
**/api/v2/*guid*/neighbours_within/*threshold*/with_quality_cutoff/*cutoff* ** specify quality cutoff; uses default output format   
**/api/v2/*guid*/neighbours_within/*threshold*/in_format/*returned_format* **  specify quality cutoff and output format    

Recover masked sequences
------------------------
**/api/v2/*guid*/sequence**  recover masked sequence for *guid*  

Multiple sequence alignments
----------------------------
**/api/v2/multiple_alignment/guids**   requires POST; see docs for details.  return multiple sequence alignment for an arbitrary set of sequences, either in json or html format.  
**/api/v2/multiple_alignment_cluster/*clustering_algorithm*/*cluster_id*/*output_format* ** return multiple sequence alignments of members of cluster cluster_id; output_format can be json, json-records, html or fasta

Clustering
----------
[/api/v2/clustering](/api/v2/clustering) List the clustering settings operating  
**/api/v2/clustering/*clustering_algorithm*/change_id**  Return the change_id, an incrementing integer which rises as changes are made to clusters  
**/api/v2/clustering/*clustering_algorithm*/clusters**  Lists clusters, and number of mixed and unmixed samples in each.   
**/api/v2/clustering/*clustering_algorithm*/guids2clusters**  Return a guid -> cluster lookup  
**/api/v2/clustering/*clustering_algorithm*/guids2clusters/after_change_id/*change_id* ** Return a guid -> cluster lookup after some particular point in time.  
**/api/v2/clustering/*clustering_algorithm*/cluster_ids**  Return unique cluster_ids for *clustering_algorithm*  
**/api/v2/clustering/*clustering_algorithm*/*cluster_id*/network** return a Cytoscape.js json string.  See ui/cytoscape_viewer1.html for example code consuming this.

Insert into server   
-------------------
**/api/v2/insert** insert into the server requires POST; see docs for details  

Server config & testing
---------------------------------
**/**  display the routes available  (this page)  
[/api/v2/server_config](/api/v2/server_config)  Describe server config.  Disabled if not in debug mode.  
[/api/v2/server_memory_usage](/api/v2/server_memory_usage) Return the server's log of internal memory and disc usage     
**/api/v2/server_memory_usage/*nrows* ** Return the [last *nrows* of the] server's log of internal memory and disc usage   
[/api/v2/server_time](/api/v2/server_time) Return server time   
[/api/v2/snpceiling](/api/v2/snpceiling) Return maximum snv stored by the server  
[/api/v2/nucleotides_excluded](/api/v2/nucleotides_excluded) List nucleotides masked (ignored) by the server in distance computations  
**/api/v2/mirror**  (requires POST)  returns the dictionary posted to the server. Can be used for testing network connectivity.    
**/api/v2/raise_error/*component*/*token*/** raises an error internally.  Can be used to test error logging.  Disabled unless in debug mode.

Notes
-----------
*italics* describe parameters which can be supplied.  
*guid* refers to the identifier of a sequence.   
All endpoints support GET unless otherwise specified.  

