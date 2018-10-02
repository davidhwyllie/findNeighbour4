Frontend
--------
**/ui/info** [server status page](/ui/info)   
**/ui/server_status/*absdelta*/*stats_type*/*nrows* **  visualise server disc and application memory usage 

Search/describe sequences in the server
-----------------------------------------------------------------------
**/api/v2/guids**  list all guids (sequence identifiers) in the server    
**/api/v2/guids_with_quality_over/*cutoff* ** list all guids with quality (proportion of Ns in the sequence) over *cutoff*    
**/api/v2/guids_and_examination_times** list all guids and their examination (i.e. insertion) time   
**/api/v2/annotations** describe annotations (e.g. quality) for all sequences   

Describe properties/neighbours of a single sequence, identified by a guid
-------------------------------------------------------------------------
**/api/v2/*guid*/exists**  test whether it exists
**/api/v2/*guid*/neighbours_within/*threshold* ** specifies threshold, uses default quality cutoff and output format   
**/api/v2/*guid*/neighbours_within/*threshold*/with_quality_cutoff/*cutoff* **  specifying quality cutoff uses default output format, as specified in MAXN_PROP_DEFAULT in config file   
**/api/v2/*guid*/neighbours_within/*threshold*/in_format/*returned_format* **  specify quality cutoff and output format    

Recover masked sequences
------------------------
**/api/v2/*guid*/sequence**  recover masked sequences for *guid*  

Multiple sequence alignments
----------------------------
**/api/v2/multiple_alignment/guids**   requires POST; see docs for details.  return multiple sequence alignment for an arbitrary set of sequences, either in json or html format.
**/api/v2/multiple_alignment/*clustering_algorithm*/*cluster_id*/*output_format* ** return multiple sequence alignments of members of cluster  

Clustering
----------
**/api/v2/clustering** List the clustering settings operating  
**/api/v2/clustering/*clustering_algorithm*/change_id**  Return the change_id, an incrementing integer which rises are changes are made to clusters  
**/api/v2/clustering/*clustering_algorithm*/guids2clusters**  Return a guid -> cluster lookup  
**/api/v2/clustering/*clustering_algorithm*/guids2clusters/after_change_id/*change_id* ** Return a guid -> cluster lookup after some particular point in time.

Insert into server   
-------------------
**/api/v2/insert** insert into the server requires POST; see docs for details  

Server config & testing
---------------------------------
**/**  display the routes available  (this page)  
**/api/v2/server_config**  Describe server config.  Disabled if not in debug mode.  
**/api/v2/server_memory_usage** Return the server's log of internal memory and disc usage     
**/api/v2/server_memory_usage/*nrows* ** Return the [last *nrows* of the] server's log of internal memory and disc usage   
**/api/v2/server_time** Return server time   
**/api/v2/nucleotides_excluded** List nucleotides masked (ignored) by the server in distance computations  
**/api/v2/mirror**  (requires POST)  returns the dictionary posted to the server. Can be used for testing network connectivity. 
**/api/v2/raise_error/*component*/*token*/** raises an error internally.  Can be used to test error logging.  Disabled unless in debug mode.

Notes
-----------
*italics* describe parameters which can be supplied.  
*guid* refers to the identifier of a sequence.   
All endpoints support GET unless otherwise specified.  

