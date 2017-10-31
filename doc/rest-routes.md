Recommended REST endpoints
==========================
Describing server configuration
---------------------------------
/v2/server_config  
/v2/server_memory_usage  
/v2/server_time  
/v2/nucleotides_excluded  

Insert into server   
-------------------
/v2/insert [POST; see code for details]  

Search/describe all sequences in the server, each identified by a guid
-----------------------------------------------------------------------
/v2/guids  
/v2/guids_with_quality_over/*cutoff*  
/v2/guids_and_examination_times  
/v2/annotations  
/v2/neighbours_within/*threshold*  
/v2/nucleotides_excluded

Describe properties/neighbours of a single sequence, identified by a guid
-------------------------------------------------------------------------
* test whether it exists  
/v2/*guid*/exists

* specifies threshold, uses default quality cutoff and output format  
/v2/*guid*/neighbours_within/*threshold*

* specifying quality cutoff  
uses default output format, as specified in MAXN_PROP_DEFAULT in config file  
/v2/*guid*/neighbours_within/*threshold*/with_quality_cutoff/*cutoff*

* specify quality cutoff and output format  
/v2/*guid*/neighbours_within/*threshold*/in_format/*returned_format*

Compare two sequences in detail
-------------------------------
* includes positions of variation between sequences  
/v2/*guid1*/*guid2*/detailed_comparison

Legacy endpoints  
================  
/sample/guids/*reference*  
/sample/guids_cutoff/*reference*/<float:cutoff>		  
/sample/guids_and_time/*reference*  
/sample/annotation/*reference*  
/sample/walks/processed/*guid*/*reference*/*method*  
/sample/walks/snp/*guid*/*reference*/*threshold*/*method*  
/sample/findneighbour/snp/*guid*/*reference*/*threshold*/*method*/*cutoff* 
/sample/neighbours/*guid*/*reference*/*threshold*/*method*  
