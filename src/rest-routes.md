Recommended REST endpoints
==========================
* Describing server configuration   
/v2/server_config  
/v2/server_memory_usage  
/v2/server_time  
/v2/nucleotides_excluded  

* Insert into server   
/v2/insert [POST; see code for details]  

* Search/describe all sequences in the server, each identified by a guid   
/v2/guids  
/v2/guids_with_quality_over/<float:cutoff>  
/v2/guids_with_quality_over/<int:cutoff>  
/v2/guids_and_examination_times  
/v2/annotations  
/v2/neighbours_within/<int:threshold>  

* Describe properties/neighbours of a single sequence, identified by a guid   
/v2/<string:guid>/exists  
/v2/<string:guid>/neighbours_within<int:threshold>  
/v2/<string:guid>/<int:threshold>/neighbours_within/<float:cutoff>  
/v2/<string:guid>/<int:threshold>/neighbours_within/<float:cutoff>/<int:returned_format>  
/v2/<string:guid>/<int:threshold>/neighbours_within/<int:cutoff>  
/v2/<string:guid>/<int:threshold>/neighbours_within/<int:cutoff>/<int:returned_format>  
/v2/<string:guid1>/<string:guid2>/detailed_comparison  



Legacy endpoints
================
/sample/guids/<string:reference>
/sample/guids_cutoff/<string:reference>/<float:cutoff>		
/sample/guids_cutoff/<string:reference>/<int:cutoff>		
/sample/guids_and_time/<string:reference>
/sample/annotation/<string:reference>
/sample/walks/processed/<string:guid>/<string:reference>/<string:method>
/sample/walks/snp/<string:guid>/<string:reference>/<int:threshold>/<string:method>
/sample/findneighbour/snp/<string:guid>/<string:reference>/<int:threshold>/<string:method>/<float:cutoff>
/sample/findneighbour/snp/<string:guid>/<string:reference>/<int:threshold>/<string:method>/<int:cutoff>
/sample/neighbours/<string:guid>/<string:reference>/<int:threshold>/<string:method>