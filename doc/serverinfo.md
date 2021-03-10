Server status information
=========================
When operating, the server will report various statistics about its operation.  
To do so, it keeps a log, recording disc and memory usage before and after relevant events, including  
- sample insertion  
- memory repacking  
- database repacking  
- clustering

#### Detailed report
[Show detailed report](/api/v2/monitor)  (may renders slowly initially)

#### Memory usage now
[Show current server memory usage](/api/v2/server_memory_usage/1/html) 

### Key database stats now
[Show key database stats](/api/v2/server_database_usage/1)   
Ratio of sequences to rows holding sequence-sequence pairs (storage_ratio) is a key statistic.   
Values > 30 are associated with a 'write optimised' format and are characterised by large databases and slow recovery.

