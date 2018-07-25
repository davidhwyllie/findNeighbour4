# Demonstrations of the server using real and simulated data

* AC587 - a collection of 43 mapped samples containing TB data, as well as 38 control TB samples.  The latter are added before the 43 related samples, as they are used by the server to estimate expected N frequencies in real data.
To run the demo:
- make sure mongodb is running
- from the src directory  
-- start the server  
``` python findNeighbour3-server.py ../demos/AC587/config/config.json ```  
-- run the software adding samples to the server  
``` python demo_ac587.py ```

TODO: depict combinations of clusters over time with both algorithms
check that mixed samples go into both.  ** IMPORTANT ** haven't observed this happening, except in unit tests
a gif with the clusters on the L and the MSA on the R might work