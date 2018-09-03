Integration testing
-------------------

Detecting mixtures
--------------------
Generates a series of phylogenies and sequences; simulates mixtures between members of this phylogeny.
First, a server must be running:  
```python findNeighbour3-server.py ../demos/simulation/config/config.json ```

Data is simulated    
```python make_simulation.py 50 10000 0.02 0.01 0.05 0.1 20 1 ../output/simulation_set_1  ```

and the simulated data is fed to the findNeighbour3 server.  Responses are recorded.    
```python run_simulation.py  ../output/simulation_set_1  ```

Outputs are written to disc as excel files in ../output/simulation_set_1 in the above example.
An R script to summarise outputs is available:
``` STILL TO DO ```
