# Demonstrations of server functioning


## Using real data
Due to their size, the real data is not included in github.  It can at present be downloaded from the below links:

[AC587 (a small dataset, < 100 samples)](https://www.dropbox.com/s/3ohmf475aa9enxn/AC587.zip?dl=0)  
[AA041 (a larger dataset, ~ 1000 samples)](https://www.dropbox.com/s/2jk0vx2gax1q6nt/AA041.zip?dl=0)  
[England201618 (a >7000 sample dataset described in the publication *M. tuberculosis microvariation is common and is associated with transmission:  analysis of three years prospective universal sequencing in England*)](https://www.dropbox.com/s/5y5d0m7qpbhbzxl/england201618.tar?dl=0)

Download the zip files, unzip them, and place them in /demos.
  
The input for the relatedness server are fasta files derived from reference mapped, consensus base called data.
The sample data provided here is from the Public Health England bioinformatics pipeline used for TB processing,
which is freely available at https://github.com/oxfordmmm/CompassCompact.


### AC587  
a collection of 43 mapped samples containing related TB isolates, as well as unrelated controls TB samples.
The latter are added before the 43 related samples, as they are used by the server to estimate expected N frequencies in real data.
To run the demo:
- make sure mongodb is running
- from the src directory  
-- start the server  
``` pipenv run python findNeighbour4-server.py ../demos/AC587/config/config.json ```  
-- run the software adding samples to the server  
``` pipenv run python demo_ac587.py ```

### AA041
a larger collection of ~ 1000 mapped samples containing related TB data.
To run the demo:
- make sure mongodb is running
- from the src directory  
-- start the server  
``` pipenv run python findNeighbour4-server.py ../demos/AA041/config/config.json ```  
-- run the software adding samples to the server  
``` pipenv run python demo_aa041.py ```

### England201618
a collection of over 7,000 TB samples from England 2016 to 2018.  This data set includes  
i) .csv files containing the positions of [where different bases map with high quality to the same site across each genome](https://github.com/davidhwyllie/VCFMIX)  
ii) fasta files in which these positions are marked with [IUPAC codes](https://www.bioinformatics.org/sms/iupac.html).  
These fasta files can be loaded by scripts similar to those above.
-- the script demo_pheall.py is suitable for doing this.

-- if not already running start the server with 
```./fn4_startup.sh ../demos/phe_dev/config_phe_dev_nocl.json```		
Note: (1) this config file does not do clustering; (2) it will run processes in the background using pipenv.  

-- edit demo_pheall.py to change the fastadir to the directory in which they are contained  
-- then submit the files to the server     
``` pipenv run python3 demo_pheall.py ```  

 



