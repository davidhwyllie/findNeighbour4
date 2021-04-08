""" reads the covid genome mask here:
 https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf
 https://github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/problematic_sites_sarsCov2.vcf

 and outputs in in the zero-indexed format required by findNeighbour4

https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473/14

 Excludes the 'mask' and 'caution' bases
"""

import os


inputfile = os.path.join("..","reference",'covid_exclusion.vcf')
outputfile = os.path.join('..','reference','covid-exclude.txt')

with open(inputfile,'rt') as f:
    lines = f.readlines()

nMasked =0 
with open(outputfile, 'wt') as f:
    for line in lines:
        if not line[0] == '#':
            pos = int(line.split('\t')[1])-1
            nMasked +=1
            f.write('{0}\n'.format(pos))
print("Complete.  Wrote {0} sites to {1}".format(nMasked, outputfile))



