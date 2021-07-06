# FAST SITE REMOVAL
#
# Given a rate file generated in IQ-Tree (-wsr option) and an alignment, 
# remove the X proportion of fastest (highest rate) sites in the alignment.
#
# Usage: python fast_site_removal.py [alignment] [ratefile] [percent_of_fastest_sites_to_remove: 0 - 1]


import sys
import operator

# load in rate file
rf = open(sys.argv[2],'r').readlines()

rate_d = {}
for line in rf:
    if (line.startswith('#') == False) and (line.startswith('Site') == False):
        rate_d[line.split('\t')[0]] = float(line.split('\t')[1])

# sort the dictionary

sorted_rate_d = sorted(list(rate_d.items()), key=operator.itemgetter(1))

# load in the alignment

aln = open(sys.argv[1],'r').read().split('>')

aln_d = {}

for seq in aln[1:]:
    head = seq.split('\n')[0]
    sequence = []
    for part in seq.split('\n')[1:]:
        for base in part.strip():
            sequence.append(base)
    aln_d[head] = sequence

# how many sites do you want to see?

perc_sites = 1-float(sys.argv[3])

# output the filtered alignment

total_sites = len(rf)-1
number_sites_to_extract = round(total_sites * perc_sites)
sites_to_extract = []

n = 0
while n < number_sites_to_extract:
    sites_to_extract.append(sorted_rate_d[n][0])
    n += 1
    
sites_to_extract = [int(x) for x in sites_to_extract]
sites_to_extract.sort()
sites_to_extract = [str(x) for x in sites_to_extract]

# write out the filtered alignment
out = open(sys.argv[2]+'.'+str(float(sys.argv[3])*100)+'fsr','w')
for seq in list(aln_d.keys()):
    out.write('>'+seq+'\n')
    for s in sites_to_extract:
        out.write(aln_d[seq][int(s)-1])
    out.write('\n')
out.close()
