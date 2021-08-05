# concatenation.py
'''
Concatenate alignments together for phylogenomic analyses.
Sequence headers must start with the species identifier followed by
a period (e.g. >Homo_sapiens.PROTEINID).

Output is a concatenated alignment file (concatenation.DATE.fasta) and
an overview of gene/site presence per species (.species_stat.tab file).

Usage: python concatenation.py '*fasta.aln'
'''

from glob import glob
import sys
from datetime import datetime

# record all species
species = []
for fname in glob(sys.argv[1]):
    fasta = open(fname,'r').read().split('>')[1:]
    species += [s.split('.')[0] for s in fasta]

species_d = {}
counts_d = {}

for sp in list(set(species)):
    species_d[sp] = ''
    counts_d[sp] = 0

# add sequence information
for aln in glob(sys.argv[1]):
    aln_d = {}
    aln_file = open(aln,'r').read().split('>')[1:]
    length = len(aln_file[0].split('\n',1)[1].replace('\n',''))
    for seq in aln_file:
        aln_d[seq.split('.')[0]] = seq.split('\n',1)[1].replace('\n','')
    for sp in species_d:
        try:
            species_d[sp] += aln_d[sp]
            counts_d[sp] += 1
        except:
            species_d[sp] += ('-'*length)

# output concatenated alignment
lengths = []
for sp in species_d:
    lengths.append(len(species_d[sp]))
if len(set(lengths)) > 1:
    print('\nError: there can only be one sequence per species in each alignment\n')
else:
    out = open('concatenation.'+str(datetime.now()).split(' ')[0].replace('-','_')+'.fasta','w')
    for sp in species_d:
        out.write('>'+sp+'\n'+species_d[sp]+'\n')
    out.close()

# output gene presence/absence statistics
total_genes = len(glob(sys.argv[1]))
total_length = len(species_d[sp])

out = open('concatenation.'+str(datetime.now()).split(' ')[0].replace('-','_')+'.species_stats.tab','w')
out.write('species\t%genes\t%sites\n')
for sp in species_d:
    out.write(sp+'\t'+str(round(100*(float(counts_d[sp])/total_genes),2))+'\t'+str(round(100*((total_length - float(species_d[sp].count('-')))/total_length),2))+'\n')
out.close()
    
