# exonerate_protein_prediction.py 

'''
Predict proteins from a genome using exonerate and protein predictions from a transcriptome.

Results are in results/ (.proteins.fasta)

usage: python exonerate_protein_prediction.py [proteins.fasta] [genome.fasta] [threads] [genetic code (int)]

dependencies:
# exonerate
conda install -c bioconda exonerate
# blast
conda install -c bioconda blast
# parallel
conda install -c conda-forge parallel
# cd-hit
conda install -c bioconda cd-hit
'''

import subprocess
from glob import glob
import sys
from statistics import median
from datetime import datetime

print('\nPredicting proteins for genome: ' + sys.argv[2] + '\nUsing proteins from: ' + sys.argv[1])
print('Start time: ' + str(datetime.now()))
subprocess.call('mkdir temp/', shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

# map proteins to genome using tblastn
print('\nMapping proteins to genome using tBLASTn...')
# Make blast databases
subprocess.call('makeblastdb -in ' + sys.argv[2] + ' -out temp/' + sys.argv[2] + '.db -dbtype nucl', shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
# find homologs using tblastn
subprocess.call('tblastn -query ' + sys.argv[1] + ' -db temp/' + sys.argv[2] + '.db -max_target_seqs 100 -num_threads ' + sys.argv[3] + ' -evalue 1e-5 -outfmt 6 -out temp/genome_mapping.tblastn -db_gencode ' + sys.argv[4], shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


# go through blast outputs and identify homology regions
print('Parsing blast output...')

'''
Record query and each hit, joining the hsps into a full coding region.
Note each hit for each query and if there are mutliple hits on the 
same scaffold treat them as different hits if they are far away
from the previous hsp (2000 bp - max gap size is 2000 bp).
'''

# adjust maximum intron size here
max_gap_size = 2000

# load in blast files
blast = open('temp/genome_mapping.tblastn','r').readlines()
blast_d = {} # query_m: [scaffold, start, stop, strand]

def strand_check(blast_line):
    if int(blast_line.split('\t')[8]) > int(blast_line.split('\t')[9]):
        strand = '-'
    else:
        strand = '+'
    return strand

def start_stop(blast_line, strand):
    if strand == '+':
        start = int(blast_line.split('\t')[8])
        stop = int(blast_line.split('\t')[9])
    else:
        start = int(blast_line.split('\t')[9])
        stop = int(blast_line.split('\t')[8])
    return [start, stop]

# record initial hit        
q = blast[0].split('\t')[0]
h = blast[0].split('\t')[1]
prev_strand = strand_check(blast[0])
ss = start_stop(blast[0],prev_strand)
prev_start, prev_stop = ss[0], ss[1]
m = 1
blast_d[q+'_'+str(m)] = [h,prev_start,prev_stop,prev_strand]

# go through blast output file
for line in blast[1:]:
    strand = strand_check(line)
    ss = start_stop(line,strand)
    start, stop = ss[0], ss[1]
    if (line.split('\t')[0] == q) and (line.split('\t')[1] == h) and (prev_strand == strand):
        if (start > prev_stop) and (start-prev_stop < max_gap_size):
            blast_d[q+'_'+str(m)][2] = stop
        elif (stop < prev_start) and (prev_start-stop < max_gap_size):
            blast_d[q+'_'+str(m)][1] = start
        else:
            m += 1
            q, h, prev_start, prev_stop, prev_strand = line.split('\t')[0], line.split('\t')[1], start, stop, strand
            blast_d[q+'_'+str(m)] = [line.split('\t')[1],prev_start,prev_stop,prev_strand]
    elif (line.split('\t')[0] == q) and (line.split('\t')[1] == h) and (prev_strand != strand):
        q, h, prev_start, prev_stop, prev_strand = line.split('\t')[0], line.split('\t')[1], start, stop, strand
        blast_d[q+'_'+str(m)] = [line.split('\t')[1],prev_start,prev_stop,prev_strand]    
    elif (line.split('\t')[0] == q) and (line.split('\t')[1] != h):
        q, h, prev_start, prev_stop, prev_strand = line.split('\t')[0], line.split('\t')[1], start, stop, strand
        blast_d[q+'_'+str(m)] = [line.split('\t')[1],prev_start,prev_stop,prev_strand]
    elif (line.split('\t')[0] != q):
        m = 1
        q, h, prev_start, prev_stop, prev_strand = line.split('\t')[0], line.split('\t')[1], start, stop, strand
        blast_d[q+'_'+str(m)] = [line.split('\t')[1],prev_start,prev_stop,prev_strand]

# write out all sequences for exonerate
subprocess.call('mkdir temp/sequences/', shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

# proteins
print('Extracting sequences...')
proteins = open(sys.argv[1],'r').read().split('>')[1:]
protein_d = {}
for p in proteins:
    protein_d[p.split('\n')[0]] = p.split('\n',1)[1].replace('\n','')
for p in blast_d:
    out = open('temp/sequences/'+p+'.fasta','w')
    out.write('>'+p+'\n'+protein_d[p.rsplit('_',1)[0]]+'\n')
    out.close()

# extract regions mapping to genome scaffolds
def reverse_complement(seq):
    complement = ''
    for b in seq.upper()[::-1]:
        if b == 'A':
            complement += 'T'
        elif b == 'T':
            complement += 'A'
        elif b == 'C':
            complement += 'G'
        elif b == 'G':
            complement += 'C'
    return complement

genome = open(sys.argv[2],'r').read().split('>')[1:]
genome_d = {}
homology_d = {}
for s in genome:
    genome_d[s.split('\n')[0]] = s.split('\n',1)[1].replace('\n','')
for i in blast_d:
    out = open('temp/sequences/'+i+'.homology.fasta','w')
    if blast_d[i][-1] == '+':
        out.write('>'+i+'_'+blast_d[i][0]+'\n'+genome_d[blast_d[i][0]][blast_d[i][1]-1:blast_d[i][2]]+'\n')
        homology_d[i] = genome_d[blast_d[i][0]][blast_d[i][1]-1:blast_d[i][2]]
    else: # if on the negative strand output the reverse complement
        out.write('>'+i+'_'+blast_d[i][0]+'\n'+reverse_complement(genome_d[blast_d[i][0]][blast_d[i][1]-1:blast_d[i][2]])+'\n')
        homology_d[i] = reverse_complement(genome_d[blast_d[i][0]][blast_d[i][1]-1:blast_d[i][2]])
    out.close()

# run exonerate
print('Using Exonerate to identify exons...')
subprocess.call("find temp/sequences/ -type f -name '*.homology.fasta' | awk -F '.homology.fasta' '{print $1}' | cut -f 3 -d '/' | parallel -j " + sys.argv[3] + " 'exonerate -m p2g  --geneticcode " + sys.argv[4] + " --maxintron 20000 --score 100 --showtargetgff -q temp/sequences/{}.fasta -t temp/sequences/{}.homology.fasta | egrep -w exon > temp/sequences/{}.exon.gff'", shell = True)

# extract coding regions and combine gff files into one final file
print('Extracting coding regions...')
subprocess.call('mkdir results/', shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
out = open('results/'+sys.argv[2]+'.cds.fasta','w')

n = 1
for seq in blast_d:
    gff = open('temp/sequences/'+seq+'.exon.gff','r').readlines()
    if len(gff) == 0:
        continue
    else:
        sequence = ''
        for hit in gff:
            sequence += homology_d[seq][int(hit.split('\t')[3])-1:int(hit.split('\t')[4])]
        out.write('>seq'+str(n)+'.'+seq+'_'+blast_d[seq][0]+'\n'+sequence+'\n')
        n += 1
out.close()

# translate proteins using fastatranslate
print('Translating coding regions...')
subprocess.call('fastatranslate -f results/'+sys.argv[2]+'.cds.fasta --geneticcode ' + sys.argv[4] + ' > temp/translation.fasta',shell = True,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
# clean headers
translations = open('temp/translation.fasta','r').read()
out = open('temp/translation.fasta','w')
out.write(translations.replace(' ','_').replace(':','_').replace('[','').replace(']','').replace('(','').replace(')','')) 
out.close()

# blast proteins to initial dataset
print('Comparing predicted proteins to original proteome...')
subprocess.call('diamond makedb --in ' + sys.argv[1] + ' --db temp/proteins.db', shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
subprocess.call('diamond blastp -c1 --query temp/translation.fasta --max-target-seqs 1 --db temp/proteins.db --outfmt 6 --sensitive --evalue 1e-5 --out temp/translation.blastp --threads ' + sys.argv[3], shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
# record blast results
blast = open('temp/translation.blastp', 'r').readlines()
protein_blast_d = {}
for line in blast:
    protein_blast_d[line.split('\t')[0]] = [line.split('\t')[1], float(line.split('\t')[-2])]

# identify the best translation per protein (longest and with best blast hit)
translation = open('temp/translation.fasta','r').read().split('>')[1:]
translation_d = {}
for seq in translation:
    translation_d[seq.split('\n')[0]] = seq.split('\n',1)[1].replace('\n','')
    
seq_d = {}
for p in translation_d:
    try:
        protein_blast_d[p]
        if (seq_d[p.split('.')[0]][1] == 'NA'):
            seq_d[p.split('.')[0]] = [translation_d[p],protein_blast_d[p][0],protein_blast_d[p][1],len(translation_d[p].replace('*',''))]
        elif (seq_d[p.split('.')[0]][2] > protein_blast_d[p][1]):
            seq_d[p.split('.')[0]] = [translation_d[p],protein_blast_d[p],protein_blast_d[p][1],len(translation_d[p].replace('*',''))]
    except:
        try:
            if (seq_d[p.split('.')[0]][1] != 'NA'):
                continue
            elif (seq_d[p.split('.')[0]][3] < translation[p][1]):
                seq_d[p.split('.')[0]] = [translation_d[p],protein_blast_d[p][0],protein_blast_d[p][1],len(translation_d[p].replace('*',''))]
        except:
            seq_d[p.split('.')[0]] = [translation_d[p],'NA','NA',len(translation_d[p].replace('*',''))]
    
# output the predicted proteins and cluster at 100%
print('Writing finished proteins...')
out = open('results/'+sys.argv[2]+'.proteins.fasta','w')
for seq in seq_d:
    sequence = max(seq_d[seq][0].split('*')) # split at stop codons - take the longest uninterupted sequence
    if (len(sequence) > 50) and (seq_d[seq][1] != 'NA'):
        out.write('>'+seq+'\n'+sequence+'\n')
out.close()

# cluster genes at 99%
print('Clustering resulting proteins at 99%...')
subprocess.call('cd-hit -i results/'+sys.argv[2]+'.proteins.fasta -c 0.99 -o temp/cluster -M 4000 -T ' + sys.argv[3], shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
subprocess.call('mv temp/cluster results/'+sys.argv[2]+'.proteins.fasta',shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

subprocess.call('rm -r temp/',shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
print('\nFinish time: ' + str(datetime.now()))
