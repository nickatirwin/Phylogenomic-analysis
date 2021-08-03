# Phylogenomic-analysis

A series of scripts that are useful for analyzing phylogenomic data (e.g., working with phylogenomic trees and comparative genomics).

## Scripts

### exonerate_protein_prediction.py

Predict proteins from a genome using a given proteome for reference (e.g. proteome of a closely related species, transcriptome, or a set of genes). The gene models can be rough but may provide a useful starting place and is useful for estimating genome completeness or identifying genes for phylogenomic analyses. 

In brief, proteins are mapped to the genome using tBLASTn (e < 1e-5), mapped regions are extracted, query proteins are used as a model for Exonerate to identify exons, exons are combined into coding regions, and coding regions are translated. The resulting protein predictions are compared against the original proteome and protein models with a hit back to the original dataset are retained and clustered at 99% to reduce redundancy. The pipeline is fairly quick and tends to run in 10-30 minutes with 60 threads for the datasets I've tested.

Required dependencies:
```
# exonerate
conda install -c bioconda exonerate
# blast
conda install -c bioconda blast
# parallel
conda install -c conda-forge parallel
# cd-hit
conda install -c bioconda cd-hit
```
Usage:
```
python exonerate_protein_prediction.py [proteome.fasta] [genome.fasta] [threads] [genetic code (integer)]
```
### fast_site_removal.py

This script will remove a given proportion of the fastest evolving sites in an alignment. Site rates should be estimated using IQ-Tree (https://github.com/Cibiv/IQ-TREE) using the -wsr option (produces a .rate file). Removal of fast evolving sites can be useful for assessing long branch attraction and assessing the phylogenetic support for a given topology.

Usage:
```
python fast_site_removal.py [alignment file] [.rate file] [proportion of sites to remove]
```
