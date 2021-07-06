# Phylogenomic-analysis

Useful scripts for conducting phylogenomic analyses.

## Overview

A series of scripts that are useful for analyzing phylogenomic data (e.g., working with phylogenomic trees and comparative genomics).

## Scripts

### fast_site_removal.py

This script will remove a given proportion of the fastest evolving sites in an alignment. Site rates should be estimated using IQ-Tree (https://github.com/Cibiv/IQ-TREE) using the -wsr option (produces a .rate file). Removal of fast evolving sites can be useful for assessing long branch attraction and assessing the phylogenetic support for a given topology.

```
# remove the top 25% fastest evolving sites in an alignment
python fast_site_removal.py ProteinA.fasta.aln ProteinA.fasta.aln.rate 0.25
```

