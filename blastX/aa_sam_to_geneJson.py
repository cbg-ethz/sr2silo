# Script to convert SAM to json gene, sequence insertions
from __future__ import annotations

## read in each line in the sam file

### if the read_id is not in the dictionary, add it

###### make a dictionary for AA_sequences for the gene

###### make a dictionary for insertions for the gene

### get the gene name from the reference name

### get the cleartext sequence from the aligned query sequence and the cigar string

### get the insertions from the aligned query sequence and the cigar string

### add the gene name: sequence, and insertions to the dictionaries

## write the dictionaries to a json file



