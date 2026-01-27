# Processing Module

Core functions for BAM processing, translation, and alignment.

## Data Structures

### AlignedRead

::: sr2silo.process.AlignedRead

### Gene

::: sr2silo.process.Gene

### GeneSet

::: sr2silo.process.GeneSet

### NucInsertion

::: sr2silo.process.NucInsertion

### AAInsertion

::: sr2silo.process.AAInsertion

### AAInsertionSet

::: sr2silo.process.AAInsertionSet

## BAM Conversion

::: sr2silo.process.bam_to_sam

::: sr2silo.process.sam_to_bam

::: sr2silo.process.bam_to_fasta_query

::: sr2silo.process.sort_bam_file

::: sr2silo.process.sort_and_index_bam

::: sr2silo.process.sort_sam_by_qname

::: sr2silo.process.get_gene_set_from_ref

## Translation and Alignment

::: sr2silo.process.nuc_to_aa_alignment

::: sr2silo.process.parse_translate_align

::: sr2silo.process.parse_translate_align_in_batches

::: sr2silo.process.curry_read_with_metadata

## Read Merging

::: sr2silo.process.paired_end_read_merger

## Exceptions

::: sr2silo.process.ZeroFilteredReadsError
