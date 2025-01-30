#!/bin/bash

# ==== Database ====
diamond makedb --in ref/sequences_assembley.fasta -d ref/hxb_pol_db


# ==== Alignment ====
diamond blastx -d ref/hxb_pol_db \
		-q output.fastq \
		-o diamond_blastx.sam \
		--evalue 1 \
		--gapopen 6 \
		--gapextend 2 \
		--outfmt 101 \
		--matrix BLOSUM62 \
		--unal 1 \
		--max-hsps 1 \
		--more-sensitive
