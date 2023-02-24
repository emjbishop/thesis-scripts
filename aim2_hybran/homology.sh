#!/bin/bash

# User provided number of cores for BLAST
CORES=$1

# Project directory
PROJ="/home/ebishop/thesis_projects/aim2_hybran/homology"

# Annotation directory
ANOT="/grp/valafar/data/annotation/hybran"

# BLAST database paths
DBX="/grp/valafar/data/depot/hybran.pub/blast_nr_db/nr"
DBN="/grp/valafar/data/depot/hybran.pub/blast_nt_db/nt"

# 1. Get all the Hybran GBK file paths, excluding H37Rv
ls $ANOT/*.gbk | grep -v H37Rv > $PROJ/gbk_fofn

# 2. Get all the novel gene names from all the isolates (they're of the format ORFXXXX)
# Exclude standard Rv duplicates (ORF0001 - ORF0011)
echo "Getting all the novel gene names"
grep -u 'gene="ORF\+' $ANOT/*.gbk | awk -F '"' '{print $2}' | sort -u | awk '{if (NR>11) print}' > $PROJ/orfs_list

# 3. Generate a multifasta of the novel gene nucleotide sequences
echo "Creating the novel gene multifasta"
python make_multifasta.py

# 4. Run blastx and blastn (adjust num_threads if needed), keeping only the top first matches
# Output is tab separated: novel gene id, match id, e-value, bit score
echo "Running BLASTX"
blastx -db $DBX -query $PROJ/novel_orfs.fasta -out $PROJ/blastx.out -outfmt "6 qseqid sseqid evalue bitscore" -max_target_seqs 1 -num_threads $CORES

echo "Running BLASTN"
blastn -db $DBN -query $PROJ/novel_orfs.fasta -out $PROJ/blastn.out -outfmt "6 qseqid sseqid evalue bitscore" -max_target_seqs 1 -num_threads $CORES
