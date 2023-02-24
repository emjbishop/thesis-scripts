#!/bin/bash

IN="/home/ebishop/thesis_projects/aim2_hybran/rnaseq/rnaseq_gitlab/data_raw"

# Run this script within the 'tpm-calc' conda environment

# 0. Move into output directory
cd /home/ebishop/thesis_projects/aim2_hybran/rnaseq/rnaseq_gitlab/data_intermed

# 1. Trim and filter the raw reads
echo "Trimming and filering raw reads"
trimmomatic PE $IN/ERR3039930_1.fastq.gz $IN/ERR3039930_2.fastq.gz 1_paired1.fq.gz 1_unpaired1.fq.gz 1_paired2.fq.gz 1_unpaired2.fq.gz ILLUMINACLIP:/home/ebishop/anaconda/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 

# 2. Build a GMAP index (will be used by GSNAP)
echo "Building GMAP index"
gmap_build -D . -d 2_gmap_index_H37Rv $IN/H37Rv.fasta.gz -c NC_000962.3 -g

# 3. Align the reads
echo "Aligning reads"
gsnap -A sam -o 3_gsnap.sam --gunzip -D . -d 2_gmap_index_H37Rv 1_paired1.fq.gz 1_paired2.fq.gz 

# 4. Sort the SAM file
echo "Sorting SAM file"
samtools sort -o 4_sorted.sam 3_gsnap.sam 

# 5. Use AGAT to convert Hybran's output GFF3 to GTF
echo "Converting GFF to GTF"
agat_convert_sp_gff2gtf.pl --gff $IN/H37Rv-NCBI_nofasta.gff -o 5_H37Rv.gtf
# Change chromosome name to match SAM file
sed 's/L_1/1/g' 5_H37Rv.gtf > 5_H37Rv_chr1.gtf 

# 6. Convert input SAM file to BAM
echo "Converting SAM to BAM"
samtools view -b 4_sorted.sam > 6_sorted.bam

# 7. Use TPMCalculator to get transcripts per million (TPM) scores for each feature
echo "Running TPMCalculator"
TPMCalculator -g 5_H37Rv_chr1.gtf -b 6_sorted.bam -k locus_tag

