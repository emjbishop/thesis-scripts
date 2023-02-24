#!/bin/bash

# run  like this:
# >cat isolates.txt | create_merged_vcf.sh
# where isolates.txt is a list of isolate IDs of desired isolates with vcfs to create the merged.vcf from
conda activate gwa
mkdir vcfs
isolates="-"
# Edit public-genomes/ VCFs
for isolate in $(cat $isolates)
do
	if [ -f $GROUPHOME/data/variants/public-genomes/$isolate.vcf.gz ]
	then
		vcfgz=$isolate.vcf.gz
		vcf=$isolate.vcf
		zcat $GROUPHOME/data/variants/public-genomes/$vcfgz | \
		sed 's/INFO\t/INFO\tFORMAT\t/; s#|$#|\tGT\t1#; s#CSQ=$#CSQ=\tGT\t1#' > vcfs/$vcf
		bgzip -f $vcf
		tabix -p vcf $vcfgz
	elif [ -f $GROUPHOME/data/variants/pbhoover/$isolate.vcf.gz ]
	then
		vcfgz=$isolate.vcf.gz
		vcf=$isolate.vcf
		zcat $GROUPHOME/data/variants/pbhoover/$vcfgz | \
		sed 's/INFO\t/INFO\tFORMAT\t/; s#|$#|\tGT\t1#;s#CSQ=$#CSQ=\tGT\t1#; s/HETERO/PASS/; 10i ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' > vcfs/$vcf
		bgzip -f $vcf
		tabix -p vcf $vcfgz
	fi
done
bcftools merge -m none -0 -O z vcfs/*.vcf.gz > merged.vcf.gz 
rm -rf vcfs
