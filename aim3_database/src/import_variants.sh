#!/bin/bash

file_path=$GROUPHOME/data/variants/pbhoover

for f in $file_path/*.vcf.gz; do
   [ -f "$f" ] || continue
   python vep_vcf_parser.py $f
done
