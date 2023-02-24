#!/bin/sh

fp=$GROUPHOME/data/depot/database/intermediate_files

awk 'BEGIN{OFS=","} NR==FNR{Arr[$0]++;next}{if($1 in Arr){print $1, $2, $3, $4, $5, $6, $7, $8, $10}}' $fp/isolates.csv $GROUPHOME/metadata/dst.txt > $fp/isolate_dst.csv
