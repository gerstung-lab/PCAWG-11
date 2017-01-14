#!/bin/bash
#BSUB -q normal
#BSUB -J vcfAnnotateFinal[1-2703]
#BSUB -o log/vcfAnnotateFinal-%J-%I.out
#BSUB -e log/vcfAnnotateFinal-%J-%I.err
#BSUB -R "span[hosts=1] select[mem>4800] rusage[mem=4800]"
#BSUB -M 4800
#BSUB -n 1

INPUT_FOLDER="../final/final_consensus_12oct_passonly/snv_mnv"
OUTPUT_FOLDER="../final/annotated_002/snv_mnv"
OVERWRITE=true

FILES=(`ls $INPUT_FOLDER/*.vcf.gz`)
INPUT=${FILES[(($LSB_JOBINDEX-1))]}
echo $INPUT
STEM=`basename $INPUT | sed s/.vcf.gz//g`
ID=`echo $STEM | awk '{x=$0; gsub("\\\..*","",x); print x}'`
DPFILE="../dp/20161213_vanloo_wedge_consSNV_prelimConsCNAallStar/2_subclones/$ID_subclonal_structure.txt.gz"
if [ ! -f "$DPFILE" ]; then
echo "No DP file. Exit."
exit 0
fi
OUTPUT="$OUTPUT_FOLDER/$STEM"
if [ ! -f "$OUTPUT_FOLDER/$STEM.complete_annotation.vcf.bgz" ] || [ "$OVERWRITE" = true ]; then
Rscript vcfAnnotateFinal.R $INPUT $OUTPUT_FOLDER/$STEM.complete_annotation.vcf
else
echo "$STEM.output exists. skipping."
fi

