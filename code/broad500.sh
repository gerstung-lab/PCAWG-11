#!/bin/bash
#BSUB -q normal
#BSUB -J broad500[1-500]
#BSUB -o log/broad500-%J-%I.out
#BSUB -e log/broad500-%J-%I.err
#BSUB -R "span[hosts=1] select[mem>4800] rusage[mem=4800]"
#BSUB -M 4800
#BSUB -n 1

XMIN=3
RHO=0

INPUT_FOLDER="../broad500/VCF"
OUTPUT_FOLDER="../broad500/annotated_vcf_rho$RHO_xmin$XMIN"
OVERWRITE=true

mkdir -p $OUTPUT_FOLDER

FILES=(`ls $INPUT_FOLDER/*.vcf`)
INPUT=${FILES[(($LSB_JOBINDEX-1))]}
echo $INPUT
STEM=`basename $INPUT | sed s/.vcf//g`
OUTPUT="$OUTPUT_FOLDER/$STEM"
if [ ! -f "$OUTPUT_FOLDER/$STEM.complete_annotation.vcf.bgz" ] || [ "$OVERWRITE" = true ]; then
Rscript broad500.R $INPUT $OUTPUT_FOLDER/$STEM.complete_annotation.vcf $XMIN $RHO
else
echo "$STEM.output exists. skipping."
fi

