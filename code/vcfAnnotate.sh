#!/bin/bash
#BSUB -q normal
#BSUB -J vcfAnnotate[1-2266]
#BSUB -o log/vcfAnnotate-%J-%I.out
#BSUB -e log/vcfAnnotate-%J-%I.err
#BSUB -R "span[hosts=1] select[mem>4800] rusage[mem=4800]"
#BSUB -M 4800
#BSUB -n 1

INPUT_FOLDER="../subs/2016-03/sanger"
OUTPUT_FOLDER="../subs/2016-03/annotated"
OVERWRITE=false

FILES=(`ls $INPUT_FOLDER/*mnv.vcf`)
INPUT=${FILES[(($LSB_JOBINDEX-1))]}
echo $INPUT
STEM=`basename $INPUT | sed s/.vcf//g`
OUTPUT="$OUTPUT_FOLDER/$STEM"
if [ ! -f "$OUTPUT.vagrent.vcf" ]; then
source /lustre/scratch112/sanger/cgppipe/PanCancerFinal/final.bash.setup
AnnotateVcf.pl -i $INPUT -o $OUTPUT.vagrent.vcf -c /lustre/scratch112/sanger/kr2/PanCancerFinal/ref/vagrent/e74/Homo_sapiens.GRCh37.74.vagrent.cache.gz -sp Homo_sapiens -as GRCh37
fi
if [ ! -f "$OUTPUT_FOLDER/$STEM.complete_annotation.vcf.bgz" ] || [ "$OVERWRITE" = true ]; then
Rscript vcfAnnotate.R $OUTPUT.vagrent.vcf
else
echo "$STEM.output exists. skipping."
fi

