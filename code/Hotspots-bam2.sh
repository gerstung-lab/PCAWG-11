#!/bin/bash
#BSUB -q normal
#BSUB -J hotspot[1-2961]
#BSUB -o log/hotspot-%J-%I.out
#BSUB -e log/hotspot-%J-%I.err
#BSUB -R "span[hosts=1] select[mem>4800] rusage[mem=4800]"
#BSUB -M 4800
#BSUB -n 1

/software/R-3.3.0/bin/Rscript Hotspots-bam2.R $LSB_JOBINDEX


