#!/bin/bash
#BSUB -q normal
#BSUB -J trinuc[1-269]
#BSUB -o log/trinuc-%J-%I.out
#BSUB -e log/trinuc-%J-%I.err
#BSUB -R "span[hosts=1] select[mem>2400] rusage[mem=2400]"
#BSUB -M 2400
#BSUB -n 1


Rscript trinuc2.R $LSB_JOBINDEX

