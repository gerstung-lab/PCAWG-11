# Code accompanying _The evolutionary history of 2,658 cancers_

**Moritz Gerstung and Santiago Gonzalez, on behalf of the PCAWG-11 Evolution and Heterogeneity Working Group.**

This repository contains analysis code for the aforementioned publication, especially, Figures 1b,e-h, 2, 4, 5 & Extended Data Figures 3, 6, 8 and 9. The output of the code 
is an `rmarkdown` html report, hosted at <https://gerstung-lab.github.io/PCAWG-11>.

Note that some referenced files in `ref` and `final` folders are protected access and are note part of this repository.

## Contents

The `code` folder contains a number of scripts to run `MutationTime.R` to estimate the timing of point mutations and copy number gains in individual samples. This step can be parallelised across individual samples.

* `code/PCAWG-functions.R`
A collection of utility functions. This is sourced before running. This functions also sources a number of working group output files, 
from `ref` and `final` folder, which are not part of the repository at the moment.
* `code/VCF-annotate.R`
An R script to annotate point mutations and estimate the timing of copy number gains. It is usually envoked via 
``` 
$ bsub < VCF-annotate.sh
```` 
which submits an array job to an LSF cluster.
Input to this script are (assumed in `final`, not provided):
	* consensus copy number segments
	* consensus subclonal reconstruction 
	* VCF files for SNVs and indels
	* inferred gender

	The run time ranges from 5 minutes to about 8 hours, depending on the number of variants per sample (median about 10 minutes).
	
* The main code file is `code/PCAWG-final.R`; this is run on R-3.3.3 using
```{R}
> rmarkdown::render("code/PCAWG-final.R", clean=FALSE)
```
The output is a html report with mixed code, results and figures, which can be found at `run/PCAWG-final.html`. Runtime is approximately 4hr, with most time spent on loading all data.
It requires about 50Gb of memory per core. 
* Experimental code can be found in `PCAWG-sandbox.R`.
* Code used for timing simulations is in  `PCAWG-simulate.R`.
* `modules/MutationTime.R` contain the MutationTimeR submodule. Don't forget to run `git submodule init`.

## Dependencies

The main analysis workflow is run in

```
R-3.3.3
```

Additionally it requires the following R libraries
```
Rsamtools
VariantAnnotation
MutationTimeR <https://github.com/gerstung-lab/MutationTimeR>
rstan
rmarkdown
```

## Output
Main output is `run/PCAWG-final.html`, hosted at <https://gerstung-lab.github.io/PCAWG-11/> in addition to this repository.

