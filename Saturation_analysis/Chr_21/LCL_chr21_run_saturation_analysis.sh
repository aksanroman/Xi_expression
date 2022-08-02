#!/bin/bash

#Run LCL saturation analysis for Chr21

bsub Rscript LCL_chr21_SaturationAnalysis.R 45 1
bsub Rscript LCL_chr21_SaturationAnalysis.R 40 100
bsub Rscript LCL_chr21_SaturationAnalysis.R 35 100
bsub Rscript LCL_chr21_SaturationAnalysis.R 30 100
bsub Rscript LCL_chr21_SaturationAnalysis.R 25 100
bsub Rscript LCL_chr21_SaturationAnalysis.R 20 100
bsub Rscript LCL_chr21_SaturationAnalysis.R 15 100
bsub Rscript LCL_chr21_SaturationAnalysis.R 10 100