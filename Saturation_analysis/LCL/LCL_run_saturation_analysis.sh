#!/bin/bash

#Run LCL saturation analysis

bsub Rscript LCL_NPX_PAR_SaturationAnalysis.R 106 1
bsub Rscript LCL_NPX_PAR_SaturationAnalysis.R 100 100
bsub Rscript LCL_NPX_PAR_SaturationAnalysis.R 90 100
bsub Rscript LCL_NPX_PAR_SaturationAnalysis.R 80 100
bsub Rscript LCL_NPX_PAR_SaturationAnalysis.R 70 100
bsub Rscript LCL_NPX_PAR_SaturationAnalysis.R 60 100
bsub Rscript LCL_NPX_PAR_SaturationAnalysis.R 50 100
bsub Rscript LCL_NPX_PAR_SaturationAnalysis.R 40 100
bsub Rscript LCL_NPX_PAR_SaturationAnalysis.R 30 100
bsub Rscript LCL_NPX_PAR_SaturationAnalysis.R 20 100
bsub Rscript LCL_NPX_PAR_SaturationAnalysis.R 10 100

bsub Rscript LCL_NPY_SaturationAnalysis.R 59 1
bsub Rscript LCL_NPY_SaturationAnalysis.R 55 100
bsub Rscript LCL_NPY_SaturationAnalysis.R 50 100
bsub Rscript LCL_NPY_SaturationAnalysis.R 45 100
bsub Rscript LCL_NPY_SaturationAnalysis.R 40 100
bsub Rscript LCL_NPY_SaturationAnalysis.R 35 100
bsub Rscript LCL_NPY_SaturationAnalysis.R 30 100
bsub Rscript LCL_NPY_SaturationAnalysis.R 25 100
bsub Rscript LCL_NPY_SaturationAnalysis.R 20 100
bsub Rscript LCL_NPY_SaturationAnalysis.R 15 100
bsub Rscript LCL_NPY_SaturationAnalysis.R 10 100
