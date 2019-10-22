#!/bin/bash

bsub -N -o ~/hpc_log/cluster_FULL.out -n 3 -W 48:00 -R "rusage[mem=10000]" "R --vanilla --slave < ~/eval_pmodel/rscript_calib_FULL.R >~/hpc_log/rscript_calib_FULL.Rout"
