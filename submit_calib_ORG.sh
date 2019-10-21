#!/bin/bash

bsub -N -o ~/hpc_log/cluster_ORG.out  -n 3 -W 18:00 -R "rusage[mem=10000]" "R --vanilla --slave < ~/eval_pmodel/rscript_calib_ORG.R  >~/hpc_log/rscript_calib_ORG.Rout"
