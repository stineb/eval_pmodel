#!/bin/bash

bsub -N -o ~/hpc_log/cluster_BRC.out  -n 3 -W 12:00 -R "rusage[mem=10000]" "R --vanilla --slave < ~/eval_pmodel/rscript_calib_BRC.R  >~/hpc_log/rscript_calib_BRC.Rout"
