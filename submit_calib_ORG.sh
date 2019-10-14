#!/bin/bash

bsub -N -o ~/hpc_log/grid.out -n 3 -W 10:00 -R "rusage[mem=10000]" "R --vanilla --slave < eval_pmodel/rscript_calib_ORG.R >hpc_log/rscript_calib_ORG.Rout"


# bsub -N -o ~/hpc_log/grid.out -n 3 -W 1:00 -R "rusage[mem=10000]" "R --vanilla --slave < eval_pmodel/rscript_calib_BRC.R >hpc_log/rscript_calib_BRC.Rout"
