# P-model evaluation

## Workflow

1. Calibrate and evaluate for site-level run setups
	- `submit_calib_*.sh` runs `rscript_calib_*.R`

2. Get statistics and create figures with `eval_pmodel2.Rmd`
 	- Loads outputs of `rscript_calib_*.R`

3. Global simulations for ORG, BRC, and FULL using calibrate parameters from 1.
	- Run `./linkdirs_sofun_ORG_GLOBAL.py`