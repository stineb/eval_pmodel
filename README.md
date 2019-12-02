# P-model evaluation

## Workflow

1. Calibrate and evaluate for site-level run setups
  - Single calibration/evaluation for all setups done with `./rscript_calib_<SETUP>.R` locally. Calibration and evaluation results are in `./calib_results/`.
  - Out-of-bag calibration/evaluation done on Euler (see also ./rscript_calib_<SETUP>.R). Job submissions with `submit_calib_<SETUP>.sh` which runs `rscript_calib_<SETUP>.R`. Calibration and evaluation results are in `./calib_results/oob_<SETUP>/`.

2. Get statistics and create all figures for site-scale evaluations of GPP with `eval_pmodel2.Rmd`, and for LUE with `eval_lue.Rmd`.

3. Global simulations for ORG, BRC, and FULL using calibrate parameters from 1.
	- Run `./linkdirs_sofun_ORG_GLOBAL.py`
	- Run simulations and process output
	- Run `eval_pmodel_global.Rmd`
	
4. Create sites table with `plot_table_siteselection.Rmd`

	
## Global simulations	

!----------------                          !----------!----------                                                                                    !
! Simulation name                          ! Setup    ! fAPAR forcing file                                                                           !
!----------------                          !----------!----------                                                                                    !
! global_ORG_MODIS-C006_MOD15A2_2000_2016  ! ORG      ! MODIS-C006_MOD15A2__LAI_FPAR__LPDAAC__GLOBAL_0.5degree__UHAM-ICDC__2000_2018__MON__fv0.02.nc !        
! global_BRC_MODIS-C006_MOD15A2_2000_2016  ! BRC      ! MODIS-C006_MOD15A2__LAI_FPAR__LPDAAC__GLOBAL_0.5degree__UHAM-ICDC__2000_2018__MON__fv0.02.nc !        
! global_FULL_MODIS-C006_MOD15A2_2000_2016 ! FULL     ! MODIS-C006_MOD15A2__LAI_FPAR__LPDAAC__GLOBAL_0.5degree__UHAM-ICDC__2000_2018__MON__fv0.02.nc !         
! global_ORG_fAPAR3g_v2_2000_2016          ! ORG      ! fAPAR3g_v2_1982_2016_FILLED.nc                                                               !
! global_BRC_fAPAR3g_v2_2000_2016          ! BRC      ! fAPAR3g_v2_1982_2016_FILLED.nc                                                               !
! global_FULL_fAPAR3g_v2_2000_2016         ! FULL     ! fAPAR3g_v2_1982_2016_FILLED.nc                                                               !
!----------------                          !----------!----------                                                                                    !     
