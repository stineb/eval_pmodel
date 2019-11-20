# P-model evaluation

## Workflow

1. Calibrate and evaluate for site-level run setups
	- `submit_calib_*.sh` runs `rscript_calib_*.R`
	- Outputs for calibration and evaluation done with the same set (complete) of sites is in `calib_results/out_eval_<setup>.Rdata`, `calib_results/params_opt_<setup>.csv`, and `out_gensa_<setup>.Rdata`
	- Outputs for calibration and evaluation done with the leave-one-site-out calibration is in subdirectories `org`, `brc`, and `full`; consodilated outputs from leave-one-site-out calibration is in `data/out_oob_<setup>.Rdata`

2. Get statistics and create figures with `eval_pmodel2.Rmd`
 	- Loads outputs of `rscript_calib_*.R`

3. Global simulations for ORG, BRC, and FULL using calibrate parameters from 1.
	- Run `./linkdirs_sofun_ORG_GLOBAL.py`
	- Run simulations and process output
	- Run `eval_pmodel_global.Rmd`
	
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
