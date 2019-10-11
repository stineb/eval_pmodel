##------------------------------------------
## Environment
##------------------------------------------
library(rsofun)
load_dependencies_rsofun()
systr <- "''"    # for Mac
# systr <- ""      # for Linux

##------------------------------------------
## Simulation settings
##------------------------------------------
# source("filter_days.R")

path_siteinfo <- "~/eval_pmodel/siteinfo_pet_fluxnet2015.csv"
siteinfo <- rsofun::metainfo_Tier1_sites_kgclimate_fluxnet2015 %>% 
  dplyr::filter(!(sitename %in% c("DE-Akm", "IT-Ro1"))) %>%  # excluded because fapar data could not be downloaded (WEIRD)
  dplyr::filter(!(sitename %in% c("AU-ASM", "AU-Wom"))) %>%  # excluded because no GPP data was found in FLUXNET file
  dplyr::filter(sitename != "FI-Sod") %>%  # excluded because some temperature data is missing
  dplyr::filter( c4 %in% c(FALSE, NA) & classid != "CRO" & classid != "WET" ) %>% 
  
  ## test
  dplyr::filter(sitename %in% c("FR-Pue", "FR-LBr", "IT-Noe")) %>% 
  
  write_csv(path = path_siteinfo)

settings_sims <- list(
  siteinfo        = path_siteinfo,
  ensemble        = TRUE,
  setup           = "site",
  name            = "fluxnet2015",
  dir_sofun       = "~/sofun/",
  path_output     = "~/sofun_outputs/output_fluxnet2015_sofun/202/",
  path_output_nc  = "~/sofun_outputs/output_nc_fluxnet2015_sofun/s202/",
  path_input      = "~/sofun_inputs/input_fluxnet2015_sofun/",
  grid            = NA,
  implementation  = "fortran",
  in_ppfd         = TRUE,
  in_netrad       = FALSE,
  recycle         = 1,
  spinupyears     = 10,
  calibvars       = c("gpp"),
  soilmstress     = TRUE,
  tempstress      = TRUE,
  loutdgpp        = TRUE,
  loutdwcont      = TRUE,
  loutdaet        = TRUE,
  loutdpet        = TRUE,
  loutdalpha      = TRUE
  )

##------------------------------------------
## Input settings
##------------------------------------------
settings_input <-  list(
  data                     = NA,
  temperature              = "fluxnet2015",
  precipitation            = "fluxnet2015",
  vpd                      = "fluxnet2015",
  ppfd                     = "fluxnet2015",
  netrad                   = "fluxnet2015",  #  c("fluxnet2015", "watch_wfdei"),
  patm                     = "fluxnet2015",
  netrad                   = NA,
  cloudcover               = "cru",
  path_watch_wfdei         = "~/data/watch_wfdei/",
  path_cru                 = "~/data/cru/ts_4.01/",
  path_MODIS_FPAR_MCD15A3H = "~/data/fluxnet_subsets/fapar_MODIS_FPAR_MCD15A3H_gee_MCD15A3H_fluxnet2015_gee_subset/",
  path_MODIS_EVI_MOD13Q1   = "~/data/fluxnet_subsets/fapar_MODIS_EVI_MOD13Q1_gee_MOD13Q1_fluxnet2015_gee_subset/",
  path_co2                 = "~/data/co2/cCO2_rcp85_const850-1765.dat",
  path_fluxnet2015         = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1d/original/unpacked/",
  path_fluxnet2015_hh      = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_0.5h/original/unpacked/",
  get_from_remote          = FALSE,
  settings_gee             = get_settings_gee( 
    bundle = "fpar", 
    python_path = "/Users/benjaminstocker/Library/Enthought/Canopy_64bit/User/bin/python",
    gee_path = "~/gee_subset/gee_subset/"
    )
  )


##------------------------------------------
## Model setup
##------------------------------------------
setup_sofun <- list(
  model      = "pmodel",
  dir        = "~/sofun",
  do_compile = FALSE,
  simsuite   = FALSE
  )


##------------------------------------------
## Calibration settings
##------------------------------------------
## Use only sites for calibration for which ANN method by Stocker et al. (2018) worked fine,
## and exclude sites where C4 vegetation is present.
flue_sites <- readr::read_csv( "~/data/flue/flue_stocker18nphyt.csv" ) %>%
              dplyr::filter( !is.na(cluster) ) %>% 
              distinct(site) %>% 
              pull(site)

calibsites <- rsofun::metainfo_Tier1_sites_kgclimate_fluxnet2015 %>% 
  dplyr::filter(!(sitename %in% c("DE-Akm", "IT-Ro1"))) %>%  # excluded because fapar data could not be downloaded (WEIRD)
  dplyr::filter(!(sitename %in% c("AU-ASM", "AU-Wom"))) %>%  # excluded because no GPP data was found in FLUXNET file
  dplyr::filter(sitename != "FI-Sod") %>%  # excluded because some temperature data is missing
  dplyr::filter( c4 %in% c(FALSE, NA) & classid != "CRO" & classid != "WET" ) %>%
  dplyr::filter( sitename %in% flue_sites ) %>%
  pull(sitename)

# ## Calibration sites for TerraP, excluding CA-Obs
# calibsites <- c("AU-Tum", "CA-NS3", "CA-NS6", "DE-Geb", "DE-Hai", "DE-Kli", "FI-Hyy", "FR-Fon", "FR-LBr", "FR-Pue", "IT-Cpz", "NL-Loo", "US-Ha1", "US-MMS", "US-UMB", "US-WCr")
calibsites <- c("FR-Pue", "FR-LBr", "IT-Noe")

## Define calibration settings common for all setups
settings_calib <- list(
  method           = "gensa",
  targetvars       = c("gpp"),
  timescale        = list( gpp = "d" ),
  path_fluxnet2015 = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1d/original/unpacked/",
  path_fluxnet2015_hh= "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_0.5h/original/unpacked/",
  path_gepisat     = "~/data/gepisat/v3_fluxnet2015/daily_gpp/",
  maxit            = 4, # (5 for gensa) (30 for optimr)    #
  sitenames        = calibsites,  # XXX debug dplyr::filter(siteinfo$light, c4 %in% c(FALSE, NA) )$mysitename, # calibrate for non-C4 sites
  filter_temp_max  = 35.0,
  filter_drought   = FALSE,
  metric           = "rmse",
  dir_results      = "~/eval_pmodel/calib_results"
 )


##//////////////////////////////////////////
## ORG
##------------------------------------------
## Prepare specific files and settings
##------------------------------------------
## Define calibration setup-specific simulation parameters
settings_sims$loutdgpp    = TRUE
settings_sims$soilmstress = FALSE
settings_sims$tempstress  = FALSE
settings_sims$loutdrd     = FALSE
settings_sims$loutdtransp = FALSE
settings_sims$loutdwcont  = FALSE
settings_sims$loutdaet    = FALSE
settings_sims$loutdpet    = FALSE
settings_sims$loutdalpha  = FALSE

## Prepare the model setup for this calibration set
settings_sims <- prepare_setup_sofun( 
  settings = settings_sims,
  setup = setup_sofun,
  write_paramfils = TRUE 
  )

## Define fAPAR input data and re-write input files for SOFUN
settings_input$fapar = "MODIS_FPAR_MCD15A3H"
settings_input$splined_fapar = TRUE

# inputdata <- prepare_input_sofun(
#   settings_input = settings_input,
#   settings_sims = settings_sims,
#   return_data = FALSE,
#   overwrite_csv_climate = TRUE,
#   overwrite_climate = TRUE,
#   overwrite_csv_fapar = TRUE,
#   overwrite_fapar = TRUE,
#   verbose = TRUE
#   )

## Additional setup-specific calibration-settings
## Specify data source for observations to which model is calibrated
settings_calib_ORG <- settings_calib
settings_calib_ORG$name = "ORG"
settings_calib_ORG$par = list( kphio = list( lower=0.01, upper=0.12, init=0.05 ) )
settings_calib_ORG$datasource = list( gpp = "fluxnet2015_NT" )
settings_calib_ORG$filter_temp_min = NA
settings_calib_ORG$filter_soilm_min = NA

### Evaluation settings
mylist <- readr::read_csv("~/eval_pmodel/myselect_fluxnet2015.csv") %>% 
  dplyr::filter( use==1 ) %>% 
  dplyr::pull( Site )

settings_eval <- list(
  sitenames = settings_sims$sitenames,
  sitenames_siteplots = mylist,
  agg = 8,
  path_fluxnet2015_d = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1d/original/unpacked/",
  path_fluxnet2015_w = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_7d/original/unpacked/",
  path_fluxnet2015_m = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1m/original/unpacked/",
  path_fluxnet2015_y = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1y/original/unpacked/",
  path_gepisat_d     = "~/data/gepisat/v3_fluxnet2015/daily_gpp/",
  benchmark = list( gpp = c("fluxnet2015_NT") ),
  remove_premodis = TRUE
  )

  
##------------------------------------------
### Out of bag calibration for ORG
##------------------------------------------
ddf_obs_calib <- get_obs_calib( 
  settings_calib = settings_calib_ORG, 
  settings_sims, 
  settings_input 
  )

ddf_obs_eval  <- get_obs_eval( 
  settings_eval = settings_eval, 
  settings_sims = settings_sims, 
  overwrite = TRUE, 
  light = TRUE 
  )

out_oob <- oob_calib_eval_sofun(
  setup = setup_sofun, 
  settings_calib = settings_calib_ORG, 
  settings_eval = settings_eval, 
  settings_sims = settings_sims, 
  settings_input = settings_input, 
  ddf_obs_calib = ddf_obs_calib,
  ddf_obs_eval = ddf_obs_eval
  )

save(out_oob, file = "~/eval_pmodel/data/out_oob_ORG.Rdata")


# ## Setup `BRC`

# The setup `BRC` is different from `ORG` in that it implements the temperature dependence of the quantum yield efficiency parameter after Bernacchi et al., 2003 PCE (included by setting `settings_sims$tempstress  = TRUE`)

# ```{r, eval=TRUE, message=FALSE, warning=FALSE}
# ## Define calibration setup-specific simulation parameters
# settings_sims$soilmstress = FALSE
# settings_sims$tempstress  = TRUE

# ## Prepare the model setup for this calibration set
# settings_sims <- prepare_setup_sofun( 
#   settings = settings_sims, 
#   setup = setup_sofun,
#   write_paramfils = TRUE 
#   )

# ## Define fAPAR input data and re-write input files for SOFUN
# settings_input$fapar = "MODIS_FPAR_MCD15A3H"
# settings_input$splined_fapar = TRUE

# inputdata <- prepare_input_sofun( 
#   settings_input = settings_input, 
#   settings_sims = settings_sims, 
#   return_data = FALSE, 
#   overwrite_csv_fapar = TRUE, 
#   overwrite_fapar = TRUE, 
#   verbose = TRUE
#   )

# ## Additional setup-specific calibration-settings
# ## Specify data source for observations to which model is calibrated
# settings_calib_BRC <- settings_calib
# settings_calib_BRC$name = "BRC"
# settings_calib_BRC$par = list( kphio = list( lower=0.01, upper=0.2, init=0.05 ) )
# settings_calib_BRC$datasource = list( gpp = "fluxnet2015_NT" )
# settings_calib_BRC$filter_temp_min = NA
# settings_calib_BRC$filter_soilm_min = NA

# ## Get observational data (GPP based on NT method) used as target for calibration
# if (!exists("ddf_obs_NT")){
#   ddf_obs_NT <- get_obs_calib( settings_calib_BRC, settings_sims, settings_input )
#   save( ddf_obs_NT, file = "./data/ddf_obs_NT.Rdata" )
# }

# set.seed(1982)
# settings_calib_BRC <- calib_sofun(
#   setup          = setup_sofun,
#   settings_calib = settings_calib_BRC,
#   settings_sims  = settings_sims,
#   settings_input = settings_input,
#   ddf_obs        = ddf_obs_NT
#   )
# ```


# ### Setup `FULL`

# Define settings, prepare setup and input data and do the calibration for the setup `FULL`. This is based on all data (including dry and cold days), with empirical soil moisture and low temperature stress functions and calibrates multiple parameters at once (`kphio` and parameters related to the soil moisture and temperature stress functions.

# ```{r, eval=TRUE, message=FALSE, warning=FALSE}
# ## Define calibration setup-specific simulation parameters
# settings_sims$soilmstress = TRUE
# settings_sims$tempstress  = TRUE

# ## Prepare the model setup for this calibration set
# settings_sims <- prepare_setup_sofun( 
#   settings = settings_sims, 
#   setup = setup_sofun,
#   write_paramfils = TRUE 
#   )

# ## Same fAPAR data used as above - nothing to be done 
# ## Define fAPAR input data and re-write input files for SOFUN
# settings_input$fapar = "MODIS_FPAR_MCD15A3H"
# settings_input$splined_fapar = TRUE

# inputdata <- prepare_input_sofun(
#   settings_input = settings_input,
#   settings_sims = settings_sims,
#   return_data = FALSE,
#   overwrite_csv_fapar = TRUE, # WARNING: ASSUMING "LINEAR" KNITTING 
#   overwrite_fapar = TRUE, # WARNING: ASSUMING "LINEAR" KNITTING 
#   verbose = TRUE
#   )

# ## Additional setup-specific calibration-settings
# ## Specify data source for observations to which model is calibrated
# settings_calib_FULL <- settings_calib
# settings_calib_FULL$name = "FULL"
# settings_calib_FULL$par = list( kphio       = list( lower=0.01, upper=0.4, init=0.1 ),
#                                 soilm_par_a = list( lower=0.0,  upper=1.0, init=0.2 ),
#                                 soilm_par_b = list( lower=0.0,  upper=2.0, init=0.2 ) )
# settings_calib_FULL$datasource = list( gpp = "fluxnet2015_NT" )
# settings_calib_FULL$filter_temp_min = NA
# settings_calib_FULL$filter_soilm_min = NA

# ## load obs data if not available
# if (!exists("ddf_obs_NT")){
#   ddf_obs_NT <- get_obs_calib( settings_calib_FULL, settings_sims, settings_input )
#   save( ddf_obs_NT, file = "./data/ddf_obs_NT.Rdata" )
# }


# ## Same observational data as for setup 'RED' (filtering is done afterwards, inside calib_sofun())
# set.seed(1982)
# settings_calib_FULL <- calib_sofun(
#   setup          = setup_sofun,
#   settings_calib = settings_calib_FULL,
#   settings_sims  = settings_sims,
#   settings_input = settings_input,
#   ddf_obs        = ddf_obs_NT
#   )
# ```

# <!-- #### Additional test for `FULL`
# Test the sensitivity of calibrated parameters to the observational data: repeat calibration of each site separately.
# ```{r}
# # ## test 
# # settings_calib_FULL_1site <- calib_sofun(
# #   setup          = setup_sofun,
# #   settings_calib = settings_calib_FULL,
# #   settings_sims  = settings_sims,
# #   settings_input = settings_input,
# #   ddf_obs        = dplyr::filter( ddf_obs_NT, sitename=="FR-Pue" ),
# #   sitename       = "FR-Pue"
# #   )

# # list_calib_bysite <- purrr::map( as.list(calibsites), 
# #                                  ~calib_sofun(
# #                                               setup          = setup_sofun,
# #                                               settings_calib = settings_calib_FULL,
# #                                               settings_sims  = settings_sims,
# #                                               settings_input = settings_input,
# #                                               ddf_obs        = dplyr::filter( ddf_obs_NT, sitename==. ),
# #                                               sitename       = .
# #                                               ) )
# # save( list_calib_bysite, file = "data/list_calib_bysite.Rdata" )
# ```
#  -->

# ### Sensitivity to fAPAR forcing data with `FULL_FPARitp` and `FULL_EVI`

# Do the same as described under `FULL`, but with alternative fAPAR data. Calibrated parameters are written to `"params_opt_FULL_FPARitp.csv"` and `"params_opt_FULL_EVI.csv"`.

# #### Setup `FULL_FPARitp`
# ```{r, eval=TRUE, message=FALSE, warning=FALSE}
# ## Define calibration setup-specific simulation parameters
# settings_sims$soilmstress = TRUE
# settings_sims$tempstress  = TRUE

# ## Prepare the model setup for this calibration set
# settings_sims <- prepare_setup_sofun( 
#   settings = settings_sims, 
#   setup = setup_sofun,
#   write_paramfils = TRUE 
#   )

# ## Define fAPAR input data and re-write input files for SOFUN
# settings_input$fapar = "MODIS_FPAR_MCD15A3H"
# settings_input$splined_fapar = FALSE

# inputdata <- prepare_input_sofun( 
#   settings_input = settings_input, 
#   settings_sims = settings_sims, 
#   return_data = FALSE, 
#   overwrite_csv_fapar = TRUE,
#   overwrite_fapar = TRUE,
#   verbose = TRUE
#   )

# ## Additional setup-specific calibration-settings
# ## Specify data source for observations to which model is calibrated
# settings_calib_FULL_FPARitp <- settings_calib
# settings_calib_FULL_FPARitp$name = "FULL_FPARitp"
# settings_calib_FULL_FPARitp$par = list( kphio       = list( lower=0.01, upper=0.4, init=0.1 ),
#                                         soilm_par_a = list( lower=0.0, upper=1.0, init=0.2 ),
#                                         soilm_par_b = list( lower=0.0, upper=2.0, init=0.2 ) )
# settings_calib_FULL_FPARitp$datasource = list( gpp = "fluxnet2015_NT" )
# settings_calib_FULL_FPARitp$filter_temp_min = NA
# settings_calib_FULL_FPARitp$filter_soilm_min = NA

# if (!exists("ddf_obs_NT")) load("./data/ddf_obs_NT.Rdata")

# set.seed(1982)
# settings_calib_FULL_FPARitp <- calib_sofun(
#   setup          = setup_sofun,
#   settings_calib = settings_calib_FULL_FPARitp,
#   settings_sims  = settings_sims,
#   settings_input = settings_input,
#   ddf_obs        = ddf_obs_NT
#   )
# ```


# #### Setup `FULL_EVI`
# ```{r, eval=TRUE, message=FALSE, warning=FALSE}
# ## Define calibration setup-specific simulation parameters
# settings_sims$soilmstress = TRUE
# settings_sims$tempstress  = TRUE

# ## Prepare the model setup for this calibration set
# settings_sims <- prepare_setup_sofun( 
#   settings = settings_sims, 
#   setup = setup_sofun,
#   write_paramfils = FALSE 
#   )
# ## Define fAPAR input data and re-write input files for SOFUN
# settings_input$fapar = "MODIS_EVI_MOD13Q1"
# settings_input$splined_fapar = FALSE

# inputdata <- prepare_input_sofun( 
#   settings_input = settings_input, 
#   settings_sims = settings_sims, 
#   return_data = FALSE, 
#   overwrite_csv_fapar = TRUE,
#   overwrite_fapar = TRUE,
#   verbose = TRUE
#   )

# ## Additional setup-specific calibration-settings
# ## Specify data source for observations to which model is calibrated
# settings_calib_FULL_EVI <- settings_calib
# settings_calib_FULL_EVI$name = "FULL_EVI"
# settings_calib_FULL_EVI$par = list( kphio       = list( lower=0.01, upper=0.4, init=0.1 ),
#                                     soilm_par_a = list( lower=0.0,  upper=1.0, init=0.2 ),
#                                     soilm_par_b = list( lower=0.0,  upper=2.0, init=0.2 ) )
# settings_calib_FULL_EVI$datasource = list( gpp = "fluxnet2015_NT" )
# settings_calib_FULL_EVI$filter_temp_min = NA
# settings_calib_FULL_EVI$filter_soilm_min = NA

# if (!exists("ddf_obs_NT")) load("./data/ddf_obs_NT.Rdata")

# settings_calib_FULL_EVI <- calib_sofun(
#   setup          = setup_sofun,
#   settings_calib = settings_calib_FULL_EVI,
#   settings_sims  = settings_sims,
#   settings_input = settings_input,
#   ddf_obs        = ddf_obs_NT
#   )
# ```


# ### Sensitivity to target GPP data with `FULL_DT` and `FULL_Ty` and `FULL_NTsub`

# Do the same as described under `ORG`, but with alternative calibration target data: GPP based on the daytime decomposition method and GPP on GePiSaT method. Note that GPP data decomposed based on the GePiSaT method is much sparser than the FLUXNET 2015 data (NT and DT decomposition). To have a fair comparison, we use days that have valid data in all datasets GePiSaT, NT, and DT. This requires the `FULL` setup to be repeated with only a subset of the data (`FULL_NTsub`). Calibrated parameters are written to `"params_opt_FULL_DT.csv"` and `"params_opt_FULL_Ty.csv"`.

# #### Setup `FULL_DT`

# ```{r, eval=TRUE, message=FALSE, warning=FALSE}
# ## Define calibration setup-specific simulation parameters
# settings_sims$soilmstress = TRUE
# settings_sims$tempstress  = TRUE

# ## Prepare the model setup for this calibration set
# settings_sims <- prepare_setup_sofun( 
#   settings = settings_sims, 
#   setup = setup_sofun,
#   write_paramfils = FALSE 
#   )

# ## Define fAPAR input data and re-write input files for SOFUN
# settings_input$fapar = "MODIS_FPAR_MCD15A3H"
# settings_input$splined_fapar = TRUE

# inputdata <- prepare_input_sofun( 
#   settings_input = settings_input, 
#   settings_sims = settings_sims, 
#   return_data = FALSE, 
#   overwrite_csv_fapar = TRUE, 
#   overwrite_fapar = TRUE, 
#   verbose = TRUE
#   )

# ## Additional setup-specific calibration-settings
# ## Specify data source for observations to which model is calibrated
# source("filter_days.R")
# settings_calib_FULL_DT <- settings_calib
# settings_calib_FULL_DT$name <- "FULL_DT"
# settings_calib_FULL_DT$par <- list( kphio       = list( lower=0.01, upper=0.4, init=0.1 ),
#                                     soilm_par_a = list( lower=0.0,  upper=1.0, init=0.2 ),
#                                     soilm_par_b = list( lower=0.0,  upper=2.0, init=0.2 ) )
# settings_calib_FULL_DT$datasource = list( gpp = "fluxnet2015_DT" )
# settings_calib_FULL_DT$filter_days = c("fluxnet2015_Ty", "fluxnet2015_DT", "fluxnet2015_NT")
# settings_calib_FULL_DT$filter_temp_min = NA
# settings_calib_FULL_DT$filter_soilm_min = NA

# ## Get observational data (GPP based on NT method) used as target for calibration
# if (!exists("ddf_obs_DT")){
#   ddf_obs_DT <- get_obs_calib( settings_calib_FULL_DT, settings_sims, settings_input ) %>% 
#     filter_days( settings_calib_FULL_DT$filter_days, settings_calib_FULL_DT$path_gepisat )
#   save(ddf_obs_DT, file = "data/ddf_obs_DT.Rdata")
# } else {
#   load("data/ddf_obs_DT.Rdata")
# }
# print(paste0("Total number of calibration target data points: ", sum(!is.na(ddf_obs_DT$gpp_obs))))

# settings_calib_FULL_DT <- calib_sofun(
#   setup          = setup_sofun,
#   settings_calib = settings_calib_FULL_DT,
#   settings_sims  = settings_sims,
#   settings_input = settings_input,
#   ddf_obs        = ddf_obs_DT
#   )
# ```

# #### Setup `FULL_NTsub`

# ```{r, eval=TRUE, message=FALSE, warning=FALSE}
# ## Define calibration setup-specific simulation parameters
# settings_sims$soilmstress = TRUE
# settings_sims$tempstress  = TRUE

# ## Prepare the model setup for this calibration set
# settings_sims <- prepare_setup_sofun( 
#   settings = settings_sims, 
#   setup = setup_sofun,
#   write_paramfils = FALSE 
#   )

# ## Define fAPAR input data and re-write input files for SOFUN
# settings_input$fapar = "MODIS_FPAR_MCD15A3H"
# settings_input$splined_fapar = TRUE

# inputdata <- prepare_input_sofun( 
#   settings_input = settings_input, 
#   settings_sims = settings_sims, 
#   return_data = FALSE, 
#   overwrite_csv_fapar = TRUE, 
#   overwrite_fapar = TRUE, 
#   verbose = TRUE
#   )

# ## Additional setup-specific calibration-settings
# ## Specify data source for observations to which model is calibrated
# settings_calib_FULL_NTsub <- settings_calib
# settings_calib_FULL_NTsub$name <- "FULL_NTsub"
# settings_calib_FULL_NTsub$par <- list(  kphio       = list( lower=0.01, upper=0.4, init=0.1 ),
#                                         soilm_par_a = list( lower=0.0,  upper=1.0, init=0.2 ),
#                                         soilm_par_b = list( lower=0.0,  upper=2.0, init=0.2 ) )
# settings_calib_FULL_NTsub$datasource = list( gpp = "fluxnet2015_NT" )
# settings_calib_FULL_NTsub$filter_days = c("fluxnet2015_Ty", "fluxnet2015_DT", "fluxnet2015_NT")
# settings_calib_FULL_NTsub$filter_temp_min = NA
# settings_calib_FULL_NTsub$filter_soilm_min = NA

# ## Get observational data (GPP based on NT method) used as target for calibration
# if (!exists("ddf_obs_NTsub")){
#   ddf_obs_NTsub <- get_obs_calib( settings_calib_FULL_NTsub, settings_sims, settings_input ) %>% 
#     filter_days( settings_calib_FULL_NTsub$filter_days, settings_calib_FULL_NTsub$path_gepisat )
#   save( ddf_obs_NTsub, file = "data/ddf_obs_NTsub.Rdata")
# } else {
#   load("data/ddf_obs_NTsub.Rdata")
# }
# print(paste0("Total number of calibration target data points: ", sum(!is.na(ddf_obs_NTsub$gpp_obs))))

# settings_calib_FULL_NTsub <- calib_sofun(
#   setup          = setup_sofun,
#   settings_calib = settings_calib_FULL_NTsub,
#   settings_sims  = settings_sims,
#   settings_input = settings_input,
#   ddf_obs        = ddf_obs_NTsub
#   )
# ```


# #### Setup `FULL_Ty`

# ```{r, eval=TRUE, message=FALSE, warning=FALSE}
# ## Define calibration setup-specific simulation parameters
# settings_sims$soilmstress = TRUE
# settings_sims$tempstress  = TRUE

# ## Prepare the model setup for this calibration set
# settings_sims <- prepare_setup_sofun( 
#   settings = settings_sims, 
#   setup = setup_sofun,
#   write_paramfils = FALSE 
#   )

# ## Define fAPAR input data and re-write input files for SOFUN
# settings_input$fapar = "MODIS_FPAR_MCD15A3H"
# settings_input$splined_fapar = TRUE

# inputdata <- prepare_input_sofun( 
#   settings_input = settings_input, 
#   settings_sims = settings_sims, 
#   return_data = FALSE, 
#   overwrite_climate = FALSE, 
#   overwrite_csv_fapar = TRUE, 
#   overwrite_fapar = TRUE, 
#   verbose = TRUE
#   )

# ## Additional setup-specific calibration-settings
# ## Specify data source for observations to which model is calibrated
# settings_calib_FULL_Ty <- settings_calib
# settings_calib_FULL_Ty$name = "FULL_Ty"
# settings_calib_FULL_Ty$par = list( kphio          = list( lower=0.01, upper=0.4, init=0.1 ),
#                                    soilm_par_a    = list( lower=0.0, upper=1.0, init=0.2 ),
#                                    soilm_par_b    = list( lower=0.0, upper=2.0, init=0.2 ) )
# settings_calib_FULL_Ty$datasource = list( gpp = "fluxnet2015_Ty" )
# settings_calib_FULL_Ty$filter_days = c("fluxnet2015_Ty", "fluxnet2015_DT", "fluxnet2015_NT")
# settings_calib_FULL_Ty$filter_temp_min = NA
# settings_calib_FULL_Ty$filter_soilm_min = NA

# ## Get observational data (GPP based on NT method) used as target for calibration
# if (!exists("ddf_obs_Ty")){
#   ddf_obs_Ty <- get_obs_calib( settings_calib_FULL_Ty, settings_sims, settings_input ) %>% 
#     filter_days( settings_calib_FULL_Ty$filter_days, settings_calib_FULL_Ty$path_gepisat )
#   save( ddf_obs_Ty, file = "data/ddf_obs_Ty.Rdata" )
# } else {
#   load("data/ddf_obs_Ty.Rdata")
# }
# print(paste0("Total number of calibration target data points: ", sum(!is.na(ddf_obs_Ty$gpp_obs))))

# settings_calib_FULL_Ty <- calib_sofun(
#   setup          = setup_sofun,
#   settings_calib = settings_calib_FULL_Ty,
#   settings_sims  = settings_sims,
#   settings_input = settings_input,
#   ddf_obs        = ddf_obs_Ty
#   )
# ```

