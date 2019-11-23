##------------------------------------------
## Environment
##------------------------------------------
library(rsofun)
load_dependencies_rsofun()
systr <- "''"    # for Mac
# systr <- ""      # for Linux
overwrite <- TRUE

##------------------------------------------
## Simulation settings
##------------------------------------------
# source("filter_days.R")

path_siteinfo <- "~/eval_pmodel/siteinfo_pet_fluxnet2015.csv"
siteinfo <- rsofun::metainfo_Tier1_sites_kgclimate_fluxnet2015 %>% 
  dplyr::filter(!(sitename %in% c("DE-Akm", "IT-Ro1"))) %>%  # excluded because fapar data could not be downloaded (WEIRD)
  # dplyr::filter(!(sitename %in% c("AU-ASM", "AU-Wom"))) %>%  # excluded because no GPP data was found in FLUXNET file
  dplyr::filter(sitename != "FI-Sod") %>%  # excluded because some temperature data is missing
  dplyr::filter( c4 %in% c(FALSE, NA) & classid != "CRO" & classid != "WET" ) %>% 
  
  # ## test
  # dplyr::filter(sitename %in% c("FR-Pue", "FR-LBr", "IT-Noe")) %>% 
  
  write_csv(path = path_siteinfo)

settings_sims <- list(
  siteinfo       = path_siteinfo,
  ensemble       = TRUE,
  setup          = "site",
  name           = "fluxnet2015",
  dir_sofun      = "~/sofun/",
  path_output    = "~/sofun_outputs/output_fluxnet2015_sofun/202/",
  path_output_nc = "~/sofun_outputs/output_nc_fluxnet2015_sofun/s202/",
  path_input     = "~/sofun_inputs/input_fluxnet2015_sofun/",
  grid           = NA,
  implementation = "fortran",
  in_ppfd        = TRUE,
  in_netrad      = FALSE,
  recycle        = 1,
  spinupyears    = 10,
  calibvars      = c("gpp"),
  soilmstress    = TRUE,
  tempstress     = TRUE,
  loutdgpp       = TRUE,
  loutdwcont     = FALSE,
  loutdaet       = FALSE,
  loutdpet       = FALSE,
  loutdalpha     = FALSE,
  loutdgpp       = TRUE,
  loutdrd        = FALSE,
  loutdtransp    = FALSE,
  loutdwcont     = FALSE,
  loutdaet       = FALSE,
  loutdpet       = FALSE,
  loutdalpha     = FALSE
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
    bundle = "evi", 
    python_path = "/Users/benjaminstocker/Library/Enthought/Canopy_64bit/User/bin/python",
    gee_path = "~/gee_subset/gee_subset/"
    ),
  fapar = "MODIS_EVI_MOD13Q1",
  splined_fapar = FALSE
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
## Prepare the model setup for this calibration set
##------------------------------------------
settings_sims <- prepare_setup_sofun( 
  settings = settings_sims,
  setup = setup_sofun,
  write_paramfils = FALSE 
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
  # dplyr::filter(!(sitename %in% c("AU-Wom"))) %>%  # excluded because no GPP data was found in FLUXNET file
  dplyr::filter(sitename != "FI-Sod") %>%  # excluded because some temperature data is missing
  dplyr::filter( c4 %in% c(FALSE, NA) & classid != "CRO" & classid != "WET" ) %>%
  dplyr::filter( sitename %in% flue_sites ) %>%
  pull(sitename)

# ## Calibration sites for TerraP, excluding CA-Obs
# calibsites <- c("AU-Tum", "CA-NS3", "CA-NS6", "DE-Geb", "DE-Hai", "DE-Kli", "FI-Hyy", "FR-Fon", "FR-LBr", "FR-Pue", "IT-Cpz", "NL-Loo", "US-Ha1", "US-MMS", "US-UMB", "US-WCr")
# calibsites <- c("FR-Pue", "FR-LBr", "IT-Noe")

## Define calibration settings common for all setups
settings_calib <- list(
  method           = "gensa",
  targetvars       = c("gpp"),
  timescale        = list( gpp = "d" ),
  path_fluxnet2015 = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1d/original/unpacked/",
  path_fluxnet2015_hh= "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_0.5h/original/unpacked/",
  path_gepisat     = "~/data/gepisat/v3_fluxnet2015/daily_gpp/",
  maxit            = 10, # (5 for gensa) (30 for optimr)    #
  sitenames        = calibsites,
  filter_temp_max  = 35.0,
  filter_drought   = FALSE,
  metric           = "rmse",
  dir_results      = "~/eval_pmodel/calib_results",
  name = "FULL_EVI",
  par = list( kphio = list( lower=0.01, upper=0.4, init=0.1 ),
                            soilm_par_a = list( lower=0.0,  upper=1.0, init=0.2 ),
                            soilm_par_b = list( lower=0.0,  upper=2.0, init=0.2 ) ),
  datasource = list( gpp = "fluxnet2015_NT" ),
  filter_temp_min = NA,
  filter_soilm_min = NA
 )

##------------------------------------------
## Evaluation settings
##------------------------------------------
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


# 
# ##//////////////////////////////////////////
# ## FULL_EVI
# ##------------------------------------------
# ## Prepare input files
# ##------------------------------------------
inputdata <- prepare_input_sofun(
  settings_input        = settings_input,
  settings_sims         = settings_sims,
  return_data           = FALSE,
  overwrite_csv_climate = FALSE,
  overwrite_climate     = FALSE,
  overwrite_csv_fapar   = TRUE,
  overwrite_fapar       = TRUE,
  verbose               = TRUE
  )


##------------------------------------------
### Out of bag calibration for FULL_EVI
##------------------------------------------
filn <- "~/eval_pmodel/data/ddf_obs_calib_NT.Rdata"
if (file.exists(filn)){
  load(filn)
} else {
  ddf_obs_calib <- get_obs_calib( 
    settings_calib = settings_calib, 
    settings_sims, 
    settings_input 
  )
  save(ddf_obs_calib, file = filn)  
}

filn <- "~/eval_pmodel/data/obs_eval_NT.Rdata"
if (file.exists(filn)){
  load(filn)
} else {
  obs_eval  <- get_obs_eval( 
    settings_eval = settings_eval, 
    settings_sims = settings_sims, 
    overwrite     = TRUE, 
    light         = TRUE,
    add_forcing   = FALSE
  )
  save(obs_eval, file = filn)
}  


##------------------------------------------
## Single calibration and evaluation for FULL_EVI
## Using 75% of data for training and 25% for testing
##------------------------------------------
set.seed(1982)
settings_calib <- calib_sofun(
  setup          = setup_sofun,
  settings_calib = settings_calib,
  settings_sims  = settings_sims,
  settings_input = settings_input,
  ddf_obs        = ddf_obs_calib
)

filn <- "./data/mod_FULL_EVI.Rdata"
if (file.exists(filn)){
  load(filn)
} else {
  ## Update parameters
  param_filn <- paste0( settings_calib$dir_results, "/params_opt_", settings_calib$name, ".csv")
  params_opt <- readr::read_csv( param_filn )
  nothing <- update_params( params_opt, settings_sims$dir_sofun )

  ## run at evaluation sites
  mod <- runread_sofun(
    settings = settings_sims, 
    setup = setup_sofun
  )
  save(mod, file = filn)
} 


## evaluate at calib sites only (for comparison)
# settings_eval$sitenames <- settings_calib$sitenames
out_eval_FULL_EVI <- eval_sofun(
  mod,
  settings_eval,
  settings_sims,
  obs_eval = obs_eval,
  overwrite = TRUE,
  light = FALSE
  )

## write to file
save(out_eval_FULL_EVI, file = paste0(settings_calib$dir_results, "/out_eval_FULL_EVI.Rdata"))

save(settings_eval, file = "./data/settings_eval_FULL_EVI.Rdata")
save(settings_sims, file = "./data/settings_sims_FULL_EVI.Rdata")
save(settings_calib, file = "./data/settings_calib_FULL_EVI.Rdata")

