##------------------------------------------
## Environment
##------------------------------------------
library(rsofun)
load_dependencies_rsofun()
# systr <- "''"    # for Mac
systr <- ""      # for Linux

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
  loutdwcont      = FALSE,
  loutdaet        = FALSE,
  loutdpet        = FALSE,
  loutdalpha      = FALSE
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
  maxit            = 4, # (5 for gensa) (30 for optimr)    #
  sitenames        = calibsites,
  filter_temp_max  = 35.0,
  filter_drought   = FALSE,
  metric           = "rmse",
  dir_results      = "~/eval_pmodel/calib_results"
 )


##//////////////////////////////////////////
## FULL
##------------------------------------------
## Prepare specific files and settings
##------------------------------------------
## Define calibration setup-specific simulation parameters
settings_sims$loutdgpp    = TRUE
settings_sims$soilmstress = TRUE
settings_sims$tempstress  = TRUE
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
  write_paramfils = FALSE 
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
settings_calib_FULL <- settings_calib
settings_calib_FULL$name = "FULL"
settings_calib_FULL$par = list( kphio       = list( lower=0.01, upper=0.4, init=0.1 ),
                                soilm_par_a = list( lower=0.0,  upper=1.0, init=0.2 ),
                                soilm_par_b = list( lower=0.0,  upper=2.0, init=0.2 ) )
settings_calib_FULL$datasource = list( gpp = "fluxnet2015_NT" )
settings_calib_FULL$filter_temp_min = NA
settings_calib_FULL$filter_soilm_min = NA

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
filn <- "~/eval_pmodel/data/ddf_obs_calib_NT.Rdata"
if (file.exists(filn)){
  load(filn)
} else {
  ddf_obs_calib <- get_obs_calib( 
    settings_calib = settings_calib_FULL, 
    settings_sims, 
    settings_input 
  )
  save(ddf_obs_calib, file = filn)  
}

filn <- "~/eval_pmodel/data/ddf_obs_eval_NT.Rdata"
if (file.exists(filn)){
  load(filn)
} else {
  ddf_obs_eval  <- get_obs_eval( 
    settings_eval = settings_eval, 
    settings_sims = settings_sims, 
    overwrite = TRUE, 
    light = TRUE 
  )
  save(ddf_obs_eval, file = filn)
}  

out_oob <- oob_calib_eval_sofun(
  setup = setup_sofun, 
  settings_calib = settings_calib_FULL, 
  settings_eval = settings_eval, 
  settings_sims = settings_sims, 
  settings_input = settings_input, 
  ddf_obs_calib = ddf_obs_calib,
  ddf_obs_eval = ddf_obs_eval
  )

save(out_oob, file = "~/eval_pmodel/data/out_oob_FULL.Rdata")
