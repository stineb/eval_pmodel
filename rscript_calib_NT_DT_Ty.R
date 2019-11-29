##------------------------------------------
## Environment
##------------------------------------------
library(rsofun)
load_dependencies_rsofun()
systr <- "''"    # for Mac
# systr <- ""      # for Linux
overwrite <- TRUE
source("filter_days.R")

##------------------------------------------
## Simulation settings
##------------------------------------------
path_siteinfo <- "~/eval_pmodel/siteinfo_pet_fluxnet2015.csv"
siteinfo <- rsofun::metainfo_Tier1_sites_kgclimate_fluxnet2015 %>% 
  dplyr::filter(!(sitename %in% c("DE-Akm", "IT-Ro1"))) %>%  # excluded because fapar data could not be downloaded (WEIRD)
  # dplyr::filter(!(sitename %in% c("AU-ASM", "AU-Wom"))) %>%  # excluded because no GPP data was found in FLUXNET file
  dplyr::filter(sitename != "FI-Sod") %>%  # excluded because some temperature data is missing
  dplyr::filter( c4 %in% c(FALSE, NA) & classid != "CRO" & classid != "WET" ) %>% 
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
  threshold_GPP            = 0.5,
  get_from_remote          = FALSE,
  settings_gee             = get_settings_gee( 
    bundle = "fpar", 
    python_path = "/Users/benjaminstocker/Library/Enthought/Canopy_64bit/User/bin/python",
    gee_path = "~/gee_subset/gee_subset/"
    ),
  fapar = "MODIS_FPAR_MCD15A3H",
  splined_fapar = TRUE
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
## Prepare input files
##------------------------------------------
# inputdata <- prepare_input_sofun(
#   settings_input        = settings_input,
#   settings_sims         = settings_sims,
#   return_data           = FALSE,
#   overwrite_csv_climate = FALSE,
#   overwrite_climate     = FALSE,
#   overwrite_csv_fapar   = TRUE,
#   overwrite_fapar       = TRUE,
#   verbose               = TRUE
#   )


##//////////////////////////////////////////
## DT
##------------------------------------------

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

  ## Define calibration settings common for all setups
  settings_calib_DT <- list(
    method           = "gensa",
    targetvars       = c("gpp"),
    timescale        = list( gpp = "d" ),
    path_fluxnet2015 = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1d/original/unpacked/",
    path_fluxnet2015_hh= "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_0.5h/original/unpacked/",
    threshold_GPP            = 0.5,
    path_gepisat     = "~/data/gepisat/v3_fluxnet2015/daily_gpp/",
    maxit            = 10, # (5 for gensa) (30 for optimr)    #
    sitenames        = calibsites,
    filter_temp_max  = 35.0,
    filter_drought   = FALSE,
    metric           = "rmse",
    dir_results      = "~/eval_pmodel/calib_results",
    name = "FULL_DT",
    par = list( kphio       = list( lower=0.06, upper=0.1, init=0.085 ),
                soilm_par_a = list( lower=0.0,  upper=1.0, init=0.0 ),
                soilm_par_b = list( lower=0.0,  upper=1.5, init=0.6 ) ),
    datasource = list( gpp = "fluxnet2015_DT" ),
    filter_temp_min = NA,
    filter_soilm_min = NA,
    filter_days = c("fluxnet2015_Ty", "fluxnet2015_DT", "fluxnet2015_NT")
    )


  ##------------------------------------------
  ## Evaluation settings
  ##------------------------------------------
  mylist <- readr::read_csv("~/eval_pmodel/myselect_fluxnet2015.csv") %>% 
    dplyr::filter( use==1 ) %>% 
    dplyr::pull( Site )

  settings_eval_DT <- list(
    sitenames = settings_sims$sitenames,
    sitenames_siteplots = mylist,
    agg = 8,
    path_fluxnet2015_d = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1d/original/unpacked/",
    path_fluxnet2015_w = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_7d/original/unpacked/",
    path_fluxnet2015_m = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1m/original/unpacked/",
    path_fluxnet2015_y = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1y/original/unpacked/",
    path_gepisat_d     = "~/data/gepisat/v3_fluxnet2015/daily_gpp/",
    benchmark = list( gpp = c("fluxnet2015_DT") ),
    remove_premodis = TRUE,
    filter_days = c("fluxnet2015_Ty", "fluxnet2015_DT", "fluxnet2015_NT")
    )


  ##------------------------------------------
  ## Collect caliration target data
  ##------------------------------------------
  filn <- "~/eval_pmodel/data/obs_calib_DT.Rdata"
  if (file.exists(filn)){
    load(filn)
  } else {
    obs_calib <- get_obs_calib( 
      settings_calib = settings_calib_DT, 
      settings_sims, 
      settings_input 
    )
    save(obs_calib, file = filn)  
  }


  ##------------------------------------------
  ## Single calibration
  ##------------------------------------------
  # set.seed(1982)
  # settings_calib_DT <- calib_sofun(
  #   setup          = setup_sofun,
  #   settings_calib = settings_calib_DT,
  #   settings_sims  = settings_sims,
  #   settings_input = settings_input,
  #   ddf_obs        = obs_calib
  #   )


  ##------------------------------------------
  ## Run model
  ##------------------------------------------
  ## Update parameters
  filn <- paste0( settings_calib_DT$dir_results, "/params_opt_", settings_calib_DT$name, ".csv")
  params_opt <- readr::read_csv( filn )
  nothing <- update_params( params_opt, settings_sims$dir_sofun )

  ## run at evaluation sites
  filn <- "./data/mod_DT.Rdata"
  if (file.exists(filn)){
    load(filn)
  } else {
    mod_DT <- runread_sofun(
      settings = settings_sims, 
      setup = setup_sofun
    )
    save(mod_DT, file = filn)
  }


  ##------------------------------------------
  ## Get evaluation benchmark data
  ##------------------------------------------
  filn <- "./data/obs_eval_DT.Rdata"
  if (file.exists(filn)){
    load(filn)
  } else {
    obs_eval_DT  <- get_obs_eval( 
      settings_eval = settings_eval_DT, 
      settings_sims = settings_sims, 
      overwrite     = TRUE, 
      light         = TRUE,
      add_forcing   = FALSE
    )
    save(obs_eval_DT, file = filn)
  }  


##//////////////////////////////////////////
## NT
##------------------------------------------

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

  ## Define calibration settings common for all setups
  settings_calib_NT <- list(
    method           = "gensa",
    targetvars       = c("gpp"),
    timescale        = list( gpp = "d" ),
    path_fluxnet2015 = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1d/original/unpacked/",
    path_fluxnet2015_hh= "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_0.5h/original/unpacked/",
    threshold_GPP            = 0.5,
    path_gepisat     = "~/data/gepisat/v3_fluxnet2015/daily_gpp/",
    maxit            = 10, # (5 for gensa) (30 for optimr)    #
    sitenames        = calibsites,
    filter_temp_max  = 35.0,
    filter_drought   = FALSE,
    metric           = "rmse",
    dir_results      = "~/eval_pmodel/calib_results",
    name = "FULL_NT",
    par = list( kphio       = list( lower=0.06, upper=0.1, init=0.085 ),
                soilm_par_a = list( lower=0.0,  upper=1.0, init=0.0 ),
                soilm_par_b = list( lower=0.0,  upper=1.5, init=0.6 ) ),
    datasource = list( gpp = "fluxnet2015_NT" ),
    filter_temp_min = NA,
    filter_soilm_min = NA,
    filter_days = c("fluxnet2015_Ty", "fluxnet2015_DT", "fluxnet2015_NT")
    )


  ##------------------------------------------
  ## Evaluation settings
  ##------------------------------------------
  mylist <- readr::read_csv("~/eval_pmodel/myselect_fluxnet2015.csv") %>% 
    dplyr::filter( use==1 ) %>% 
    dplyr::pull( Site )

  settings_eval_NTsub <- list(
    sitenames = settings_sims$sitenames,
    sitenames_siteplots = mylist,
    agg = 8,
    path_fluxnet2015_d = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1d/original/unpacked/",
    path_fluxnet2015_w = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_7d/original/unpacked/",
    path_fluxnet2015_m = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1m/original/unpacked/",
    path_fluxnet2015_y = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1y/original/unpacked/",
    path_gepisat_d     = "~/data/gepisat/v3_fluxnet2015/daily_gpp/",
    benchmark = list( gpp = c("fluxnet2015_NT") ),
    remove_premodis = TRUE,
    filter_days = c("fluxnet2015_Ty", "fluxnet2015_DT", "fluxnet2015_NT")
    )


  ##------------------------------------------
  ## Collect caliration target data
  ##------------------------------------------
  filn <- "~/eval_pmodel/data/obs_calib_NTsub.Rdata"
  if (file.exists(filn)){
    load(filn)
  } else {
    obs_calib <- get_obs_calib( 
      settings_calib = settings_calib_NT, 
      settings_sims, 
      settings_input 
    )
    save(obs_calib, file = filn)  
  }


  ##------------------------------------------
  ## Single calibration
  ##------------------------------------------
  # set.seed(1982)
  # settings_calib_NT <- calib_sofun(
  #   setup          = setup_sofun,
  #   settings_calib = settings_calib_NT,
  #   settings_sims  = settings_sims,
  #   settings_input = settings_input,
  #   ddf_obs        = obs_calib
  #   )


  ##------------------------------------------
  ## Run model
  ##------------------------------------------
  ## Update parameters
  filn <- paste0( settings_calib_NT$dir_results, "/params_opt_", settings_calib_NT$name, ".csv")
  params_opt <- readr::read_csv( filn )
  nothing <- update_params( params_opt, settings_sims$dir_sofun )

  ## run at evaluation sites
  filn <- "./data/mod_NT.Rdata"
  if (file.exists(filn)){
    load(filn)
  } else {
    mod_NT <- runread_sofun(
      settings = settings_sims, 
      setup = setup_sofun
    )
    save(mod_NT, file = filn)
  }


  ##------------------------------------------
  ## Get evaluation benchmark data
  ##------------------------------------------
  filn <- "./data/obs_eval_NTsub.Rdata"
  if (file.exists(filn)){
    load(filn)
  } else {
    obs_eval_NTsub  <- get_obs_eval( 
      settings_eval = settings_eval_NTsub, 
      settings_sims = settings_sims, 
      overwrite     = TRUE, 
      light         = TRUE,
      add_forcing   = FALSE
    )
    save(obs_eval_NTsub, file = filn)
  }  

##//////////////////////////////////////////
## Ty
##------------------------------------------

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

  ## Define calibration settings common for all setups
  settings_calib_Ty <- list(
    method           = "gensa",
    targetvars       = c("gpp"),
    timescale        = list( gpp = "d" ),
    path_fluxnet2015 = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1d/original/unpacked/",
    path_fluxnet2015_hh= "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_0.5h/original/unpacked/",
    threshold_GPP            = 0.5,
    path_gepisat     = "~/data/gepisat/v3_fluxnet2015/daily_gpp/",
    maxit            = 10, # (5 for gensa) (30 for optimr)    #
    sitenames        = calibsites,
    filter_temp_max  = 35.0,
    filter_drought   = FALSE,
    metric           = "rmse",
    dir_results      = "~/eval_pmodel/calib_results",
    name = "FULL_Ty",
    par = list( kphio       = list( lower=0.06, upper=0.1, init=0.085 ),
                soilm_par_a = list( lower=0.0,  upper=1.0, init=0.0 ),
                soilm_par_b = list( lower=0.0,  upper=1.5, init=0.6 ) ),
    datasource = list( gpp = "fluxnet2015_Ty" ),
    filter_temp_min = NA,
    filter_soilm_min = NA,
    filter_days = c("fluxnet2015_Ty", "fluxnet2015_DT", "fluxnet2015_NT")
    )


  ##------------------------------------------
  ## Evaluation settings
  ##------------------------------------------
  mylist <- readr::read_csv("~/eval_pmodel/myselect_fluxnet2015.csv") %>% 
    dplyr::filter( use==1 ) %>% 
    dplyr::pull( Site )

  settings_eval_Ty <- list(
    sitenames = settings_sims$sitenames,
    sitenames_siteplots = mylist,
    agg = 8,
    path_fluxnet2015_d = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1d/original/unpacked/",
    path_fluxnet2015_w = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_7d/original/unpacked/",
    path_fluxnet2015_m = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1m/original/unpacked/",
    path_fluxnet2015_y = "~/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1y/original/unpacked/",
    path_gepisat_d     = "~/data/gepisat/v3_fluxnet2015/daily_gpp/",
    benchmark = list( gpp = c("fluxnet2015_Ty") ),
    remove_premodis = TRUE,
    filter_days = c("fluxnet2015_Ty", "fluxnet2015_DT", "fluxnet2015_NT")
    )


  ##------------------------------------------
  ## Collect caliration target data
  ##------------------------------------------
  filn <- "~/eval_pmodel/data/obs_calib_Ty.Rdata"
  if (file.exists(filn)){
    load(filn)
  } else {
    obs_calib <- get_obs_calib( 
      settings_calib = settings_calib_Ty, 
      settings_sims, 
      settings_input 
    )
    save(obs_calib, file = filn)  
  }


  ##------------------------------------------
  ## Single calibration
  ##------------------------------------------
  # set.seed(1982)
  # settings_calib_Ty <- calib_sofun(
  #   setup          = setup_sofun,
  #   settings_calib = settings_calib_Ty,
  #   settings_sims  = settings_sims,
  #   settings_input = settings_input,
  #   ddf_obs        = obs_calib
  #   )


  ##------------------------------------------
  ## Run model
  ##------------------------------------------
  ## Update parameters
  filn <- paste0( settings_calib_Ty$dir_results, "/params_opt_", settings_calib_Ty$name, ".csv")
  params_opt <- readr::read_csv( filn )
  nothing <- update_params( params_opt, settings_sims$dir_sofun )

  ## run at evaluation sites
  filn <- "./data/mod_Ty.Rdata"
  if (file.exists(filn)){
    load(filn)
  } else {
    mod_Ty <- runread_sofun(
      settings = settings_sims, 
      setup = setup_sofun
    )
    save(mod_Ty, file = filn)
  }


  ##------------------------------------------
  ## Get evaluation benchmark data
  ##------------------------------------------
  filn <- "./data/obs_eval_Ty.Rdata"
  if (file.exists(filn)){
    load(filn)
  } else {
    obs_eval_Ty  <- get_obs_eval( 
      settings_eval = settings_eval_Ty, 
      settings_sims = settings_sims, 
      overwrite     = TRUE, 
      light         = TRUE,
      add_forcing   = FALSE
    )
    save(obs_eval_Ty, file = filn)
  }  


##//////////////////////////////////////////
## Combine obs_eval data for DT, NT, and Ty
##------------------------------------------
merged <- select(    obs_eval_DT$ddf,    sitename, date, gpp_DT = gpp ) %>% 
  left_join( select( obs_eval_NTsub$ddf, sitename, date, gpp_NT = gpp ), by = c("sitename", "date") ) %>% 
  left_join( select( obs_eval_Ty$ddf,    sitename, date, gpp_Ty = gpp ), by = c("sitename", "date") ) %>% 
  mutate( use = ifelse( !is.na(gpp_DT) & !is.na(gpp_NT) & !is.na(gpp_Ty), TRUE, FALSE ) ) %>% 
  mutate( gpp_DT = ifelse(use, gpp_DT, NA), gpp_NT = ifelse(use, gpp_NT, NA), gpp_Ty = ifelse(use, gpp_Ty, NA) )
sum(!is.na(merged$gpp_DT))
sum(!is.na(merged$gpp_NT))
sum(!is.na(merged$gpp_Ty))
obs_eval_DT$ddf$gpp    <- merged$gpp_DT
obs_eval_NTsub$ddf$gpp <- merged$gpp_NT
obs_eval_Ty$ddf <- obs_eval_Ty$ddf %>% 
  dplyr::select(-gpp) %>% 
  left_join(dplyr::select(merged, sitename, date, gpp = gpp_Ty), by = c("sitename", "date"))


# x-daily
merged <- select(    obs_eval_DT$xdf,    sitename, inbin, gpp_DT = gpp ) %>% 
  left_join( select( obs_eval_NTsub$xdf, sitename, inbin, gpp_NT = gpp ), by = c("sitename", "inbin") ) %>% 
  left_join( select( obs_eval_Ty$xdf,    sitename, inbin, gpp_Ty = gpp ), by = c("sitename", "inbin") ) %>% 
  mutate( use = ifelse( !is.na(gpp_DT) & !is.na(gpp_NT) & !is.na(gpp_Ty), TRUE, FALSE ) ) %>% 
  mutate( gpp_DT = ifelse(use, gpp_DT, NA), gpp_NT = ifelse(use, gpp_NT, NA), gpp_Ty = ifelse(use, gpp_Ty, NA) )
sum(!is.na(merged$gpp_DT))
sum(!is.na(merged$gpp_NT))
sum(!is.na(merged$gpp_Ty))
obs_eval_DT$xdf$gpp    <- merged$gpp_DT
obs_eval_NTsub$xdf$gpp <- merged$gpp_NT
obs_eval_Ty$xdf <- obs_eval_Ty$xdf %>% 
  dplyr::select(-gpp) %>% 
  left_join(dplyr::select(merged, sitename, inbin, gpp = gpp_Ty), by = c("sitename", "inbin"))


##//////////////////////////////////////////
## Evaluate all
##------------------------------------------
## DT
##------------------------------------------
# settings_eval_DT$sitenames <- settings_calib_DT$sitenames
out_eval_DT <- eval_sofun(
  mod_DT,
  settings_eval_DT,
  settings_sims,
  obs_eval = obs_eval_DT,
  overwrite = TRUE,
  light = FALSE
  )

## write to files
save(out_eval_DT, file = paste0(settings_calib_DT$dir_results, "/out_eval_FULL_DT.Rdata"))
save(settings_eval_DT,  file = "./data/settings_eval_FULL_DT.Rdata")
save(settings_calib_DT, file = "./data/settings_calib_FULL_DT.Rdata")

##------------------------------------------
## NT
##------------------------------------------
# settings_eval_NTsub$sitenames <- settings_calib_NT$sitenames
out_eval_NTsub <- eval_sofun(
  mod_NT,
  settings_eval_NTsub,
  settings_sims,
  obs_eval = obs_eval_NTsub,
  overwrite = TRUE,
  light = FALSE
  )

## write to files
save(out_eval_NTsub, file = paste0(settings_calib_NT$dir_results, "/out_eval_FULL_NTsub.Rdata"))
save(settings_eval_NTsub,  file = "./data/settings_eval_FULL_NTsub.Rdata")
save(settings_calib_NT, file = "./data/settings_calib_FULL_NTsub.Rdata")

##------------------------------------------
## Ty
##------------------------------------------
# settings_eval_Ty$sitenames <- settings_calib_Ty$sitenames
out_eval_Ty <- eval_sofun(
  mod_Ty,
  settings_eval_Ty,
  settings_sims,
  obs_eval = obs_eval_Ty,
  overwrite = TRUE,
  light = FALSE
  )

## write to files
save(out_eval_Ty, file = paste0(settings_calib_Ty$dir_results, "/out_eval_FULL_Ty.Rdata"))
save(settings_eval_Ty,  file = "./data/settings_eval_FULL_Ty.Rdata")
save(settings_calib_Ty, file = "./data/settings_calib_Ty.Rdata")

save(settings_sims,     file = "./data/settings_sims_FULL.Rdata")

