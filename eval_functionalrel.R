#' Evaluates the functional relationships.
#'
#' Evaluates the functional relationships in the data using General Additive Models (using function \code{gam} from the \code{mgcv} package) and the P-model using its function.
#' 
#' @param df A data frame containing daily observational data as output from function \code{get_obs_eval()}.
#' @param nam_target A character string specifying the column name of the target variable for which the functional relationships are to be determined.
#' @param filnam_params A character string specifying the path of the parameter file name for SOFUN.
#' @param overwrite if \code{TRUE}, the GAM model object is overwritten
#' @param ndays_agg An integer specifying the level of data aggregation by number of days. If \code{ndays_agg=NULL}, no aggregation is done.
#'
#' @return A data frame with aggregated data including GAM predictions and P-model results as columns \code{lue_gam} and \code{lue_mod}, respectively.
#' @export
#'
#' @examples eval_response( df, overwrite = TRUE, ndays_agg = 10 )
#' 
eval_functionalrel <- function( df, nam_target, predictors, overwrite = FALSE, ... ){

  # ## xxx debug
  # df <- out_eval_RED$data$ddf
  
  # ## rename (should go outside this function)
  # df <- df %>%  mutate( lue_obs = gpp_obs / (fapar * ppfd_fluxnet2015), lue_mod = gpp_mod / (fapar * ppfd_fluxnet2015) ) %>%
  #               mutate( lue_obs = ifelse( is.nan(lue_obs), NA, lue_obs ), lue_mod = ifelse( is.nan(lue_mod), NA, lue_mod ) ) %>%
  #               mutate( lue_obs = remove_outliers(lue_obs), lue_mod = remove_outliers(lue_mod) ) %>%
  #               dplyr::rename( vpd = vpd_fluxnet2015, ppfd = ppfd_fluxnet2015, soilm = soilm_obs_mean )
  
  # ##------------------------------------------------------------
  # ## Aggregate to multi-day periods
  # ## periods should start with the 1st of January each year, otherwise can't compute mean seasonal cycle
  # ##------------------------------------------------------------
  # if (!is.null(ndays_agg)){

  #   # ## Generate vector of starting dates of X-day periods, making sure the 1st of Jan is always the start of a new period
  #   # listyears <- seq( ymd("1990-01-01"), ymd("2018-01-01"), by = "year" )                  
  #   # breaks <- purrr::map( as.list(listyears), ~seq( from=., by=paste0( ndays_agg, " days"), length.out = ceiling(365 / ndays_agg)) ) %>% Reduce(c,.)
  
  #   # ## take mean across periods
  #   # df_agg <- df %>%  mutate( inbin = cut( date, breaks = breaks, right = FALSE ) ) %>%
  #   #   group_by( sitename, inbin ) %>%
  #   #   summarise_at( vars(one_of(c(nam_target, predictors))), funs(mean(., na.rm=TRUE) )) %>%
  #   #   mutate_at( vars(one_of(c(nam_target, predictors))), funs(ifelse(is.nan(.), NA, .)))  %>%
  #   #   mutate( date = ymd(as.character(inbin)) ) %>%
  #   #   dplyr::select( -inbin ) %>%
  #   #   ungroup()

  # } else {
  #   df_agg <- df
  # }

  ## filter out data if any of the variables is NA. Different for gpp_mod and gpp_obs
  df_training <- df %>% dplyr::filter_at( vars(one_of(c(nam_target, predictors))), all_vars(!is.na(.)) ) %>% 
  
    ## filter days with temperature below zero
    dplyr::filter( temp > 0.0 )

  ## train the neural network at observed daily GPP
  if (!file.exists("data/gam.Rdata")||overwrite){
    set.seed(1982)
    
    ## create formula with splines for each predictor "s(predictor)"
    forml  <- as.formula(  paste0( nam_target, " ~ s(", paste( predictors, collapse=") + s(" ), ")" ) )
    
    ## train model
    gam <- mgcv::gam( forml, data = df_training, method = "REML" )

    # gam_mod <- mgcv::gam( lue_mod ~ s(temp) + s(vpd) + s(soilm), data = df_training, method = "REML" )
    # gam.check(gam)
    # summary(gam)
    # plot(gam)
    save( gam,     file = "data/gam.Rdata" )
    # save( gam_mod, file = "data/gam_mod.Rdata" )
  } else {
    load( gam,     file = "data/gam.Rdata" )
    # load( gam_mod, file = "data/gam_mod.Rdata" )
 }

  ## predict values with GAM
  predicted     <- predict( gam,     df_training )
  # predicted_mod <- predict( gam_mod, df_training )

  ## calculate values with P-model
  # params_opt <- readr::read_csv( filnam_params )
  # df_training <- df_training %>%  
  #   mutate( lue_gam = predicted, lue_mod_gam = predicted_mod ) %>%
  #   group_by( date, sitename ) %>%
  #   nest() %>%
  #   mutate( out_pmodel = purrr::map( data, ~rpmodel(  tc = .$temp, 
  #                                                     vpd = .$vpd, 
  #                                                     co2 = 300, 
  #                                                     elv = 300, 
  #                                                     kphio = params_opt$kphio, 
  #                                                     fapar = NA, 
  #                                                     ppfd = NA, 
  #                                                     method="full"
  #                                                     ) ) ) %>%
  #   mutate( lue_mod = purrr::map_dbl( out_pmodel, varnam_pmodel ) ) %>%
  #   unnest( data )
              
  ## xxx test
  # with( filter(df, sitename=="FR-LBr"), plot( date, lue_obs, type="l", ylim=c(0,1)))
  # with( filter(df_agg, sitename=="FR-LBr"), lines( date, lue_obs, col="red"))
  # with( filter(df_training, sitename=="FR-LBr"), lines( date, lue_gam, col="green"))
  # with( filter(df_training, sitename=="FR-LBr"), lines( date, lue_mod, col="cyan"))

  # ## evaluate performance of GAM and P-model predictions
  # stats_gam <- with( df_training,  analyse_modobs( lue_gam, lue_obs, heat = TRUE ) )
  # stats_mod <- with( df_training,  analyse_modobs( lue_mod, lue_obs, heat = TRUE ) )
  
  ##-------------------------------------
  ## Evaluate GAM and P-model
  ##-------------------------------------
  df_eval <- purrr::map( as.list(predictors), ~eval_response_byvar(., df_training, gam, all_predictors = predictors, nam_target = nam_target, nsample = 12 ) )
  names(df_eval) <- predictors

  # ## temperature
  # eval_response_byvar( df_training, gam, evalvar = "temp", predictors = c("temp", "vpd", "soilm"), kphio = params_opt$kphio, varmin = 0, varmax = 40, nsample = 12, ylim=c(0,0.5), ... ) 

  # ## vpd
  # eval_response_byvar( df_training, gam, evalvar = "vpd", predictors = c("temp", "vpd", "soilm"), kphio = params_opt$kphio, varmin = 0, varmax = 3000, nsample = 12, ylim=c(0,0.5), ... )

  # ## soilm
  # eval_response_byvar( df_training, gam, evalvar = "soilm", predictors = c("temp", "vpd", "soilm"), kphio = params_opt$kphio, varmin = 0, varmax = 1.0, nsample = 12, ylim=c(0,0.5), ... )

  # ## Return aggregated data with GAM predictions and P-model results
  # df_agg <- df_agg %>% left_join( dplyr::select( df_training, sitename, date, lue_gam, lue_mod, lue_mod_gam ), by = c("sitename", "date") )

  return(df_eval)

}

eval_response_byvar <- function( evalvar, df, gam, all_predictors, nam_target, nsample, makepdf=FALSE, ... ){

  ## Get evaluation data
  evaldata <- get_eval_data( evalvar, df, all_predictors, nsample )

  ## predict for GAM with evaluation data
  predicted_eval <- predict( gam, evaldata )
  evaldata <- evaldata %>% mutate( var_gam = predicted_eval )

  # ## predict for P-model with evaluation data
  # evaldata <- evaldata %>%  mutate( id = 1:n() ) %>%
  #                           group_by( id ) %>%
  #                           nest() %>%
  #                           mutate( out_pmodel = purrr::map( data, ~rpmodel(  tc = .$temp, 
  #                                                                             vpd = .$vpd, 
  #                                                                             co2 = 300, 
  #                                                                             elv = 300, 
  #                                                                             kphio = kphio, 
  #                                                                             fapar = NA, 
  #                                                                             ppfd = NA, 
  #                                                                             method="full"
  #                                                                             ) ) ) %>%
  #                           mutate( lue_mod = purrr::map_dbl( out_pmodel, varnam_pmodel ) ) %>%
  #                           unnest( data )

  ## summarise by temperature steps
  eval_sum <- evaldata %>% get_quantiles( evalvar, nam_target = "var_gam" )
  names(eval_sum) <- paste0( names(eval_sum), "_gam" )

  # eval_sum <- evaldata %>% group_by_( evalvar ) %>%
  #                          summarise( 
  #                                     median_gam = median( var_gam ),  
  #                                     q33_gam = quantile(  var_gam, 0.33 ),
  #                                     q66_gam = quantile(  var_gam, 0.66 ),
  #                                     q25_gam = quantile(  var_gam, 0.25 ),
  #                                     q75_gam = quantile(  var_gam, 0.75 )

  #                                     # median_mod = median( lue_mod ),  
  #                                     # q33_mod = quantile( lue_mod, 0.33 ),
  #                                     # q66_mod = quantile( lue_mod, 0.66 ),
  #                                     # q25_mod = quantile( lue_mod, 0.25 ),
  #                                     # q75_mod = quantile( lue_mod, 0.75 )
  #                                     ) %>%
  #                          mutate( evalvar=evalvar )

#   ## plot response in observational and simulated data
#   if (makepdf) filn <- paste0("fig/gam_response_", evalvar, ".pdf")
#   if (makepdf) print( paste( "Creating plot", filn ))
#   if (makepdf) pdf( filn )
# 
#     par(las=0)
# 	  plot( eval_sum[[evalvar]], eval_sum$median_gam, type = "l", col="black", xlab = evalvar, ... )
# 	  polygon( c(eval_sum[[evalvar]], rev(eval_sum[[evalvar]])), c(eval_sum$q33_gam, rev(eval_sum$q66_gam)), col=rgb(0,0,0,0.2), border = NA )
# 
#     # lines( eval_sum[[evalvar]], eval_sum$median_mod, col="tomato" )
#     # polygon( c(eval_sum[[evalvar]], rev(eval_sum[[evalvar]])), c(eval_sum$q33_mod, rev(eval_sum$q66_mod)), col=add_alpha("tomato", 0.2), border = NA )
# 
#   if (makepdf) dev.off()

  return( eval_sum )

}


## aggregates to multi-day periods
aggregate_mean <- function( ddf, ndays_agg, dovars, year_start, year_end ){

  if (ndays_agg==1){

    return(ddf)

  } else {

    ## Generate periods: vector of starting dates of X-day periods, making sure the 1st of Jan is always the start of a new period
    listyears <- seq( ymd( paste0( year_start, "-01-01" ) ), ymd( paste0( year_end, "-01-01" ) ), by = "year" )                  
    breaks <- purrr::map( as.list(listyears), ~seq( from=., by=paste0( ndays_agg, " days"), length.out = ceiling(365 / ndays_agg)) ) %>% Reduce(c,.)
    
    ## take mean across periods and repalce NaN by NA
    ddf_agg <- ddf %>%  mutate( inbin = cut( date, breaks = breaks, right = FALSE ) ) %>%
      group_by( sitename, inbin ) %>%
      summarise_at( vars(one_of(dovars)), funs(mean(., na.rm=TRUE) )) %>%
      mutate_at( vars(one_of(dovars)), funs(ifelse(is.nan(.), NA, .)))  %>%
      mutate( date = ymd(as.character(inbin)) ) %>%
      dplyr::select( -inbin ) %>%
      ungroup()
    return(ddf_agg)
  
  }  
}


## applies P-model over each row in the data frame where temp, and vpd are separate columns
apply_pmodel <- function( df, varnam_pmodel, kphio ){

  df <- df %>% 
    mutate( id=1:n() ) %>%
    group_by( id ) %>%
    nest() %>%
    mutate( out_pmodel = purrr::map( data, ~rsofun::rpmodel(  tc = .$temp, 
                                                              vpd = .$vpd, 
                                                              co2 = 300, 
                                                              elv = 300, 
                                                              kphio = kphio, 
                                                              fapar = NA, 
                                                              ppfd = NA, 
                                                              method="full"
                                                              ) ) ) %>%
    mutate( var_pmodel = purrr::map_dbl( out_pmodel, varnam_pmodel ) ) %>%
    unnest( data ) %>%
    rename( varnam_pmodel = var_pmodel )

  return(df)
}

get_eval_data <- function( evalvar, df, all_predictors, nsample ){

  if (evalvar %in% all_predictors){
    predictors <- all_predictors[-which(all_predictors==evalvar)]
  } else {
    predictors <- all_predictors
  }

  varmin <- min( df[[evalvar]], na.rm=TRUE )
  varmax <- max( df[[evalvar]], na.rm=TRUE )

  ## create synthetic data for predictors, sampling from the observed values independently for each predictor
  evaldata <- expand.grid(  seq( varmin, varmax, length.out = 30 ),
                            sample( unlist( df[ predictors[1] ] ), nsample, replace = TRUE ), 
                            sample( unlist( df[ predictors[2] ] ), nsample, replace = TRUE )
                            ) %>% 
              as_tibble() %>%
              setNames( c( evalvar, predictors ) )

  return(evaldata)
}

get_quantiles <- function( df, evalvar, nam_target ){

  ## summarise by temperature steps
  df_agg <- df %>% group_by_( evalvar ) %>%
                   summarise( 
                              median = median(    eval(as.name(nam_target)) ),  
                              q33    = quantile(  eval(as.name(nam_target)), 0.33 ),
                              q66    = quantile(  eval(as.name(nam_target)), 0.66 ),
                              q25    = quantile(  eval(as.name(nam_target)), 0.25 ),
                              q75    = quantile(  eval(as.name(nam_target)), 0.75 )
                              ) %>%
                   mutate( evalvar=evalvar )
  return(df_agg)
}

prepare_data_functionalrel <- function( df, predictors, ... ){

  df <- df %>%  mutate( lue_obs = gpp_obs / (fapar * ppfd_fluxnet2015) ) %>%
                mutate( lue_obs = ifelse( is.nan(lue_obs), NA, lue_obs ) ) %>%
                mutate( lue_obs = remove_outliers(lue_obs) ) %>%
                dplyr::rename( vpd = vpd_fluxnet2015, ppfd = ppfd_fluxnet2015, soilm = soilm_obs_mean ) %>%
                aggregate_mean( ... )
  return(df)  
}

