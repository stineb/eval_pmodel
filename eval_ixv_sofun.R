#' Evaluates SOFUN model outputs.
#'
#' Calculates a set of perfomance metrics for model outputs, compared against observational data.
#' Currently only evaluations of GPP model outputs, compared agains FLUXNET 2015 data, are implemented.
#'
#' @param agg An integer specifying the number of days over which to aggregate
#' @param settings A named list of data frames containing containing model outputs. 
#' The names of list elements corresponds to site names.
#' @param settings_eval A list specifying evaluation settings 
#' (see vignette eval_sofun.pdf for more information and examples)
#' @param settings_sims A named list containing the simulation settings 
#' (see vignette_rsofun.pdf for more information and examples)
#' @param ddf_obs (Optional) A named list of data frames containing observational data for each sites. 
#' The names of list elements corresponds to site names. Defaults to \code{NA} 
#'
#' @return A list containing data frames of modelled and observed values aggregated to several temporal scales 
#' (ddf for daily, xdf for X-daily, mdf for monthly, adf for annual), data frames of respective performance metrics,
#' and functions for plotting respective data.
#' @export
#'
#' @examples out_eval <- eval_sofun( mod, settings_eval, settings_sims )
#' 
eval_ixv_sofun <- function( agg, mod, settings_eval, settings_sims, ddf_obs ){


  ## get sites for which no model output is available and overwrite settings_eval$sitenames
  missing_mod <- purrr::map_lgl( mod, ~identical(., NA ) ) %>% which() %>% names()
  settings_eval$sitenames <- settings_eval$sitenames[which(!(settings_eval$sitenames %in% missing_mod))]
  
  ##------------------------------------------------------------
  ## Get daily model output
  ##------------------------------------------------------------
  # missing_mod <- purrr::map_lgl( mod$daily, ~identical(., NA ) ) %>% which() %>% names()
  # missing_mod <- purrr::map_lgl( mod$daily, ~identical(., NULL ) ) %>% which() %>% names()
  
  ddf_mod <- lapply( as.list(settings_eval$sitenames),  function(x) dplyr::select( mod$daily[[x]], date, gpp_mod = gpp ) %>% mutate( sitename = x ) ) %>%
    bind_rows()

  ##------------------------------------------------------------
  ## Get observations for evaluation
  ##------------------------------------------------------------
  out_agg <- agg_ddf( agg, ddf_obs, settings_eval, settings_sims )
  
  ## detach
  xdf <- out_agg$xdf

  ##------------------------------------------------------------
  ## Aggregate model output data to annual/monthly/weekly, only for selected sites,
  ## and merge into respective observational data frame 
  ##------------------------------------------------------------
  print("Aggregating model outputs...")

  ## mean across multi-day period
  xdf <- ddf_mod %>% 
  	# mutate( year = year(date), week = week(date) ) %>%
  	mutate( year = year(date), inbin = cut( date, breaks = out_agg$breaks_xdf, right = FALSE ) ) %>%
    group_by( sitename, inbin ) %>%
    summarise( gpp_mod_mean = mean( gpp_mod, na.rm = TRUE ), gpp_mod_min = min( gpp_mod, na.rm = TRUE ), gpp_mod_max = max( gpp_mod, na.rm = TRUE ), n_mod = sum(!is.na(gpp_mod)) ) %>%
    dplyr::rename( gpp_mod = gpp_mod_mean ) %>%
  	right_join( xdf, by = c("sitename", "inbin") )

  ##------------------------------------------------------------
  ## Get mean seasonal cycle (by week (or X-day period) of year)
  ##------------------------------------------------------------
  print("Evaluate mean seasonal cycle by X-day periods...")
  meanxoydf <- xdf %>%  mutate( xoy = yday(inbin) ) %>%
                group_by( sitename, xoy ) %>% 
                summarise( obs_mean = mean( gpp_obs, na.rm=TRUE ), obs_min = min( gpp_obs, na.rm=TRUE ), obs_max = max( gpp_obs, na.rm=TRUE ),
                           mod_mean = mean( gpp_mod, na.rm=TRUE ), mod_min = min( gpp_mod, na.rm=TRUE ), mod_max = max( gpp_mod, na.rm=TRUE )
                           ) %>%
                mutate( obs_min = ifelse( is.infinite(obs_min), NA, obs_min ), obs_max = ifelse( is.infinite(obs_max), NA, obs_max ) ) %>%
                mutate( obs_mean = interpol_lin(obs_mean), obs_min = interpol_lin(obs_min), obs_max = interpol_lin( obs_max ), site=sitename )

  ##------------------------------------------------------------
  ## Get IXV (inter-day variability) as daily value minus mean by site and DOY
  ##------------------------------------------------------------
	print("Evaluate inter-X-day variability...")
	ixvdf <- xdf %>%  mutate( xoy = yday(inbin) ) %>%
            left_join( dplyr::rename( meanxoydf, gpp_mod_mean = mod_mean, gpp_obs_mean = obs_mean ), by = c("sitename", "xoy") ) %>%
						mutate( gpp_mod = gpp_mod - gpp_mod_mean, gpp_obs = gpp_obs - gpp_obs_mean ) %>%
						dplyr::select( -gpp_obs_mean, -gpp_mod_mean, -obs_min, -obs_max, -mod_min, -mod_max )
	
	stats_bysite <- ixvdf %>% 
	  group_by( sitename ) %>%
	  nest() %>%
	  mutate( nxdays_obs = purrr::map( data, ~sum(!is.na( .$gpp_obs )  ) ), nxdays_mod = purrr::map( data, ~sum(!is.na( .$gpp_mod )  ) ) ) %>%
	  unnest( nxdays_obs, nxdays_mod ) %>%
	  filter( nxdays_obs > 2 & nxdays_mod > 2 ) %>%
	  mutate( linmod = purrr::map( data, ~lm( gpp_obs ~ gpp_mod, data = . ) ),
	          stats  = purrr::map( data, ~get_stats( .$gpp_mod, .$gpp_obs ) ) ) %>%
	  mutate( data   = purrr::map( data, ~add_fitted(.) ) ) %>%
	  unnest( stats )
	
	stats <- with( ixvdf, get_stats( gpp_mod, gpp_obs ) )

	print("Done with eval_ixv_sofun().")

  ##------------------------------------------------------------
  ## Mod. vs. obs. of m of IXV correlation: x_(x,i) - mean_x( x_(x,i) )
  ##------------------------------------------------------------
  plot_modobs_anomalies_xdaily <- function( makepdf = FALSE ){  # using ixvdf, ixvdf_stats
    # source("analyse_modobs.R")
    if (makepdf && !dir.exists(settings_eval$dir_figs)) system( paste0( "mkdir -p ", settings_eval$dir_figs))
    if (makepdf) pdf( paste0( settings_eval$dir_figs, "/modobs_anomalies_xdaily.pdf" ) )
      modobs_anomalies_xdaily <- with( ixvdf, analyse_modobs(
        gpp_mod, 
        gpp_obs, 
        col=rgb(0,0,0,0.05), 
        ylab = expression( paste("observed GPP (gC m"^-2, "d"^-1, ")" ) ), 
        xlab = expression( paste("simulated GPP (gC m"^-2, "d"^-1, ")" ) )
        ))
      out <- ixvdf_stats %>%  mutate( purrr::map( data, ~lines( fitted ~ gpp_mod, data = ., col=rgb(0,0,1,0.1) ) ) )  # to have it sorted: %>% mutate( data = purrr::map( data, ~arrange( ., gpp_mod ) ) )
      title( "IXV correlation" )
    if (makepdf) dev.off()


    ## histogram of X-daily anomalies from mean seasonal cycle based on XOY
    ##------------------------------------------------------------
    if (makepdf && !dir.exists(settings_eval$dir_figs)) system( paste0( "mkdir -p ", settings_eval$dir_figs))
    if (makepdf) pdf( paste0( settings_eval$dir_figs, "/hist_anomalies_xdaily.pdf" ) )
      par(las=1)
      ## do not plot, to get density
      outhist1 <- with( ixvdf, hist( gpp_obs, breaks = 20, plot = FALSE ) )
      outhist2 <- with( ixvdf, hist( gpp_mod, breaks = outhist1$breaks, freq = FALSE, plot = FALSE ) )

      ## plot with proper y-axis
      plot(outhist1, freq = FALSE, col = rgb(0,0,0,0.3), main = "Anomalies in X-day periods", xlab = expression( paste("GPP anomaly (gC m"^-2, "d"^-1, ")" ) ), ylim = c(0,max(outhist1$density, outhist2$density)))
      plot(outhist2, freq = FALSE, add=TRUE, col = rgb(1,0,0,0.3))

      mtext( bquote( sigma[obs] == .(format( sd(ixvdf$gpp_obs, na.rm = TRUE), digits = 3)) ), side=3, adj=0, line=0 ) 
      mtext( bquote( sigma[mod] == .(format( sd(ixvdf$gpp_mod, na.rm = TRUE), digits = 3)) ), side=3, adj=0, line=-1 )  
      legend("topright", c("observed", "modelled"), fill = c(rgb(0,0,0,0.3), rgb(1,0,0,0.3)), bty = "n")
    if (makepdf) dev.off()
    return(modobs_anomalies_xdaily)
  }

  plot <- list( modobs_anomalies_xdaily = plot_modobs_anomalies_xdaily )
  
	return( list( stats=stats, stats_bysite=stats_bysite ) )
}


#' Get observational data for model evaluation
#'
#' Gets observational data for model evaluation, given the benchmarking variable, and data source used 
#' for the evaluation. This information is specified in the evaluation settings (argument \code{settings_eval}).
#'
#' @param agg An integer specifying the number of days over which to aggregate
#' @param settings_eval A list specifying evaluation settings (see vignette eval_sofun.pdf for more information and examples)
#' @param settings_sims A named list containing the simulation settings (see vignette_rsofun.pdf for more information and examples)
#' @param overwrite (Optional) A logical specifying whether temporary data stored in \code{./tmpdir} should be overwritten. Defaults to \code{TRUE}.
#'
#' @return A list containing data frames of observed values aggregated to several temporal scales (ddf for daily, xdf for X-daily, mdf for monthly, adf for annual).
#' @export
#'
#' @examples obs <- get_obs_eval( settings_eval, settings_sims )
#' 
agg_ddf <- function( agg, ddf, obs_eval, settings_eval, settings_sims ){

  ##------------------------------------------------------------
  ## Aggregate to multi-day periods
  ## periods should start with the 1st of January each year, otherwise can't compute mean seasonal cycle
  ##------------------------------------------------------------
  ## Generate vector of starting dates of X-day periods, making sure the 1st of Jan is always the start of a new period
  listyears <- seq( ymd("1990-01-01"), ymd("2018-01-01"), by = "year" )                  
  breaks <- purrr::map( as.list(listyears), ~seq( from=., by=paste0( agg, " days"), length.out = ceiling(365 / agg)) ) %>% Reduce(c,.)

  ## take mean across periods
  xdf <- ddf %>% mutate( inbin = cut( date, breaks = breaks, right = FALSE ) ) %>%
                 group_by( sitename, inbin ) %>%
                 summarise( gpp_obs_mean = mean( gpp_obs, na.rm = TRUE ), gpp_obs_min = min( gpp_obs, na.rm = TRUE ), gpp_obs_max = max( gpp_obs, na.rm = TRUE ), n_obs = sum(!is.na(gpp_obs)) ) %>%
                 dplyr::rename( gpp_obs = gpp_obs_mean ) %>%
                 mutate( gpp_obs = ifelse(is.nan(gpp_obs), NA, gpp_obs ), gpp_obs_min = ifelse(is.infinite(gpp_obs_min), NA, gpp_obs_min ), gpp_obs_max = ifelse(is.infinite(gpp_obs_max), NA, gpp_obs_max ) )

  return( list( xdf=xdf, breaks_xdf = breaks ) )
}

add_fitted <- function( data ){
  linmod <- lm( gpp_obs ~ gpp_mod, data = data, na.action = "na.exclude" )
  data$fitted <- fitted( linmod )
  return(data)  
}

add_fitted_alternativenames <- function( data ){
  linmod <- lm( obs_mean ~ mod_mean, data = data, na.action = "na.exclude" )
  data$fitted <- fitted( linmod )
  return(data)  
}

interpol_lin <- function(vec){
  out <- try( approx( seq(length(vec)), vec, xout = seq(length(vec)) )$y )
  if (class(out)=="try-error") out <- vec * NA
  return(out)
}
