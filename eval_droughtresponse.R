#' Evaluate the model's drought response
#'
#' xxx
#' 
#' @param df A data frame containing all data for all sites, modelled and observed
#' @param obsvarnam The variable name of the observational data, given as a column in \code{df}
#' @param modvarnam The variable name of the modelled data, given as a column in \code{df}
#'
#' @return A data frame of data aligned by drought events and aggregated across events (for different quantiles)
#' @export
#'
#' @examples xxx
#' 
eval_droughtresponse <- function( mod, obs, path_flue, before, after, leng_threshold, nbins, do_norm=FALSE ){

  ## Get fLUE Stocker et al., 2018 publicly available data
  df_flue <- read_csv( path_flue ) %>% 
    dplyr::select(-year, -doy, -cluster) %>% 
    dplyr::rename( isevent = is_flue_drought )

  ## do some weird cleaning and do some weird cleaning
  mod <- mod %>% 
    # na.omit.list() %>% 
    bind_rows() %>% 
    dplyr::select( site=sitename, date, gpp ) %>% 
    dplyr::rename(gpp_mod = gpp)

  ## Get observational data
  obs <- obs %>% 
    dplyr::select( site=sitename, date, gpp ) %>% 
    dplyr::rename(gpp_obs = gpp)

  ## combine data frames into one
  df_modobs <- obs %>% 
    left_join( mod, by=c("site", "date") ) %>% 
    mutate( bias_gpp = gpp_mod - gpp_obs )

  ## Rearrange data. Function returns list( df_dday, df_dday_aggbydday )
  dovars <- colnames( dplyr::select( df_modobs, -date, -site ) )

  out_align <- align_events( df_modobs, df_flue, dovars, leng_threshold, before, after, nbins, do_norm=do_norm )

  return( out_align$df_dday_aggbydday )
}

na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
