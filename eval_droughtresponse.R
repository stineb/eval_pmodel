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
eval_droughtresponse <- function( df_modobs, df_isdrought, obsvarnam, modvarnam ){

  ## Get information which dates are classified as 'soil moisture drought' based on fLUE (Stocker et al., 2018)

  ## Rearrange data. Function returns list( df_dday, df_dday_aggbydday )
  out_align <- align_events( df_modobs, df_isevent, dovars, before=20, after=100, do_norm=FALSE, normbin=2 )

}