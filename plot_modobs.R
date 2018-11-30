##------------------------------------------------------------
## Mod. vs. obs for daily values (absolute)
##------------------------------------------------------------
## observed vs. modelled
plot_modobs_daily <- function( out_eval, subtitle = "", label = "", makepdf = FALSE, ... ){ 
  dir_figs <- "./fig/"
  if (makepdf && !dir.exists(dir_figs)) system( paste0( "mkdir -p ", dir_figs))
  if (makepdf) filn <- paste0( dir_figs, "/modobs_daily_", label, ".pdf" )
  if (makepdf) print( paste( "Plotting to file:", filn ) )
  if (makepdf) pdf( filn )
  modobs_ddf <- with( out_eval$data$ddf, 
    analyse_modobs( 
      gpp_mod, 
      gpp_obs, 
      heat=TRUE, 
      ylab = expression( paste("observed GPP (gC m"^-2, "d"^-1, ")" ) ), 
      xlab = expression( paste("simulated GPP (gC m"^-2, "d"^-1, ")" ) ),
      plot.title = label,
      plot.subtitle = subtitle,
      ...
    ) )
  if (makepdf) dev.off()
  return( modobs_ddf )
}


##------------------------------------------------------------
## Mod. vs. obs. for ggregated values (absolute) aggregated to X-day periods
##------------------------------------------------------------
## observed vs. modelled
plot_modobs_xdaily <- function( out_eval, subtitle = "", label = "", makepdf = FALSE, ... ){    # using xdf
  dir_figs <- "./fig/"
  if (makepdf && !dir.exists(dir_figs)) system( paste0( "mkdir -p ", dir_figs))
  if (makepdf) filn <- paste0( dir_figs, "/modobs_xdaily_", label, ".pdf" )
  if (makepdf) print( paste( "Plotting to file:", filn ) )
  if (makepdf) pdf( filn )
  modobs_xdf <- with( out_eval$data$xdf, 
    analyse_modobs( 
      gpp_mod, 
      gpp_obs, 
      heat=TRUE, 
      ylab = expression( paste("observed GPP (gC m"^-2, "d"^-1, ")" ) ), 
      xlab = expression( paste("simulated GPP (gC m"^-2, "d"^-1, ")" ) ),
      plot.title = label,
      ...
      ) )
  if (makepdf) dev.off()
  return(modobs_xdf)
}


##------------------------------------------------------------
## Mod. vs. obs for monthly values (absolute)
##------------------------------------------------------------
## observed vs. modelled
plot_modobs_monthly <- function( out_eval, subtitle = "", label = "", makepdf = FALSE, ... ){  # using mdf
  dir_figs <- "./fig/"
  if (makepdf && !dir.exists(dir_figs)) system( paste0( "mkdir -p ", dir_figs))
  if (makepdf) filn <- paste0( dir_figs, "/modobs_monthly_", label, ".pdf" )
  if (makepdf) print( paste( "Plotting to file:", filn ) )
  if (makepdf) pdf( filn )
  modobs_mdf <- with( out_eval$data$mdf, 
    analyse_modobs( 
      gpp_mod, 
      gpp_obs, 
      heat = TRUE, 
      ylab = expression( paste("observed GPP (gC m"^-2, "d"^-1, ")" ) ), 
      xlab = expression( paste("simulated GPP (gC m"^-2, "d"^-1, ")" ) ),
      plot.title = label,
      ...
    ) )
  if (makepdf) dev.off()
  return( modobs_mdf )
}
