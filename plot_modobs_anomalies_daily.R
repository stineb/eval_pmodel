##------------------------------------------------------------
## Mod. vs. obs. of IDV (interday variability) correlation: x_(d,i) - mean_d( x_(d,i) )
##------------------------------------------------------------
plot_modobs_anomalies_daily <- function( out_eval, label="", makepdf = FALSE ){
  dir_figs <- "./fig/"
  if (makepdf && !dir.exists(dir_figs)) system( paste0( "mkdir -p ", dir_figs))
  if (makepdf) filn <- paste0( dir_figs, "/modobs_anomalies_daily_", label, ".pdf" )
  if (makepdf) print( paste( "Plotting to file:", filn ) )
  if (makepdf) pdf( filn )
    modobs_anomalies_daily <- with( out_eval$data$idvdf, rsofun::analyse_modobs(
      gpp_mod,
      gpp_obs,
      col=rgb(0,0,0,0.05),
      ylab = expression( paste("observed GPP (gC m"^-2, "d"^-1, ")" ) ),
      xlab = expression( paste("simulated GPP (gC m"^-2, "d"^-1, ")" ) )
      ))
    out <- out_eval$data$idvdf_stats %>%  mutate( purrr::map( data, ~lines( fitted ~ gpp_mod, data = ., col=rgb(0,0,1,0.05) ) ) )  # to have it sorted: %>% mutate( data = purrr::map( data, ~arrange( ., gpp_mod ) ) )
    title( label )
  if (makepdf) dev.off()
  
  ## histogram of daily anomalies from mean seasonal cycle based on DOY
  ##------------------------------------------------------------
  if (makepdf && !dir.exists(dir_figs)) system( paste0( "mkdir -p ", dir_figs))
  if (makepdf) filn <- paste0( dir_figs, "/hist_anomalies_daily_", label, ".pdf" )
  if (makepdf) print( paste( "Plotting to file:", filn ) )
  if (makepdf) pdf( filn )
    par(las=1)
    out <- with( out_eval$data$idvdf, hist( gpp_obs, breaks = 50, col = rgb(0,0,0,0.3), freq = FALSE, main = label, ylim = c(0,0.6), xlab = expression( paste("GPP anomaly (gC m"^-2, "d"^-1, ")" ) ) ) )
    with( out_eval$data$idvdf, hist( gpp_mod, breaks = out$breaks, col = rgb(1,0,0,0.3), freq = FALSE, add = TRUE ) )
    mtext( bquote( sigma[obs] == .(format( sd(out_eval$data$idvdf$gpp_obs, na.rm = TRUE), digits = 3)) ), side=3, adj=0, line=0 ) 
    mtext( bquote( sigma[mod] == .(format( sd(out_eval$data$idvdf$gpp_mod, na.rm = TRUE), digits = 3)) ), side=3, adj=0, line=-1 )  
    legend("topright", c("observed", "modelled"), fill = c(rgb(0,0,0,0.3), rgb(1,0,0,0.3)), bty = "n")
  if (makepdf) dev.off()

  return(modobs_anomalies_daily)

}   
