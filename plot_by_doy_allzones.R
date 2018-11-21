##------------------------------------------------------------
## Wrapper for mean seasonality by site (daily) for all climate zones
##------------------------------------------------------------
plot_by_doy_allzones <- function( out_eval, dashed1 = NA, dashed2 = NA, makepdf = FALSE, pattern="" ){
  system( "mkdir -p fig/meandoy_byzone" )

  # tmp <- purrr::map( out_eval$data$meandoydf_byclim_stats$data, ~plot_by_doy_byzone(., dashed1 = dashed1, makepdf = makepdf, pattern = pattern ) )
  tmp <- purrr::map( 
  	as.list(seq(nrow(out_eval$data$meandoydf_byclim_stats))) , 
  	~plot_by_doy_byzone( 
  		out_eval$data$meandoydf_byclim_stats$data[[.]], 
  		dashed1 = NA, #dashed1$data$meandoydf_byclim_stats$data[[.]], 
  		dashed2 = NA, #dashed2$data$meandoydf_byclim_stats$data[[.]], 
  		makepdf = makepdf, 
  		pattern = pattern ) 
  	)
}


##------------------------------------------------------------
## Mean seasonality with aggregated data from one climate zone 
##------------------------------------------------------------
plot_by_doy_byzone <- function( df, dashed1 = NA, dashed2 = NA, makepdf = FALSE, pattern="" ){

  if (df$nsites[1]>4){
  	dir_figs <- "./fig/meandoy_byzone"
  	if (makepdf && !dir.exists(dir_figs)) system( paste0( "mkdir -p ", dir_figs))
  	if (makepdf) filn <- paste0( dir_figs, "/meandoy_byzone_", df$climatezone[1], "_", pattern, ".pdf" )
  	if (makepdf) print( paste( "Plotting to file:", filn ) )
  	if (makepdf) pdf( filn )
      par(las=1)
      yrange <- range( df$mod_min, df$mod_max, df$obs_min, df$obs_max, na.rm = TRUE )
      plot(  df$doy, df$obs_mean, type="l", ylim = yrange, ylab = expression( paste("simulated GPP (gC m"^-2, "d"^-1, ")" ) ), xlab = "DOY" )
      polygon( c(df$doy, rev(df$doy)), c(df$obs_min, rev(df$obs_max)), border = NA, col = rgb(0,0,0,0.3)  )
      lines( df$doy, df$mod_mean, col="red", lwd=1.75 )
      if (!is.na(dashed1)) lines( dashed1$doy, dashed1$mod_mean, col="red", lwd=0.75, lty=1 )
      if (!is.na(dashed2)) lines( dashed2$doy, dashed2$mod_mean, col="springgreen3", lwd=0.75, lty=1 )
      polygon( c(df$doy, rev(df$doy)), c(df$mod_min, rev(df$mod_max)), border = NA, col = rgb(1,0,0,0.3)  )
      title( df$climatezone[1] )
      mtext( bquote( italic(N) == .( df$nsites[1])), side=3, line=1, cex=1.0, adj=1.0 )
    if (makepdf) dev.off()
  } else {
    rlang::warn( paste0("plot_by_doy_byzone(): Number of sites below 5 for climate zone ", df$climatezone[1]) )
  }

}