##------------------------------------------------------------
## Wrapper for mean seasonality by site (daily) for all climate zones
##------------------------------------------------------------
plot_by_doy_allzones <- function( out_eval, out_eval1 = NA, out_eval2 = NA, makepdf = FALSE, label="", col2="royalblue", col3="darkgoldenrod", lab1="NA", lab2="NA", lab3="NA" ){
  system( "mkdir -p fig/meandoy_byzone" )

  # tmp <- purrr::map( out_eval$data$meandoydf_byclim_stats$data, ~plot_by_doy_byzone(., out_eval1 = out_eval1, makepdf = makepdf, pattern = pattern ) )
  tmp <- purrr::map( 
  	as.list(seq(nrow(out_eval$data$meandoydf_byclim_stats))) , 
  	~plot_by_doy_byzone( 
  		out_eval$data$meandoydf_byclim_stats$data[[.]], 
  		df2 = out_eval1$data$meandoydf_byclim_stats$data[[.]], 
  		df3 = out_eval2$data$meandoydf_byclim_stats$data[[.]], 
  		makepdf = makepdf, label=label, col2=col2, col3=col3, lab1=lab1, lab2=lab2, lab3=lab3 
  		) 
  	)
}


##------------------------------------------------------------
## Mean seasonality with aggregated data from one climate zone 
##------------------------------------------------------------
plot_by_doy_byzone <- function( df, df2 = NA, df3 = NA, makepdf = FALSE, label="", col2="royalblue", col3="darkgoldenrod", lab1="NA", lab2="NA", lab3="NA" ){

  if (df$nsites[1]>4){
  	dir_figs <- "./fig/meandoy_byzone"
  	if (makepdf && !dir.exists(dir_figs)) system( paste0( "mkdir -p ", dir_figs))
  	if (makepdf) filn <- paste0( dir_figs, "/meandoy_byzone_", stringr::str_replace( df$climatezone[1], " ", "_" ), "_", label, ".pdf" )
  	if (makepdf) print( paste( "Plotting to file:", filn ) )
  	if (makepdf) pdf( filn )
      par(las=1)
      yrange <- range( df$mod_min, df$mod_max, df$obs_min, df$obs_max, na.rm = TRUE )
      plot(  df$doy, df$obs_mean, type="l", ylim = yrange, ylab = expression( paste("simulated GPP (gC m"^-2, "d"^-1, ")" ) ), xlab = "DOY" )
      polygon( c(df$doy, rev(df$doy)), c(df$obs_min, rev(df$obs_max)), border = NA, col = rgb(0,0,0,0.3)  )
      lines( df$doy, df$mod_mean, col="red", lwd=1.75 )
      if (!identical(NA,df2)) lines( df2$doy, df2$mod_mean, col=col2, lwd=1.75, lty=1 )
      if (!identical(NA,df3)) lines( df3$doy, df3$mod_mean, col=col3, lwd=1.75, lty=1 )
      polygon( c(df$doy, rev(df$doy)), c(df$mod_min, rev(df$mod_max)), border = NA, col = rgb(1,0,0,0.3)  )
      title( df$climatezone[1] )
      mtext( bquote( italic(N) == .( df$nsites[1])), side=3, line=1, cex=1.0, adj=1.0 )
      legend("topright", c(lab1, lab2, lab3), col=c("red", col2, col3), lwd = c(1.75,1.75,1.75), bty = "n" )
    if (makepdf) dev.off()
  } else {
    rlang::warn( paste0("plot_by_doy_byzone(): Number of sites below 5 for climate zone ", df$climatezone[1] ) )
  }

}