plot_modobs_spatial_annual <- function( out_eval, annual_pooled_stats = NA, spatial_stats = NA, label = "", makepdf = FALSE, ... ){
  
  if (!identical(out_eval$data$linmod_meandf, NA)){
    
    dir_figs <- "./fig/"
    if (makepdf && !dir.exists(dir_figs)) system( paste0( "mkdir -p ", dir_figs))
    if (makepdf) filn <- paste0( dir_figs, "/modobs_spatial_annual_", label, ".pdf" )
    if (makepdf) print( paste( "Plotting to file:", filn ) )
    if (makepdf) pdf( filn, width = 5, height = 5 )
    
    par(las=1, mar=c(4,4.5,4,1))
    
    ## set up plotting and add linear regression line for means by site
    with( out_eval$data$meandf, 
        plot( gpp_mod, gpp_obs, 
            pch=16, col=rgb(0,0,0,0.5), type = "n", 
            ylab = expression( paste("observed GPP (gC m"^-2, "yr"^-1, ")" ) ), 
            xlab = expression( paste("simulated GPP (gC m"^-2, "yr"^-1, ")" ) ),
            main = label,
            ... ) )
    abline( out_eval$data$linmod_meandf, col="red")
    lines(c(-9999,9999), c(-9999,9999), lty=3)
    
    ## plot black regression lines of annual values within sites
    out <- out_eval$data$adf_stats %>% mutate( purrr::map( data, ~lines( fitted ~ gpp_mod, data = . ) ) )  # to have it sorted: %>% mutate( data = purrr::map( data, ~arrange( ., gpp_mod ) ) )
        
    ## Add annotations for statistics of annual values (pooled)
    if (!identical(NA, out_eval$metrics$gpp$fluxnet2015$annual_pooled)) mtext( bquote( italic(R)^2 == .(format( out_eval$metrics$gpp$fluxnet2015$annual_pooled$rsq,  digits = 2) ) ), adj = 1, cex = 0.8, line=2 )
    if (!identical(NA, out_eval$metrics$gpp$fluxnet2015$annual_pooled)) mtext( paste0( "RMSE = ",       format( out_eval$metrics$gpp$fluxnet2015$annual_pooled$rmse, digits = 3 ) ),  adj = 1, cex = 0.8, line=1 )
    
    ## Add annotations for statistics of means by site (~spatial)
    if (!identical(NA, out_eval$metrics$gpp$fluxnet2015$spatial))       mtext( bquote( italic(R)^2 == .(format( out_eval$metrics$gpp$fluxnet2015$spatial$rsq,       digits = 2) ) ), adj = 0, cex = 0.8, line=2, col="red" )
    if (!identical(NA, out_eval$metrics$gpp$fluxnet2015$spatial))       mtext( paste0( "RMSE = ",       format( out_eval$metrics$gpp$fluxnet2015$spatial$rmse,      digits = 3 ) ),  adj = 0, cex = 0.8, line=1, col="red" )
    # if (!identical(NA, out_eval$metrics$gpp$fluxnet2015$spatial))       mtext( paste0( "slope = ",      format( out_eval$metrics$gpp$fluxnet2015$spatial$meanslope, digits = 3 ) ),  adj = 0, cex = 0.8, col="red" )
    
    if (makepdf) dev.off()


    ## And for each site individually

    # ## Histogram of slopes
    # ##------------------------------------------------------------
    # ## (Uncomment to plot as inset in spatial-IAV plot) 
    # # u <- par("usr")
    # # v <- c(
    # #   grconvertX(u[1:2], "user", "ndc"),
    # #   grconvertY(u[3:4], "user", "ndc")
    # # )
    # # v_orig <- v
    # # v <- c( v[1]+0.03, v[1]+0.2*v[2], v[3]+0.50*v[4], v[3]+0.72*v[4] )
    # # par( fig=v, new=TRUE, mar=c(0,0,0,0), mgp=c(3,0.5,0) )
    # if (makepdf && !dir.exists(dir_figs)) system( paste0( "mkdir -p ", dir_figs))
    # if (makepdf) pdf( paste0( dir_figs, "/hist_slopes_anomalies_annual.pdf" ) )
    #   hist( out_eval$data$annual_bysite_stats$slope, xlim=c(-5,5), cex.axis=0.7, axes=FALSE, col="grey70", main="", breaks = 50, xlab="slope" )
    #   abline( v=1.0, col="red" )
    #   axis( 1, cex.axis=1.0, xlab="slope" )
    #   title( "Slopes of annual regressions" )
    # if (makepdf) dev.off()
    # 
    # ## Histogram of R2
    # ##------------------------------------------------------------
    # ## (Uncomment to plot as inset in spatial-IAV plot) 
    # # u <- par("usr")
    # # v <- c(
    # #   grconvertX(u[1:2], "user", "ndc"),
    # #   grconvertY(u[3:4], "user", "ndc")
    # # )
    # # v_orig <- v
    # # v <- c( v[1]+0.03, v[1]+0.2*v[2], v[3]+0.50*v[4], v[3]+0.72*v[4] )
    # # par( fig=v, new=TRUE, mar=c(0,0,0,0), mgp=c(3,0.5,0) )
    # if (makepdf && !dir.exists(dir_figs)) system( paste0( "mkdir -p ", dir_figs))
    # if (makepdf) pdf( paste0( dir_figs, "/hist_r2_anomalies_annual.pdf" ) )
    #   hist( out_eval$data$annual_bysite_stats$rsq, xlim=c(-1,1), cex.axis=0.7, axes=FALSE, col="grey70", main="", breaks = 12, xlab= bquote( italic(R)^2 ) )
    #   abline( v=1.0, col="red" )
    #   axis( 1, cex.axis=1.0, xlab = bquote( italic(R)^2 ) )
    #   title( bquote( bold(Slopes ~ of ~ italic(R)^2) ) )
    # if (makepdf) dev.off()

  } else {
    print("plot_modobs_spatial_annual(): Not enough sites to get spatial correlations.")
  }
  
  
}