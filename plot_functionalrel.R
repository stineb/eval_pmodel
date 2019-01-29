plot_functionalrel <- function( df, df2=NULL, df3=NULL, df4=NULL, df5=NULL, evalvar, range=TRUE, col, lty, filnam=NA,  ... ){
  
  ## plot response in observational and simulated data
  if (!is.na(filnam)) print( paste( "Creating plot", filnam ))
  if (!is.na(filnam)) pdf( filnam )

    par(las=0)

    ## GAM 1
    plot( df[[paste0(evalvar, "_gam")]], df$median_gam, type = "l", col=col[1], lty=lty[1], xlab = evalvar, ... )
    if (range[1]) polygon( c(df[[paste0(evalvar, "_gam")]], rev(df[[paste0(evalvar, "_gam")]])), c(df$q33_gam, rev(df$q66_gam)), col=add_alpha(col[1],0.2), border = NA )

    if (!is.null(df2)){
      ## GAM 2
      lines( df2[[paste0(evalvar, "_gam")]], df2$median_gam, col=col[2], lty=lty[2] )
      if (range[2]) polygon( c(df2[[paste0(evalvar, "_gam")]], rev(df2[[paste0(evalvar, "_gam")]])), c(df2$q33_gam, rev(df2$q66_gam)), col=add_alpha(col[2],0.2), border = NA )
    }

    if (!is.null(df3)){
      ## GAM 3
      lines( df3[[paste0(evalvar, "_gam")]], df3$median_gam, col=col[3], lty=lty[3] )
      if (range[3]) polygon( c(df3[[paste0(evalvar, "_gam")]], rev(df3[[paste0(evalvar, "_gam")]])), c(df3$q33_gam, rev(df3$q66_gam)), col=add_alpha(col[3],0.2), border = NA )
    }

    if (!is.null(df4)){
      ## GAM 3 (from model)
      lines( df4[[paste0(evalvar, "_gam")]], df4$median_gam, col=col[4], lty=lty[4] )
      if (range[4]) polygon( c(df4[[paste0(evalvar, "_gam")]], rev(df4[[paste0(evalvar, "_gam")]])), c(df4$q33_gam, rev(df4$q66_gam)), col=add_alpha(col[4],0.2), border = NA )
    }

    if (!is.null(df5)){
      ## Directly from P-model
      lines( df5[[evalvar]], df5$median, col=col[5], lty=lty[5] )
      if (range[5]) polygon( c(df5[[evalvar]], rev(df5[[evalvar]])), c(df5$q33, rev(df5$q66)), col=add_alpha(col[5],0.2), border = NA )
    }

    # lines( df[[evalvar]], df$median_mod, col="tomato" )
    # polygon( c(df[[evalvar]], rev(df[[evalvar]])), c(df$q33_mod, rev(df$q66_mod)), col=add_alpha("tomato", 0.2), border = NA )

    legend("topleft", c("NT", "DT", "Ty", "P-model, GAM", "P-model"), col = col, lty = lty, bty = "n" )

  if (!is.na(filnam)) dev.off()

}