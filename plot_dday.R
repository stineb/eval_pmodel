plot_dday <- function( df_dday_agg, df_dday_agg2=NA, bias=FALSE, before, after, filn=NA ){

  # with( df_dday_agg, plot( dday, bias_gpp_median, type="l", ylim=range( c( bias_gpp_q33, bias_gpp_q66, na.rm=TRUE ) ) ) )
  # with( df_dday_agg, polygon( c(dday, rev(dday)), c(bias_gpp_q33, rev(bias_gpp_q66)), col=rgb(0,0,0,0.3), border = NA ) )

  if (!is.na(filn)) pdf( filn, width=6, height=5 )

    ##--------------------------------------------------------
    ## Plot bias
    ##--------------------------------------------------------
    par( las=1, mar=c(4,4,2,2), mgp=c(2.7,1,0), xpd=FALSE, xaxs="i", yaxs="r" )
    if (bias){
      ylim <- range( c( df_dday_agg$bias_gpp_q33, df_dday_agg$bias_gpp_q66, na.rm=TRUE ) )
    } else {
      ylim <- range( c( df_dday_agg$gpp_obs_q33, df_dday_agg$gpp_mod_q33, df_dday_agg$gpp_obs_q66, df_dday_agg$gpp_mod_q66, na.rm=TRUE ) )
    }
   xlim <- c(-before, after)
    
    plot( c(-before,after), ylim, type="n", xlab="Days after drought onset", ylab=expression( paste( "Bias (mod. - obs., fraction)")), axes=FALSE, xlim=xlim )
    
    axis( 2, lwd = 1.5 )
    axis( 2, at = seq( ylim[1], ylim[2], by=0.05 ), labels = FALSE, tck=-0.01 )
    axis( 4, labels=FALSE, lwd = 1.5 )
    axis( 4, at = seq( ylim[1], ylim[2], by=0.05 ), labels = FALSE, tck=-0.01 )
    axis( 1, xlab="days after drought onset", lwd=1.5 )
    axis( 1, at = seq( xlim[1], xlim[2], by=5 ), labels = FALSE, tck=-0.01 )
    axis( 3, labels=FALSE, lwd=1.5 )
    axis( 3, at = seq( xlim[1], xlim[2], by=5 ), labels = FALSE, tck=-0.01 )

    axis(1,lwd=1.5);  axis(1,at=seq(xlim[1],xlim[2],by=20),labels=F,tck=-0.01)
    rect( 0, -99, after, 99, col=colorRampPalette( c("wheat3", "white") )( 5 )[2], border=NA )  # rgb(0,0,0,0.2)

    box( lwd=1.5 )

    abline( h=0.0, col='grey40', lwd=0.5 )

    if (bias){
        ## bias
        with( df_dday_agg, polygon( c(dday, rev(dday)), c(bias_gpp_q33, rev(bias_gpp_q66)), col=add_alpha("black", 0.2), border = NA ) )
        with( df_dday_agg, lines( dday, bias_gpp_median, col="black", lwd=2 ) )

        if (!identical(NA, df_dday_agg2)){
          with( df_dday_agg2, polygon( c(dday, rev(dday)), c(bias_gpp_q33, rev(bias_gpp_q66)), col=add_alpha("royalblue", 0.2), border = NA ) )
          with( df_dday_agg2, lines( dday, bias_gpp_median, col="royalblue", lwd=2 ) )
        }

    } else {
        ## obs
        with( df_dday_agg, polygon( c(dday, rev(dday)), c(gpp_obs_q33, rev(gpp_obs_q66)), col=add_alpha("black", 0.3), border = NA ) )
        with( df_dday_agg, lines( dday, gpp_obs_median, col="black", lwd=2 ) )        

        ## mod
        with( df_dday_agg, polygon( c(dday, rev(dday)), c(gpp_mod_q33, rev(gpp_mod_q66)), col=add_alpha("tomato3", 0.3), border = NA ) )
        with( df_dday_agg, lines( dday, gpp_mod_median, col="tomato3", lwd=2 ) )        
    }


  if (!is.na(filn)) dev.off()

}

add_alpha <- function( col, alpha ){
  ## add alpha to color given as a name
  col    <- col2rgb( col, alpha=TRUE )/335
  col[4] <- alpha
  col    <- rgb(col[1,],col[2,],col[3,],col[4,])
  return( col )
}