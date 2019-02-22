plot_linmod <- function( df, df2=NULL, df3=NULL, sitename, fitted, lab_legend=NA, label="", ... ){
  # print(sitename)
  filn <- paste0("fig/spatial_annual_bysite/plot_spatial_annual_", sitename, "_", label, ".pdf")
  print(paste("creating file", filn))
  pdf( filn, width=5, height=5 )
  par(las=1)
  xlim <- range(c(df$gpp_mod,df$gpp_obs), na.rm=TRUE)
  ylim <- xlim
  plot( df$gpp_mod, df$gpp_obs, 
        pch=16, 
        xlab = expression( paste("simulated GPP (gC m"^-2, "yr"^-1, ")" ) ),
        ylab = expression( paste("observed GPP (gC m"^-2, "yr"^-1, ")" ) ), 
        col=add_alpha("royalblue3", 0.7),
        ... )
  if (!is.null(df2)) points( df2$gpp_mod, df2$gpp_obs, col=add_alpha("springgreen3",0.7), pch=16 )
  if (!is.null(df3)) points( df3$gpp_mod, df3$gpp_obs, col="darkgoldenrod", pch=16 )
  
  # linmod <- lm( df$gpp_obs ~ df$gpp_mod )
  df_fitted <- tibble( gpp_mod=df$gpp_mod, fitted=df$fitted ) %>%
    arrange( gpp_mod )
  lines( df_fitted$gpp_mod, df_fitted$fitted, col="royalblue3" )

  df_fitted <- tibble( gpp_mod=df2$gpp_mod, fitted=df2$fitted ) %>%
    arrange( gpp_mod )
  lines( df_fitted$gpp_mod, df_fitted$fitted, col="springgreen3" )
  
  lines(c(-9999,9999), c(-9999,9999), lty=3)
  title(sitename)
  
  if (!is.na(lab_legend)) legend( "bottomright", lab_legend, col=c("royalblue3", "springgreen3"), bty="n", lty=1 )
  
  dev.off()
}