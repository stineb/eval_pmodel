plot_modobs_spatial_annual2 <- function( out_eval ){
  
  get_start_end <- function(df){
    df_start <- df %>% 
      arrange(mod) %>% 
      drop_na(mod, fitted) %>% 
      slice(1)
    df_end <- df %>% 
      arrange(desc(mod)) %>% 
      drop_na(mod, fitted) %>% 
      slice(1)
    out <- tibble(
      xmin = df_start$mod, 
      xmax = df_end$mod,
      ymin = df_start$fitted,
      ymax = df_end$fitted )
    return(out)
  }
  df <- out_eval$gpp$fluxnet2015$data$adf_stats %>% 
    mutate(start_end = purrr::map(data, ~get_start_end(.))) %>% 
    tidyr::unnest(start_end)

  rsq_lab_annual <-  format(out_eval$gpp$fluxnet2015$metrics$annual_pooled$rsq, digits = 2)
  rmse_lab_annual <- format(out_eval$gpp$fluxnet2015$metrics$annual_pooled$rmse, digits = 3)
  
  rsq_lab_spatial <-  format(out_eval$gpp$fluxnet2015$metrics$spatial$rsq, digits = 2)
  rmse_lab_spatial <- format(out_eval$gpp$fluxnet2015$metrics$spatial$rmse, digits = 3)
  
  gg <- df %>% 
    ggplot() +
    geom_segment(aes(x=xmin, y=ymin, xend=xmax, yend=ymax)) +
    geom_line(data = fortify(out_eval$gpp$fluxnet2015$data$linmod_meandf), aes(x = mod, y = .fitted), color="red") +
    geom_abline(intercept=0, slope=1, linetype="dotted") +
    theme_classic() +
    xlim(0,NA) +
    ylim(0,NA) +
    labs(
      subtitle = bquote( bold("Annual:") ~ italic(R)^2 == .(rsq_lab_annual) ~~
                           RMSE == .(rmse_lab_annual) ~ "\n" ~
                         bold("Spatial:") ~ italic(R)^2 == .(rsq_lab_spatial) ~~
                           RMSE == .(rmse_lab_spatial) ) ) #+
    #geom_text_repel(aes(x = xmin, y = ymin, label = sitename), size = 2, segment.colour = "white")
  

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
    #   hist( out_eval$gpp$fluxnet2015$data$annual_bysite_stats$slope, xlim=c(-5,5), cex.axis=0.7, axes=FALSE, col="grey70", main="", breaks = 50, xlab="slope" )
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
    #   hist( out_eval$gpp$fluxnet2015$data$annual_bysite_stats$rsq, xlim=c(-1,1), cex.axis=0.7, axes=FALSE, col="grey70", main="", breaks = 12, xlab= bquote( italic(R)^2 ) )
    #   abline( v=1.0, col="red" )
    #   axis( 1, cex.axis=1.0, xlab = bquote( italic(R)^2 ) )
    #   title( bquote( bold(Slopes ~ of ~ italic(R)^2) ) )
    # if (makepdf) dev.off()
# 
#   } else {
#     print("plot_modobs_spatial_annual(): Not enough sites to get spatial correlations.")
#   }
  
  return(gg)
}
