# ##------------------------------------------------------------
# ## Wrapper for mean seasonality by site (daily) for selected sites only sites
# ##------------------------------------------------------------
# plot_by_doy_allsites <- function( out_eval, out_eval1 = NA, out_eval2 = NA, out_eval3 = NA, 
#   makepdf = FALSE, label="", 
#   lab1="", lab2="", lab3="", lab4="", ... 
#   ){   # using meandoydf_stats
#   system( "mkdir -p fig/meandoy_bysite" )
#   # mylist <- readr::read_csv("myselect_fluxnet2015.csv") %>% filter( use==1 ) %>% dplyr::select( -use ) %>% unlist()
#   mylist <- c("AU-Tum", "CA-NS3", "CA-NS6", "CA-Obs", "DE-Geb", "DE-Hai", "DE-Kli", "FI-Hyy", "FR-Fon", "FR-LBr", "FR-Pue", "IT-Cpz", "NL-Loo", "US-Ha1", "US-MMS", "US-UMB", "US-WCr")
  
#   tmp <- purrr::map( 
#     filter( meandoydf_stats, sitename %in% mylist )$data, 
#     ~plot_by_doy_bysite(., makepdf = makepdf) 
#     )


# }

plot_by_doy_allsites <- function( out_eval, out_eval2 = NA, out_eval3 = NA, out_eval4 = NA, 
  makepdf = FALSE, label="", 
  lab1="", lab2="", lab3="", lab4="", ... ){
  system( "mkdir -p fig/meandoy_bysite" )

  ## get info about climate zones of sites
  df_sites <- rsofun::metainfo_Tier1_sites_kgclimate_fluxnet2015 %>%
    mutate( hemisphere = ifelse( lat>0, "north", "south" ) ) %>%
    mutate( kghcode = paste0(koeppen_code, "_", hemisphere))

  df_zones <- df_sites %>% 
    group_by( koeppen_code, hemisphere ) %>% 
    summarise( nsites = n() ) %>% 
    left_join( dplyr::rename( koeppen_legend, koeppen_code = Code ), by = "koeppen_code" ) %>% 
    arrange( -nsites ) %>%
    dplyr::filter( nsites>4 ) %>%
    left_join( dplyr::select(df_sites, kghcode, koeppen_code, hemisphere), by=c("koeppen_code", "hemisphere") )

  # sites_per_zone <- purrr::map( 
  #   as.list(nrow(df_zones)), 
  #   ~dplyr::filter( df_sites, koeppen_code==df_zones$koeppen_code[.] ) & hemisphere==df_zones$hemisphere[.] )
  # names(sites_per_zone) <- df_zones$koeppen_code

  ## subset sites depending
  df <- out_eval$data$meandoydf_stats %>% 
    left_join( dplyr::select(df_sites, sitename, kghcode), by="sitename" ) %>%
    dplyr::filter( kghcode %in% df_zones$kghcode )
  
  df2 <- out_eval2$data$meandoydf_stats %>% 
    left_join( dplyr::select(df_sites, sitename, kghcode), by="sitename" ) %>%
    dplyr::filter( kghcode %in% df_zones$kghcode )
  
  df3 <- out_eval3$data$meandoydf_stats %>% 
    left_join( dplyr::select(df_sites, sitename, kghcode), by="sitename" ) %>%
    dplyr::filter( kghcode %in% df_zones$kghcode )
  
  # df4 <- out_eval4$data$meandoydf_stats %>% 
  #   left_join( dplyr::select(df_sites, sitename, kghcode), by="sitename" ) %>%
  #   dplyr::filter( kghcode %in% df_zones$kghcode )

  list_ndata <- purrr::map(as.list(seq(nrow(df))), ~nrow(filter(out_eval$data$ddf, sitename==df$sitename[[.]])) ) 
  
  # tmp <- purrr::map( out_eval$data$meandoydf_byclim_stats$data, ~plot_by_doy_byzone(., out_eval2 = out_eval2, makepdf = makepdf, pattern = pattern ) )
  tmp <- purrr::map( 
    as.list(seq(nrow(df))), 
    ~plot_by_doy_bysite( 
      df$data[[.]], 
      df2 = df2$data[[.]], 
      df3 = df3$data[[.]],
      # df4 = df4$data[[.]],
      zone = dplyr::filter( rsofun::metainfo_Tier1_sites_kgclimate_fluxnet2015, sitename==df$sitename[[.]] ) %>%
             dplyr::select( koeppen_code ),
      classid = dplyr::filter( rsofun::metainfo_Tier1_sites_kgclimate_fluxnet2015, sitename==df$sitename[[.]] ) %>%
                dplyr::select( classid ),
      ndata = list_ndata[[.]],
      lab1=lab1, lab2=lab2, lab3=lab3, lab4=lab4,
      makepdf=makepdf, label=label, 
      ...
      ) 
    )
}

##------------------------------------------------------------
## Mean seasonality (daily) for one site
##------------------------------------------------------------
plot_by_doy_bysite <- function( df, df2 = NA, df3 = NA, df4 = NA, 
  zone = "", classid = "", ndata = NA, makepdf = FALSE , label="", 
  col2="royalblue", col3="darkgoldenrod", col4=rbeni::add_alpha("springgreen3", 0.5), 
  lab1="", lab2="", lab3="", lab4="" 
  ){

  col2 <- "royalblue"
  col3 <- "darkgoldenrod"
  col4 <- rbeni::add_alpha("springgreen3", 0.5)

  dir_figs <- paste0( "./fig/meandoy_bysite", label, "/", zone, "/" )
  if (makepdf && !dir.exists(dir_figs)) system( paste0( "mkdir -p ", dir_figs ) )
  if (makepdf) pdf( paste0( dir_figs, "/meandoy_bysite_", df$site[1], ".pdf" ) )
    par(las=1)
    yrange <- range( c(df$mod_min, df$mod_max, df$obs_min, df$obs_max), na.rm = TRUE )
    
    # obs
    plot(  df$doy, df$obs_mean, type="l", ylim = yrange, ylab = expression( paste("simulated GPP (gC m"^-2, "d"^-1, ")" ) ), xlab = "DOY" )
    # polygon( c(df$doy, rev(df$doy)), c(df$obs_min, rev(df$obs_max)), border = NA, col = rgb(0,0,0,0.3)  )
    
    # mod: FULL
    lines( df$doy, df$mod_mean, col="red", lwd=1.75 )
    polygon( c(df$doy, rev(df$doy)), c(df$mod_min, rev(df$mod_max)), border = NA, col = rgb(1,0,0,0.3)  )

    # mod 2
    if (!identical(NA,df2)) lines( df2$doy, df2$mod_mean, col=col2, lwd=1.75, lty=1 )

    # mod 3
    if (!identical(NA,df3)) lines( df3$doy, df3$mod_mean, col=col3, lwd=1.75, lty=1 )

    # # mod4
    # if (!identical(NA,df4)) lines( df4$doy, df4$mod_mean, col=col4, lwd=1.75, lty=1 )

    title( df$climatezone[1] )
    
    mtext( bquote( italic(N) == .(ndata) ), side=3, line=1, cex=1.0, adj=1.0 )
    
    # legend("topright", c(lab1, lab2, lab3, lab4), col=c("red", col2, col3, col4), lwd = 1.75, bty = "n" )
    legend("topright", c(lab1, lab2, lab3), col=c("red", col2, col3), lwd = 1.75, bty = "n" )
    
    title( paste( df$site[1], " (", zone, classid, ")" ) )

  if (makepdf) dev.off()
}