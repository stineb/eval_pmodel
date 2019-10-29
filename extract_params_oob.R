# Extract mean R2 across left-out evaluations
extract_params_oob <- function(param, path, sitenames){
  
  extract_params_oob_bysite <- function(path, sitename){
    df_params <- read_csv(paste0(path, "/params_opt_leftout_", sitename, ".csv")) %>% 
      dplyr::mutate(sitename = sitename)
    return(df_params)
  }
  df_params <- purrr::map_dfr(
    as.list(sitenames),
    ~extract_params_oob_bysite(path, .)
  )
  return(df_params)  
}