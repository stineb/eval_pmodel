# Extract mean R2 across left-out evaluations
extract_mean_rsq <- function(out_oob){
  extract_rsq <- function(mylist){
    out <- mylist$gpp$fluxnet2015$metrics$xdaily_pooled$rsq
    if (!is.na(out) && !is.null(out)) return(out)
  }
  na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
  null.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.null(x)))]) }
  out_oob <- na.omit.list(out_oob)
  list_rsq <- purrr::map(as.list(settings_calib$sitenames),
                         ~extract_rsq(out_oob[[.]]))
  list_rsq <- null.omit.list(list_rsq)
  mean_rsq <- list_rsq %>% unlist() %>% mean()
  return(mean_rsq)  
}