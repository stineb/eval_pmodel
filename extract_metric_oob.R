# Extract mean R2 across left-out evaluations
extract_metric_oob <- function(out_oob, metric="rsq"){
  sitenames <- names(out_oob)[-which(names(out_oob)=="AALL")]
  extract_metric <- function(mylist, metric){
    out <- mylist$gpp$fluxnet2015$metrics$xdaily_pooled[[metric]]
    if (!is.na(out) && !is.null(out)) return(out)
  }
  na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
  null.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.null(x)))]) }
  out_oob <- na.omit.list(out_oob)
  list_metric <- purrr::map(as.list(sitenames),
                         ~extract_metric(out_oob[[.]], metric))
  names(list_metric) <- sitenames
  list_metric <- null.omit.list(list_metric)
  df_metric <- tibble(sitename = names(list_metric),  metric = unlist(list_metric)) %>% 
    setNames(c("sitename", metric))
  return(df_metric)  
}