#' @title organ_subtype
#' @description Devide the dataframe by subtypes
#' @param in_seq dataframe contain the subtype and segment information
#' @param subtype_col The col contain the segment information
#' @return list contain the internal gene, HA and NA dataframe
#' @details DETAILS
#' @rdname organ_subtype
#' @export 
#' @import dplyr

organ_subtype <- function(in_seq)  {
  in_seq$sub <- in_seq$segment
  # colnames(in_seq)[subtype_col] <- 'segment'
  A <- group_split(in_seq, segment) %>% as.list() %>% lapply(as.data.frame)
  names(A) <- lapply(A, function(x) {x[1,'segment']})
  
  # %>% group_split(segment) %>% as.list() %>% lapply(as.data.frame)
  ## NA subtype
  if(any(names(A) %in% 'NA')) {
    B <- A[['NA']] 
    B$sub <- stringr::str_extract(B$Subtype, pattern = 'N.*')
    B <- B[!(is.na(B$sub)), ]
    na_sub <- B %>% group_split(sub) %>% as.list() %>% lapply(as.data.frame)
    names(na_sub) <- lapply(na_sub, function(x) {x[1,'sub']})
    A <- A[-which(names(A) %in% c('NA'))]
  }else(na_sub='')
  
  
  if(any(names(A) %in% 'HA')) {
    B <- A[['HA']]
    B$sub <- stringr::str_extract(B$Subtype, pattern = 'H[0-9]+')
    B <- B[!(is.na(B$sub)), ] 
    
    ha_sub <- B %>% group_split(sub) %>% as.list() %>% lapply(as.data.frame)
    names(ha_sub) <- lapply(ha_sub, function(x) {x[1,'sub']})
    A <- A[-which(names(A) %in% c('HA'))]
  }else(ha_sub='')
  
  test <- c(A, na_sub, ha_sub)
  
  return(test[!test == ''])
}
