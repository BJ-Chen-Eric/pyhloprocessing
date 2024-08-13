#' @title str_split2df
#' @description split the strings into dataframe without autofill
#' @param header_vector string vector
#' @param pattern split pattern, Default: '/'
#' @return a splited dataframe
#' @details DETAILS
#' @examples 
#' header_cleaning(c('a_b_c'), pattern='_')
#' @rdname header_cleaning
#' @export 
#' @import dplyr
#' @importFrom stringr str_split
header_cleaning <- function(header_vector, pattern='/') {
  test <- stringr::str_split(header_vector, pattern = pattern) %>% do.call(what=rbind) %>% as.data.frame()
  A <- lapply(stringr::str_split(header_vector, pattern = pattern), length) %>% unlist() %>% as.data.frame() %>% 
    mutate(index=seq(1, nrow(test)))
  for(i in unique(A$.)) {
    if(identical(i, ncol(test))) {next}
    test[A[A$. %in% i, 'index'], (i+1):ncol(test)] <- ''
  }
  return(test)
}
