#' @title table_DF
#' @description Build the frequency table into dataframe
#' @param x data vector
#' @param percent whether to percentage, Default: F
#' @return frequency dataframe
#' @examples 
#' table_DF(c(1, 1, 1, 2, 2))
#' @rdname table_DF
#' @import dplyr
#' @import rlang
#' @export 

table_DF <- function(x, percent=F) {
  A <- table(x) %>% as.data.frame(stringsAsFactors=F)
  if(isTRUE(percent))  {
    A <- table(x) %>% prop.table() %>% as.data.frame(stringsAsFactors=F) %>% mutate(Freq=100*.data$Freq)
  }
  return(A)
}
