#' @title read_as_list
#' @description Read multiple files
#' @param path The direction of the files
#' @param prefix The prefix of the file which would like to import together, Default: ''
#' @param sep The csv/txt separation symbol, Default: 'auto'
#' @param header Whether the colnames or not, Default: 'auto'
#' @param file_type Which type of files would like to import together, Default: 'txt'
#' @return A list containing the imported files
#' @details DETAILS
#' @rdname read_as_list
#' @export 
#' @import dplyr
read_as_list <- function(path, prefix='', sep='auto', header='auto', file_type='txt') {
  out_list <- list()
  file_path <- paste(path ,list.files(path = path, pattern = prefix), sep = '')
  # file_type <- sub(str_extract(file_path[1], pattern = '\\..*'), pattern = '.', replacement = '')
  file_path <- fs::as_fs_path(file_path)
  if(file_type == 'tree')  {
    for(i in file_path) {
      out_list[[i]] <- ape::read.nexus(i)
    }
  }
  if(file_type == 'fasta')  {
    out_list <- file_path %>% purrr::map(.f = function(path){seqinr::read.fasta(path, as.string = F, whole.header = T)})
  }
  if(file_type == 'RData')  {
    out_list <- file_path %>% purrr::map(.f = function(path){readRDS(path)})
  }
  # if(file_type == 'csv')  {
  #   out_list <- file_path %>% purrr::map(.f = function(path){read.csv(path, header = T, sep = ',')})
  # }
  if(file_type %in% c('txt', 'csv'))  {
    for(i in file_path) {
      out_list[[i]] <- data.table::fread(i, header = header, stringsAsFactors = F, sep = sep) %>% 
        as.data.frame()
    }
  }
  names(out_list) <- file_path
  pattern <- c('.*\\/', '\\.[a-z\\.]+', '\\.[a-z]+')
  for(i in 1:2)  {
    names(out_list) <- sub(names(out_list), pattern = pattern[i], replacement = '')
  }
  return(out_list)
}
