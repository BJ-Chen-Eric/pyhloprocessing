#' @title write_fasta
#' @description Generate Fasta file local file
#' @param sequence sequences vector
#' @param header header vector
#' @param file_out generated fasta file direction
#' @rdname write_fasta
#' @export 
#' @import dplyr
#' @importFrom utils write.csv write.table
write_fasta <- 
function(sequence, header, file_out) {
  out <- matrix(0,0,0) %>% as.data.frame()
  if(identical(stringr::str_extract(header[1], pattern = '>'), '>')) {
    header <- header
  }else{header <- paste('>', header, sep = '')}
  
  out[seq(1, length(sequence)*2, 2), 1] <- header
  out[seq(2, length(sequence)*2, 2), 1] <- sequence
  
  write.table(out, file = file_out, quote = F, col.names = F, row.names = F)
}


