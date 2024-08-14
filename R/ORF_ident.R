#' @title ORF_ident
#' @description identiifed Open reading frame
#' @param fasta_path direction contain the aligned fasta file
#' @param coor_gap_thres coordinate gap proportion up limitation, Default: 86
#' @param out_fasta_path identified ORF write direction, Default: working direction
#' @param out_fasta_prefix identified ORF file prefix, Default: ORF
#' @return write the open reading frame out
#' @details DETAILS
#' @rdname ORF_ident
#' @export 
#' @importFrom data.table fread 
#' @import here
#' @import rlang
#' @import stringr
ORF_ident <- 
function(fasta_path, coor_gap_thres=86, out_fasta_path=paste(getwd(), '/', sep = ''), out_fasta_prefix='ORF') {
  aligned <- read_as_list(path = fasta_path, prefix = 'fa|fasta', file_type = 'fasta')
  outlier <- data.table::fread(here::here('outlier', 'outlier.csv')) %>% as.data.frame()
  
  # j="gisaid_IRD_merged_H5_aligned"
  # j="all_clean_H15_aligned" 
  out_list <- list()
  # ratio_list <- list()
  # non_nt_list <- list()
  stop_codon_list <- list()
  # cons_dist_list <- list()
  # gap_cum_table <- list()
  # gap_cum_plot <- list()
  
  for(j in names(aligned)) {
    cat(paste('Start processing ', j, '\n', sep=''))
    A <- aligned[[j]] %>% do.call(what=rbind) %>% as.data.frame() 
    # seg <- stringr::str_extract(j, pattern = '[A-Z]+[0-9]+|[A-Z]{2}_') %>% 
    #   stringr::str_remove(pattern = '_')
    seg <- stringr::str_extract(j, pattern = 'PB2|PB1|PA|NP|MP|NS|N1|N2|N3|N4|N5|N6|N7|N8|N9|H1|H2|H3|H4|H5|H6|H7|H8|H9|H10|H11|H12|H13|H14|H15|H16')
    
    
    if(grepl(pattern = '(H7|H5)', j))  { #grepl(pattern = '(H7|H5|H9)', j)
      cat(paste('Performing Gap fill procedure', '\n', sep = ''))
      # calculate coordinate gap
      gap <- apply(A[, 1:ncol(A)], MARGIN = 2, FUN = function(x) {grep(x, pattern='-') %>% length()}) %>% as.data.frame()
      gap$index <- seq(1, nrow(gap))
      
      # build consensus, the greatest nucleic acid
      cons <- c()
      for(i in 1:ncol(A))  {
        test <- table(A[, i])
        cons <- c(cons, which.max(test) %>% names())
      }
      cons <- paste(cons, collapse = '')
      
      # identified Start codon
      atg1 <- gregexpr(cons, pattern = 'atg')[[1]][1]
      
      # limit the searching region, 1. start+500~last-200, 2. coordinates have more than 95% gap
      fill <- gap[((atg1)+500):(nrow(gap)-200), ]
      fill <- fill[fill$. >= nrow(A)*0.95, ]
      
      # identified the high gap proportion continuous resgion
      fill$in2 <- c(fill$index[-1], 0)
      fill$in2 <- fill$in2-fill$index
      fill[fill$in2 <= 10, 'in2'] <- 1
      rownames(fill) <- seq(1, nrow(fill),1)
      v <- (fill[!(fill$in2 %in% 1), ] %>% rownames())
      
      
      if(length(v) == 1)  {
        v <- ifelse(identical(v, character(0)), yes = nrow(fill), no = v)
      }
      
      # identified the continuous regions
      row <- c(1, v)
      sub_fill_list <- list()
      x=0
      for(i in 2:(length(row)))  {
        B <- fill[row[i-1]:row[i], ]
        if(nrow(B) > 100) {
          x=x+1
          sub_fill_list[[x]] <- B
        }else(next)
      }
      
      for(i in seq_along(sub_fill_list))  {
        sub_fill <- sub_fill_list[[i]]
        cat(paste('Processing region: ', sub_fill[2,2], ' to ', sub_fill[nrow(sub_fill),2], '\n', sep = ''))
        test <- A[, (sub_fill[2,2]-5):(sub_fill[nrow(sub_fill),2]+5)]
        l <- ncol(test)
        
        test$fill <- apply(test[1:nrow(test), ], MARGIN = 1, FUN = function(x) {paste(x, collapse = '') %>% str_remove_all(pattern = '-')}) #
        test$sub_l <- l-nchar(test$fill)
        for(k in 1:nrow(test))  {
          test[k, 'fill'] <- paste(test[k, 'fill'], rep('-', test[k, 'sub_l']) %>% paste(collapse = ''), sep = '')
        }
        test <- stringr::str_split(test$fill, pattern = '') %>% do.call(what=rbind) %>% as.data.frame()
        A[, (sub_fill[2,2]-5):(sub_fill[nrow(sub_fill),2]+5)] <- test
      }
    }
    
    
    # test <- apply(A, MARGIN = 1, FUN = function(x) {paste(x, collapse = '')})
    # out <- matrix(0,0,0) %>% as.data.frame()
    # out[seq(1, length(test)*2, 2), 1] <- paste('>', names(test), sep = '')
    # out[seq(2, length(test)*2, 2), 1] <- test
    # write.table(out, file = paste(args$out_fasta_dir, j, '_fill.fa', sep = ''), sep = '\n'
    #             , quote = F, col.names = F, row.names = F)
    
    
    gap <- apply(A[, 1:ncol(A)], MARGIN = 2, FUN = function(x) {grep(x, pattern='-') %>% length()}) %>% 
      as.data.frame() %>% mutate(index=row_number(), prop=./nrow(A)*100)
    
    gap <- gap[gap$prop < coor_gap_thres, ]
    A <- A[, gap$index]
    
    
    cons <- c()
    for(i in 1:ncol(A))  {
      test <- table(A[, i])
      cons <- c(cons, which.max(test) %>% names())
    }
    atg1 <- (paste(cons, collapse = '') %>% gregexpr(pattern = 'atg') %>% unlist())[1]
    
    test <- apply(A[1:nrow(A), atg1:ncol(A)], MARGIN = 1, FUN = function(x) {paste(x, collapse = '')})
    
    # cut into codon
    s <- seq(1, nchar(test[1]), 3)
    e <- seq(3, nchar(test[1]), 3)
    
    e <- e[1:min(length(s), length(e))]
    s <- s[1:min(length(s), length(e))]
    
    B <- sapply(test, function(x) {substring(x, s, e)}) %>% t()
    
    
    if(!seg %in% 'NS') {
      B <- apply(B, MARGIN = 2, function(x) {grep(x, pattern = '(taa|tag|tga)') %>% length()}) %>% as.data.frame()
      stop1 <- which.max(B$.)
      test <- A[, atg1:(stop1*3+atg1-1-3)]
      B <- sapply(apply(test[1:nrow(test), ], MARGIN = 1, FUN = function(x) {paste(x, collapse = '')})
                  , function(x) {substring(x, seq(1, nchar(x), 3), seq(3, nchar(x), 3))}) %>% t()
      B <- apply(B, MARGIN = 2, function(x) {grep(x, pattern = '(taa|tag|tga)') %>% length()}) %>% as.data.frame()
      stop_codon_list[[j]] <- B %>% sum()
    }
    
    if(seg %in% 'NS') {
      for(i in seq_len(nrow(B))) {
        stop <- B[i, ] %in% c('taa', 'tag', 'tga') %>% which()
        stop <- stop[stop >190] %>% min()
        if(!is.infinite(stop)) {
          B[i, (stop):(ncol(B))] <- '---'
        }
      }
      stop_codon_list[[j]] <- apply(B, MARGIN = 2, function(x) {grep(x, pattern = '(taa|tag|tga)') %>% length()}) %>% as.data.frame() %>% sum()
      
      test <- apply(B %>% as.data.frame(), MARGIN = 1, FUN = function(x) {paste(x, collapse = '')})
      test <- sapply(test, function(x) {substring(x, seq(1, nchar(x), 1), seq(1, nchar(x), 1))}) %>% t() %>% as.data.frame()
    }
    
    
    test$result <- apply(test[1:nrow(test), ], MARGIN = 1, FUN = function(x) {paste(x, collapse = '')})
    if(grepl(rownames(test)[1], pattern='EPI')) {
      test$epi <- str_extract(rownames(test), pattern = 'EPI_ISL_[0-9]+|IRD_[0-9]+')
    }
    
    
    # out_list[[j]] <- test
    cat(paste('Finish ', j, ' processing', '\n', sep=''))
    
    write_fasta(sequence = test$result, header = rownames(test), 
                file_out = paste(out_fasta_path, out_fasta_prefix, '_', seg, '.fa', sep = ''))
    
    
  }
  
}




