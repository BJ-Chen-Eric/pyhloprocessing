#' @title run_mafft
#' @description perform alignment through MAFFT
#' @param input_path the direction contain the fasta files for alignment, the file have to be .fa or .fasta
#' @return aligned fasta file in the offered direction
#' @details DETAILS
#' @rdname run_mafft
#' @export 
run_mafft <- function(input_path) {
     # Try to execute the 'mafft --version' command
    result <- tryCatch({
      system("mafft --version > /dev/null 2>&1")
    }, error = function(e) {
      return(FALSE)
    })
    
    # Check if the result indicates MAFFT is available
    if (result == 0) {
      message("MAFFT is installed and ready to use.")
    } else {
      message("MAFFT is not installed on your system.")
      message("Please install MAFFT by following the instructions here: https://mafft.cbrc.jp/alignment/software/")
    }
    
    # List of FASTA files to align
    file_list <- list.files(path = input_path, pattern = "\\.fasta$|\\.fa$", full.names = TRUE)
    
    # Loop through each file, align with MAFFT, and save the output
    for (file in file_list) {
      output_file <- sub("\\.fasta$|\\.fa$", "_aligned.fasta", file)  # Create an output filename
      system(paste("mafft --auto --anysymbol --thread 8", file, ">", output_file))  # Run MAFFT
    }
  }




