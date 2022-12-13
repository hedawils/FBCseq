#' Calculate depmap score mean for depmap cell lines, and N per gene
#'
#'For a subset of cell lines
#'
#' @param depmap depmap matrix (see example)
#' @param depmap_ids vector of depmap_ids for which score should be calculated
#' @param col_names name of returned pvalue column
#'
#' @return dataframe of mean per gene
#'
#' @export
#'
#' @example
#' depmap <- AuronR::read_from_Auron_db("depmap_crispr")
#' depmap <- AuronR::from_json(depmap, json_col = "data", identifier_col = "depmap_id")
#' depmap_ids <- AuronR::read_from_Auron_db("depmap_metadata", columns = "depmap_id") %>% pull(depmap_id)
#' depmap_mean_and_n(depmap, depmap_ids)


process_fastq <- function(control_fastq_files, experimental_fastq_files, embedded_primer = "GCCCTTGGTGGAGGC"){
  
  library_names <- c(control_fastq_files, experimental_fastq_files)
  
  list.files(pattern = ".fastq")
  libs <- list()
  for (f in 1:length(library_names)){
    dada2::fastqFilter(library_names[f], "trim.fastq.gz", truncLen = (120+nchar(embedded_primer)), trimLeft = nchar(embedded_primer), primer.fwd = embedded_primer, n =1e+05,compress = TRUE)
    lib <- dada2::derepFastq("trim.fastq.gz", n = 1e+05)
    libs[[f]] <- data.frame("seq" = names(lib$uniques), "count" = lib$uniques)
    unlink("trim.fastq.gz")
  }
  
  rbind.data.frame()
}

#
# library(tidyverse)
# list(x, y, z) %>% reduce(left_join, by = "i")
