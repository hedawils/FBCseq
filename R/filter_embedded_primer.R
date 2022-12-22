#' Calculate depmap score mean for depmap cell lines, and N per gene
#'
#'For a subset of cell lines
#'
#' @param seq sequence
#' @param read_length Miseq read length - default = 250
#' @param primer sequence of reverse primer
#' 
#'
#' @export
#'
#' @example



filter_embedded_primer <- function(DNA, primer) {
  library(Biostrings)
  read_length <- nchar(DNA[1] %>% as.character)
  primer_pattern <- paste0(primer, "|", primer %>% DNAString %>% reverseComplement %>% as.character)
  ### filter sequences that have primer sequence (either orientation) at 3' end
  DNA <- DNA[(DNA %>% substr(read_length-nchar(primer), read_length)) %>% str_detect(primer_pattern)] 
  ### trim off primer from 3' end, and trim read lenght to multiple of 3
  DNA <- substr(DNA, 1 + read_length %% 3, read_length - nchar(primer))
  return(DNA)}
