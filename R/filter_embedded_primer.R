#' Filter embedded CH1 primer from reads
#'
#' Identify reads containing correct reverse primer, filter them, then trim reads
#' such that HCDR3 is in correct reading frame
#'
#' @param seq sequence
#' @param read_length read length of single end read
#' @param primer sequence of reverse primer
#' 
#'
#' @export
#'




filter_embedded_primer <- function(DNA, primer, read_length) {
  
  #primer_pattern <- paste0(primer, "|", primer %>% DNAString %>% reverseComplement %>% as.character)
  ### filter sequences that have primer sequence (either orientation) at 3' end
  DNA <- DNA[(DNA %>% substr(read_length-nchar(primer), read_length)) %>% 
               str_detect(primer %>% DNAString %>% reverseComplement %>% as.character)]  %>% as.character
  
  #DNA <- substr(DNA, 1 + read_length %% 3, read_length - nchar(primer))
  return(DNA)}
