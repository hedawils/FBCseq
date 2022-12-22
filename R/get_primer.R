#' Get reverse primer sequences for subset of HCDR3s
#'
#' Trims reverse complemented HCDR3 DNA sequence to appropriate primer length, based upon annealing temperature calculation
#'
#' @param sequence sequence to trim from 3' end to obtain primer. reverse complement of HCDR3 template
#' @param annealing annealing temperature of primers, default = 62
#'
#'
#' @export
#'
#' @example
#' get_primer("CCAGATGTCGGCATTATTAGGCCATCTTCTGGCACAGAAATAGGTGGCCGTGTCCGAGGCTGTCAGACTGTTCATTTTCAGATCCATCGTGGTCGAGGTTTTGGA", annealing = 62)
#' ### "CCAGATGTCGGCATTATTAGGCCATCTTCT"

get_primer <- function(sequence, annealing = 62){
  success <- FALSE
  while (!success) {
    sequence <- substr(sequence,1,nchar(sequence)-1)
    DNA <- alphabetFrequency(DNAString(sequence))
    tm <- 64.9 + 41 * (DNA[2] + DNA[3] - 16.4)/(DNA[1] + DNA[2] + DNA[3] + DNA[4])
    success <- tm < annealing
  }
  return(sequence)
}
