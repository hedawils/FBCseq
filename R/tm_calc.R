#' Calculate tm of DNA oligo sequence
#'
#'
#' @param sequence primer sequence, as string
#'
#'
#' @export
#'
#' @example
#' tm_calc("GCCCTTGGTGGAGGC")
#' #### 46.15714



tm_calc <- function(sequence)
{sequence <- substr(sequence,1,nchar(sequence)-1) 
DNA <- Biostrings::alphabetFrequency(DNAString(sequence))
tm <- 64.9 + 41 * (DNA[2] + DNA[3] - 16.4)/(DNA[1] + DNA[2] + DNA[3] + DNA[4])
return(tm)
}
