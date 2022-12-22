#' Extract and clean rabbit HCDR3s
#'
#' Splices out HCDR3s, then filters HCDR3s that are in frame and do not havce stop codons.
#'
#' @param aa AA string set of translated sequences
#' 
#'
#' @export
#'
#' @example



extract_HCDR3s <- function(aa) {
  HCDR3 <- splice_HCDR3s(aa)
  HCDR3 <- HCDR3[!is.na(HCDR3)] #This code cleans up HCDR3 vector
  HCDR3 <- HCDR3[!(grepl("VT.SS", HCDR3))] #Remove sequences without proper HCDR3 extraction
  HCDR3 <- HCDR3[!(grepl("*", HCDR3, fixed = TRUE))] #Remove sequences without proper HCDR3 extraction
  return(HCDR3)}
