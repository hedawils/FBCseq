#' Extract and clean rabbit HCDR3s
#'
#'For a subset of cell lines
#'
#' @param aa AA string set of translated sequences
#' 
#'
#' @export
#'
#' @example



extract_HCDR3s <- function(aa) {
  HCDR3 <- sapply(strsplit(sapply(strsplit(as.character(aa), "WG.G"), "[", 1),
                           ".TYFC|A.YFC|AT.FC|ATY.C"), "[", 2)
  HCDR3 <- HCDR3[!is.na(HCDR3)] #This code cleans up HCDR3 vector
  HCDR3 <- HCDR3[!(grepl("VT.SS", HCDR3))] #Remove sequences without proper HCDR3 extraction
  HCDR3 <- HCDR3[!(grepl("*", HCDR3, fixed = TRUE))] #Remove sequences without proper HCDR3 extraction
  return(HCDR3)}
