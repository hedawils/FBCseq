#' Extract and clean rabbit HCDR3s
#'
#' Splices out HCDR3s, then filters HCDR3s that are in frame and do not havce stop codons.
#'
#' @param aa AA string set of translated sequences
#' @param tail_pattern fixed sequence in FW4 of rabbit VH domains. sequences without this are removed
#' @param upstream_junction protein sequence REGEX upstream of rabbit HCDR3
#' @param downstream_junction protein sequence REGEX downstream of rabbit HCDR3
#'
#'
#' @export
#'
#' @example



extract_HCDR3s <- function(aa, tail_pattern = "VT.SS", upstream_junction = ".TYFC|A.YFC|AT.FC|ATY.C", downstream_junction = "WG.G") {
  HCDR3 <- splice_rabbit_HCDR3s(aa, upstream_junction, downstream_junction)
  HCDR3 <- HCDR3[!is.na(HCDR3)] #This code cleans up HCDR3 vector
  HCDR3 <- HCDR3[!(grepl("VT.SS", HCDR3))] #Remove sequences without proper HCDR3 extraction
  HCDR3 <- HCDR3[!(grepl("*", HCDR3, fixed = TRUE))] #Remove sequences without proper HCDR3 extraction
  return(HCDR3)}
