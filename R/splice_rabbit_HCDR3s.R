#' Splice out rabbit HCDR3 sequence
#'
#' Splices using amino acid regular expressions that recognize rabbit HCDR3 junctional sequences
#'
#' @param aa AA string set of translated sequences
#' @param upstream_junction protein sequence REGEX upstream of rabbit HCDR3
#' @param downstream_junction protein sequence REGEX downstream of rabbit HCDR3
#' @export
#'
#' @example


splice_rabbit_HCDR3s <- function(aa, upstream_junction = ".TYFC|A.YFC|AT.FC|ATY.C", downstream_junction = "WG.G") {
  HCDR3 <- sapply(strsplit(sapply(strsplit(as.character(aa), downstream_junction), "[", 1),
                           upstream_junction), "[", 2)
  return(HCDR3)}
