#' Splice out rabbit HCDR3 sequence
#'
#' Splices using amino acid regular expressions that recognize rabbit HCDR3 junctional sequences
#'
#' @param aa AA string set of translated sequences
#' 
#'
#' @export
#'
#' @example


splice_rabbit_HCDR3s <- function(aa) {
  HCDR3 <- sapply(strsplit(sapply(strsplit(as.character(aa), "WG.G"), "[", 1),
                           ".TYFC|A.YFC|AT.FC|ATY.C"), "[", 2)
  return(HCDR3)}
