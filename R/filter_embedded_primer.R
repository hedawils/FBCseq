#' Calculate depmap score mean for depmap cell lines, and N per gene
#'
#'For a subset of cell lines
#'
#' @param seq sequence
#' @param read_length Miseq read length - default = 250
#'
#'
#' @export
#'
#' @example



filter_embedded_primer <- function(seq, read_length) seq[grepl(primer, substr(seq, read_length -
                                                              nchar(primer), read_length)) ]
