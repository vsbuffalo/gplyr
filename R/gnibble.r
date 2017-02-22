# gnibble.r -- simple wrappers for making gnibbles from things

#' Create a gnibble -- a genome tibble
#'
#' @param ... Arguments passed to \code{tibble()}.
#' @param  chrom_lengths Dataframe of chromosome lengths, with columns \code{chrom}, \code{lengths}.
#'
#' @export
gnibble <- function(..., chrom_lengths) {
  out <- tibble::tibble(...)
  class(out) <- union('gnibble', class(out))
  if (!missing(chrom_lengths))
    chrom_lengths(out) <- chrom_lengths
  out
}




