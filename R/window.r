# window.r -- helper functions

#' @export
sim_ranges <- function(n, chrom_lengths, max_len=10e3) {
  chroms_i <- sample(seq_along(chrom_lengths$chrom), n, replace=TRUE)  
  chrom <- chrom_lengths$chrom[chroms_i]
  start <- sapply(chrom_lengths$length[chroms_i],
                function(end) floor(runif(1, 1, end-max_len-1)))
  end <- start + floor(runif(n, 1, max_len))
  gnibble(chrom, start=start, end=end, chrom_lengths=chrom_lengths)
}

#' Add Chromosome Lengths to Dataframe attributes
#'
#' @param .data Dataframe containing genomic range data.
#' @param chrom_lengths Dataframe of chromosome lengths (with columns
#' \code{chrom} and \code{length}.
#'
#' @export
`chrom_lengths<-` <- function(.data, value) {
  attr(.data, 'chrom_lengths') <- value
  .data
}

has_chrom_lengths <- function(x) 'chrom_lengths' %in% names(attributes(x))


#' Get Chromosome Lengths from Dataframe attributes
#'
#' @param .data Dataframe containing genomic range data.
#'
#' @export
chrom_lengths <- function(.data) {
  if (!has_chrom_lengths(.data)) stop("attribute 'chrom_lengths' is not set")
  attr(.data, 'chrom_lengths') 
}


#' Bin a Dataframe of Genomic Ranges into Windows
#'
#' Simple strand-ignoring binning procedure for genomic ranges stored in a
#' dataframe. The dataframe must have columns \code{chrom}, \code{start},
#' \code{end}. 
#'
#' @param .data Dataframe with columns \code{chrom}, \code{start}, \code{end}
#' @param chrom_lengths Dataframe of chromosome names and lengths in columns
#' \code{chrom} and \code{length}.
#' @param width Width of window in base pairs.
#' @param chroms Optional vector of chromosome order. 
#'
#' @export
create_windows <- function(.data, width) {
  .data <- dplyr::arrange_(.data, 'chrom', 'start')
  chrom_lengths <- chrom_lengths(.data)
  tiles <- unlist(GenomicRanges::tile(with(chrom_lengths, 
                   GenomicRanges::GRanges(chrom, IRanges::IRanges(0, length))),
                                      width=width))
  windows <- gnibble(chrom=factor(as.vector(GenomicRanges::seqnames(tiles)),
                                 levels=chrom_lengths$chrom),
                   start=GenomicRanges::start(tiles),
                   end=GenomicRanges::end(tiles),
                   width=end-start,
                   key=factor(sprintf("%s:%d-%s", chrom, start, end)))
  gr <- with(.data, GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end)))
  # find overlaps (same order as original data tibble)
  olaps <- GenomicRanges::findOverlaps(gr, tiles, select="first")
  .data$window <- factor(windows$key[olaps], levels=windows$key)
  chrom_lengths(.data) <- chrom_lengths
  class(.data) <- union('gnibble', class(.data))
  attr(.data, 'windows') <- windows
  .data
}

has_windows <- function(.data) {
  atts <- attributes(.data)
  return('windows' %in% names(atts))
}

#' @export
windows <- function(.data) {
  if (!has_windows(.data))
    stop("no 'windows' attribute found")
  atts <- attributes(.data)
  atts$windows
}

#' @export
`windows<-` <- function(.data, value) {
  attr(.data, 'windows') <- value
  .data
}

#' Separate the Window Column into its Metadata columns \code{chrom},
#' \code{start}, \code{end}
#' 
#' @param .data Dataframe with a column \code{window} (and \codes{windows}
#'        attribute). 
#' @param remove If \code{TRUE}, remove columns \code{chrom}, \code{start}, 
#'        and \code{end}.
#
#' @export
separate_window <- function(.data, remove=TRUE) {
  windows <- windows(.data)
  i <- match(.data$window, windows$key)
  out <- .data
  out$chrom <- factor(windows$chrom[i], levels=levels(windows$chrom))
  out$wstart <- windows$start[i]
  out$wend <- windows$end[i]
  added_cols <- c('chrom', 'wstart', 'wend')
  # reorder columns so chrom, start end first.
  new_cols <- c(added_cols, setdiff(colnames(out), added_cols))
  if (remove)
    new_cols <- setdiff(new_cols, 'window')
  out <- out[, new_cols]
  class(out) <- union('gnibble', class(.data))
  chrom_lengths(out) <- chrom_lengths(.data)
  windows(out) <- windows(.data)
  return(out)
}

#' @export
append_wcenter <- function(.data, remove=TRUE) {
  dplyr::mutate(.data, wcenter=(wstart+wend)/2)
}

make_key <- function(chrom, start, end) sprintf("%s:%d-%s", chrom, start, end)

#' Unite Window Columns into Single Window Key Column
#' 
#' @param .data Dataframe with columns \code{chrom}, \code{start}, and 
#'        \code{end}, (and \codes{windows} attribute). 
#' @param remove If \code{TRUE}, remove columns \code{chrom}, \code{start}, 
#'        and \code{end}.
#'
#' @export
unite_window <- function(.data, remove=TRUE) {
  windows <- windows(.data)
  keys <- make_key(.data$chrom, .data$wstart, .data$wend)
  if (!all(unique(keys) %in% windows$key)) {
    stop("some window keys are not in windows attributes dataframe!")
  }
  out <- .data
  # use existing factors in windows key
  out$window <- windows$key[match(keys, windows$key)]
  added_cols <- remove_cols <- c('window')
  if (remove)
    remove_cols <- c('chrom', 'start', 'end') 
  out <- out[, c('window', setdiff(colnames(out), remove_cols))]
  chrom_lengths(out) <- chrom_lengths(.data)
  windows(out) <- windows(.data)
  return(out)
}

#' @export
append_wcumpos <- function(.data) {
   mutate(append_wcenter(.data), wcumpos=cumsum(wcenter))
}


