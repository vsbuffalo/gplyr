## window.r -- helper functions

sim_ranges <- function(n, chrom_lengths, max_len=10e3) {
  chroms_i <- sample(seq_along(chrom_lengths$chrom), n, replace=TRUE)  
  chrom <- chrom_lengths$chrom[chroms_i]
  start <- sapply(chrom_lengths$length[chroms_i],
                function(end) floor(runif(1, 1, end-max_len-1)))
  end <- start + floor(runif(n, 1, max_len))
  out <- tibble::tibble(chrom, start=start, end=end)
  attr(out, 'chrom_lengths') <- chrom_lengths
  out
}

bin_start <- function(x) as.numeric(gsub("\\(([^,]+),([^]]+)]", "\\1", x))
bin_end <- function(x) as.numeric(gsub("\\(([^,]+),([^]]+)]", "\\2", x))

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
  windows <- tibble::tibble(chrom=
                            factor(as.vector(GenomicRanges::seqnames(tiles)),
                                 levels=chrom_lengths$chrom),
                   start=GenomicRanges::start(tiles),
                   end=GenomicRanges::end(tiles),
                   width=end-start,
                   key=factor(sprintf("%s:%d-%s", chrom, start, end)))
  gr <- with(.data, GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end)))
  # find overlaps (same order as original data tibble)
  olaps <- GenomicRanges::findOverlaps(gr, tiles, select="first")
  .data$window <- factor(windows$key[olaps], levels=windows$key)
  attr(.data, 'chrom_lengths') <- chrom_lengths
  attr(.data, 'windows') <- windows
  .data
}

get_windows <- function(.data) {
  atts <- attributes(.data)
  if (!('windows' %in% names(atts)))
    stop("no 'windows' attribute found")
  atts$windows
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
  window_cols <- c('chrom', 'start', 'end')
  windows <- get_windows(.data)
  i <- match(.data$window, windows$key)
  out <- .data
  out$chrom <- factor(windows$chrom[i], levels=levels(windows$chrom))
  out$wstart <- windows$start[i]
  out$wend <- windows$end[i]
  added_cols <- c('chrom', 'start', 'end', 'wstart', 'wend')
  # reorder columns so chrom, start end first.
  new_cols <- c(added_cols, setdiff(colnames(out), added_cols))
  if (remove)
    new_cols <- setdiff(new_cols, 'window')
  out <- out[, new_cols]
  chrom_lengths(out) <- chrom_lengths(.data)
  return(out)
}

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
  windows <- get_windows(.data)
  keys <- make_key(.data$chrom, .data$start, .data$end)
  if (!all(unique(keys) %in% windows$key)) {
    browser()
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
  return(out)
}

append_wcumpos <- function(.data) {
  if (!('wcenter' %in% colnames(.data)))
    .data <- .data %>% append_wcenter()
  arrange(.data, 'chrom', 'win_start') %>% mutate(wcumpos=cumsum(wcenter))
}


