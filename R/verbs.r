# verbs.r -- custom verbs for gnibbles

# group_by <- function (.data, ..., add = FALSE) {
#   UseMethod('group_by', .data)
# }

# #' @export
# group_by.gnibble <- function(.data, ..., add = FALSE) {
#   class(.data) <- setdiff(class(.data), 'gnibble')
#   out <- group_by_(.data, .dots = lazyeval::lazy_dots(...), add = add)
#   class(out) <- union(c('grouped_df', 'gnibble'), class(out))
#   out
# }

#' Mutate function for gnibble
#' @param .data A gnibble dataframe.
#' @param ... Name-value pairs of expressions. Use \code{NULL} to drop a variable.
#' @param .dots: Used to work around non-standard evaluation. 
#' @importFrom dplyr mutate_
#' @export
mutate_.gnibble <- function(.data, ..., .dots) {
  dots <- lazyeval::all_dots(.dots, ..., all_named = TRUE)
  class(.data) <- setdiff(class(.data), 'gnibble')
  out <- dplyr::mutate_(dplyr::tbl_df(.data), .dots = dots)
  chrom_lengths(out) <- chrom_lengths(.data)
  if (has_windows(.data))
    attr(out, 'windows') <- windows(.data)
  class(out) <- union('gnibble', class(out))
  out
}

#' Filter function for gnibble
#' @param .data A gnibble dataframe.
#' @param ... Name-value pairs of expressions. Use \code{NULL} to drop a variable.
#' @param .dots: Used to work around non-standard evaluation. 
#' @importFrom dplyr filter_
#' @export
filter_.gnibble <- function(.data, ..., .dots) {
  dots <- lazyeval::all_dots(.dots, ...)
  class(.data) <- setdiff(class(.data), 'gnibble')
  out <- filter_(tbl_df(.data), .dots = dots)
  if (has_windows(.data))
    attr(out, 'windows') <- windows(.data)
  class(out) <- union('gnibble', class(out))
  out
}
