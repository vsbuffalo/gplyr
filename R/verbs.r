# verbs.r -- custom verbs for gnibbles

group_by <- function (.data, ..., add = FALSE) {
  UseMethod('group_by', .data)
}

#' @export
group_by.gnibble <- function(.data, ..., .dots, add = FALSE) {
  class(.data) <- setdiff(class(.data), 'gnibble')
  out <- group_by_(.data, .dots = lazyeval::lazy_dots(...), add = add)
  class(out) <- union(c('grouped_df', 'gnibble'), class(out))
  out
}

# not need (yet)
#' @export
# summarise_.gnibble <- function(.data, ..., .dots) {
#   dots <- lazyeval::all_dots(.dots, ..., all_named = TRUE)
#   class(.data) <- setdiff(class(.data), 'gnibble')
#   out <- dplyr::summarise_(tbl_df(.data), .dots = dots)
#   class(out) <- union('gnibble', class(out))
#   out
# }

#' @export
arrange_.gnibble <- function(.data, ..., .dots) {
  dots <- lazyeval::all_dots(.dots, ..., all_named = TRUE)
  class(.data) <- setdiff(class(.data), 'gnibble')
  out <- dplyr::arrange_(tbl_df(.data), .dots = dots)
  class(out) <- union('gnibble', class(out))
  out
}

# @export
mutate_.gnibble <- function(.data, ..., .dots) {
  dots <- lazyeval::all_dots(.dots, ..., all_named = TRUE)
  class(.data) <- setdiff(class(.data), 'gnibble')
  out <- dplyr::mutate_(dplyr::tbl_df(.data), .dots = dots)
  chrom_lengths(out) <- chrom_lengths(.data)
  if (has_windows(.data))
    attr(out, 'windows') <- get_windows(.data)
  class(out) <- union('gnibble', class(out))
  out
}
