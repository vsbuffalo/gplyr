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

# #' @export
# group_by.data.frame <- function(.data, ..., add = FALSE) {
#   dplyr::group_by(.data, ..., add=add)
# }

# @export
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
