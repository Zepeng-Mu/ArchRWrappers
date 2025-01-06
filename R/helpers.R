#' Title
#'
#' @param a
#' @param b
#'
#' @return
#' @export
#'
#' @examples
"%&%" <- function(a, b) {
  paste0(a, b)
}

imessage <- function(x) {
  cli::col_br_red(x)
}
