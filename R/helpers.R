#' String concatenation operator
#'
#' A convenient infix operator for concatenating strings without spaces.
#' This is equivalent to \code{paste0(a, b)}.
#'
#' @param a First string or vector to concatenate
#' @param b Second string or vector to concatenate
#'
#' @return A character vector of concatenated strings
#' @export
#'
#' @examples
#' "Hello" %&% " " %&% "World" # Returns "Hello World"
#' "chr" %&% 1:3 # Returns c("chr1", "chr2", "chr3")
"%&%" <- function(a, b) {
  paste0(a, b)
}

#' Print colored message to console
#'
#' This function prints a message in bright red color to the console using
#' the cli package. Useful for highlighting important messages or warnings.
#'
#' @param x Character string or object to print
#'
#' @return Invisibly returns the colored string
#' @export
#'
#' @examples
#' \dontrun{
#' imessage("This is an important message!")
#' }
imessage <- function(x) {
  cli::col_br_red(x)
}
