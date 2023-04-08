#' find_parenthesis_pairs
#' find pairs of parenthesis
find_parenthesis_pairs <- function(s) {
  stack <- list()
  pairs <- list()
  for (i in 1:nchar(s)) {
    c <- substr(s, i, i)
    if (c == "(") {
      stack <- c(stack, i)
    } else if (c == ")") {
      if (length(stack) == 0) {
        stop("Mismatched parenthesis")
      }
      open <- tail(stack, 1)
      stack <- head(stack, -1)
      pairs <- c(pairs, list(c(open, i)))
    }
  }
  if (length(stack) > 0) {
    stop("Mismatched parenthesis")
  }
  pairs
}
