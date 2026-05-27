#' Two-Sided Moving Average Smoothing
#'
#' @description
#' Applies a simple two-sided moving average smoother to a numeric vector.
#' This is used throughout the scanCP pipeline for stabilizing detector
#' statistics and spacing curves.
#'
#' @param x Numeric vector to smooth.
#' @param n Integer. Total number of neighbors used in the moving average.
#'   Must be >= 1.
#' @param circular Logical. If TRUE, the smoothing wraps around the ends.
#'
#' @return A numeric vector of the same length as `x` containing the smoothed
#'   values. Endpoints may be `NA` if `circular = FALSE`.
#'
#' @details
#' This function is a thin wrapper around `stats::filter()` with a symmetric
#' moving average kernel. It is intentionally lightweight and dependency-free.
#'
#' @examples
#' x <- rnorm(100)
#' y <- ma(x, n = 5)
#'
#' @importFrom stats filter
#' @export
ma <- function(x, n = 5, circular = FALSE) {

  # Validate input
  if (!is.numeric(x)) stop("x must be numeric.")
  if (n < 1) stop("n must be >= 1.")

  # Create equal-weight kernel
  kernel <- rep(1 / n, n)

  # Apply two-sided moving average
  stats::filter(x, kernel, sides = 2, circular = circular)
}
