#' Simulate a Piecewise-Constant Signal with Smooth Trend and Noise
#'
#' Generate a univariate signal composed of:
#' - a smooth deterministic component,
#' - a piecewise-constant step component defined by changepoint indices,
#' - additive Gaussian noise.
#'
#' This function is useful for benchmarking changepoint detection algorithms
#' and reproducing controlled simulation studies.
#'
#' @param n Integer. Length of the signal.
#' @param domain Numeric vector of length 2 giving the range of the time axis.
#' @param changepoints_idx Integer vector of changepoint locations (indices).
#'   Must lie strictly inside \eqn{(1, n)}.
#' @param shift_sizes Numeric vector giving the mean level of each segment.
#'   Must have length equal to \code{length(changepoints_idx) + 1}.
#' @param noise_sd Numeric. Standard deviation of the Gaussian noise.
#' @param smooth_fun Function. A function \eqn{f(t)} defining the smooth trend.
#' @param seed Optional integer. If supplied, sets the random seed for reproducibility.
#'
#' @return A list with components:
#' \describe{
#'   \item{t}{Time grid of length \eqn{n}.}
#'   \item{smooth}{Smooth component \eqn{f(t)}.}
#'   \item{step}{Piecewise-constant step component.}
#'   \item{z}{Final noisy signal.}
#'   \item{changepoints_idx}{Sorted changepoint indices.}
#'   \item{shift_sizes}{Segment means.}
#'   \item{params}{List of simulation parameters.}
#' }
#'
#' @examples
#' sim <- simulate_piecewise_signal_idx(
#'   n = 1000,
#'   changepoints_idx = c(300, 700),
#'   shift_sizes = c(0.5, -0.3, 0.7),
#'   noise_sd = 0.05,
#'   seed = 123
#' )
#'
#' @export
simulate_piecewise_signal_idx <- function(
    n = 1000,
    domain = c(-4, 4),
    changepoints_idx = c(300, 700),
    shift_sizes = c(0.5, -0.3, 0.7),
    noise_sd = 0.04,
    smooth_fun = function(x) 0.01 * (3*x/2 - x^3/2),
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  # Create domain
  t <- seq(domain[1], domain[2], length.out = n)

  # Validate changepoints
  if (any(changepoints_idx <= 1 | changepoints_idx >= n)) {
    stop("changepoints_idx must be inside (1, n).")
  }
  if (!is.numeric(changepoints_idx) || any(changepoints_idx %% 1 != 0)) {
    stop("changepoints_idx must be integer indices.")
  }

  # Sort indices to avoid user errors
  changepoints_idx <- sort(changepoints_idx)

  # Number of segments
  n_segments <- length(changepoints_idx) + 1

  if (length(shift_sizes) != n_segments) {
    stop(sprintf(
      "Length of 'shift_sizes' (%d) must equal number of segments (%d).",
      length(shift_sizes), n_segments
    ))
  }

  # Step function container
  g <- numeric(n)

  # Build index cuts
  cuts_idx <- c(0, changepoints_idx, n)

  for (i in seq_along(shift_sizes)) {
    idx <- (cuts_idx[i] + 1):cuts_idx[i + 1]
    g[idx] <- shift_sizes[i]
  }

  # Smooth continuous part
  f <- smooth_fun(t)

  # Final signal
  z <- f + g + rnorm(n, sd = noise_sd)

  list(
    t = t,
    smooth = f,
    step = g,
    z = z,
    changepoints_idx = changepoints_idx,
    shift_sizes = shift_sizes,
    params = list(noise_sd = noise_sd, domain = domain, n.val = n)
  )
}
