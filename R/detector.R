#' Compute the MLP-Based Change Point Detector Statistic
#'
#' @description
#' Computes a detector statistic using rolling-window MLP fits. The statistic
#' combines a *ratio* component and a *difference* component derived from the
#' residual sum of squares (RSS) of:
#'
#' - two adjacent small-window MLP fits (window size `w`)
#' - one large-window MLP fit (window size `2w`)
#'
#' The detector is used to highlight potential changepoints in the signal.
#'
#' @param y Numeric vector. Original signal.
#' @param fit_mlp_res List. Output of [fit_mlp()], containing:
#'   \itemize{
#'     \item `small` — list of small-window fits
#'     \item `large` — list of large-window fits
#'   }
#' @param w Integer. Window size used in the small MLP.
#' @param a Numeric in \[0,1\]



#' @param a Weight between ratio (0) and difference (1).
#' @param b Numeric. Small correction added to denominator of ratio statistic.
#' @param scale_01 Logical. If TRUE, rescales detector to \[0,1\]



#'
#' @return A numeric vector of detector values of length `n - 2w`.
#'
#' @details
#' For each index `i`, the detector uses:
#'
#' - `rss1`: RSS of small-window MLP on segment `[i, i+w-1]`
#' - `rss2`: RSS of small-window MLP on segment `[i+w, i+2w-1]`
#' - `rss.tot`: RSS of large-window MLP on `[i, i+2w-1]`
#'
#' The ratio component is:
#' \deqn{ (rss.tot + b) / (rss1 + rss2 + b) }
#'
#' The difference component is:
#' \deqn{ rss.tot - rss1 - rss2 }
#'
#' @examples
#' # sim <- simulate_piecewise_signal_idx()
#' # fit <- fit_mlp(sim$z, w = 100)
#' # d <- calc_detector(sim$z, fit, w = 100)
#'
#' @export
calc_detector <- function(
    y,
    fit_mlp_res,
    w = 100,
    a = 1,
    b = 0,
    scale_01 = FALSE
) {

  # ------------------------------------------------------------
  # Validate inputs
  # ------------------------------------------------------------
  if (!is.numeric(y)) stop("y must be numeric.")
  if (!is.list(fit_mlp_res)) stop("fit_mlp_res must be a list.")
  if (!all(c("small", "large") %in% names(fit_mlp_res)))
    stop("fit_mlp_res must contain elements 'small' and 'large'.")

  res.small <- fit_mlp_res$small
  res.large <- fit_mlp_res$large

  n.val <- length(y)
  n.det <- n.val - 2 * w
  if (n.det <= 0) stop("Signal too short for window size w.")

  # ------------------------------------------------------------
  # Preallocate detector components
  # ------------------------------------------------------------
  diff1 <- numeric(n.det)   # ratio component
  diff2 <- numeric(n.det)   # difference component

  # ------------------------------------------------------------
  # Compute detector components
  # ------------------------------------------------------------
  for (i in seq_len(n.det)) {

    # RSS for first half-window
    rss1 <- sum((y[i:(i + w - 1)] - res.small[[i]][, 2])^2)

    # RSS for second half-window
    rss2 <- sum((y[(i + w):(i + 2 * w - 1)] - res.small[[i + w]][, 2])^2)

    # RSS for full window
    rss.tot <- sum((y[i:(i + 2 * w - 1)] - res.large[[i]][, 2])^2)

    # Ratio component
    diff1[i] <- (rss.tot + b) / (rss1 + rss2 + b)

    # Difference component
    diff2[i] <- rss.tot - rss1 - rss2
  }

  # ------------------------------------------------------------
  # Optional scaling to [0,1]
  # ------------------------------------------------------------
  if (scale_01) {

    # Normalize each component
    diff1 <- (diff1 - min(diff1)) / (max(diff1) - min(diff1))
    diff2 <- (diff2 - min(diff2)) / (max(diff2) - min(diff2))

    # Combine
    det <- (1 - a) * diff1 + a * diff2

    # Final normalization
    det <- (det - min(det)) / (max(det) - min(det))

  } else {
    det <- (1 - a) * diff1 + a * diff2
  }

  det
}
