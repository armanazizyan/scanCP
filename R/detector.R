#' Compute the MLP-Based Changepoint Detector Statistic
#'
#' @description
#' Computes the hybrid detector statistic used in scanCP. The detector
#' combines a *ratio* component and a *difference* component derived from
#' the residual sum of squares (RSS) of:
#'
#' - two adjacent small-window MLP fits (window size `w`)
#' - one large-window MLP fit (window size `2w`)
#'
#' The detector highlights locations where the large-window fit performs
#' substantially worse than the two small-window fits, indicating a
#' potential changepoint.
#'
#' @param y Numeric vector. Original signal.
#' @param fit_mlp_res List. Output of \code{\link{fit_mlp}}, containing:
#'   \itemize{
#'     \item \code{small} — list of small-window MLP predictions
#'     \item \code{large} — list of large-window MLP predictions
#'   }
#' @param w Integer. Window size used for the small-window MLP.
#' @param ma_window Integer. Moving-average smoothing window applied to the
#'   detector. Defaults to \code{w}.
#' @param use_abs_det Logical. Whether to take \code{abs(detector)} before
#'   smoothing. Defaults to \code{TRUE}.
#'
#' @return A numeric vector of detector values of length \code{n - 2w}.
#'
#' @details
#' For each index \code{i}, the detector uses:
#'
#' - \code{rss1}: RSS of the small-window MLP on \code{[i, i+w-1]}
#' - \code{rss2}: RSS of the small-window MLP on \code{[i+w, i+2w-1]}
#' - \code{rss.tot}: RSS of the large-window MLP on \code{[i, i+2w-1]}
#'
#' The ratio component is:
#' \deqn{ (rss.tot + b) / (rss1 + rss2 + b) }
#'
#' The difference component is:
#' \deqn{ rss.tot - rss1 - rss2 }
#'
#' The final detector is a convex combination:
#' \deqn{ (1 - a) \cdot \text{ratio} + a \cdot \text{difference} }
#'
#' Internal parameters \code{a}, \code{b}, and \code{scale_01} are kept
#' fixed at their defaults for simplicity and stability.
#'
#' @examples
#' \donttest{
#'   y <- rnorm(300)
#'   fit <- fit_mlp(y, w = 100)
#'   det <- calc_detector(y, fit, w = 100)
#' }
#'
#' @usage calc_detector(y, fit_mlp_res, w = 100,
#'               ma_window = w, use_abs_det = TRUE)
#'
#' @export
calc_detector <- function(
    y,
    fit_mlp_res,
    w = 100,
    ma_window = w,
    use_abs_det = TRUE
) {

  scale_01 = FALSE
  circular = FALSE
  a = 1
  b = 0

  # Validate inputs
  if (!is.numeric(y)) stop("y must be numeric.")
  if (!is.list(fit_mlp_res)) stop("fit_mlp_res must be a list.")
  if (!all(c("small", "large") %in% names(fit_mlp_res)))
    stop("fit_mlp_res must contain elements 'small' and 'large'.")

  res.small <- fit_mlp_res$small
  res.large <- fit_mlp_res$large

  n.val <- length(y)
  n.det <- n.val - 2 * w
  if (n.det <= 0) stop("Signal too short for window size w.")

  diff1 <- numeric(n.det)   # ratio component
  diff2 <- numeric(n.det)   # difference component

  for (i in seq_len(n.det)) {

    rss1 <- sum((y[i:(i + w - 1)] - res.small[[i]][, 2])^2)
    rss2 <- sum((y[(i + w):(i + 2 * w - 1)] - res.small[[i + w]][, 2])^2)
    rss.tot <- sum((y[i:(i + 2 * w - 1)] - res.large[[i]][, 2])^2)

    diff1[i] <- (rss.tot + b) / (rss1 + rss2 + b)
    diff2[i] <- rss.tot - rss1 - rss2
  }

  # Combine components
  det <- (1 - a) * diff1 + a * diff2

  # Optional scaling
  if (scale_01) {
    det <- (det - min(det)) / (max(det) - min(det))
  }

  # Optional absolute value
  if (use_abs_det)
    det <- abs(det)

  # Smooth detector
  det <- as.numeric(ma(det, n = ma_window, circular = circular))
  det[is.na(det)] <- 0

  det
}
