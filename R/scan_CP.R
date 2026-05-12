#' Deep Learning–Based Changepoint Detection (Core Pipeline)
#'
#' @description
#' Implements the core scanCP changepoint detection pipeline:
#' 1. Rolling-window MLP fitting
#' 2. Detector construction and smoothing
#' 3. ECDF-based changepoint extraction
#' 4. Local refinement and piecewise correction
#'
#' This function performs the full univariate changepoint analysis and
#' returns the detected changepoints at the top level for convenience.
#' The optional global MLP smoother is intentionally excluded and can be
#' added via \code{\link{add_global_mlp}}.
#'
#' @param y Numeric vector. The raw input signal.
#' @param w Integer. Window size used for rolling MLPs and detector
#'   construction. Defaults to 100.
#' @param ma_window Integer. Moving-average smoothing window for the
#'   detector. Defaults to \code{w}.
#' @param threshold Either \code{"auto"} (default) or a numeric ECDF
#'   threshold for selecting significant peaks.
#' @param threshold_tails Numeric vector of length 2 giving the
#'   \code{c(left_tail_cutoff, right_tail_cutoff)} used in automatic
#'   threshold selection. Defaults to \code{c(0.6, 0.95)}.
#' @param min_cp_distance Integer. Minimum distance between detected
#'   changepoints (peak separation). Defaults to \code{2 * w}.
#' @param margin Integer. Local refinement margin around each detected
#'   changepoint. Defaults to \code{floor(w / 2)}.
#' @param use_abs_det Logical. Whether to use \code{abs(detector)} before
#'   peak detection. Defaults to \code{TRUE}.
#' @param mlp_control List of parameters passed to \code{\link{fit_mlp}}
#'   (e.g., hidden layer sizes, learning rates, iteration limits).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{changepoints}{Vector of detected changepoint locations.}
#'   \item{detector}{The smoothed detector statistic.}
#'   \item{fit_mlp}{Rolling-window MLP fit results.}
#'   \item{decomposition}{Full output from
#'     \code{\link{decompose_signal_core}} including corrected signal,
#'     piecewise constant component, ECDF values, spacing curve, and
#'     refinement details.}
#' }
#'
#' @details
#' This function is the primary user-facing entry point for univariate
#' changepoint detection. It exposes only the most important parameters
#' while keeping advanced MLP hyperparameters inside \code{mlp_control}.
#'
#' The global smoother is intentionally not included here; use
#' \code{\link{add_global_mlp}} to attach a global MLP trend estimate.
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   y <- c(rnorm(200, 0, 1),
#'          rnorm(200, 3, 1),
#'          rnorm(200, -2, 1))
#'
#'   res <- scan_cp(y, w = 100)
#'   res$changepoints
#' }
#'
#' @export
scan_cp <- function(
    y,
    w = 100,
    ma_window = w,
    threshold = "auto",
    threshold_tails = c(0.2, 0.95),
    min_cp_distance = 2 * w,
    margin = floor(w / 2),
    use_abs_det = TRUE,
    mlp_control = list()
) {

  # 1. Fit rolling MLPs
  fit_args <- c(list(vec = y, w = w), mlp_control)
  fit_res <- do.call(fit_mlp, fit_args)

  # 2. Compute detector
  det <- calc_detector(
    y = y,
    fit_mlp_res = fit_res,
    w = w,
    ma_window = ma_window,
    use_abs_det = use_abs_det
  )

  # 3. Structural decomposition
  core <- decompose_signal_core(
    y = y,
    detector = det,
    w = w,
    ma_window = ma_window,
    threshold = threshold,
    right_tail_cutoff = threshold_tails[2],
    left_tail_cutoff  = threshold_tails[1],
    minpeakdistance = min_cp_distance,
    margin = margin,
    use_abs_det = use_abs_det
  )

  # 4. Extract CPs at top level
  cps <- core$changepoints

  # 5. Print CPs for user convenience
  if (length(cps) == 0) {
    message("No changepoints detected.")
  } else {
    message("Detected changepoints at: ", paste(cps, collapse = ", "))
  }

  # 6. Return structured output
  list(
    changepoints   = cps,
    detector       = det,
    fit_mlp        = fit_res,
    decomposition  = core
  )
}
