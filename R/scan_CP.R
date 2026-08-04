#' Full Changepoint Detection Pipeline Using Rolling MLPs
#'
#' @description
#' Runs the complete changepoint detection pipeline on a univariate signal.
#' This includes rolling-window MLP fitting, detector construction,
#' thresholding, structural decomposition, and changepoint extraction.
#'
#' @param y Numeric vector. The input signal.
#' @param w Integer. Window size for the small-window MLP used in
#'   \code{\link{fit_mlp}}. The large-window MLP automatically uses
#'   window size \code{2w}.
#' @param smoothing_window Integer. Window size for the moving-average smoothing
#'   applied to the detector statistic and in the decomposition. Defaults to
#'   \code{w}. This parameter replaces \code{ma_window} for simplicity.
#' @param threshold Character or numeric. If \code{"auto"}, the threshold is
#'   selected automatically using spacing-curve analysis. If \code{"auc"},
#'   changepoints are selected via net AUC with percentile-based filtering.
#'   Otherwise, a numeric threshold may be supplied directly.
#' @param threshold_tails Numeric vector of length 2. Left and right tail
#'   cutoffs used when estimating the automatic threshold. Defaults to
#'   \code{c(0.2, 0.95)}.
#' @param auc_energy_threshold Numeric or NULL. When \code{threshold = "auc"}, peaks
#'   with net AUC values less than \code{auc_energy_threshold * sum(all AUC values)}
#'   are filtered out. If NULL (default), computed as \code{1 / k} where \code{k}
#'   is the number of detected peaks, providing an adaptive threshold.
#' @param min_cp_distance Integer. Minimum separation (in indices) between
#'   detected changepoints. Defaults to \code{2 * w}.
#' @param margin Integer. Margin used during local refinement of changepoint
#'   locations. Defaults to \code{floor(w / 2)}.
#' @param use_abs_det Logical. Whether to use the absolute value of the
#'   detector statistic when identifying changepoints.
#' @param mlp_control A named list of hyperparameters passed directly to
#'   \code{\link{fit_mlp}}. Any subset may be supplied; unspecified values
#'   fall back to defaults defined in \code{fit_mlp()}.
#' @param parallel Logical. If \code{TRUE}, rolling MLP fitting is
#'   performed in parallel via \code{\link{fit_mlp}}. If \code{FALSE}, all
#'   computation is performed serially. Setting \code{parallel = FALSE} is
#'   recommended for CRAN checks and for systems without parallel support.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{changepoints}}{Estimated changepoint locations.}
#'   \item{\code{detector}}{The detector statistic computed by
#'         \code{\link{calc_detector}}.}
#'   \item{\code{fit_mlp}}{Rolling MLP fit results returned by
#'         \code{\link{fit_mlp}}.}
#'   \item{\code{decomposition}}{Full structural decomposition returned by
#'         \code{\link{decompose_signal_core}}.}
#' }
#'
#' @details
#' The pipeline proceeds in three main stages:
#' \enumerate{
#'   \item Rolling-window MLP fitting via \code{\link{fit_mlp}}.
#'   \item Detector construction via \code{\link{calc_detector}}.
#'   \item Structural decomposition and changepoint extraction via
#'         \code{\link{decompose_signal_core}}.
#' }
#'
#' Detected changepoints are printed for user convenience and returned as part
#' of the output list. Parallel computation is optional and controlled by the
#' \code{parallel} argument.
#'
#'
#' @examples
#' # Minimal example
#' set.seed(1)
#' y <- c(rnorm(200, 0), rnorm(200, 3))
#'
#' \donttest{
#'   # Full pipeline (parallel disabled for CRAN)
#'   res <- scan_cp(y, w = 20, parallel = FALSE)
#' }
#'
#'
#' @export
scan_cp <- function(
    y,
    w = 100,
    smoothing_window = NULL,
    threshold = "auto",
    threshold_tails = c(0.2, 0.95),
    auc_energy_threshold = NULL,
    min_cp_distance = 2 * w,
    margin = floor(w / 2),
    use_abs_det = TRUE,
    mlp_control = list(),
    parallel = FALSE
) {

  # Set smoothing_window to w if not provided
  if (is.null(smoothing_window))
    smoothing_window <- w

  # 1. Fit rolling MLPs
  fit_args <- c(list(vec = y, w = w, parallel = parallel), mlp_control)
  fit_res <- do.call(fit_mlp, fit_args)

  # 2. Compute detector
  det <- calc_detector(
    y = y,
    fit_mlp_res = fit_res,
    w = w,
    ma_window = smoothing_window,
    use_abs_det = use_abs_det
  )

  # 3. Structural decomposition
  core <- decompose_signal_core(
    y = y,
    detector = det,
    w = w,
    ma_window = smoothing_window,
    threshold = threshold,
    auc_energy_threshold = auc_energy_threshold,
    right_tail_cutoff = threshold_tails[2],
    left_tail_cutoff  = threshold_tails[1],
    min_cp_distance = min_cp_distance,
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

