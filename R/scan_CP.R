#' Full Deep Learning–Based Changepoint Detection Pipeline
#'
#' @description
#' Runs the complete scanCP pipeline:
#' 1. Rolling-window MLP fitting
#' 2. Detector computation
#' 3. Structural decomposition (changepoints + corrected signal)
#' 4. Optional global MLP smoothing
#'
#' All internal parameters can be passed through structured lists.
#'
#' @param y Numeric vector. The raw signal.
#' @param w Integer. Window size for rolling MLPs and detector.
#' @param mlp_params List of parameters passed to [fit_mlp()].
#' @param detector_params List of parameters passed to [calc_detector()].
#' @param decomp_params List of parameters passed to [decompose_signal_core()].
#' @param global_mlp Logical. Whether to run the optional global smoother.
#' @param global_params List of parameters passed to [fit_global_mlp()].
#'
#' @return A list containing:
#' \describe{
#'   \item{fit_mlp}{Rolling MLP fit results.}
#'   \item{detector}{Detector statistic.}
#'   \item{decomposition}{Structural decomposition results.}
#'   \item{global}{Optional global MLP smoothing results.}
#' }
#'
#' @export
scan_cp <- function(
    y,
    w = 100,
    mlp_params      = list(),
    detector_params = list(),
    decomp_params   = list(),
    global_mlp      = FALSE,
    global_params   = list()
) {

  # ------------------------------------------------------------
  # 1. Fit rolling MLPs
  # ------------------------------------------------------------
  fit_args <- c(list(vec = y, w = w), mlp_params)
  fit_res <- do.call(fit_mlp, fit_args)

  # ------------------------------------------------------------
  # 2. Compute detector
  # ------------------------------------------------------------
  det_args <- c(list(y = y, fit_mlp_res = fit_res, w = w), detector_params)
  det <- do.call(calc_detector, det_args)

  # ------------------------------------------------------------
  # 3. Structural decomposition
  # ------------------------------------------------------------
  decomp_args <- c(list(y = y, diff = det, w = w), decomp_params)
  core <- do.call(decompose_signal_core, decomp_args)

  # ------------------------------------------------------------
  # 4. Optional global smoother
  # ------------------------------------------------------------
  if (global_mlp) {
    global_args <- c(list(corrected_signal = core$corrected_signal), global_params)
    global <- do.call(fit_global_mlp, global_args)
  } else {
    global <- NULL
  }

  # ------------------------------------------------------------
  # 5. Return full pipeline output
  # ------------------------------------------------------------
  list(
    fit_mlp      = fit_res,
    detector     = det,
    decomposition = core,
    global       = global
  )
}
