#' Core Signal Decomposition (Changepoints + Piecewise Correction)
#'
#' @description
#' Extracts changepoints from a detector statistic using ECDF thresholding,
#' spacing-curve analysis, and local refinement. Produces the corrected
#' signal and its piecewise-constant component. This function does *not*
#' perform global MLP smoothing.
#'
#' @param y Numeric vector. Original signal.
#' @param detector Numeric vector. Detector statistic produced by
#'   \code{\link{calc_detector}}.
#' @param w Integer. Window size used in detector construction.
#' @param ma_window Integer. Moving-average smoothing window for the detector
#'   and spacing curve. Defaults to \code{w}.
#' @param right_tail_cutoff Numeric. Upper ECDF cutoff for automatic
#'   threshold selection. Defaults to \code{0.95}.
#' @param left_tail_cutoff Numeric. Lower ECDF cutoff for automatic
#'   threshold selection. Defaults to \code{0.6}.
#' @param threshold Either \code{"auto"} (default) or a numeric ECDF
#'   threshold for selecting significant peaks.
#' @param use_abs_det Logical. Whether to use \code{abs(detector)} before
#'   peak detection. Defaults to \code{TRUE}.
#' @param min_cp_distance Integer. Minimum distance between detected peaks.
#'   Defaults to \code{2 * w}.
#' @param margin Integer. Local refinement margin around each changepoint.
#'   Defaults to \code{floor(w / 2)}.
#' @param circular Logical. Whether moving-average smoothing wraps around.
#'   Defaults to \code{FALSE}. Typically left unchanged.
#'
#' @return A list containing:
#' \describe{
#'   \item{corrected_signal}{The refined, piecewise-corrected signal.}
#'   \item{piecewise_constant}{Estimated piecewise-constant component.}
#'   \item{raw_correction}{Cumulative correction vector.}
#'   \item{changepoints}{Refined changepoint locations.}
#'   \item{changepoint_ecdf}{ECDF values of selected peaks.}
#'   \item{shift_values}{Estimated shifts at each changepoint.}
#'   \item{smoothed_detector}{Smoothed detector statistic.}
#'   \item{ecdf_spacing}{Spacing-curve values (if \code{threshold = "auto"}).}
#'   \item{threshold_used}{Final threshold applied.}
#'   \item{local_maxima}{Matrix of detected local maxima.}
#' }
#'
#' @details
#' The function identifies local maxima of the smoothed detector, evaluates
#' their significance using the ECDF of the detector, optionally computes a
#' spacing curve for automatic threshold selection, and refines each
#' changepoint using a two-cluster k-means split within a local window.
#'
#' @importFrom pracma findpeaks
#' @importFrom stats ecdf median kmeans rnorm
#' @export
decompose_signal_core <- function(
    y,
    detector,
    w = 200,
    ma_window = 100,
    right_tail_cutoff = 0.95,
    left_tail_cutoff  = 0.6,
    threshold = "auto",
    use_abs_det = TRUE,
    min_cp_distance = NULL,
    circular = FALSE,
    margin = NULL
) {

  # --- 0. Preprocessing -------------------------------------------------------

  if (use_abs_det)
    detector <- abs(detector)

  n <- length(y)

  if (is.null(min_cp_distance))
    min_cp_distance <- 2 * w

  if (is.null(margin))
    margin <- floor(w / 2)

  # --- 1. Local maxima of raw detector ---------------------------------------

  # (Removed raw_peaks — unused)

  # --- 2. Smooth detector -----------------------------------------------------

  sm.det <- as.numeric(ma(detector, n = ma_window, circular = circular))
  sm.det[is.na(sm.det)] <- 0

  # --- 3. Local maxima of smoothed detector ----------------------------------

  l.max <- pracma::findpeaks(sm.det, minpeakdistance = min_cp_distance)
  if (is.null(l.max) || nrow(l.max) == 0) {
    return(list(
      corrected_signal   = y,
      piecewise_constant = rep(0, n),
      raw_correction     = rep(0, n),
      changepoints       = integer(0),
      changepoint_ecdf   = numeric(0),
      shift_values       = numeric(0),
      smoothed_detector  = sm.det,
      ecdf_spacing       = NULL,
      threshold_used     = NA,
      local_maxima       = NULL
    ))
  }

  # --- 4. ECDF significance ---------------------------------------------------

  f_det <- stats::ecdf(sm.det)
  ecdf_vals <- f_det(l.max[, 1])
  l.max <- cbind(l.max, ecdf_vals)  # col 5 = ECDF

  # --- 5. Spacing curve (only if needed) -------------------------------------

  if (is.character(threshold) && threshold == "auto") {
    x.det <- sort(sm.det)
    dx <- diff(x.det)
    s <- as.numeric(ma(dx, n = ma_window, circular = circular))
    s[is.na(s)] <- 0

    thr <- select_best_spike(
      s,
      right_tail_cutoff = right_tail_cutoff,
      left_tail_cutoff  = left_tail_cutoff
    )
  } else if (is.numeric(threshold)) {
    thr <- threshold
    s <- NULL
  } else {
    stop("threshold must be 'auto' or numeric.")
  }

  # --- 6. Changepoints (raw) --------------------------------------------------

  cps_raw <- l.max[l.max[, 5] >= thr, , drop = FALSE]

  if (nrow(cps_raw) == 0) {
    return(list(
      corrected_signal   = y,
      piecewise_constant = rep(0, n),
      raw_correction     = rep(0, n),
      changepoints       = integer(0),
      changepoint_ecdf   = numeric(0),
      shift_values       = numeric(0),
      smoothed_detector  = sm.det,
      ecdf_spacing       = s,
      threshold_used     = thr,
      local_maxima       = l.max
    ))
  }

  cps <- cps_raw[, 2] + w

  # --- 7. Local refinement ----------------------------------------------------

  refined <- c()
  shifts <- c()

  for (cp in sort(cps)) {

    lo <- max(1, cp - margin)
    hi <- min(n, cp + margin)

    segment <- y[lo:hi]

    b <- kmeans(segment, centers = 2)
    bsf <- best_split_free(b$cluster)

    corrected_cp <- cp + (margin - bsf$index)
    corrected_cp <- min(max(corrected_cp, 1), n)
    refined <- c(refined, corrected_cp)

    left_idx  <- max(1, corrected_cp - margin) : corrected_cp
    right_idx <- corrected_cp : min(n, corrected_cp + margin)

    shifts <- c(
      shifts,
      stats::median(y[left_idx]) - stats::median(y[right_idx])
    )
  }

  # --- 8. Build correction vector --------------------------------------------

  cps <- refined
  starts <- cps + 1
  ends <- c(cps[-1], n)

  cor.vec <- numeric(n)
  cum.shifts <- cumsum(shifts)

  for (i in seq_along(starts)) {
    cor.vec[starts[i]:ends[i]] <- cum.shifts[i]
  }

  new.y <- y + cor.vec

  # --- 9. Piecewise constant component ---------------------------------------

  idx_start <- c(1, cps + 1)
  idx_end   <- c(cps, n)

  step_mean <- numeric(n)
  for (i in seq_along(idx_start)) {
    step_mean[idx_start[i]:idx_end[i]] <-
      mean((y - new.y)[idx_start[i]:idx_end[i]])
  }

  # --- 10. Return -------------------------------------------------------------

  list(
    corrected_signal   = new.y,
    piecewise_constant = step_mean,
    raw_correction     = cor.vec,
    changepoints       = cps,
    changepoint_ecdf   = cps_raw[, 5],
    shift_values       = shifts,
    smoothed_detector  = sm.det,
    ecdf_spacing       = s,
    threshold_used     = thr,
    local_maxima       = l.max
  )
}




#' Fit Global MLP and Compute Residual Diagnostics
#'
#' @description
#' Fits a global MLP model to a corrected signal and computes:
#' - smooth trend estimate
#' - residuals
#'
#' This function does NOT compute a piecewise-constant component.
#' That is handled entirely in `decompose_signal_core()`.
#'
#' @param corrected_signal Numeric vector. Output of [decompose_signal_core()].
#'
#' @return A list containing:
#' \describe{
#'   \item{smooth_curve}{Global MLP smooth fit.}
#'   \item{residual}{Residuals from the smooth fit.}
#'   \item{mlp_model}{The fitted RSNNS MLP model.}
#' }
#'
#' @importFrom RSNNS mlp
#' @export
fit_global_mlp <- function(corrected_signal) {

  n <- length(corrected_signal)
  x <- 1:n

  # Fit global MLP
  sub_x <- as.matrix(scale(x))
  sub_y <- as.matrix(corrected_signal)

  net <- RSNNS::mlp(
    sub_x, sub_y,
    size = 64,
    learnFuncParams = c(0.01, 0.001),
    maxit = 2000,
    linOut = TRUE,
    hiddenActFunc = "Act_Logistic"
  )

  smooth <- net$fitted.values
  resid <- corrected_signal - smooth

  list(
    smooth_curve = smooth,
    residual     = resid,
    mlp_model    = net
  )
}
