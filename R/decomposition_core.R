#' Core Signal Decomposition (Changepoints + Piecewise Correction)
#'
#' @description
#' Performs the changepoint extraction, ECDF thresholding, spacing-curve
#' spike selection, and local refinement. Produces the corrected signal
#' and piecewise-constant component, but does NOT fit the global MLP.
#'
#' @param y Numeric vector. Original signal.
#' @param diff Numeric vector. Detector statistic.
#' @param w Integer. Window size used in detector.
#' @param ma_window Integer. Moving-average window.
#' @param right_tail_cutoff Numeric. ECDF upper cutoff.
#' @param left_tail_cutoff Numeric. ECDF lower cutoff.
#' @param threshold "auto" or numeric ECDF threshold.
#'
#' @return A list containing:
#' \describe{
#'   \item{corrected_signal}{Signal after piecewise correction.}
#'   \item{piecewise_constant}{Piecewise constant component.}
#'   \item{raw_correction}{Correction vector.}
#'   \item{changepoints}{Corrected changepoint indices.}
#'   \item{changepoint_ecdf}{ECDF values of selected changepoints.}
#'   \item{shift_values}{Estimated shifts.}
#'   \item{p_values}{Wilcoxon p-values.}
#'   \item{smoothed_detector}{Smoothed detector.}
#'   \item{ecdf_spacing}{Spacing curve.}
#'   \item{threshold_used}{Final threshold.}
#'   \item{local_maxima}{Local maxima matrix.}
#' }
#'
#' @importFrom pracma findpeaks
#' @importFrom stats ecdf median
#' @export
decompose_signal_core <- function(
    y,
    diff,
    w = 200,
    ma_window = 100,
    right_tail_cutoff = 0.95,
    left_tail_cutoff  = 0.6,
    threshold = "auto"
) {

  n <- length(y)

  # 1. Local maxima of raw detector
  pre.sm.local <- pracma::findpeaks(diff, minpeakdistance = w)

  # 2. Smooth detector
  ma.diff <- as.numeric(ma(diff, ma_window))
  ma.diff[is.na(ma.diff)] <- 0

  # 3. Local maxima of smoothed detector
  l.max <- pracma::findpeaks(ma.diff, minpeakdistance = w)

  # 4. ECDF significance
  f_diff <- stats::ecdf(ma.diff)
  ecdf_vals <- f_diff(l.max[, 1])
  l.max <- cbind(l.max, ecdf_vals)

  # 5. ECDF spacing curve
  x.diff <- sort(ma.diff)
  dx <- diff(x.diff)
  s <- as.numeric(ma(dx, w, circular = FALSE))
  s[is.na(s)] <- 0

  # 6. Threshold selection
  if (is.character(threshold) && threshold == "auto") {
    thr <- select_best_spike(
      s,
      right_tail_cutoff = right_tail_cutoff,
      left_tail_cutoff  = left_tail_cutoff
    )
  } else if (is.numeric(threshold)) {
    thr <- threshold
  } else {
    stop("threshold must be 'auto' or numeric.")
  }

  # 7. Changepoints
  cps_raw <- l.max[l.max[, 5] >= thr, , drop = FALSE]
  cps <- cps_raw[, 2] + w

  # 8. Local correction
  marg <- w
  cor.cps <- c()
  shift.vals <- c()
  wc.p <- c()

  for (cp in sort(cps)) {

    lo <- max(1, cp - marg)
    hi <- min(n, cp + marg)

    segment <- y[lo:hi]

    b <- kmeans(segment, centers = 2)
    bsf <- best_split_free(b$cluster)

    corrected_cp <- cp + (marg - bsf$index)
    corrected_cp <- min(max(corrected_cp, 1), n)
    cor.cps <- c(cor.cps, corrected_cp)

    w1_start <- max(1, cp - bsf$index)
    w1_end   <- min(n, cp + marg - bsf$index)

    w2_start <- max(1, cp + marg - bsf$index)
    w2_end   <- min(n, cp + 2 * marg - bsf$index)

    shift.vals <- c(
      shift.vals,
      stats::median(y[w1_start:w1_end]) - stats::median(y[w2_start:w2_end])
    )

    wc.p <- c(
      wc.p,
      wilcox.test(
        y[w1_start:w1_end],
        y[w2_start:w2_end],
        exact = FALSE
      )$p.value
    )
  }

  # 9. Correction vector
  cps <- cor.cps
  starts <- cps + 1
  ends <- c(cps[-1], n)

  cor.vec <- numeric(n)
  cumulative.shift.vals <- cumsum(shift.vals)

  for (i in seq_along(starts)) {
    cor.vec[starts[i]:ends[i]] <- cumulative.shift.vals[i]
  }

  # 10. Apply correction
  new.y <- y + cor.vec

  # 11. Piecewise constant component
  idx_start <- c(1, cps + 1)
  idx_end <- c(cps, n)

  step_mean <- numeric(n)
  for (i in seq_along(idx_start)) {
    step_mean[idx_start[i]:idx_end[i]] <-
      mean((y - new.y)[idx_start[i]:idx_end[i]])
  }

  list(
    corrected_signal   = new.y,
    piecewise_constant = step_mean,
    raw_correction     = cor.vec,
    changepoints       = cps,
    changepoint_ecdf   = cps_raw[, 5],
    shift_values       = shift.vals,
    p_values           = wc.p,
    smoothed_detector  = ma.diff,
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
