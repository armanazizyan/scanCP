#' Combine Multiple Detector Statistics into a Joint Detector
#'
#' @description
#' Combines several univariate detector vectors into a single joint detector
#' using L1, L2, or max aggregation. Optionally applies moving-average
#' smoothing and computes contribution weights.
#'
#' @param det.lst List of numeric detector vectors (all same length).
#' @param method One of `"L1"`, `"L2"`, `"max"`.
#' @param ma_window Optional integer. If provided, each detector is smoothed
#'   using `ma()` with this window size.
#' @param circular Logical. Whether smoothing wraps around.
#' @param scale_contributions Logical. Whether to compute contribution weights.
#' @param scale_joint Logical. Whether to scale the joint detector to \[0,1\].
#'
#' @return A list containing:
#' \describe{
#'   \item{joint}{Joint detector vector.}
#'   \item{s.mat}{Scaled detector matrix.}
#'   \item{det.mat}{Raw detector matrix.}
#'   \item{contrib}{Contribution weights (or NULL).}
#' }
#'
#' @export
combine_detectors <- function(
    det.lst,
    method = c("L1", "L2", "max"),
    ma_window = NULL,
    circular = FALSE,
    scale_contributions = TRUE,
    scale_joint = TRUE
) {

  method <- match.arg(method)

  # 1. Optional smoothing
  if (!is.null(ma_window)) {
    det.lst <- lapply(det.lst, ma, n = ma_window, circular = circular)
  }

  # 2. Convert to matrix
  det.mat <- do.call(cbind, lapply(det.lst, as.vector))

  # 3. Scaling (currently identity)
  scale.mat <- det.mat

  # 4. Joint detector
  joint.det <- switch(
    method,
    L1  = rowSums(abs(scale.mat)),
    L2  = sqrt(rowSums(scale.mat^2)),
    max = apply(scale.mat, 1, max)
  )

  # 5. Scale joint detector
  if (scale_joint) {
    rng <- range(joint.det, na.rm = TRUE)
    joint.det <- (joint.det - rng[1]) / (rng[2] - rng[1])
  }

  # 6. Contribution weights
  if (scale_contributions) {
    contrib <- switch(
      method,
      L1  = abs(scale.mat),
      L2  = scale.mat^2,
      max = pmax(scale.mat, 0)
    )
    contrib <- contrib / rowSums(contrib)
  } else {
    contrib <- NULL
  }

  list(
    joint   = as.numeric(joint.det),
    s.mat   = scale.mat,
    det.mat = det.mat,
    contrib = contrib
  )
}

#' ECDF-Based Changepoint Detection for a Detector Statistic
#'
#' @description
#' Applies the ECDF + spacing-curve thresholding method to a detector vector
#' and returns raw changepoints (no correction).
#'
#' @param diff Numeric detector vector.
#' @param w Integer. Window size used in detector construction.
#' @param ma_window Integer. Moving-average window for smoothing.
#' @param right_tail_cutoff Numeric. ECDF upper cutoff.
#' @param left_tail_cutoff Numeric. ECDF lower cutoff.
#' @param threshold "auto" or numeric ECDF threshold.
#' @param circular Logical. Whether smoothing wraps around.
#'
#' @return A list containing:
#' \describe{
#'   \item{changepoints}{Raw changepoint indices.}
#'   \item{ecdf_values}{ECDF values of selected maxima.}
#'   \item{threshold_used}{Final threshold.}
#'   \item{local_maxima}{Matrix of local maxima.}
#'   \item{smoothed_detector}{Smoothed detector.}
#'   \item{ecdf_spacing}{Spacing curve.}
#' }
#'
#' @importFrom pracma findpeaks
#' @importFrom stats ecdf
#' @export
detect_cp_ecdf <- function(
    diff,
    w = 200,
    ma_window = 100,
    right_tail_cutoff = 0.95,
    left_tail_cutoff  = 0.6,
    threshold = "auto",
    circular = FALSE
) {

  # 1. Smooth detector
  ma.diff <- as.numeric(ma(diff, n = ma_window, circular = circular))

  # 2. Local maxima
  l.max <- pracma::findpeaks(replace_na(ma.diff, 0), minpeakdistance = 2 * w)

  if (is.null(l.max) || nrow(l.max) == 0) {
    return(list(
      changepoints      = integer(0),
      ecdf_values       = numeric(0),
      threshold_used    = NA,
      local_maxima      = NULL,
      smoothed_detector = ma.diff
    ))
  }

  # 3. ECDF significance
  f_diff <- stats::ecdf(ma.diff)
  ecdf_vals <- f_diff(l.max[, 1])
  l.max <- cbind(l.max, ecdf_vals)

  # 4. ECDF spacing curve
  x.diff <- sort(ma.diff)
  dx <- diff(x.diff)
  s <- as.numeric(ma(dx, n = w, circular = FALSE))
  s[is.na(s)] <- 0

  # 5. Threshold selection
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

  # 6. Raw changepoints
  cps_raw <- l.max[l.max[, 5] >= thr, , drop = FALSE]
  cps <- cps_raw[, 2] + w

  list(
    changepoints      = cps,
    ecdf_values       = cps_raw[, 5],
    threshold_used    = thr,
    local_maxima      = l.max,
    smoothed_detector = ma.diff,
    ecdf_spacing      = s
  )
}

#' Multivariate Synchronized Changepoint Detection
#'
#' @description
#' Implements the synchronized multivariate changepoint pipeline:
#' 1. Fit univariate MLP-based detectors for each dimension.
#' 2. Combine detectors into a joint detector using L1/L2/max.
#' 3. Detect changepoints on the joint detector.
#' 4. Compute per-dimension contributions at each detected CP.
#'
#' @param Y Numeric matrix (n × p). Each column is a signal dimension.
#' @param w Integer. Window size for rolling MLPs and detector.
#' @param method Character. Combination method: "L1", "L2", or "max".
#' @param mlp_params List of parameters passed to `fit_mlp()`.
#' @param detector_params List of parameters passed to `calc_detector()`.
#' @param combine_params List of parameters passed to `combine_detectors()`.
#' @param ecdf_params List of parameters passed to `detect_cp_ecdf()`.
#'
#' @return A list containing:
#' \describe{
#'   \item{detectors}{List of univariate detectors (one per dimension).}
#'   \item{joint_detector}{Joint detector vector.}
#'   \item{contributions}{Matrix of per-dimension contributions.}
#'   \item{changepoints}{Detected synchronized changepoints.}
#'   \item{cp_contributions}{Contribution of each dimension at each CP.}
#'   \item{ecdf_output}{Full output of `detect_cp_ecdf()`.}
#' }
#'
#' @export
scan_cp_multi_sync <- function(
    Y,
    w = 100,
    method = c("L1", "L2", "max"),
    mlp_params      = list(),
    detector_params = list(),
    combine_params  = list(),
    ecdf_params     = list()
) {

  method <- match.arg(method)

  if (!is.matrix(Y))
    stop("Y must be a numeric matrix.")

  p <- ncol(Y)
  n <- nrow(Y)

  # ------------------------------------------------------------
  # 1. Fit univariate detectors for each dimension
  # ------------------------------------------------------------
  det.lst <- vector("list", p)
  names(det.lst) <- colnames(Y) %||% paste0("dim", seq_len(p))

  for (j in seq_len(p)) {

    yj <- Y[, j]

    # Fit MLPs
    fit_args <- c(list(vec = yj, w = w), mlp_params)
    fit_res <- do.call(fit_mlp, fit_args)

    # Compute detector
    det_args <- c(list(y = yj, fit_mlp_res = fit_res, w = w), detector_params)
    det.lst[[j]] <- do.call(calc_detector, det_args)
  }

  # ------------------------------------------------------------
  # 2. Combine detectors into joint detector
  # ------------------------------------------------------------
  comb_args <- c(list(det.lst = det.lst, method = method), combine_params)
  comb <- do.call(combine_detectors, comb_args)

  joint.det <- comb$joint
  contrib.mat <- comb$contrib   # n × p matrix

  # ------------------------------------------------------------
  # 3. Detect changepoints on joint detector
  # ------------------------------------------------------------
  ecdf_args <- c(list(diff = joint.det, w = w), ecdf_params)
  ecdf_out <- do.call(detect_cp_ecdf, ecdf_args)

  cps <- ecdf_out$changepoints

  # ------------------------------------------------------------
  # 4. Contribution of each dimension at each CP
  # ------------------------------------------------------------
  if (length(cps) > 0) {
    cp_contrib <- contrib.mat[cps, , drop = FALSE]
    rownames(cp_contrib) <- paste0("cp_", cps)
  } else {
    cp_contrib <- matrix(0, nrow = 0, ncol = p)
  }

  # ------------------------------------------------------------
  # 5. Return full synchronized output
  # ------------------------------------------------------------
  list(
    detectors        = det.lst,
    joint_detector   = joint.det,
    contributions    = contrib.mat,
    changepoints     = cps,
    cp_contributions = cp_contrib,
    ecdf_output      = ecdf_out
  )
}
