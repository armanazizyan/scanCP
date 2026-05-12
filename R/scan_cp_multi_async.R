#' Multivariate Asynchronous Changepoint Detection
#'
#' @description
#' Applies the full univariate scanCP pipeline independently to each
#' dimension (column) of a multivariate signal matrix. This corresponds
#' to the *asynchronous* changepoint setting, where each dimension may
#' have changepoints at different locations.
#'
#' @param Y Numeric matrix of size n × p. Each column is a separate signal.
#' @param w Integer. Window size for rolling MLPs and detector.
#' @param mlp_params List of parameters passed to [fit_mlp()].
#' @param detector_params List of parameters passed to [calc_detector()].
#' @param decomp_params List of parameters passed to [decompose_signal_core()].
#' @param global_mlp Logical. Whether to run the optional global smoother.
#' @param global_params List of parameters passed to [fit_global_mlp()].
#'
#' @return A list of length p, where each element contains:
#' \describe{
#'   \item{fit_mlp}{Rolling MLP fit results for that dimension.}
#'   \item{detector}{Detector statistic for that dimension.}
#'   \item{decomposition}{Structural decomposition results.}
#'   \item{global}{Optional global MLP smoothing results.}
#' }
#'
#' @details
#' This function simply loops over columns and calls [scan_cp()] on each.
#' No cross‑dimension information is used. This is appropriate when
#' changepoints are *not* synchronized across dimensions.
#'
#' @export
scan_cp_multi_async <- function(
    Y,
    w = 100,
    mlp_params      = list(),
    detector_params = list(),
    decomp_params   = list(),
    global_mlp      = FALSE,
    global_params   = list()
) {

  if (!is.matrix(Y))
    stop("Y must be a numeric matrix.")

  p <- ncol(Y)
  out <- vector("list", p)
  names(out) <- colnames(Y) %||% paste0("dim", seq_len(p))

  for (j in seq_len(p)) {

    yj <- Y[, j]

    out[[j]] <- scan_cp(
      y               = yj,
      w               = w,
      mlp_params      = mlp_params,
      detector_params = detector_params,
      decomp_params   = decomp_params,
      global_mlp      = global_mlp,
      global_params   = global_params
    )
  }

  out
}
