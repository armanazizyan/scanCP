#' Asynchronous Multivariate Changepoint Detection
#'
#' @description
#' Applies the univariate \code{\link{scan_cp}} pipeline independently to
#' each column of a multivariate signal. Each dimension is processed
#' separately (asynchronously), producing a list of univariate results.
#'
#' @param Y Numeric matrix (n × p). Each column is a signal dimension.
#' @param w Integer. Window size for rolling MLPs and detector construction.
#'   Defaults to 100.
#' @param ma_window Integer. Moving-average smoothing window for the detector.
#'   Defaults to \code{w}.
#' @param threshold Either \code{"auto"} (default) or a numeric ECDF threshold.
#' @param threshold_tails Numeric vector of length 2 giving
#'   \code{c(left_tail_cutoff, right_tail_cutoff)} for automatic thresholding.
#'   Defaults to \code{c(0.6, 0.95)}.
#' @param min_cp_distance Integer. Minimum distance between detected
#'   changepoints. Defaults to \code{2 * w}.
#' @param margin Integer. Local refinement margin. Defaults to \code{floor(w/2)}.
#' @param use_abs_det Logical. Whether to use \code{abs(detector)} before
#'   smoothing. Defaults to \code{TRUE}.
#' @param mlp_control List of parameters passed to \code{\link{fit_mlp}}.
#'
#' @return A named list of length \code{p}, where each element is the output
#'   of \code{\link{scan_cp}} applied to the corresponding column of \code{Y}.
#'
#' @examples
#' \dontrun{
#'   Y <- cbind(
#'     x1 = c(rnorm(200), rnorm(200, 3)),
#'     x2 = c(rnorm(200), rnorm(200, -2))
#'   )
#'
#'   res <- scan_cp_multi_async(Y, w = 100)
#'   res$x1$changepoints
#'   res$x2$changepoints
#' }
#'
#' @export
scan_cp_multi_async <- function(
    Y,
    w = 100,
    ma_window = w,
    threshold = "auto",
    threshold_tails = c(0.6, 0.95),
    min_cp_distance = 2 * w,
    margin = floor(w / 2),
    use_abs_det = TRUE,
    mlp_control = list()
) {

  if (!is.matrix(Y))
    stop("Y must be a numeric matrix.")

  p <- ncol(Y)
  out <- vector("list", p)
  names(out) <- colnames(Y) %||% paste0("dim", seq_len(p))

  for (j in seq_len(p)) {
    out[[j]] <- scan_cp(
      y = Y[, j],
      w = w,
      ma_window = ma_window,
      threshold = threshold,
      threshold_tails = threshold_tails,
      min_cp_distance = min_cp_distance,
      margin = margin,
      use_abs_det = use_abs_det,
      mlp_control = mlp_control
    )
  }

  out
}
