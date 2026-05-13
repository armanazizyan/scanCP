#' Fit Rolling MLP Models for Change Point Detection
#'
#' @description
#' Fits two multilayer perceptron (MLP) models over rolling windows of a
#' univariate signal. A small-window MLP (window size \code{w}) captures
#' local structure, while a large-window MLP (window size \code{2w})
#' captures broader trends. Their residual behavior is used by
#' \code{\link{calc_detector}} to construct the changepoint detector.
#'
#' All MLP hyperparameters are supplied through a unified
#' \code{mlp_control} list for consistency with \code{\link{scan_cp}}.
#'
#' @param vec Numeric vector. The input signal.
#' @param w Integer. Window size for the small-window MLP. The large-window
#'   MLP automatically uses window size \code{2w}.
#' @param mlp_control A named list of MLP hyperparameters. Any subset may be
#'   supplied; unspecified values fall back to defaults. Supported fields:
#'   \describe{
#'     \item{\code{act1}}{Activation function for the small-window MLP
#'       (default: \code{"Act_Logistic"}).}
#'     \item{\code{act2}}{Activation function for the large-window MLP
#'       (default: \code{"Act_Logistic"}).}
#'     \item{\code{lr1}}{Learning rate for the small-window MLP
#'       (default: \code{0.001}).}
#'     \item{\code{lr2}}{Learning learning rate for the large-window MLP
#'       (default: \code{0.001}).}
#'     \item{\code{hl1}}{Number of hidden units for the small-window MLP
#'       (default: \code{6}).}
#'     \item{\code{hl2}}{Number of hidden units for the large-window MLP
#'       (default: \code{8}).}
#'     \item{\code{ep1}}{Training epochs for the small-window MLP
#'       (default: \code{1000}).}
#'     \item{\code{ep2}}{Training epochs for the large-window MLP
#'       (default: \code{2000}).}
#'   }
#'
#' @param parallel Logical. If \code{TRUE} (default), rolling-window MLPs are
#'   fitted using parallel computation via \pkg{parallel}, \pkg{foreach},
#'   and \pkg{doSNOW}. If \code{FALSE}, all computation is performed
#'   serially. Setting \code{parallel = FALSE} is recommended for CRAN checks
#'   and for systems where parallel backends are unavailable.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{small}}{List of fitted values for each rolling window of size \code{w}.}
#'   \item{\code{large}}{List of fitted values for each rolling window of size \code{2w}.}
#' }
#'
#' @details
#' Each rolling window is standardized before fitting. The returned matrices
#' contain three columns:
#' \enumerate{
#'   \item standardized input \code{x}
#'   \item fitted values \code{y_hat}
#'   \item original indices
#' }
#'
#' When \code{parallel = TRUE}, a cluster is created using
#' \code{parallel::makeCluster()} and progress bars are displayed via
#' \pkg{doSNOW}. When \code{parallel = FALSE}, the function falls back to
#' serial execution using \code{\%do\%}.
#'
#' @import RSNNS
#' @import foreach
#' @import doSNOW
#' @import parallel
#' @importFrom utils modifyList txtProgressBar setTxtProgressBar
#' @export
fit_mlp <- function(
    vec,
    w = 100,
    mlp_control = list(),
    parallel = TRUE
) {

  # Defaults
  default_mlp_control <- list(
    act1 = "Act_Logistic",
    act2 = "Act_Logistic",
    lr1  = 0.001,
    lr2  = 0.001,
    hl1  = 6,
    hl2  = 8,
    ep1  = 300,
    ep2  = 600
  )

  mlp <- modifyList(default_mlp_control, mlp_control)

  # Extract
  act1 <- mlp$act1; act2 <- mlp$act2
  lr1  <- mlp$lr1;  lr2  <- mlp$lr2
  hl1  <- mlp$hl1;  hl2  <- mlp$hl2
  ep1  <- mlp$ep1;  ep2  <- mlp$ep2

  n.val <- length(vec)
  x <- 1:n.val
  y <- vec

  # -----------------------------
  # Parallel or serial backend
  # -----------------------------
  if (parallel) {
    num_cores <- max(1, parallel::detectCores() - 1)
    cl <- parallel::makeCluster(num_cores)
    doSNOW::registerDoSNOW(cl)

    pb <- txtProgressBar(min = 0, max = n.val - w, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    `%dopar_or_do%` <- `%dopar%`
  } else {
    # Serial mode: no cluster, no progress bar
    opts <- list()
    `%dopar_or_do%` <- `%do%`
  }

  # -----------------------------
  # Small-window MLP
  # -----------------------------
  if (parallel) message("Fitting small-window MLPs...")

  res.small <- foreach::foreach(
    idx = 1:(n.val - w + 1),
    .packages = "RSNNS",
    .options.snow = opts
  ) %dopar_or_do% {

    i <- idx
    sub_x <- as.matrix(scale(x[i:(i + w - 1)]))
    sub_y <- as.matrix(y[i:(i + w - 1)])

    net <- RSNNS::mlp(
      sub_x, sub_y,
      size = hl1,
      learnFuncParams = lr1,
      maxit = ep1,
      linOut = TRUE,
      hiddenActFunc = act1
    )

    cbind(sub_x, net$fitted.values, i:(i + w - 1))
  }

  if (parallel) close(pb)

  # -----------------------------
  # Large-window MLP
  # -----------------------------
  if (parallel) {
    pb <- txtProgressBar(min = 0, max = n.val - 2 * w, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }

  if (parallel) message("Fitting large-window MLPs...")

  res.large <- foreach::foreach(
    idx = 1:(n.val - 2 * w + 1),
    .packages = "RSNNS",
    .options.snow = opts
  ) %dopar_or_do% {

    i <- idx
    sub_x <- as.matrix(scale(x[i:(i + 2 * w - 1)]))
    sub_y <- as.matrix(y[i:(i + 2 * w - 1)])

    net <- RSNNS::mlp(
      sub_x, sub_y,
      size = hl2,
      learnFuncParams = lr2,
      maxit = ep2,
      linOut = TRUE,
      hiddenActFunc = act2
    )

    cbind(sub_x, net$fitted.values, i:(i + 2 * w - 1))
  }

  if (parallel) {
    close(pb)
    parallel::stopCluster(cl)
  }

  list(
    small = res.small,
    large = res.large
  )
}
