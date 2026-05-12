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
#'     \item{\code{lr2}}{Learning rate for the large-window MLP
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
#' Parallel computation is used to accelerate rolling-window fitting.
#'
#' @import RSNNS
#' @import foreach
#' @import doSNOW
#' @import parallel
#' @export
fit_mlp <- function(
    vec,
    w = 100,
    mlp_control = list()
) {

  # ------------------------------------------------------------
  # Default hyperparameters
  # ------------------------------------------------------------
  default_mlp_control <- list(
    act1 = "Act_Logistic",
    act2 = "Act_Logistic",
    lr1  = 0.001,
    lr2  = 0.001,
    hl1  = 6,
    hl2  = 8,
    ep1  = 1000,
    ep2  = 2000
  )

  # Merge defaults with user overrides
  mlp <- modifyList(default_mlp_control, mlp_control)

  # Extract parameters
  act1 <- mlp$act1; act2 <- mlp$act2
  lr1  <- mlp$lr1;  lr2  <- mlp$lr2
  hl1  <- mlp$hl1;  hl2  <- mlp$hl2
  ep1  <- mlp$ep1;  ep2  <- mlp$ep2

  # ------------------------------------------------------------
  # Setup
  # ------------------------------------------------------------
  n.val <- length(vec)
  x <- 1:n.val
  y <- vec

  num_cores <- max(1, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(num_cores)
  doSNOW::registerDoSNOW(cl)

  # ------------------------------------------------------------
  # Small-window MLP
  # ------------------------------------------------------------
  pb <- txtProgressBar(min = 0, max = n.val - w, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  message("Fitting small-window MLPs...")

  res.small <- foreach::foreach(
    idx = 1:(n.val - w + 1),
    .packages = c("RSNNS"),
    .options.snow = opts
  ) %dopar% {

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

    y_hat <- net$fitted.values
    cbind(sub_x, y_hat, i:(i + w - 1))
  }

  close(pb)

  # ------------------------------------------------------------
  # Large-window MLP
  # ------------------------------------------------------------
  pb <- txtProgressBar(min = 0, max = n.val - 2 * w, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  message("Fitting large-window MLPs...")

  res.large <- foreach::foreach(
    idx = 1:(n.val - 2 * w + 1),
    .packages = c("RSNNS"),
    .options.snow = opts
  ) %dopar% {

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

    y_hat <- net$fitted.values
    cbind(sub_x, y_hat, i:(i + 2 * w - 1))
  }

  close(pb)
  parallel::stopCluster(cl)

  list(
    small = res.small,
    large = res.large
  )
}
