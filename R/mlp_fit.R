#' Fit Rolling MLP Models for Change Point Detection
#'
#' @description
#' Fits two multilayer perceptron (MLP) models over rolling windows of a
#' univariate signal. The first model uses a window of size `w`, and the
#' second uses a window of size `2w`. The double-window model captures
#' uncertainty across adjacent segments and is used in the detector statistic.
#'
#' Parallel computation is used to accelerate rolling-window fitting.
#'
#' @param vec Numeric vector. The original signal.
#' @param w Integer. Window size for the smaller MLP.
#' @param act1,act2 Character strings. Activation functions for the smaller
#'   and larger MLPs (RSNNS activation names).
#' @param lr1,lr2 Numeric. Learning rates for the smaller and larger MLPs.
#' @param hl1,hl2 Integer. Number of hidden units for the smaller and larger MLPs.
#' @param ep1,ep2 Integer. Number of training epochs for the smaller and larger MLPs.
#'
#' @return A list of two elements:
#' \describe{
#'   \item{small}{List of fitted values for window size `w`.}
#'   \item{large}{List of fitted values for window size `2w`.}
#' }
#'
#' @details
#' Each rolling window is standardized before fitting. The returned matrices
#' contain three columns:
#' \itemize{
#'   \item scaled input `x`
#'   \item fitted values `y_hat`
#'   \item original indices
#' }
#'
#' @import RSNNS
#' @import foreach
#' @import doSNOW
#' @import parallel
#' @export
fit_mlp <- function(
    vec,
    w = 100,
    act1 = "Act_Logistic",
    act2 = "Act_Logistic",
    lr1 = 0.001, lr2 = 0.001,
    hl1 = 6, hl2 = 8,
    ep1 = 1000, ep2 = 2000
) {

  # ------------------------------------------------------------
  # Setup
  # ------------------------------------------------------------
  n.val <- length(vec)
  x <- 1:n.val
  y <- vec

  # Number of cores
  num_cores <- max(1, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(num_cores)
  doSNOW::registerDoSNOW(cl)

  # ------------------------------------------------------------
  # Progress bar for small-window MLP
  # ------------------------------------------------------------
  pb <- txtProgressBar(min = 0, max = n.val - w, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  message("Fitting small-window MLPs...")

  # ------------------------------------------------------------
  # Small-window rolling MLP
  # ------------------------------------------------------------
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
  # Progress bar for large-window MLP
  # ------------------------------------------------------------
  pb <- txtProgressBar(min = 0, max = n.val - 2 * w, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  message("Fitting large-window MLPs...")

  # ------------------------------------------------------------
  # Large-window rolling MLP
  # ------------------------------------------------------------
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

  # ------------------------------------------------------------
  # Cleanup
  # ------------------------------------------------------------
  parallel::stopCluster(cl)

  # ------------------------------------------------------------
  # Return
  # ------------------------------------------------------------
  list(
    small = res.small,
    large = res.large
  )
}
