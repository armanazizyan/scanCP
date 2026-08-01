#' Plot Signal Decomposition
#'
#' @description
#' Creates a single-column visualization of signal decomposition with four panels:
#' 1. Original signal with detected changepoints
#' 2. Piecewise constant component (estimated level shifts)
#' 3. Smooth component (trend curve from global MLP)
#' 4. Residuals from MLP fit
#'
#' @param y Numeric vector. Original signal.
#' @param decomp List. Output from \code{\link{decompose_signal_core}}.
#' @param mlp_fit List. Output from \code{\link{fit_global_mlp}}. If NULL,
#'   the smooth curve and residuals panels will be skipped.
#' @param title Character. Main title for the plot. Defaults to
#'   \code{"Signal Decomposition"}.
#'
#' @return Invisibly returns a ggplot object (or list of ggplot objects if
#'   using gridExtra).
#'
#' @details
#' The four-panel decomposition shows:
#' \itemize{
#'   \item Panel 1 (Original Signal): Raw data with vertical dashed lines
#'         marking detected changepoints.
#'   \item Panel 2 (Piecewise Constant): Estimated step function capturing
#'         level shifts at changepoints.
#'   \item Panel 3 (Smooth Curve): Global MLP trend fit to the corrected signal.
#'   \item Panel 4 (Residuals): Residuals from the MLP fit (should resemble
#'         white noise if decomposition is successful).
#' }
#'
#' Changepoint locations are marked with red vertical dashed lines across all
#' panels for easy reference.
#'
#' @examples
#' \donttest{
#'   set.seed(123)
#'   y <- c(rnorm(100, 0), rnorm(100, 3), rnorm(100, 1))
#'   scp.res <- scan_cp(y, w = 20, threshold = "auc")
#'   decomp <- scp.res$decomposition
#'   mlp_fit <- fit_global_mlp(decomp$corrected_signal)
#'   plot_decomposition(y, decomp, mlp_fit)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_vline theme_minimal
#'   theme element_text labs element_blank
#' @importFrom gridExtra grid.arrange
#' @export
plot_decomposition <- function(
    y,
    decomp,
    mlp_fit = NULL,
    title = "Signal Decomposition"
) {

  n <- length(y)
  x <- 1:n
  changepoints <- decomp$changepoints
  piecewise_const <- decomp$piecewise_constant

  # --- 1. Original Signal Panel ---

  df_original <- data.frame(x = x, y = y)
  p1 <- ggplot(df_original, aes(x = x, y = y)) +
    geom_line(color = "#2C3E50", linewidth = 0.4, alpha = 0.8) +
    {if (length(changepoints) > 0) {
      geom_vline(
        xintercept = changepoints,
        linetype = "dashed",
        color = "#E74C3C",
        linewidth = 0.7,
        alpha = 0.8
      )
    }} +
    labs(
      x = NULL,
      y = "Original Signal",
      title = "1. Original Signal"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

  # --- 2. Piecewise Constant Panel ---

  df_pwc <- data.frame(x = x, y = piecewise_const)
  p2 <- ggplot(df_pwc, aes(x = x, y = y)) +
    geom_line(color = "#F39C12", linewidth = 0.5, alpha = 0.8) +
    {if (length(changepoints) > 0) {
      geom_vline(
        xintercept = changepoints,
        linetype = "dashed",
        color = "#E74C3C",
        linewidth = 0.7,
        alpha = 0.8
      )
    }} +
    labs(
      x = NULL,
      y = "Piecewise Constant",
      title = "2. Piecewise Constant Component"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

  # --- 3. Smooth Curve Panel (if mlp_fit provided) ---

  if (!is.null(mlp_fit) && !is.null(mlp_fit$smooth_curve)) {
    smooth_curve <- mlp_fit$smooth_curve
    df_smooth <- data.frame(x = x, y = smooth_curve)
    p3 <- ggplot(df_smooth, aes(x = x, y = y)) +
      geom_line(color = "#27AE60", linewidth = 0.5, alpha = 0.8) +
      {if (length(changepoints) > 0) {
        geom_vline(
          xintercept = changepoints,
          linetype = "dashed",
          color = "#E74C3C",
          linewidth = 0.7,
          alpha = 0.8
        )
      }} +
      labs(
        x = NULL,
        y = "Smooth Curve",
        title = "3. Global MLP Smooth Component"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 11, face = "bold", hjust = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
  } else {
    p3 <- NULL
  }

  # --- 4. Residuals Panel (if mlp_fit provided) ---

  if (!is.null(mlp_fit) && !is.null(mlp_fit$residual)) {
    residuals <- mlp_fit$residual
    df_resid <- data.frame(x = x, y = residuals)
    p4 <- ggplot(df_resid, aes(x = x, y = y)) +
      geom_line(color = "#9B59B6", linewidth = 0.4, alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "solid", color = "gray50", alpha = 0.5) +
      {if (length(changepoints) > 0) {
        geom_vline(
          xintercept = changepoints,
          linetype = "dashed",
          color = "#E74C3C",
          linewidth = 0.7,
          alpha = 0.6
        )
      }} +
      labs(
        x = "Index",
        y = "Residuals",
        title = "4. MLP Residuals"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 11, face = "bold", hjust = 0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
  } else {
    p4 <- NULL
  }

  # --- Combine all panels ---

  plot_list <- list(p1, p2, p3, p4)
  plot_list <- plot_list[!sapply(plot_list, is.null)]

  if (length(plot_list) == 4) {
    combined <- gridExtra::grid.arrange(
      plot_list[[1]],
      plot_list[[2]],
      plot_list[[3]],
      plot_list[[4]],
      ncol = 1,
      heights = c(1, 1, 1, 1),
      top = title
    )
  } else if (length(plot_list) == 2) {
    combined <- gridExtra::grid.arrange(
      plot_list[[1]],
      plot_list[[2]],
      ncol = 1,
      top = title
    )
  } else {
    combined <- plot_list[[1]]
  }

  return(invisible(combined))
}
