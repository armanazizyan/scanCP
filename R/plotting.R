#' Interactive Plot of Rolling MLP Fits
#'
#' @description
#' Creates an interactive Plotly visualization showing multiple rolling-window
#' MLP fits over a univariate signal. A slider allows the user to move through
#' different window positions and compare model fits dynamically.
#'
#' @param y Numeric vector. The original signal.
#' @param fit_mlp_res A list containing rolling MLP fit results. Expected to be
#'   a list of two lists, each containing matrices with fitted values.
#' @param t.chp.ind Optional numeric vector of changepoint indices to draw as
#'   vertical dashed lines.
#' @param step Integer. Step size for the slider increments.
#' @param start Integer. Starting index for the first model window.
#' @param w Integer. Window size used in the rolling MLP fitting.
#'
#' @return A Plotly object.
#'
#' @importFrom plotly plot_ly add_lines layout
#' @export
plot_mlp_fits_interactive <- function(
    y,
    fit_mlp_res,
    t.chp.ind = NA,
    step = 50,
    start = 75,
    w = 100
) {

  # ------------------------------------------------------------
  # Basic setup
  # ------------------------------------------------------------
  n <- length(y)

  # fit_mlp_res structure:
  #   fit_mlp_res[[1]] = model 2 + model 3
  #   fit_mlp_res[[2]] = model 1
  res.list1 <- fit_mlp_res[[2]]               # Model 1 fits
  res.list2 <- fit_mlp_res[[1]]               # Model 2 fits
  res.list3 <- fit_mlp_res[[1]][w + 1:n]      # Model 3 fits (shifted)

  # Indices where models were fitted
  idx_list <- seq(start, n - 2 * w, by = step)
  n_models <- length(idx_list)

  # ------------------------------------------------------------
  # Changepoint vertical lines (optional)
  # ------------------------------------------------------------
  if (!is.na(t.chp.ind)) {
    cp_segments <- lapply(t.chp.ind, function(cp) {
      list(
        type = "line",
        x0 = cp, x1 = cp,
        y0 = min(y), y1 = max(y),
        line = list(color = "blue", dash = "dash")
      )
    })
  } else {
    cp_segments <- NULL
  }

  # ------------------------------------------------------------
  # Initial plot (frame 1)
  # ------------------------------------------------------------
  i0 <- idx_list[1]

  p <- plotly::plot_ly() %>%
    # Original signal
    plotly::add_lines(
      x = 1:n, y = y,
      name = "Original Signal",
      line = list(color = "rgba(0,0,0,0.3)")
    ) %>%

    # Model 1
    plotly::add_lines(
      x = res.list1[[i0]][, 3],
      y = res.list1[[i0]][, 2],
      name = "Model 1",
      line = list(color = "red", width = 3)
    ) %>%

    # Model 2
    plotly::add_lines(
      x = res.list2[[i0]][, 3],
      y = res.list2[[i0]][, 2],
      name = "Model 2",
      line = list(color = "green", width = 3)
    ) %>%

    # Model 3
    plotly::add_lines(
      x = res.list3[[i0]][, 3],
      y = res.list3[[i0]][, 2],
      name = "Model 3",
      line = list(color = "purple", width = 3)
    ) %>%

    # Layout with slider placeholder
    plotly::layout(
      title = "Multiple MLP Fits – Interactive Slider",
      xaxis = list(title = "Index"),
      yaxis = list(title = "Value"),
      shapes = cp_segments,
      sliders = list(list(
        active = 0,
        currentvalue = list(prefix = "Model index: "),
        steps = lapply(seq_len(n_models), function(k) {
          list(
            label = as.character(idx_list[k]),
            method = "animate",
            args = list(
              list(paste0("frame", k)),
              list(
                mode = "immediate",
                frame = list(duration = 0, redraw = FALSE),
                transition = list(duration = 0)
              )
            )
          )
        })
      ))
    )

  # ------------------------------------------------------------
  # Animation frames
  # ------------------------------------------------------------
  frames <- lapply(seq_len(n_models), function(k) {
    i <- idx_list[k]

    list(
      name = paste0("frame", k),
      data = list(
        NULL,  # placeholder for original signal (unchanged)
        list(  # Model 1
          x = res.list1[[i]][, 3],
          y = res.list1[[i]][, 2],
          mode = "lines",
          line = list(color = "red", width = 3)
        ),
        list(  # Model 2
          x = res.list2[[i]][, 3],
          y = res.list2[[i]][, 2],
          mode = "lines",
          line = list(color = "green", width = 3)
        ),
        list(  # Model 3
          x = res.list3[[i]][, 3],
          y = res.list3[[i]][, 2],
          mode = "lines",
          line = list(color = "purple", width = 3)
        )
      )
    )
  })

  # Attach frames to the plotly object
  p$x$frames <- frames

  p
}
