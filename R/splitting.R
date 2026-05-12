#' Best Split for Binary Labels (Free Assignment)
#'
#' @description
#' Given a sequence of labels in \{1, 2\}, this function finds the best split
#' point that maximizes the number of correctly assigned labels on each side.
#' The left and right segments are each assigned their majority class.
#'
#' @param labels Integer vector containing only 1 and 2.
#'
#' @return A list with:
#' \describe{
#'   \item{index}{The optimal split location (between index and index + 1).}
#'   \item{total_correct}{Maximum number of correctly assigned labels.}
#'   \item{left_majority}{Majority class on the left segment.}
#'   \item{right_majority}{Majority class on the right segment.}
#'   \item{accuracy_if_assigning_majorities}{Proportion correctly assigned.}
#' }
#'
#' @details
#' For each possible split point \(i\), the function computes:
#' \itemize{
#'   \item left correct = max(#1_left, #2_left)
#'   \item right correct = max(#1_right, #2_right)
#' }
#' and selects the split maximizing their sum.
#'
#' @examples
#' labels <- c(1,1,1,2,2,2)
#' best_split_free(labels)
#'
#' @export
best_split_free <- function(labels) {

  # Validate
  if (length(labels) < 2)
    stop("Need at least 2 labels to split.")

  if (any(is.na(labels)))
    stop("Labels contain NA; please remove or impute.")

  u <- sort(unique(labels))
  if (!all(u %in% c(1, 2)))
    stop("Labels must be in {1, 2}.")

  n <- length(labels)

  # Cumulative counts
  cum1 <- cumsum(labels == 1)
  cum2 <- cumsum(labels == 2)

  total1 <- cum1[n]
  total2 <- cum2[n]

  # Left and right correct counts
  left_correct  <- pmax(cum1, cum2)
  right_correct <- pmax(total1 - cum1, total2 - cum2)

  # Only splits 1..(n-1)
  total_correct <- left_correct[1:(n - 1)] + right_correct[1:(n - 1)]

  # Best split
  best_i <- which.max(total_correct)

  list(
    index = best_i,
    total_correct = total_correct[best_i],
    left_majority  = if (cum1[best_i] >= cum2[best_i]) 1 else 2,
    right_majority = if ((total1 - cum1[best_i]) >= (total2 - cum2[best_i])) 1 else 2,
    accuracy_if_assigning_majorities = total_correct[best_i] / n
  )
}


#' Select the Most Significant Spike in an ECDF Spacing Curve
#'
#' @description
#' Identifies the most prominent spike in a spacing curve derived from the
#' ECDF of a smoothed detector statistic. This is used to automatically
#' determine a significance threshold for changepoint selection.
#'
#' @param s Numeric vector. Smoothed spacing curve.
#' @param right_tail_cutoff Numeric in (0,1). Exclude spikes with ECDF
#'   probability above this value (typically near 1).
#' @param left_tail_cutoff Numeric in (0,1). Exclude spikes with ECDF
#'   probability below this value (to avoid trivial early spikes).
#'
#' @return A numeric value in (0,1) representing the selected ECDF threshold,
#'   or `NA` if no valid spike is found.
#'
#' @details
#' Uses `pracma::findpeaks()` to identify local maxima. Prominence is computed
#' as:
#' \deqn{ peak\_height - max(left\_base, right\_base) }
#'
#' The spike with the largest prominence within the allowed ECDF range is
#' selected.
#'
#' @importFrom pracma findpeaks
#'
#' @examples
#' s <- runif(100)
#' select_best_spike(s)
#'
#' @export
select_best_spike <- function(
    s,
    right_tail_cutoff = 0.95,
    left_tail_cutoff  = 0.6
) {

  n <- length(s)
  p <- (1:n) / (n + 1)   # ECDF probabilities

  # Find peaks
  peaks <- pracma::findpeaks(s, sortstr = FALSE)
  if (is.null(peaks))
    return(NA)

  peak_heights <- peaks[, 1]
  peak_idx     <- peaks[, 2]
  left_base    <- peaks[, 3]
  right_base   <- peaks[, 4]

  # Prominence
  baseline <- pmax(s[left_base], s[right_base])
  prominence <- peak_heights - baseline

  # Exclude tail regions
  valid <- (p[peak_idx] < right_tail_cutoff) &
    (p[peak_idx] > left_tail_cutoff)

  if (!any(valid))
    return(NA)

  # Select spike with largest prominence
  best_idx <- peak_idx[valid][which.max(prominence[valid])]
  best_p   <- p[best_idx]

  best_p
}
