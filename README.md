# scanCP: Multilayer Perceptron–Based Changepoint Detection

`scanCP` provides a fast, robust, and fully automated changepoint detection
pipeline based on rolling multilayer perceptron (MLP) models. The method
combines local and global structure learning, RSS-based detector statistics,
ECDF thresholding, and local refinement to identify structural breaks in
univariate signals.

The package is designed for signals with abrupt level shifts, smooth trends,
and heterogeneous noise. It includes tools for detector construction,
changepoint extraction, signal correction, and visualization.

---

## Features

- Rolling-window MLP fits (window sizes `w` and `2w`)
- RSS-based hybrid detector (ratio + difference)
- Automatic threshold selection via spacing-curve analysis
- Local refinement using k-means segmentation
- Piecewise-constant correction of the signal
- Optional global MLP smoother
- Parallel computation for speed
- Clean, modular API

---

## Installation

### From CRAN (when available)

```r
install.packages("scanCP")
```

### From GitHub (development version)

```r
# install.packages("devtools")
devtools::install_github("armanazizyan/scanCP")
```

---

## Quick Start

```r
library(scanCP)

# Simulate a simple piecewise-constant signal
set.seed(123)
y <- c(
  rnorm(200, 0, 1),
  rnorm(200, 3, 1),
  rnorm(200, -2, 1),
  rnorm(200, 1, 1)
)

# Run the full changepoint pipeline
res <- scan_cp(y, w = 100)

# Extract changepoints
res$changepoints

# Plot the corrected signal
plot(y, type = "l", col = "gray70")
lines(res$corrected_signal, col = "blue", lwd = 2)
abline(v = res$changepoints, col = "red", lwd = 2)
```

---

## Method Overview

`scanCP` detects changepoints using a hybrid detector statistic derived from
rolling MLP fits:

- A small-window MLP (size `w`) models local structure.
- A large-window MLP (size `2w`) models broader structure.
- For each position, the detector combines:
  - a ratio of RSS values
  - a difference of RSS values

The detector is smoothed, thresholded using the ECDF or spacing-curve
analysis, and refined using a local k-means split. The result is a set of
precise changepoint locations and a corrected signal.

---

## Customizing the MLP

All MLP hyperparameters can be controlled via `mlp_control`:

```r
res <- scan_cp(
  y,
  w = 100,
  mlp_control = list(
    hl1 = 10,
    hl2 = 20,
    ep1 = 1500,
    act1 = "Act_Tanh"
  )
)
```

---

## Documentation

Full function documentation is available via:

```r
help(package = "scanCP")
```

---

## License

GPL-2 License © Arman Azizyan

---

## Contributing

Issues and pull requests are welcome at:

https://github.com/armanazizyan/scanCP
