#' Jarque-Bera test for OLS residual normality
#'
#' Displays the familiar EViews-style histogram with the residual descriptive
#' statistics and Jarque-Bera test to its right. All computations use the same
#' residual vector; choosing a different number of bins changes only the plot.
#'
#' @param model An unweighted, full-rank lm object.
#' @param plot Logical; draw the histogram and statistics panel (default TRUE).
#' @param breaks Number of bins, break points, or a method understood by
#'   [graphics::hist()]. Defaults to `"Sturges"`. EViews may choose different
#'   breaks, so bar counts need not match exactly.
#' @return An evreport object containing sample moments and Jarque-Bera test,
#'   returned invisibly.
#' @export
#' @examples
#' fit <- lm(mpg ~ wt + hp, mtcars)
#' ev_normality(fit)
#' ev_normality(fit, plot = FALSE)
#' ev_normality(fit, breaks = 8)
ev_normality <- function(model, plot = TRUE, breaks = "Sturges") {
  .ev_check(model)
  if (!is.logical(plot) || length(plot) != 1L || is.na(plot))
    stop("'plot' must be TRUE or FALSE.", call. = FALSE)
  x <- as.numeric(residuals(model))
  n <- length(x)
  v <- mean((x - mean(x))^2)
  if (v <= 0) stop("Residual variance is zero.", call. = FALSE)
  s <- mean((x - mean(x))^3) / v^(3/2)
  kurt <- mean((x - mean(x))^4) / v^2
  jb <- n / 6 * (s^2 + (kurt - 3)^2 / 4)
  out <- .ev_result("Residual Normality Test: Jarque-Bera",
    c(`Mean` = mean(x), `Median` = median(x), `Maximum` = max(x),
      `Minimum` = min(x), `Std. Dev.` = sd(x), `Skewness` = s,
      `Kurtosis` = kurt, `Jarque-Bera` = jb,
      `Probability` = pchisq(jb, 2L, lower.tail = FALSE)),
    details = list(`Null hypothesis` = "Normally distributed residuals",
                   `Observations` = n))
  print(out)
  if (plot) .ev_normality_plot(x, out$statistics, breaks = breaks)
  invisible(out)
}

# Base R plot: histogram on the left, EViews-inspired statistics on the right.
.ev_normality_plot <- function(x, stats, breaks) {
  # Compute breaks before changing the graphics layout (so bad breaks fail early).
  h <- graphics::hist(x, breaks = breaks, plot = FALSE)
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::layout(1)
    graphics::par(oldpar)
  }, add = TRUE)
  graphics::layout(matrix(c(1L, 2L), nrow = 1L), widths = c(3, 1.7))

  graphics::par(mar = c(4.5, 4.2, 2.5, 0.7))
  graphics::plot(h, freq = TRUE, col = "#487EB0", border = "white",
                 main = "Histogram of Residuals", xlab = "Residuals", ylab = "Frequency")

  graphics::par(mar = c(1.5, 0.4, 1.5, 0.5))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  graphics::text(0.03, 0.96, "Series: Residuals", adj = c(0, 1), font = 2, cex = 0.9)
  graphics::text(0.03, 0.88, sprintf("Observations: %d", length(x)),
                 adj = c(0, 1), cex = 0.86)
  nm <- c("Mean", "Median", "Maximum", "Minimum", "Std. Dev.",
          "Skewness", "Kurtosis")
  ys <- seq(0.75, 0.35, length.out = length(nm))
  fmt <- function(z) formatC(z, digits = 6, format = "fg", flag = "#")
  for (i in seq_along(nm)) {
    graphics::text(0.03, ys[i], nm[i], adj = c(0, 0.5), cex = 0.81)
    graphics::text(0.97, ys[i], fmt(stats[[nm[i]]]), adj = c(1, 0.5), cex = 0.81)
  }
  graphics::rect(0.01, 0.09, 0.99, 0.29, border = "#777777")
  graphics::text(0.04, 0.235, "Jarque-Bera", adj = c(0, 0.5), cex = 0.82)
  graphics::text(0.96, 0.235, fmt(stats[["Jarque-Bera"]]), adj = c(1, 0.5), cex = 0.82)
  graphics::text(0.04, 0.145, "Probability", adj = c(0, 0.5), cex = 0.82)
  graphics::text(0.96, 0.145, fmt(stats[["Probability"]]), adj = c(1, 0.5), cex = 0.82)
  invisible(h)
}
