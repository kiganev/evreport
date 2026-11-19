.ev_check <- function(model) {
  if (!inherits(model, "lm")) stop("Expected an lm object.", call. = FALSE)
  if (!is.null(model$weights)) stop("Weighted lm models are not supported.", call. = FALSE)
  if (!is.null(model$offset)) stop("Offset models are not supported.", call. = FALSE)
  X <- model.matrix(model)
  if (model$rank != ncol(X)) stop("Rank-deficient original model is not supported.", call. = FALSE)
  if (any(!is.finite(residuals(model)))) stop("Non-finite residuals.", call. = FALSE)
  if (nrow(X) <= ncol(X) + 1L) stop("Too few observations.", call. = FALSE)
  invisible(X)
}
.ev_result <- function(title, stats, table = NULL, details = list(), auxiliary_model = NULL) {
  structure(list(title = title, statistics = stats, coefficients = table,
                 details = details, auxiliary_model = auxiliary_model), class = "evreport")
}
.ev_joint <- function(y, X, restricted_cols = 1L, uncentered = FALSE) {
  n <- length(y)
  full <- lm.fit(X, y)
  rank <- full$rank
  if (rank != ncol(X)) stop("Auxiliary regression is rank deficient; simplify specification.", call. = FALSE)
  rss <- sum(full$residuals^2)
  if (rss <= 0) stop("Auxiliary residual sum of squares is zero.", call. = FALSE)
  df1 <- ncol(X) - restricted_cols
  df2 <- n - rank
  if (df1 < 1L || df2 < 1L) stop("Not enough degrees of freedom.", call. = FALSE)
  restricted <- lm.fit(X[, seq_len(restricted_cols), drop = FALSE], y)
  rss0 <- sum(restricted$residuals^2)
  F <- max(0, (rss0 - rss) / df1) / (rss / df2)
  tss <- if (uncentered) sum(y^2) else sum((y - mean(y))^2)
  if (tss <= 0) stop("Auxiliary dependent variable has no variation.", call. = FALSE)
  r2 <- max(0, 1 - rss / tss)
  list(F = F, pF = pf(F, df1, df2, lower.tail = FALSE),
       LM = n * r2, pLM = pchisq(n * r2, df1, lower.tail = FALSE),
       df1 = df1, df2 = df2, r2 = r2)
}
.ev_aux_table <- function(fit) summary(fit)$coefficients

# Return the conventional EViews left/right equation summary layout.
.ev_equation_stats <- function(y, e, k, fstat = NA_real_, fpval = NA_real_,
                               r2 = NULL, adjr2 = NULL) {
  n <- length(e)
  rss <- sum(e^2)
  if (is.null(r2)) r2 <- 1 - rss / sum((y - mean(y))^2)
  if (is.null(adjr2)) adjr2 <- 1 - (1 - r2) * (n - 1) / (n - k)
  ll <- -n / 2 * (log(2 * pi) + 1 + log(rss / n))
  c(`R-squared` = r2, `Adjusted R-squared` = adjr2,
    `S.E. of regression` = sqrt(rss / (n - k)),
    `Sum squared resid` = rss, `Log likelihood` = ll,
    `F-statistic` = unname(fstat), `Prob(F-statistic)` = unname(fpval),
    `Mean dependent var` = mean(y), `S.D. dependent var` = sd(y),
    `Akaike info criterion` = (-2 * ll + 2 * k) / n,
    `Schwarz criterion` = (-2 * ll + log(n) * k) / n,
    `Hannan-Quinn criter.` = (-2 * ll + 2 * k * log(log(n))) / n,
    `Durbin-Watson stat` = sum(diff(e)^2) / rss)
}

.ev_print_coefficients <- function(tab) {
  cat(sprintf("%-20s %13s %13s %13s %13s\n",
              "Variable", "Coefficient", "Std. Error", "t-Statistic", "Prob."))
  cat(strrep("-", 79), "\n")
  for (i in seq_len(nrow(tab))) {
    nm <- rownames(tab)[i]
    if (nm == "(Intercept)") nm <- "C"
    cat(sprintf("%-20s %13.6f %13.6f %13.6f %13.4f\n",
                nm, tab[i, 1], tab[i, 2], tab[i, 3], tab[i, 4]))
  }
  cat(strrep("-", 79), "\n")
}

.ev_print_equation_stats <- function(s) {
  left <- c("R-squared", "Adjusted R-squared", "S.E. of regression",
            "Sum squared resid", "Log likelihood", "F-statistic", "Prob(F-statistic)")
  right <- c("Mean dependent var", "S.D. dependent var", "Akaike info criterion",
             "Schwarz criterion", "Hannan-Quinn criter.", "Durbin-Watson stat", "")
  fmt <- function(x) if (length(x) == 0L || is.na(x)) "NA" else sprintf("%.6f", x)
  for (i in seq_along(left)) {
    cat(sprintf("%-23s %13s  %-23s %13s\n",
                left[i], fmt(s[[left[i]]]),
                right[i], if (nzchar(right[i])) fmt(s[[right[i]]]) else ""))
  }
}

#' Print EViews-style results
#' @param x An evreport result.
#' @param ... Not used.
#' @export
print.evreport <- function(x, ...) {
  is_equation <- identical(x$title, "Least Squares")
  is_bg <- startsWith(x$title, "Breusch-Godfrey")
  cat("\n", x$title, if (is_bg) ":" else "", "\n", sep = "")
  cat(strrep("=", 79), "\n")

  if (is_equation) {
    cat("Dependent Variable:", x$details[["Dependent Variable"]], "\n")
    cat("Method: Least Squares\n")
    cat("Included observations:", x$details[["Included observations"]], "\n\n")
    .ev_print_coefficients(x$coefficients)
    .ev_print_equation_stats(x$statistics)
  } else if (is_bg) {
    p <- as.integer(x$details[["Lags"]])
    df2 <- as.integer(x$details[["Residual degrees of freedom"]])
    cat("Null hypothesis: No serial correlation at up to", p, "lag(s)\n\n")
    cat(sprintf("%-20s %13.6f   Prob. F(%d,%d)         %.6f\n",
                "F-statistic", unname(x$statistics[["F-statistic"]]), p, df2,
                unname(x$statistics[["Prob. F"]])))
    cat(sprintf("%-20s %13.6f   Prob. Chi-Square(%d)  %.6f\n",
                "Obs*R-squared", unname(x$statistics[["Obs*R-squared"]]), p,
                unname(x$statistics[["Prob. Chi-Square"]])))
    cat("\n", strrep("=", 79), "\n", sep = "")
    cat("Test Equation:\nDependent Variable: RESID\nMethod: Least Squares\n")
    cat("Presample missing value lagged residuals set to zero.\n")
    cat("Included observations:", x$details[["Included observations"]], "\n\n")
    .ev_print_coefficients(x$coefficients)
    if (!is.null(x$details[["Equation statistics"]]))
      .ev_print_equation_stats(x$details[["Equation statistics"]])
  } else {
    for (nm in names(x$details)) {
      if (!is.list(x$details[[nm]])) cat(nm, ": ", x$details[[nm]], "\n", sep = "")
    }
    cat("\n")
    for (nm in names(x$statistics))
      cat(sprintf("%-32s %16.6f\n", nm, unname(x$statistics[[nm]])))
    if (!is.null(x$coefficients)) {
      cat("\nTest Equation:\n")
      .ev_print_coefficients(x$coefficients)
    }
  }
  invisible(x)
}
