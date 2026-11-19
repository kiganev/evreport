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
.ev_aux_table <- function(fit) {
  summary(fit)$coefficients
}
#' Print EViews-inspired results
#' @param x An evreport result.
#' @param ... Not used.
#' @export
print.evreport <- function(x, ...) {
  cat("\n", x$title, "\n", strrep("=", 68), "\n", sep = "")
  if (length(x$details)) for (nm in names(x$details)) cat(nm, ": ", x$details[[nm]], "\n", sep = "")
  if (length(x$statistics)) {
    for (nm in names(x$statistics)) {
      val <- x$statistics[[nm]]
      cat(sprintf("%-32s %16.6f\n", nm, val))
    }
  }
  if (!is.null(x$coefficients)) {
    cat("\nTest Equation:\n")
    tab <- x$coefficients
    cat(sprintf("%-19s %12s %12s %12s %12s\n", "Variable", "Coefficient", "Std. Error", "t-Statistic", "Prob."))
    cat(strrep("-", 75), "\n")
    for (i in seq_len(nrow(tab))) cat(sprintf("%-19s %12.6f %12.6f %12.6f %12.4f\n",
       if (rownames(tab)[i] == "(Intercept)") "C" else rownames(tab)[i], tab[i,1], tab[i,2], tab[i,3], tab[i,4]))
  }
  invisible(x)
}
