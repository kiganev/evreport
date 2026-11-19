#' Breusch-Pagan-Godfrey or White heteroskedasticity test
#' @param model An unweighted, full-rank lm object with an intercept.
#' @param type Either "bpg" or "white".
#' @param cross For White: FALSE uses squared slopes only; TRUE adds levels and pairwise products.
#' @return An evreport result containing three test statistics, their probabilities,
#'   auxiliary regression coefficients and an auxiliary lm model.
#' @details The auxiliary dependent variable is the original squared OLS residual.
#'   EViews displays F, Obs*R-squared and Scaled explained SS. The last statistic
#'   uses the explained sum of squares of the auxiliary regression divided by
#'   twice the squared original regression's estimated error variance.
#'   White with cross = FALSE uses squared regressors only; cross = TRUE uses
#'   levels, squares and interactions as in the supplied EViews reference reports.
#' @export
#' @examples
#' ev_hetero(lm(mpg ~ wt + hp, data = mtcars), type = "white")
ev_hetero <- function(model, type = c("bpg", "white"), cross = FALSE) {
  X <- .ev_check(model)
  type <- match.arg(type)
  if (attr(terms(model), "intercept") != 1L)
    stop("Heteroskedasticity tests currently require an intercept.", call. = FALSE)
  if (!is.logical(cross) || length(cross) != 1L || is.na(cross))
    stop("cross must be TRUE or FALSE.", call. = FALSE)

  e <- as.numeric(residuals(model))
  y <- e^2
  n <- length(y)
  base <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  if (ncol(base) < 1L) stop("At least one slope regressor is required.", call. = FALSE)

  # EViews reference cases:
  # Without cross terms: constant, each original regressor squared.
  # With cross terms: constant, squares, pairwise products and levels.
  # Arrange full-test columns in the same order as the reference EViews report.
  if (type == "white") {
    terms <- list(`(Intercept)` = rep(1, n))
    for (i in seq_len(ncol(base))) {
      xi <- base[, i]
      name_i <- colnames(base)[i]
      terms[[paste0(name_i, "^2")]] <- xi^2
      if (cross) {
        if (i < ncol(base)) {
          for (j in seq.int(i + 1L, ncol(base))) {
            name_j <- colnames(base)[j]
            terms[[paste0(name_i, "*", name_j)]] <- xi * base[, j]
          }
        }
        terms[[name_i]] <- xi
      }
    }
    A <- do.call(cbind, terms)
  } else {
    A <- cbind(`(Intercept)` = 1, base)
  }

  # Drop collinear regressors, retaining intercept and remaining full-rank columns.
  # Pivoting is deterministic but exact retained terms may depend on QR tolerance.
  qrA <- qr(A)
  kept <- sort(qrA$pivot[seq_len(qrA$rank)])
  A <- A[, kept, drop = FALSE]
  if (ncol(A) < 2L) stop("No identifiable auxiliary slopes.", call. = FALSE)
  if (colnames(A)[1L] != "(Intercept)")
    stop("Auxiliary intercept is not identified.", call. = FALSE)

  fit <- lm(y ~ A - 1)
  sa <- summary(fit)
  tab <- sa$coefficients
  rownames(tab) <- colnames(A)
  rss <- sum(residuals(fit)^2)
  tss <- sum((y - mean(y))^2)
  ess <- max(0, tss - rss)
  if (tss <= 0 || rss <= 0) stop("Auxiliary response has insufficient variation.", call. = FALSE)

  q <- ncol(A) - 1L
  df2 <- n - ncol(A)
  if (df2 <= 0) stop("Not enough residual degrees of freedom.", call. = FALSE)
  F <- (ess / q) / (rss / df2)
  r2 <- max(0, 1 - rss / tss)
  lmstat <- n * r2

  # EViews' Scaled explained SS uses sigma^2 = original SSR / (n - k).
  original_sigma2 <- sum(e^2) / (n - model$rank)
  scaled <- ess / (2 * original_sigma2^2)

  adjr2 <- 1 - (rss / df2) / (tss / (n - 1L))
  eqstats <- .ev_equation_stats(y, as.numeric(residuals(fit)),
             ncol(A), fstat = F,
             fpval = pf(F, q, df2, lower.tail = FALSE),
             r2 = r2, adjr2 = adjr2)
  out <- .ev_result(
    if (type == "bpg") "Heteroskedasticity Test: Breusch-Pagan-Godfrey" else
      "Heteroskedasticity Test: White",
    c(`F-statistic` = F, `Prob. F` = pf(F, q, df2, lower.tail = FALSE),
      `Obs*R-squared` = lmstat,
      `Prob. Chi-Square` = pchisq(lmstat, q, lower.tail = FALSE),
      `Scaled explained SS` = scaled,
      `Prob. Scaled explained SS` = pchisq(scaled, q, lower.tail = FALSE)),
    tab,
    list(`Included observations` = n,
         `Test degrees of freedom` = q,
         `Residual degrees of freedom` = df2,
         `Cross terms included` = identical(type, "white") && cross,
         `Equation statistics` = eqstats),
    fit
  )
  print(out)
  invisible(out)
}
