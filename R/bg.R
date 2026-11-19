#' Breusch-Godfrey serial correlation LM test
#' @description Uses EViews' presample-zero convention and the uncentered
#'   auxiliary R-squared for Obs*R-squared. Observations must be consecutive
#'   and correctly ordered; missing periods are not inferred.
#' @param model An unweighted, full-rank lm object.
#' @param lags Positive integer number of lagged residuals.
#' @return An evreport object containing F/LM statistics and the auxiliary lm.
#' @export
#' @examples
#' ev_bg(lm(mpg ~ wt + hp, mtcars), lags = 2)
ev_bg <- function(model, lags = 1L) {
  X <- .ev_check(model)
  if (length(lags) != 1L || !is.numeric(lags) || is.na(lags) ||
      !is.finite(lags) || lags < 1 || lags != floor(lags))
    stop("lags must be a positive integer.", call. = FALSE)
  p <- as.integer(lags)
  e <- as.numeric(residuals(model))
  n <- length(e)
  if (p >= n || n <= ncol(X) + p) stop("Too many lags.", call. = FALSE)
  L <- matrix(0, n, p)
  for (j in seq_len(p)) L[(j+1):n, j] <- e[1:(n-j)]
  colnames(L) <- paste0("RESID(-", seq_len(p), ")")
  Z <- cbind(X, L)
  aux <- lm.fit(Z, e)
  if (aux$rank != ncol(Z)) stop("Auxiliary regression is rank deficient.", call. = FALSE)
  rss <- sum(aux$residuals^2)
  rss0 <- sum(lm.fit(X, e)$residuals^2)
  df2 <- n - ncol(Z)
  if (rss <= 0) stop("Zero auxiliary residual variance.", call. = FALSE)
  F <- max(0, (rss0-rss)/p) / (rss/df2)
  LM <- n * max(0, 1-rss/sum(e^2))
  fit <- lm(e ~ Z - 1)
  # use readable display names for test equation
  tab <- summary(fit)$coefficients
  rownames(tab) <- colnames(Z)
  # Z already contains the original regression's intercept, if present.
  # Fitted with `- 1` to retain the exact names of regressors, so
  # summary(fit)$fstatistic refers to a through-origin model and cannot
  # be used for the EViews test-equation overall F-test.
  # Compute the auxiliary equation's ordinary overall F explicitly.
  has_intercept <- attr(terms(model), "intercept") == 1L
  tss <- if (has_intercept) sum((e - mean(e))^2) else sum(e^2)
  r2_aux <- 1 - rss / tss
  df_model <- ncol(Z) - as.integer(has_intercept)
  df_total <- n - as.integer(has_intercept)
  r2_aux_adj <- 1 - (rss / df2) / (tss / df_total)
  aux_F <- if (df_model > 0L) {
    max(0, (tss - rss) / df_model) / (rss / df2)
  } else NA_real_
  aux_pF <- if (df_model > 0L) {
    pf(aux_F, df_model, df2, lower.tail = FALSE)
  } else NA_real_
  aux_stats <- .ev_equation_stats(e, residuals(fit), ncol(Z),
                    fstat = aux_F, fpval = aux_pF,
                    r2 = r2_aux, adjr2 = r2_aux_adj)
  out <- .ev_result("Breusch-Godfrey Serial Correlation LM Test",
    c(`F-statistic` = F, `Prob. F` = pf(F, p, df2, lower.tail = FALSE),
      `Obs*R-squared` = LM, `Prob. Chi-Square` = pchisq(LM, p, lower.tail = FALSE)),
    tab, list(`Lags` = p,
              `Residual degrees of freedom` = df2,
              `Included observations` = n,
              `Presample lagged residuals` = "zero",
              `Equation statistics` = aux_stats), fit)
  print(out)
  invisible(out)
}
