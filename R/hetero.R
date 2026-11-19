#' Breusch-Pagan-Godfrey or White heteroskedasticity test
#' @param model An unweighted, full-rank lm object with an intercept.
#' @param type Either "bpg" or "white".
#' @param cross Include pairwise cross products for the White test (default FALSE).
#' @return An evreport object with F and Obs*R-squared statistics.
#' @details For BPG, auxiliary response is squared OLS residuals and regressors
#'   are the original design matrix. For White, auxiliary regressors are the
#'   design columns, their squares and optionally cross products. This is the
#'   unstudentized n*R2 variant; EViews can also display a scaled-ESS test.
#' @export
#' @examples
#' ev_hetero(lm(mpg ~ wt + hp, mtcars), "white")
ev_hetero <- function(model, type = c("bpg", "white"), cross = FALSE) {
  X <- .ev_check(model)
  type <- match.arg(type)
  if (attr(terms(model), "intercept") != 1L)
    stop("Heteroskedasticity tests currently require an intercept.", call. = FALSE)
  y <- as.numeric(residuals(model)^2)
  base <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  if (ncol(base) < 1L) stop("At least one slope regressor is required.", call. = FALSE)
  A <- cbind(`(Intercept)` = 1, base)
  if (type == "white") {
    sq <- base^2
    colnames(sq) <- paste0(colnames(base), "^2")
    A <- cbind(A, sq)
    if (isTRUE(cross) && ncol(base) > 1L) {
      pairs <- combn(ncol(base), 2L)
      prod <- sapply(seq_len(ncol(pairs)), function(i) base[,pairs[1,i]]*base[,pairs[2,i]])
      if (is.null(dim(prod))) prod <- matrix(prod, ncol = 1L)
      colnames(prod) <- apply(pairs, 2L, function(z) paste(colnames(base)[z], collapse = ":"))
      A <- cbind(A, prod)
    }
  }
  # Drop redundant columns deterministically (e.g., square of dummy).
  qrA <- qr(A)
  kept <- sort(qrA$pivot[seq_len(qrA$rank)])
  A <- A[, kept, drop = FALSE]
  if (ncol(A) < 2L) stop("No identifiable auxiliary slopes.", call. = FALSE)
  joint <- .ev_joint(y, A, restricted_cols = 1L)
  fit <- lm(y ~ A - 1)
  tab <- summary(fit)$coefficients
  rownames(tab) <- colnames(A)
  out <- .ev_result(if (type == "bpg") "Heteroskedasticity Test: Breusch-Pagan-Godfrey" else
    paste0("Heteroskedasticity Test: White", if (isTRUE(cross)) " (with cross terms)" else ""),
    c(`F-statistic` = joint$F, `Prob. F` = joint$pF,
      `Obs*R-squared` = joint$LM, `Prob. Chi-Square` = joint$pLM), tab,
    list(`Null hypothesis` = "Homoskedasticity",
         `Prob. F degrees of freedom` = paste0("(", joint$df1, ", ", joint$df2, ")"),
         `Prob. Chi-Square degrees of freedom` = joint$df1), fit)
  print(out)
  invisible(out)
}
