#' Ramsey RESET specification test using powers of fitted values
#'
#' @param model An unweighted, full-rank lm model.
#' @param powers Distinct integer powers of the fitted values; default 2:3.
#' @return An evreport object containing the F and LR tests, F-test and LR-test
#'   decomposition tables, unrestricted coefficients and an auxiliary lm fit.
#' @details The unrestricted equation adds fitted-value powers to the original
#'   design matrix. The LR statistic is n * log(SSR_restricted / SSR_unrestricted)
#'   and its p-value is based on a chi-squared distribution with one degree of
#'   freedom per added power. Conventional (non-robust) OLS tests only.
#' @export
#' @examples
#' ev_reset(lm(mpg ~ wt + hp, data = mtcars))
ev_reset <- function(model, powers = 2:3) {
  X <- .ev_check(model)
  if (!is.numeric(powers) || !length(powers) || anyNA(powers) ||
      any(!is.finite(powers)) || any(powers < 2) ||
      any(powers != floor(powers)) || anyDuplicated(powers))
    stop("powers must be distinct integers >= 2.", call. = FALSE)
  powers <- as.integer(powers)
  y <- as.numeric(model.response(model.frame(model)))
  yh <- as.numeric(fitted(model))
  P <- sapply(powers, function(j) yh^j)
  if (is.null(dim(P))) P <- matrix(P, ncol = 1L)
  colnames(P) <- paste0("FITTED^", powers)
  Z <- cbind(X, P)
  n <- length(y)
  k0 <- ncol(X)
  ku <- ncol(Z)
  q <- ncol(P)
  df0 <- n - k0
  dfu <- n - ku
  if (dfu <= 0L) stop("Insufficient degrees of freedom.", call. = FALSE)

  # Using a matrix in lm() retains the original design matrix exactly,
  # including factor contrasts or transformed regressors.
  aux <- lm(y ~ Z - 1)
  if (aux$rank != ku)
    stop("RESET unrestricted regression is rank deficient.", call. = FALSE)
  rss0 <- sum(residuals(model)^2)
  rssu <- sum(residuals(aux)^2)
  if (rss0 <= 0 || rssu <= 0)
    stop("Residual sum of squares must be positive.", call. = FALSE)
  ss_test <- max(0, rss0 - rssu)
  F <- (ss_test / q) / (rssu / dfu)
  pF <- pf(F, q, dfu, lower.tail = FALSE)
  ll0 <- -n / 2 * (log(2 * pi) + 1 + log(rss0 / n))
  llu <- -n / 2 * (log(2 * pi) + 1 + log(rssu / n))
  LR <- max(0, 2 * (llu - ll0))
  pLR <- pchisq(LR, q, lower.tail = FALSE)

  tab <- summary(aux)$coefficients
  rownames(tab) <- colnames(Z)
  has_intercept <- attr(terms(model), "intercept") == 1L
  tss <- if (has_intercept) sum((y - mean(y))^2) else sum(y^2)
  if (tss <= 0) stop("Dependent variable has no variation.", call. = FALSE)
  r2 <- 1 - rssu / tss
  df_model <- ku - as.integer(has_intercept)
  adjr2 <- 1 - (rssu / dfu) / (tss / (n - as.integer(has_intercept)))
  overall_F <- if (df_model > 0L) max(0, (tss - rssu) / df_model) / (rssu / dfu) else NA_real_
  overall_pF <- if (df_model > 0L) pf(overall_F, df_model, dfu, lower.tail = FALSE) else NA_real_
  eqstats <- .ev_equation_stats(y, as.numeric(residuals(aux)), ku,
    fstat = overall_F, fpval = overall_pF, r2 = r2, adjr2 = adjr2)

  f_summary <- data.frame(
    `Sum of Sq.` = c(ss_test, rss0, rssu),
    df = c(q, df0, dfu),
    `Mean Squares` = c(ss_test / q, rss0 / df0, rssu / dfu),
    row.names = c("Test SSR", "Restricted SSR", "Unrestricted SSR"),
    check.names = FALSE)
  lr_summary <- c(`Restricted LogL` = ll0, `Unrestricted LogL` = llu)
  rhs <- colnames(X)
  rhs[rhs == "(Intercept)"] <- "C"
  out <- .ev_result("Ramsey RESET Test",
    c(`F-statistic` = F, `Prob. F` = pF,
      `Likelihood ratio` = LR, `Prob. Chi-Square` = pLR),
    tab,
    list(`Dependent Variable` = paste(deparse(formula(model)[[2]]), collapse = ""),
         `Specification` = paste(c(paste(deparse(formula(model)[[2]]), collapse = ""), rhs), collapse = " "),
         `Powers` = powers,
         `Included observations` = n,
         `Test degrees of freedom` = q,
         `Residual degrees of freedom` = dfu,
         `F-test summary` = f_summary,
         `LR test summary` = lr_summary,
         `Equation statistics` = eqstats), aux)
  print(out)
  invisible(out)
}
