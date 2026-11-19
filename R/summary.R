#' EViews-style OLS equation summary
#' @param model An unweighted, full-rank lm object.
#' @return An evreport object (invisibly printed), with statistics and coefficient table.
#' @export
#' @examples
#' fit <- lm(mpg ~ wt + hp, data = mtcars)
#' ev_summary(fit)
ev_summary <- function(model) {
  .ev_check(model)
  s <- summary(model)
  e <- residuals(model)
  y <- model.response(model.frame(model))
  n <- length(e)
  k <- model$rank
  ll <- as.numeric(logLik(model))
  crit <- c(AIC = (-2 * ll + 2 * k) / n,
            `Schwarz criterion` = (-2 * ll + log(n) * k) / n,
            `Hannan-Quinn criterion` = (-2 * ll + 2 * k * log(log(n))) / n)
  fs <- s$fstatistic
  F <- if (is.null(fs)) NA_real_ else unname(fs[1])
  pF <- if (is.null(fs)) NA_real_ else pf(fs[1], fs[2], fs[3], lower.tail = FALSE)
  stats <- c(`R-squared` = s$r.squared, `Adjusted R-squared` = s$adj.r.squared,
      `S.E. of regression` = s$sigma, `Sum squared resid` = sum(e^2),
      `Log likelihood` = ll, `F-statistic` = F, `Prob(F-statistic)` = pF,
      `Mean dependent var` = mean(y), `S.D. dependent var` = sd(y), crit,
      `Durbin-Watson stat` = sum(diff(e)^2) / sum(e^2))
  out <- .ev_result("Least Squares", stats, s$coefficients,
    list(`Dependent Variable` = paste(deparse(formula(model)[[2]]), collapse = ""),
         `Included observations` = n))
  print(out)
  invisible(out)
}
