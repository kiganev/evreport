#' Jarque-Bera test for OLS residual normality
#' @param model An unweighted, full-rank lm object.
#' @return An evreport object containing sample moments and Jarque-Bera test.
#' @export
#' @examples
#' ev_normality(lm(mpg ~ wt + hp, mtcars))
ev_normality <- function(model) {
  .ev_check(model)
  x <- as.numeric(residuals(model))
  n <- length(x)
  v <- mean((x-mean(x))^2)
  if (v <= 0) stop("Residual variance is zero.", call. = FALSE)
  s <- mean((x-mean(x))^3) / v^(3/2)
  kurt <- mean((x-mean(x))^4) / v^2
  jb <- n / 6 * (s^2 + (kurt - 3)^2 / 4)
  out <- .ev_result("Residual Normality Test: Jarque-Bera",
    c(`Mean`=mean(x), `Median`=median(x), `Maximum`=max(x),
      `Minimum`=min(x), `Std. Dev.`=sd(x), `Skewness`=s,
      `Kurtosis`=kurt, `Jarque-Bera`=jb,
      `Probability`=pchisq(jb, 2L, lower.tail=FALSE)),
    details=list(`Null hypothesis`="Normally distributed residuals", `Observations`=n))
  print(out)
  invisible(out)
}
