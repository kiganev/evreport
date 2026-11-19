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
  fs <- s$fstatistic
  F <- if (is.null(fs)) NA_real_ else unname(fs[1])
  pF <- if (is.null(fs)) NA_real_ else
    unname(pf(fs[1], fs[2], fs[3], lower.tail = FALSE))
  stats <- .ev_equation_stats(y, e, k, F, pF,
                              r2 = s$r.squared, adjr2 = s$adj.r.squared)
  out <- .ev_result("Least Squares", stats, s$coefficients,
    list(`Dependent Variable` = paste(deparse(formula(model)[[2]]), collapse = ""),
         `Included observations` = n))
  print(out)
  invisible(out)
}
