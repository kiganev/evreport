#' Ramsey RESET specification test using powers of fitted values
#' @param model An unweighted, full-rank lm model.
#' @param powers Integer vector of fitted-value powers, default 2:3.
#' @return An evreport object containing F-test and likelihood-ratio-style statistics.
#' @details Uses an ordinary nested-model F-test. It is not identical to all
#'   EViews RESET variants or covariance options.
#' @export
#' @examples
#' ev_reset(lm(mpg ~ wt + hp, mtcars))
ev_reset <- function(model, powers = 2:3) {
  X <- .ev_check(model)
  if (!is.numeric(powers) || !length(powers) || anyNA(powers) ||
      any(!is.finite(powers)) || any(powers < 2) || any(powers != floor(powers)) ||
      anyDuplicated(powers)) stop("powers must be distinct integers >= 2.", call. = FALSE)
  y <- model.response(model.frame(model))
  yh <- as.numeric(fitted(model))
  P <- sapply(powers, function(j) yh^j)
  if (is.null(dim(P))) P <- matrix(P, ncol = 1)
  colnames(P) <- paste0("FITTED^", powers)
  Z <- cbind(X, P)
  aux <- lm.fit(Z, y)
  if (aux$rank != ncol(Z)) stop("RESET auxiliary regression is rank deficient.", call. = FALSE)
  rss0 <- sum(residuals(model)^2)
  rss <- sum(aux$residuals^2)
  q <- ncol(P)
  df2 <- length(y) - ncol(Z)
  if (df2 <= 0 || rss <= 0) stop("Insufficient degrees of freedom.", call. = FALSE)
  F <- max(0, (rss0-rss)/q)/(rss/df2)
  out <- .ev_result("Ramsey RESET Test", c(`F-statistic`=F,
    `Prob. F`=pf(F, q, df2, lower.tail=FALSE)),
    details=list(`Null hypothesis`="Correct functional form",
                 `Included powers`=paste(powers, collapse=", "),
                 `Prob. F degrees of freedom`=paste0("(",q,", ",df2,")")))
  print(out)
  invisible(out)
}
