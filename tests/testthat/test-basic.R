test_that("OLS summary agrees with lm statistics", {
  fit <- lm(mpg ~ wt + hp, mtcars)
  obj <- suppressMessages(capture.output(res <- ev_summary(fit)))
  expect_equal(unname(res$statistics["R-squared"]), summary(fit)$r.squared)
  expect_equal(unname(res$statistics["Sum squared resid"]), sum(resid(fit)^2))
})
test_that("BG produces two valid statistics and lag controls", {
  fit <- lm(mpg ~ wt + hp, mtcars)
  capture.output(bg <- ev_bg(fit, 2))
  expect_true(all(is.finite(bg$statistics)))
  expect_equal(unname(bg$statistics["Obs*R-squared"]),
    length(resid(fit)) * (1 - sum(resid(bg$auxiliary_model)^2) / sum(resid(fit)^2)),
    tolerance = 1e-8)
  expect_error(ev_bg(fit, 0), "positive integer")
})
test_that("diagnostics return finite test statistics", {
  fit <- lm(mpg ~ wt + hp, mtcars)
  capture.output(bp <- ev_hetero(fit, "bpg"))
  capture.output(white <- ev_hetero(fit, "white", cross = TRUE))
  capture.output(reset <- ev_reset(fit))
  capture.output(jb <- ev_normality(fit))
  expect_true(all(is.finite(bp$statistics)))
  expect_true(all(is.finite(white$statistics)))
  expect_true(all(is.finite(reset$statistics)))
  expect_true(all(is.finite(jb$statistics)))
})

test_that("EViews print layout puts coefficients ahead of paired statistics", {
  fit <- lm(mpg ~ wt + hp, mtcars)
  s <- capture.output(ev_summary(fit))
  expect_lt(which(grepl("Coefficient", s, fixed = TRUE))[1],
            which(grepl("R-squared", s, fixed = TRUE))[1])
  expect_true(any(grepl("R-squared.*Mean dependent var", s)))
  expect_true(any(grepl("Prob\\(F-statistic\\).*", s)))
  expect_false(any(grepl("Prob(F-statistic).value", s, fixed = TRUE)))

  b <- capture.output(ev_bg(fit, lags = 2))
  expect_true(any(grepl("F-statistic.*Prob. F\\(2,27\\)", b)))
  expect_true(any(grepl("Obs\\*R-squared.*Prob. Chi-Square\\(2\\)", b)))
  expect_true(any(grepl("R-squared.*Mean dependent var", b)))
})


test_that("BG test equation prints its own overall F-statistic", {
  fit <- lm(mpg ~ wt + hp, mtcars)
  lines <- capture.output(bg <- ev_bg(fit, lags = 2))
  aux <- bg$details[["Equation statistics"]]
  expect_true(is.finite(unname(aux[["F-statistic"]])))
  expect_true(is.finite(unname(aux[["Prob(F-statistic)"]])))
  expect_equal(unname(aux[["F-statistic"]]),
               unname(aux[["R-squared"]]) / 4 /
                 ((1 - unname(aux[["R-squared"]])) / 27), tolerance = 1e-8)
  expect_false(any(grepl("F-statistic\\s+NA", lines)))
})
