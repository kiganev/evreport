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
