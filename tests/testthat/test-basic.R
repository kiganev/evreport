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
  capture.output(jb <- ev_normality(fit, plot = FALSE))
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


test_that("heteroskedasticity output has full EViews layout", {
  fit <- lm(mpg ~ wt + hp, mtcars)
  for (kind in c("bpg", "white")) {
    lines <- capture.output(result <- ev_hetero(fit, type = kind))
    expect_true(any(grepl("Scaled explained SS.*Prob. Chi-Square", lines)))
    expect_true(any(grepl("F-statistic.*Prob. F", lines)))
    expect_true(any(grepl("Obs\\*R-squared.*Prob. Chi-Square", lines)))
    expect_true(any(grepl("Dependent Variable: RESID^2", lines, fixed = TRUE)))
    expect_true(any(grepl("R-squared.*Mean dependent var", lines)))
    expect_lt(which(grepl("Coefficient", lines, fixed = TRUE))[1],
              which(grepl("R-squared.*Mean dependent var", lines))[1])
    ssr <- sum(resid(fit)^2)
    aux <- result$auxiliary_model
    y <- resid(fit)^2
    ess <- sum((y - mean(y))^2) - sum(resid(aux)^2)
    expected <- ess / (2 * (ssr / df.residual(fit))^2)
    expect_equal(unname(result$statistics[["Scaled explained SS"]]),
                 expected, tolerance = 1e-8)
    auxstats <- result$details[["Equation statistics"]]
    expect_equal(unname(auxstats[["F-statistic"]]),
                 unname(result$statistics[["F-statistic"]]))
  }
})


test_that("White variants reproduce supplied EViews mtcars results", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  capture.output(no_cross <- ev_hetero(fit, type = "white", cross = FALSE))
  capture.output(full <- ev_hetero(fit, type = "white", cross = TRUE))

  expect_identical(rownames(no_cross$coefficients), c("(Intercept)", "wt^2", "hp^2"))
  expect_identical(rownames(full$coefficients),
    c("(Intercept)", "wt^2", "wt*hp", "wt", "hp^2", "hp"))

  expect_equal(unname(no_cross$statistics[["F-statistic"]]),
               0.2266091157943661, tolerance = 1e-6)
  expect_equal(unname(no_cross$statistics[["Obs*R-squared"]]),
               0.4924074271545962, tolerance = 1e-6)
  expect_equal(unname(no_cross$statistics[["Scaled explained SS"]]),
               0.4714688501474768, tolerance = 1e-6)
  expect_equal(unname(full$statistics[["F-statistic"]]),
               1.336534710167772, tolerance = 1e-6)
  expect_equal(unname(full$statistics[["Obs*R-squared"]]),
               6.543086302110516, tolerance = 1e-6)
  expect_equal(unname(full$statistics[["Scaled explained SS"]]),
               6.264855493951078, tolerance = 1e-6)
})


test_that("RESET reproduces EViews mtcars test and full layout", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  lines <- capture.output(reset <- ev_reset(fit, powers = 2:3))
  expect_equal(unname(reset$statistics[["F-statistic"]]),
               7.238403935158949, tolerance = 1e-6)
  expect_equal(unname(reset$statistics[["Likelihood ratio"]]),
               13.73752187267368, tolerance = 1e-6)
  expect_equal(unname(reset$statistics[["Prob. F"]]),
               0.003041157619883153, tolerance = 1e-6)
  expect_equal(unname(reset$statistics[["Prob. Chi-Square"]]),
               0.001039764601316361, tolerance = 1e-6)
  ftab <- reset$details[["F-test summary"]]
  expect_equal(ftab["Test SSR", "Sum of Sq."],
               68.07825905401455, tolerance = 1e-6)
  expect_equal(ftab["Restricted SSR", "Sum of Sq."],
               195.0477547414661, tolerance = 1e-6)
  expect_equal(ftab["Unrestricted SSR", "Sum of Sq."],
               126.9694956874515, tolerance = 1e-6)
  expect_identical(rownames(reset$coefficients),
                   c("(Intercept)", "wt", "hp", "FITTED^2", "FITTED^3"))
  expect_true(any(grepl("F-test summary:", lines, fixed = TRUE)))
  expect_true(any(grepl("LR test summary:", lines, fixed = TRUE)))
  expect_true(any(grepl("Unrestricted Test Equation:", lines, fixed = TRUE)))
  expect_true(any(grepl("R-squared.*Mean dependent var", lines)))
  expect_lt(which(grepl("Coefficient", lines, fixed = TRUE))[1],
            which(grepl("R-squared.*Mean dependent var", lines))[1])
  expect_true(all(is.finite(reset$details[["Equation statistics"]])))
})


test_that("normality histogram can be disabled without changing results", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  out <- capture.output(a <- ev_normality(fit, plot = FALSE))
  expect_s3_class(a, "evreport")
  expect_equal(unname(a$statistics[["Jarque-Bera"]]),
               32 / 6 * (a$statistics[["Skewness"]]^2 +
               (a$statistics[["Kurtosis"]] - 3)^2 / 4))
  expect_true(any(grepl("Jarque-Bera", out, fixed = TRUE)))
  expect_error(ev_normality(fit, plot = NA), "plot")
})

test_that("normality histogram renders in base graphics", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path, width = 10, height = 5)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  capture.output(a <- ev_normality(fit, breaks = 8))
  expect_s3_class(a, "evreport")
  expect_true(is.finite(a$statistics[["Probability"]]))
})
