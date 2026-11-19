---
title: "evreport Reference Manual"
author: "Kaloyan Ganev"
date: "19 November 2026"
lang: en-US
fontfamily: mathpazo
fontsize: 11pt
geometry: margin=2.5cm
colorlinks: true
linkcolor: blue
header-includes:
  - \usepackage{fancyhdr}
  - \pagestyle{fancy}
  - \fancyhead[L]{evreport · Reference Manual}
  - \fancyhead[R]{v0.1.6}
  - \fancyfoot[C]{\thepage}
---

# Overview

`evreport` provides **EViews-inspired console reports** and standard residual diagnostics for ordinary least-squares models estimated with R's `lm()` function. This manual documents the actual API in the **v0.1.6 source archive**. It is not an official EViews product and does not promise exact agreement across every estimation option.

The package contains five user-facing analysis functions and an S3 printing method.

| Function | Purpose |
|:--|:--|
| `ev_summary(model)` | Full OLS equation report |
| `ev_bg(model, lags = 1L)` | Breusch–Godfrey serial-correlation LM test |
| `ev_hetero(model, type = c("bpg", "white"), cross = FALSE)` | Breusch–Pagan–Godfrey / White heteroskedasticity tests |
| `ev_reset(model, powers = 2:3)` | Ramsey RESET F and likelihood-ratio tests |
| `ev_normality(model, plot = TRUE, breaks = "Sturges")` | Jarque–Bera test and residual histogram |
| `print.evreport(x, ...)` | Print a stored result again |

## Installation

From the GitHub repository (substitute the actual GitHub account for `USERNAME`):

```r
install.packages("remotes")  # only once
remotes::install_github("USERNAME/evreport")
library(evreport)
```

For a local source folder, run from its parent directory, or use an explicit path:

```r
install.packages("path/to/evreport", repos = NULL, type = "source")
```

During package development from within the repository:

```r
devtools::load_all()
```

The repository root must contain `DESCRIPTION`, `NAMESPACE`, `R/`, and `man/`. The version documented here is 0.1.6.

## Quick start

```r
library(evreport)
fit <- lm(mpg ~ wt + hp, data = mtcars)

ev_summary(fit)
ev_bg(fit, lags = 2)  # demonstration only: mtcars is not a time series
ev_hetero(fit, type = "bpg")
ev_hetero(fit, type = "white", cross = FALSE)
ev_hetero(fit, type = "white", cross = TRUE)
ev_reset(fit, powers = 2:3)
ev_normality(fit)                 # console summary plus histogram
ev_normality(fit, plot = FALSE)   # console summary only
```

**Important:** `mtcars` is cross-sectional data. The BG call illustrates syntax and should not be interpreted as a meaningful temporal serial-correlation analysis. Use correctly ordered, consecutive time-series observations for that diagnostic.

# Function reference

## `ev_summary()` — OLS equation report

```r
ev_summary(model)
```

**Argument**  
`model`: an ordinary, unweighted, full-rank `lm` model.

**Report.** The dependent variable, estimation method, observation count, coefficient table (estimate, standard error, t-statistic, p-value), and two-column equation statistics: R-squared, adjusted R-squared, S.E. of regression, residual sum of squares, log likelihood, model F-statistic and probability, dependent-variable mean and standard deviation, Akaike, Schwarz and Hannan–Quinn criteria, and Durbin–Watson statistic.

```r
fit <- lm(mpg ~ wt + hp, data = mtcars)
s <- ev_summary(fit)
s$statistics
s$coefficients
```

The information criteria are divided by sample size and count regression coefficients, excluding the error-variance parameter. This is a deliberate EViews-style convention and may differ from the output of `AIC(fit)`.

## `ev_bg()` — Breusch–Godfrey serial correlation test

```r
ev_bg(model, lags = 1L)
```

**Arguments**  
`model`: a supported `lm` object.  
`lags`: a positive integer giving the maximum residual lag.

**Null hypothesis.** No residual serial correlation through the specified lag order.

**Report.** A test F-statistic and its F-distribution probability; `Obs*R-squared` with its chi-square probability; a complete auxiliary regression of the original residuals on original regressors and lagged residuals; and the auxiliary equation's two-column summary statistics.

The function fills unavailable presample lagged residuals with **zeros** and retains the initial sample observations. Its LM statistic uses an **uncentered** auxiliary R-squared. It assumes observations are already in the correct, consecutive chronological order; it does not infer missing periods from dates.

```r
bg <- ev_bg(fit, lags = 2)
bg$statistics
bg$auxiliary_model
```

The overall F-statistic printed at the bottom of the auxiliary equation is **not** the same hypothesis as the BG F-statistic at the top: the former is the equation-wide F-test, while the latter concerns only lagged residuals.

## `ev_hetero()` — Heteroskedasticity diagnostics

```r
ev_hetero(model, type = c("bpg", "white"), cross = FALSE)
```

**Arguments**  
`model`: supported `lm` model **with an intercept**.  
`type`: `"bpg"` or `"white"`; if omitted, the default is `"bpg"`.  
`cross`: logical; affects the White specification only.

**Null hypothesis.** Homoskedastic errors.

**Report.** The F-statistic, `Obs*R-squared`, and `Scaled explained SS`, with probabilities; the entire OLS auxiliary equation with dependent variable `RESID^2`; and a two-column summary of its fit statistics.

The auxiliary dependent variable is the squared original OLS residual. For `type = "bpg"`, regressors are the original non-intercept explanatory variables plus a constant. This is the **unstudentized** auxiliary-regression variant, so results can differ from the default studentized `lmtest::bptest()`.

For `type = "white"`:

- `cross = FALSE`: a constant and **squares of the original explanatory columns only**.
- `cross = TRUE`: a constant, those squares, their levels, and their pairwise cross-products.

The implementation removes rank-redundant auxiliary columns using QR decomposition. For models with factors, transformed terms, or strongly collinear regressors, the retained terms should be inspected.

```r
w0 <- ev_hetero(fit, type = "white", cross = FALSE)
w1 <- ev_hetero(fit, type = "white", cross = TRUE)
bp <- ev_hetero(fit, type = "bpg")
w1$coefficients
w1$details[["Equation statistics"]]
```

**Scaled explained SS** is the auxiliary explained sum of squares divided by twice the squared original residual-variance estimate, where the variance estimate is `SSR / (n - k)`. The test output is designed to reproduce the form of the supplied EViews comparison reports, not every alternative White-test convention.

## `ev_reset()` — Ramsey RESET test

```r
ev_reset(model, powers = 2:3)
```

**Arguments**  
`model`: supported `lm` object.  
`powers`: a nonempty vector of **distinct integer powers of fitted values**, each at least 2; default `2:3`.

**Null hypothesis.** The added fitted-value powers have jointly zero coefficients in the unrestricted regression.

**Report.** The F-test and likelihood-ratio statistic, their degrees of freedom and probabilities; the F-test sum-of-squares decomposition (test, restricted and unrestricted SSR); the restricted and unrestricted log likelihoods; and the full unrestricted regression and its two-column statistics.

The unrestricted equation includes the original model-matrix columns plus fitted-value powers, so columns `FITTED^2` and `FITTED^3` appear with the default options. The likelihood-ratio statistic uses `n * log(SSR_restricted / SSR_unrestricted)` with a chi-square reference distribution based on the number of added powers. The method uses conventional non-robust OLS inference.

```r
reset <- ev_reset(fit)
reset$statistics
reset$details[["F-test summary"]]
reset$details[["LR test summary"]]
summary(reset$auxiliary_model)
```

An augmented regression with insufficient degrees of freedom or rank-deficient added columns is rejected.

## `ev_normality()` — Jarque–Bera test and histogram

```r
ev_normality(model, plot = TRUE, breaks = "Sturges")
```

**Arguments**  
`model`: supported `lm` object.  
`plot`: `TRUE` (default) to draw the histogram and statistics panel, or `FALSE` for console output alone.  
`breaks`: any valid histogram breaks specification accepted by `graphics::hist()`; default `"Sturges"`.

**Null hypothesis.** Normally distributed residuals.

**Report.** Residual mean, median, maximum, minimum, standard deviation, skewness, kurtosis, Jarque–Bera statistic, and chi-square probability (two degrees of freedom). The optional base-R graph puts the histogram on the left and statistics on the right.

```r
norm <- ev_normality(fit)
no_plot <- ev_normality(fit, plot = FALSE)
ev_normality(fit, breaks = 8)
norm$statistics
```

Changing `breaks` changes the visual bins, **not** the normality statistics. EViews can choose different bin boundaries. The Jarque–Bera statistic here uses residual moments with denominators `n`, as implemented in the package.

## `print.evreport()` — Print a stored report

```r
print(x)
```

When an object returned by one of the analysis functions has class `evreport`, base R dispatches to `print.evreport(x, ...)`. It recreates the relevant console layout; the `...` arguments are currently ignored.

```r
res <- ev_reset(fit)  # prints immediately
print(res)             # print again later
```

# Working with returned results

All analysis functions **print immediately** and **return an object invisibly**, so assignment both displays and stores the result.

The result is a list with class `evreport` and these common fields:

| Field | Description |
|:--|:--|
| `$title` | Report heading |
| `$statistics` | Named numeric vector of main test or model statistics |
| `$coefficients` | Coefficient matrix, if relevant; otherwise `NULL` |
| `$details` | Named list with metadata and supplementary tables |
| `$auxiliary_model` | Auxiliary/unrestricted `lm` model, if relevant; otherwise `NULL` |

```r
bg <- ev_bg(fit, lags = 2)
bg$statistics[["Obs*R-squared"]]
bg$statistics[["Prob. Chi-Square"]]
coef(bg$auxiliary_model)

white <- ev_hetero(fit, type = "white", cross = TRUE)
white$details[["Equation statistics"]]
```

Note that `$statistics` represents the **headline** diagnostics for most tests. The auxiliary equation's regression statistics are stored in `$details[["Equation statistics"]]` when applicable.

# Statistical conventions and limitations

1. **Supported models.** Functions reject weighted, rank-deficient and offset `lm` fits. Most accept ordinary intercept or no-intercept OLS, but `ev_hetero()` requires an intercept. There must be enough observations to estimate the relevant model.
2. **Sample dependence.** Diagnostic results depend on the original estimation sample. The functions do not infer or repair gaps in a time series. A BG test on arbitrarily ordered cross-sectional observations is not substantively meaningful.
3. **Inference variants.** Tests use standard, non-robust OLS variants. Some alternatives in R, EViews or textbooks use studentization, different auxiliary regressors, or different small-sample adjustments.
4. **Information criteria.** Regression summaries use sample-normalized EViews-style AIC/Schwarz/Hannan–Quinn definitions based on the regression coefficient count. These need not equal base R's generic `AIC()`/`BIC()` values.
5. **Plots.** The histogram uses base graphics and is EViews-inspired, not pixel-identical.
6. **Validation status.** This is a prototype. Several procedures were compared against user-supplied EViews reference values, but comprehensive cross-platform testing and full numerical validation have **not** been established.
7. **Output metadata.** EViews-style date, time, workfile identifier and sample period are deliberately omitted; an `lm` object does not guarantee those fields.

# Reproducible package checks

From the local package root:

```r
devtools::load_all()
devtools::test()
devtools::check()
```

For help within R:

```r
?ev_summary
?ev_bg
?ev_hetero
?ev_reset
?ev_normality
```

The package contains standard Rd help pages in `man/`. This reference manual expands those topics with workflow examples and limitations.

# Notes on attribution

`evreport` is an **independent, unofficial package** inspired by the layout of EViews reports. The package's source declares `MIT + file LICENSE`; review the repository's actual license files and authorship metadata before a public release. The source archive inspected for this manual contains placeholder author/contact fields in `DESCRIPTION`.

For the official methodology and documentation of EViews diagnostic procedures, consult the [EViews Residual Diagnostics documentation](https://help.eviews.com/content/testing-Residual_Diagnostics.html).
