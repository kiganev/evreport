# evreport (prototype 0.1.0)

An **unofficial, EViews-inspired** console report package for ordinary unweighted R `lm` models.

This is a prototype, **not yet validated against actual EViews regression reports**. The output is inspired by EViews, not a claim of bit-for-bit numerical or visual equivalence.

## Installation

The development version of `evreport` can be installed directly from GitHub using the `remotes` package.

First, install `remotes` if it is not already installed:

```r
install.packages("remotes")
```

Then install `evreport`:

```r
remotes::install_github("USERNAME/evreport")
```

Replace `USERNAME` with the GitHub username of the repository owner.

Once installed, load the package:

```r
library(evreport)
```

To update to the latest development version, run `remotes::install_github()` again.

## Examples

```r
library(evreport)
fit <- lm(mpg ~ wt + hp, data = mtcars)

ev_summary(fit)
ev_bg(fit, lags = 2)
ev_hetero(fit, type = "bpg")
ev_hetero(fit, type = "white", cross = FALSE)
ev_hetero(fit, type = "white", cross = TRUE)
ev_reset(fit, powers = 2:3)
ev_normality(fit)

# Retrieve a test result
res <- ev_bg(fit, lags = 2)
res$statistics
res$auxiliary_model
```

## Boundaries and conventions

* Only ordinary **unweighted, full-rank, non-offset `lm`** regressions supported.
* BG uses presample-zero lagged residuals and **uncentered** `n * R²`, as specified in EViews' documentation. Observations must be in proper consecutive time order. The BG F probability is an approximation, as EViews notes.
* White default is **without cross terms**; `cross = TRUE` adds pairwise products. Redundant terms (e.g. squared binary dummies) are dropped via QR rank selection. The White/BPG helper provides F and `Obs*R-squared` but **not** EViews' additional scaled explained sum of squares statistic.
* BPG here refers to the unstudentized squared-residual auxiliary regression `n * R²`, rather than R's default studentized `lmtest::bptest`. Such variants can differ.
* RESET uses fitted values raised to user-selected powers. Only the standard nested-model F test is printed.
* Information criteria are normalized by sample size and use **k regression coefficients**, excluding sigma², following the EViews-style convention.
* No automatic detection of gaps in time-series data or date ranges. EViews tests may differ for ARMA regressions, different covariance estimators, and other advanced options.
* The helper's print format is **EViews-inspired**, not exact pixel/line reproduction.

## Development/testing

```r
install.packages("testthat")
testthat::test_local("evreport")
```

Validation still needed: compare the results with EViews on identical time-series data, including intercept/no-intercept cases, near-collinearity and lag sensitivity.

Official documentation: https://help.eviews.com/content/testing-Residual_Diagnostics.html
