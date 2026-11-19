# evreport (prototype 0.1.6)

An **unofficial, EViews-inspired** console report package for ordinary unweighted R `lm` models.

This is a prototype, **not yet validated against actual EViews regression reports**. The output is inspired by EViews, not a claim of bit-for-bit numerical or visual equivalence.

## Install

From the **directory containing `evreport/`**, run:

```r
install.packages("evreport", repos = NULL, type = "source")
```

If the package is supplied as a zip archive, **unzip it first**. The archive is a source-folder distribution, not an R Windows binary package.

## Install directly from GitHub

```r
install.packages("remotes") # once only
remotes::install_github("kiganev/evreport")
```

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
ev_normality(fit)                   # histogram and statistics panel
ev_normality(fit, plot = FALSE)     # console output only
ev_normality(fit, breaks = 8)      # manually select histogram bins

# Retrieve a test result
res <- ev_bg(fit, lags = 2)
res$statistics
res$auxiliary_model
```

## Boundaries and conventions

* Only ordinary **unweighted, full-rank, non-offset `lm`** regressions supported.
* BG uses presample-zero lagged residuals and **uncentered** `n * R²`, as specified in EViews' documentation. Observations must be in proper consecutive time order. The BG F probability is an approximation, as EViews notes.
* White default uses **squared regressors only**; `cross = TRUE` uses their squares, levels, and pairwise products. The latter includes five slopes for `mpg ~ wt + hp`, matching the EViews comparison. Redundant terms (e.g. squared binary dummies) are dropped via QR rank selection. The White/BPG output includes F, `Obs*R-squared`, and `Scaled explained SS`, followed by the auxiliary squared-residual regression and the full two-column equation statistics. Date/time and sample period are omitted.
* BPG here refers to the unstudentized squared-residual auxiliary regression `n * R²`, not necessarily the same as the default studentized `lmtest::bptest`. Such variants can differ.
* RESET uses fitted values raised to user-selected powers; it reports F and likelihood-ratio tests and the unrestricted regression.
* Information criteria are normalized by sample size and use **k regression coefficients**, excluding sigma², following the EViews-style convention.
* The residual normality test plots a histogram with descriptive statistics and the Jarque-Bera probability in a side panel using base R graphics. The default histogram bins use Sturges' rule; EViews may use different bins. Set `plot = FALSE` to suppress graphics.
* No automatic detection of gaps in time-series data or date ranges. EViews tests may differ for ARMA regressions, different covariance estimators, and other advanced options.
* The helper's print format is **EViews-inspired**, not exact pixel/line reproduction.

## Development/testing

```r
install.packages("testthat")
testthat::test_local("evreport")
```

Validation still needed: compare the results with EViews on identical time-series data, including intercept/no-intercept cases, near-collinearity and lag sensitivity.

Official documentation: https://help.eviews.com/content/testing-Residual_Diagnostics.html

### RESET output

`ev_reset(fit, powers = 2:3)` now shows the EViews-style F and likelihood-ratio tests, F-test SSR decomposition, restricted/unrestricted log likelihoods, and the complete unrestricted equation. Date/time/sample fields are omitted.
