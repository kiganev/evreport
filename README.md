# evreport

**EViews-style regression reporting and diagnostic tests for R**

`evreport` is an R package that provides EViews-style output for linear regression models and econometric diagnostic tests.

The package aims to make R's econometric output more familiar to EViews users while retaining the flexibility of R.

## Features

- Comprehensive OLS regression summaries
- Breusch–Godfrey serial correlation LM test
- Breusch–Pagan–Godfrey and White heteroskedasticity tests
- Ramsey RESET specification test
- Jarque–Bera normality test

## Example

```r
library(evreport)

model <- lm(mpg ~ wt + hp, data = mtcars)

ev_summary(model)
ev_bg(model, lags = 2)
ev_hetero(model, type = "white")
ev_reset(model)
ev_normality(model)
```

## License

MIT License.