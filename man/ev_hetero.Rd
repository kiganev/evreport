\name{ev_hetero}
\alias{ev_hetero}
\title{Breusch-Pagan-Godfrey or White heteroskedasticity test}
\description{Prints EViews-style heteroskedasticity test results, the squared-residual auxiliary regression and its summary statistics.}
\usage{ev_hetero(model, type = c("bpg", "white"), cross = FALSE)}
\arguments{
\item{model}{An unweighted, full-rank lm fit with an intercept.}
\item{type}{Either "bpg" or "white".}
\item{cross}{Include pairwise cross products for the White test. Default FALSE.}
}
\details{Returns the auxiliary-regression F statistic, Obs*R-squared and Scaled explained SS, with p-values. Scaled explained SS equals auxiliary explained sum of squares divided by twice the square of the original regression's error variance estimate. EViews' date/time and sample-period metadata are omitted.}
\value{An evreport object, including statistics, coefficients, details and auxiliary_model.}
\examples{ev_hetero(lm(mpg ~ wt + hp, data = mtcars), type = "white")}
