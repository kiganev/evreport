\name{ev_reset}
\alias{ev_reset}
\title{Ramsey RESET test with EViews-style report}
\description{Prints the F and likelihood-ratio tests, the nested-regression SSR and log-likelihood summaries, and the full unrestricted equation.}
\usage{ev_reset(model, powers = 2:3)}
\arguments{
\item{model}{An unweighted, full-rank \code{lm} model.}
\item{powers}{Distinct integer fitted-value powers, default \code{2:3}.}
}
\details{The LR statistic is \code{n * log(SSR_restricted / SSR_unrestricted)}; its p-value uses a chi-square distribution with degrees of freedom equal to the number of added powers. Assumes conventional OLS covariance and a full-rank augmented specification.}
\value{An invisibly returned \code{evreport} object containing test statistics, the unrestricted \code{lm} fit, its coefficients, the two summary tables, and unrestricted equation statistics.}
\examples{
fit <- lm(mpg ~ wt + hp, data = mtcars)
ev_reset(fit, powers = 2:3)
}
