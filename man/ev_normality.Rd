\name{ev_normality}
\alias{ev_normality}
\title{Jarque-Bera residual normality test and histogram}
\description{Produces a Jarque-Bera test and, by default, an EViews-inspired
histogram with descriptive statistics alongside it.}
\usage{ev_normality(model, plot = TRUE, breaks = "Sturges")}
\arguments{
\item{model}{An unweighted, full-rank \code{lm} fit.}
\item{plot}{Logical; whether to draw the histogram and statistics panel.}
\item{breaks}{Histogram bins or break specification accepted by \code{graphics::hist}.}
}
\value{An invisible \code{evreport} object containing descriptive statistics and
Jarque-Bera test results.}
\details{Bin choices may differ from EViews, so the histogram's bar heights and
boundaries may not match exactly. The numerical statistics do not depend on \code{breaks}.}
\examples{
fit <- lm(mpg ~ wt + hp, mtcars)
ev_normality(fit)
ev_normality(fit, plot = FALSE)
}
