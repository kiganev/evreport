\name{ev_bg}
\alias{ev_bg}
\title{Breusch-Godfrey serial correlation LM test}
\description{Breusch-Godfrey serial correlation LM test}
\usage{ev_bg(model, lags = 1L)}
\arguments{
\item{model}{An unweighted, full-rank lm fit.}
\item{lags}{Positive integer lag order.}
}
\value{An evreport list with statistics, display table, and (where available) auxiliary model.}
