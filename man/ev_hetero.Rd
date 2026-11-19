\name{ev_hetero}
\alias{ev_hetero}
\title{Breusch-Pagan-Godfrey or White heteroskedasticity test}
\description{Breusch-Pagan-Godfrey or White heteroskedasticity test}
\usage{ev_hetero(model, type = c("bpg", "white"), cross = FALSE)}
\arguments{
\item{model}{An unweighted, full-rank lm fit.}
\item{type}{Either bpg or white.}
\item{cross}{Include pairwise products for the White test.}
}
\value{An evreport list with statistics, display table, and (where available) auxiliary model.}
