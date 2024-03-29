% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbla_partial.R
\name{pbla_std}
\alias{pbla_std}
\title{PBLA for stochastic epidemic model}
\usage{
pbla_std(r, beta, gamma, m = 1, A = 1, lag = 0)
}
\arguments{
\item{r}{numeric vector: removal times}

\item{beta}{matrix of rates}

\item{gamma}{numeric vector of rates}

\item{m}{positive integer shape}

\item{A}{integer patient zeros}

\item{lag}{numeric fixed lag}
}
\value{
negative log likelihood
}
\description{
Compute pair-based likelihood approximation. Supports Erlang infectious periods.
}
