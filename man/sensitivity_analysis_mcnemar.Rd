% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sensval.R
\name{sensitivity_analysis_mcnemar}
\alias{sensitivity_analysis_mcnemar}
\title{Performs sensitivity analysis for pair-matched binary data.}
\usage{
sensitivity_analysis_mcnemar(mpdifs, gamma = gamma, alt = alt)
}
\arguments{
\item{mpdifs}{vector of length nsets, where entry k contains the treated
minus control value corresponding to the kth matched pair.}

\item{gamma}{value of the sensitivity parameter.}

\item{alt}{character denoting direction of alternative hypothesis.}
}
\value{
p-value yielded by sensitivity analysis.
}
\description{
This function is internal to other functions and performs sensitivity
analysis for pair-matched binary data using McNemar's test statistic.
}
