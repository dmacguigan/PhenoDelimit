% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_clumpp_results.R
\name{read_clumpp_results}
\alias{read_clumpp_results}
\title{Summarize CLUMPP results}
\usage{
read_clumpp_results(
  wd,
  perc.var,
  model.numbers,
  apriori = FALSE,
  clust.method = "kmeans"
)
}
\arguments{
\item{wd}{working directory to store CLUMPP files, should contain a CLUMPP_permuted directory created by the dapc_clumpp_permuted function}

\item{perc.var}{vector containing cumulative percentages of variance to retain for discriminant analyses}

\item{model.numbers}{vector containing delimitation model numbers}
}
\description{
Read in results from CLUMPP miscfiles, return data frame of H' values for each delimitation model
}
