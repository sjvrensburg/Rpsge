% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{run_psge}
\alias{run_psge}
\title{Run PSGE Evolution}
\usage{
run_psge(psge_instance, verbose = TRUE, progress_fn = NULL)
}
\arguments{
\item{psge_instance}{PSGE instance from setup_psge()}

\item{verbose}{Print progress (default TRUE)}

\item{progress_fn}{Custom progress callback function (optional)}
}
\value{
Best solution found
}
\description{
Run PSGE Evolution
}
