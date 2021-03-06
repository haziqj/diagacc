% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/study_parallel.R
\name{run_study_par}
\alias{run_study_par}
\title{Title}
\usage{
run_study_par(
  object = NULL,
  B = 4,
  n = c(250, 1000),
  tau = c(0.08, 0.4),
  miss.prop = c(0.5, 0.8, 1),
  no.cores = parallel::detectCores(),
  lc.method = c("MCMC", "EM"),
  lcre.method = c("MCMC", "EM")
)
}
\arguments{
\item{object}{diagaccSim2 object}

\item{B}{Number of replications}

\item{n}{Setting for sample size}

\item{tau}{Setting for prevalence}

\item{miss.prop}{Setting for missing proportion}

\item{no.cores}{Number of cores to use.}

\item{lc.method}{Method for fitting LC models}

\item{lcre.method}{Method for fitting LCRE models}
}
\description{
Title
}
