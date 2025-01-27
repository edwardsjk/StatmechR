% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_epi_mdl.R
\name{fit_epi_mdl}
\alias{fit_epi_mdl}
\title{A function to fit a generic epi model given an epi model and a
curve generating function. Currently uses optim with the defualt
parameters}
\usage{
fit_epi_mdl(ecs, epi_mdl_func, epi_mdl_pars, error_func, ...)
}
\arguments{
\item{ecs}{a list of epidemic curves we are trying to fit.}

\item{epi_mdl_func}{the epidemic model we are trying to fit.
Should take a set of parameters and produce a curve of the given length}

\item{epi_mdl_pars}{a data frame containing the starting parameters
for the epidemic model}

\item{error_func}{the error function we are trying to minimze}

\item{...}{additional parameters to the objective function}
}
\value{
a data frame with the epidemic model parameters and the covergence result. # nolint
takes in the parameter set and the epidemic curve
}
\description{
A function to fit a generic epi model given an epi model and a
curve generating function. Currently uses optim with the defualt
parameters
}
