% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_hr.R
\name{get_hr}
\alias{get_hr}
\title{Get the Hazard Ratios with their Standard Errors}
\usage{
get_hr(fitted_model)
}
\arguments{
\item{fitted_model}{A list returned by the function \code{fit2ts} or \code{fit1ts}.}
}
\value{
A list with the following elements:
\itemize{
\item \code{HR} A vector of hazard ratios (calculated as exp(hat_beta)).
\item \code{SE_HR} A vector of Standard Errors for the hazard ratios calculated
via the delta method.
\item \code{beta} A vector of the estimated beta coefficients.
\item \code{SE_beta} A vector of the Standard Errors for the beta coefficients.
}
}
\description{
\code{get_hr()} takes as input the results of a model with covariates
estimated by \code{fit2ts} or \code{fit1ts} and returns the estimated hazard ratios
together with their standard errors.
}
