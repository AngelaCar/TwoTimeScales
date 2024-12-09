% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/competing_risks_methods.R
\name{surv2ts}
\alias{surv2ts}
\title{Survival function with two time scales}
\usage{
surv2ts(cumhaz = list())
}
\arguments{
\item{cumhaz}{a list with all the cause-specific cumulated hazard matrices
(minimum one element needs to be supplied).
If more than one cause-specific cumulated hazard is provided,
then they should all be matrices of the same dimension.}
}
\value{
a matrix containing the values of the survival function over \code{s} and \code{u}.
}
\description{
Computes the survival matrix, that contains the probability of not
experiencing an event of any cause by time \code{s} and fixed time at entry \code{u},
from a list of cause-specific cumulated hazard matrices.
}
\examples{
# Create some fake data - the bare minimum
id <- 1:20
u <- c(5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
       4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96)
s <- c(0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
       7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00)
ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1)#'

fakedata <- as.data.frame(cbind(id, u, s, ev))
fakedata2ts <- prepare_data(u = fakedata$u,
                            s_out = fakedata$s,
                            ev = fakedata$ev,
                            ds = .5)
# Fit a fake model - not optimal smoothing
fakemod <- fit2ts(fakedata2ts,
                  optim_method = "grid_search",
                  lrho = list(seq(1 ,1.5 ,.5),
                              seq(1 ,1.5 ,.5)))

# Obtain the fake cumulated hazard
fakecumhaz2ts <- cumhaz2ts(fakemod)
# Fake survival curve
fakesurv2ts <- surv2ts(cumhaz = list(fakecumhaz2ts$CumHaz))

}
