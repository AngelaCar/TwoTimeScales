% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_cif2ts_pointwise.R
\name{predict_cif2ts_pointwise}
\alias{predict_cif2ts_pointwise}
\title{Point-wise prediction of cumulative incidence over 2 time scale}
\usage{
predict_cif2ts_pointwise(fitted_models = list(), u, s)
}
\arguments{
\item{fitted_models}{a list with cause-specific hazard models}

\item{u}{The value(s) of \code{u} where prediction is required}

\item{s}{The value(s) of \code{s} where prediction is required}
}
\value{
A data.frame with one row containing:
the values of \code{u} and \code{s}for which predictions of the overall survival
(\code{surv}) probability, and the values of the cumulative incidence functions,
one for each cause, are obtained.
}
\description{
Point-wise prediction of cumulative incidence over 2 time scale
}
\examples{
id <- 1:20
u <- c(5.43, 3.25, 8.15, 5.53, 7.28, 6.61, 5.91, 4.94, 4.25, 3.86, 4.05, 6.86,
       4.94, 4.46, 2.14, 7.56, 5.55, 7.60, 6.46, 4.96)
s <- c(0.44, 4.89, 0.92, 1.81, 2.02, 1.55, 3.16, 6.36, 0.66, 2.02, 1.22, 3.96,
       7.07, 2.91, 3.38, 2.36, 1.74, 0.06, 5.76, 3.00)
ev <- c(1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1)#'
fakedata <- as.data.frame(cbind(id, u, s, ev))
# cause 1
fakedata2ts1 <- prepare_data(u = fakedata$u,
                            s_out = fakedata$s,
                            ev = fakedata$ev,
                            ds = .5)
# Fit a fake model - not optimal smoothing for cause type 1
fakemod1 <- fit2ts(fakedata2ts1,
                  optim_method = "grid_search",
                  lrho = list(seq(1 ,1.5 ,.5),
                              seq(1 ,1.5 ,.5)))
# cause 2
fakedata2ts2 <- prepare_data(u = fakedata$u,
                              s_out = fakedata$s,
                              ev = 1-(fakedata$ev),
                              ds = .5)
# Fit a fake model - not optimal smoothing for cause 2
fakemod2 <- fit2ts(fakedata2ts2,
                    optim_method = "grid_search",
                    lrho = list(seq(1 ,1.5 ,.5),
                                seq(1 ,1.5 ,.5)))
predict_cif2ts_pointwise(fitted_models = list(fakemod1, fakemod2),
                         u = 5.3, s = 4.4)
}
