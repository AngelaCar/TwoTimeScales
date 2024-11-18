## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(TwoTimeScales)

## ----dependson="mod-paper"----------------------------------------------------
summary(mod1)

## ----grid-search, dependson="bindata_nocov", cache=TRUE, eval=F---------------
#  par(mfrow = c(1,2))
#  
#  mod3 <- fit2ts(data2ts = dt2ts,
#                 Bbases_spec = list(bdeg = 3,
#                                    nseg_s = 20,
#                                    min_s = 0,
#                                    max_s = 2730,
#                                    nseg_u = 20,
#                                    min_u = 0,
#                                    max_u = 2300),
#                 optim_method = "grid_search",
#                 optim_criterion = "aic",
#                 lrho = list(seq(-1, 3, by=.2),
#                             seq(-1, 3, by = .2)),
#                 par_gridsearch = list(
#                   plot_aic = T,
#                   plot_bic = T,
#                   mark_optimal = T,
#                   plot_contour = T
#                 ))

## ----eval = F-----------------------------------------------------------------
#  mod_cov <- fit2ts(data2ts = dt2ts_cov,
#                    Bbases_spec = list(bdeg = 3,
#                                       nseg_s = 20,
#                                       min_s = 0,
#                                       max_s = 2730,
#                                       nseg_u = 20,
#                                       min_u = 0,
#                                       max_u = 2300),
#                    pord = 2,
#                    optim_method = "ucminf",
#                    optim_criterion = "aic")

