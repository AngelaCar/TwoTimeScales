## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(TwoTimeScales)

## ----1ts-data-prep, cache=T---------------------------------------------------
dt1ts <- prepare_data(s_out = reccolon2ts$timesr,
                      events = reccolon2ts$status,
                      ds = 30)
str(dt1ts)
print(dt1ts)

## ----1ts-lt, cache=T----------------------------------------------------------
dt1ts_lt <- prepare_data(s_in = reccolon2ts$entrys,
                         s_out = reccolon2ts$timesr,
                         events = reccolon2ts$status,
                         ds = 30)


## -----------------------------------------------------------------------------
dt1ts_lt2 <- prepare_data(s_in = reccolon2ts$entrys,
                          s_out = reccolon2ts$timesr,
                          events = reccolon2ts$status,
                          ds = 90)
dt1ts_lt2$bins$ns

## -----------------------------------------------------------------------------
range(reccolon2ts$timesr)
dt1ts_2 <- prepare_data(s_in = reccolon2ts$entrys,
                        s_out = reccolon2ts$timesr,
                        events = reccolon2ts$status,
                        ds = 30, min_s = 0, max_s = 3000)

str(dt1ts_2)

## ----1ts-data-cov, cache=T----------------------------------------------------
covs <- subset(reccolon2ts, select = c("rx", "node4", "sex"))
dt1ts_cov <- prepare_data(s_in = reccolon2ts$entrys,
                          s_out = reccolon2ts$timesr,
                          events = reccolon2ts$status,
                          ds = 30,
                          individual = TRUE, 
                          covs = covs)
print(dt1ts_cov)

## ----fit1ts, dependson="1ts-data-prep"----------------------------------------
# Model 1 - Default parameters (numerical optimization of aic, default param for
#                               B-splines)
m1 <- fit1ts(data1ts = dt1ts)
# Model 2 - Single data inputs
m2 <- fit1ts(y = dt1ts$bindata$y, r = dt1ts$bindata$r, bins = dt1ts$bins)

table(m1$optimal_model$alpha == m2$optimal_model$alpha)
str(m1)

## ----fit1ts-other, eval = F---------------------------------------------------
#  # Model 3 - Change specifications of the B-splines (degree, number of segments
#  #                                                   and range)
#  m3 <- fit1ts(data1ts = dt1ts,
#               Bbases_spec = list(bdeg = 2,          # quadratic B-splines
#                                  nseg_s = 20,       # 20 segments
#                                  min_s = 0,
#                                  max_s = 2730))
#  
#  # Model 4 - As m3, but change penalty order
#  m4 <- fit1ts(data1ts = dt1ts,
#               Bbases_spec = list(bdeg = 2,
#                                  nseg_s = 20,
#                                  min_s = 0,
#                                  max_s = 2730),
#               pord = 3)                             # third-degree penalty
#  
#  # Model 5 - As m3, but change optimization method to grid_search
#  m5 <- fit1ts(data1ts = dt1ts,
#               Bbases_spec = list(bdeg = 2,
#                                  nseg_s = 20,
#                                  min_s = 0,
#                                  max_s = 2730),
#               optim_method = "grid_search")         # search for optimal smoothing                                                                       over grid of values
#  
#  # Model 6 - As m5, but optimization criterion is "bic" and include grid of
#  #           values for log_10(rho)
#  m6 <- fit1ts(data1ts = dt1ts,
#               Bbases_spec = list(bdeg = 2,
#                                  nseg_s = 20,
#                                  min_s = 0,
#                                  max_s = 2730),
#               optim_method = "grid_search",
#               optim_criterion = "bic",              # use BIC rather than AIC
#               lrho = seq(-2, 3, by=.2))             # provide grid for log_10(rho)
#  

## ----fit1ts-grid, dependson="1ts-data-prep", fig.width=6, fig.align='center'----
par(mfrow = c(1,2))
m6 <- fit1ts(data1ts = dt1ts,
             Bbases_spec = list(bdeg = 2,
                                nseg_s = 20,
                                min_s = 0,
                                max_s = 2730),
             optim_method = "grid_search",
             optim_criterion = "bic",
             lrho = seq(-1, 3, by=.2),
             par_gridsearch = list(
               plot_aic = T,
               plot_bic = T,
               mark_optimal = T
             ))
par(mfrow = c(1,1))
m6.aic <- m6$AIC
m6.bic <- m6$BIC

m6.aic[1:6]; m6.bic[1:6]

## ----fit1ts-cov, dependson="1ts-data-cov", cache = T--------------------------
m7 <- fit1ts(data1ts = dt1ts_cov,
             Bbases_spec = list(nseg_s = 15,
                                min_s = 0,
                                max_s = 2730))
betas <- m7$optimal_model$beta
betas


## ----get-haz-1d, dependson = "fit1ts-cov"-------------------------------------
basehaz <- get_hazard_1d(fitted_model = m7,
                         plot_grid = c("smin" = 0, "smax" = 2730, "ds" = 10))

str(basehaz)
range(basehaz$hazard)

## ----get-hr-1d, dependson = "fit1ts-cov"--------------------------------------
hr <- get_hr(fitted_model = m7)

hr$HR

## ----plot_haz1ts, fig.width = 5, dependson="fit1ts-cov", fig.align='center', fig.height=4----
plot(m7)


## ----plot_haz1ts-3, fig.width = 5, dependson="fit1ts-cov", fig.align='center', fig.height=4----
plot(m7,
     plot_grid = c("smin" = 0, "smax" = 2730, "ds" = 10),
     plot_options= list(
       col = "purple",
       main = "Baseline hazard",
       ylab = "hazard",
       xlab = "time since recurrence",
       cex_main = 1))


## ----plot_haz1ts-3log, fig.width = 5, dependson="fit1ts-cov", fig.align='center', fig.height=4----
plot(m7,
     plot_grid = c("smin" = 0, "smax" = 2730, "ds" = 10),
     plot_options= list(
       loghazard = TRUE,
       col = "purple",
       main = "Baseline hazard",
       ylab = "log-hazard",
       xlab = "time since recurrence",
       cex_main = 1))


## ----plot_haz1ts-4, fig.width = 5, dependson="fit1ts-cov", fig.align='center', fig.height=4----
plot(m7,
     which_plot = "covariates")


## ----CIs-lev, fig.width = 6, dependson="fit1ts-cov", fig.align='center', fig.height=4----
par(mfrow = c(1,2),
    font.main = 1)
# Option 1: symmetric CIs with delta method
plot(m7,
     which_plot = "covariates",
     plot_options = list(
       HR = T,
       symmetric_CI = T,
       ylim = c(0.8, 2.1),
       main = "HR with symmetric CIs",
       cex_main = 1
     ))
abline(h=1, lty=2, col="grey")
# Option 2: non-symmetric CIs 
plot(m7,
     which_plot = "covariates",
     plot_options = list(
       HR = T,
       symmetric_CI = F,
       ylim = c(0.8, 2.1),
       main = "HR with non-symmetric CIs",
       cex_main = 1
     ))
abline(h=1, lty=2, col="grey")
par(mfrow = c(1,1))

## ----HRs-plot, fig.width = 6.5, dependson="fit1ts-cov", fig.align='center', fig.height=4----
par(mfrow = c(1,3),
    font.main = 1)
plot(m7,
     which_plot = "covariates",
     plot_options = list(
       ylim = c(-0.22, 0.8),
       main = "95% CIs",
       cex_main = 1
     ))
abline(h=0, lty=2, col="grey")
plot(m7,
     which_plot = "covariates",
     plot_options = list(
       ylim = c(-0.22, 0.8),
       confidence = .90,
       main = "90% CIs",
       cex_main = 1
     ))
abline(h=0, lty=2, col="grey")
plot(m7,
     which_plot = "covariates",
     plot_options = list(
       ylim = c(-0.22, 0.8),
       confidence = .99,
       main = "99% CIs",
       cex_main = 1
     ))
abline(h=0, lty=2, col="grey")
par(mfrow = c(1,1))

