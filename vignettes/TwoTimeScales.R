## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(eval = TRUE)

## ----setup--------------------------------------------------------------------
library(TwoTimeScales)

## ----data---------------------------------------------------------------------
str(reccolon2ts)

## ----1tsprep, cache=T---------------------------------------------------------
dt1ts <- prepare_data(s_in = reccolon2ts$entrys,
                      s_out = reccolon2ts$timesr,
                      events = reccolon2ts$status,
                      ds = 30)
str(dt1ts)

print(dt1ts)

## ----2tsprep, cache=TRUE------------------------------------------------------
dt2ts <- prepare_data(u = reccolon2ts$timer,
                      s_in = reccolon2ts$entrys,
                      s_out = reccolon2ts$timesr,
                      events = reccolon2ts$status,
                      ds = 30)
str(dt2ts)

print(dt2ts)

## ----2tscov, cache=TRUE-------------------------------------------------------
covs <- subset(reccolon2ts, select = c("rx", "node4", "sex"))
dt2ts_cov <- prepare_data(u = reccolon2ts$timer,
                          s_in = reccolon2ts$entrys,
                          s_out = reccolon2ts$timesr,
                          events = reccolon2ts$status,
                          ds = 30,
                          individual = TRUE, 
                          covs = covs)
str(dt2ts_cov)
print(dt2ts_cov)

## ----1ts-model, dependson = "1tsprep", cache=TRUE-----------------------------
m1ts <- fit1ts(data1ts = dt1ts,
               Bbases_spec = list(bdeg = 3,
                                  nseg_s = 20,
                                  min_s = 0,
                                  max_s = 2730))
str(m1ts)

## ----2ts-mod, dependson = "2tsprep", cache=TRUE-------------------------------
m2ts <- fit2ts(data2ts = dt2ts,
               Bbases_spec = list(bdeg = 3,
                                  nseg_s = 20,
                                  min_s = 0,
                                  max_s = 2730,
                                  nseg_u = 16,
                                  min_u = 0,
                                  max_u = 2300))
summary(m2ts)

## ----1ts-plot, fig = T, dependson="1ts-mod", fig.width=5, fig.align='center', fig.height=4----
plot(x = m1ts,
     plot_grid = c("smin" = 0, "smax" = 2730, "ds" = 10),
     plot_options= list(
       col = "darkblue",
       main = "Hazard",
       ylab = "hazard",
       xlab = "time since recurrence",
       cex_main = 1))

## ----2ts-plot, cache = T, fig = T, dependson="2ts-mod", fig.width=6, fig.height=4.5, fig.align='center'----
plot(x = m2ts,
     plot_grid = list(c("umin" = 0, "umax" = 2300, "du" = 10),
                      c("smin" = 0, "smax" = 2730, "ds" = 10)),
     plot_options= list(
       main = "Bi-dimensional hazard",
       ylab = "time since recurrence",
       xlab = "time since randomization",
       cex_main = 1))

