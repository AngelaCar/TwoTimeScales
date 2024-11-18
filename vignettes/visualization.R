## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(TwoTimeScales)

## ----load-mod, cache=T, echo=F------------------------------------------------
load("modelcov.Rda")

## ----basehaz, fig = T, dependson="load-mod", fig.width=6, fig.align='center', fig.height=4.5----
plot(mod_cov)

## ----original, fig = T, dependson="load-mod", fig.width=6, fig.align='center', fig.height=9----
par(mfrow = c(2,1), 
    font.main = 1)
plot(mod_cov,
     plot_options = list(
       rectangular_grid = T,               # for grid of rectangles
       original = T,                       # for plot in (t,s)-plane
       main = "Original plane - rectangular grid",
       xlab = "Time since randomization",
       ylab = "Time since recurrence"
     ))
plot(mod_cov,
     plot_options = list(
       rectangular_grid = F,                # for grid of parallelograms
       original = T,                        # for plot in (t,s)-plane
       main = "Original plane - grid of parallelograms",
       xlab = "Time since randomization",
       ylab = "Time since recurrence"
     ))
par(mfrow = c(1,1))

## ----new-grid, fig = T, chache = T, dependson="load-mod", fig.width=6, fig.align='center', fig.height=4.5----
plot(mod_cov,
     plot_grid = list(c(umin=0, umax=2300, du=5),
                      c(umin=0, umax=2730, du=5)),
     plot_options = list(n_shades = 100,
                         main = "Denser plotting grid",
                         xlab = "Time at recurrence",
                         ylab = "Time since recurrence"))

## ----loghaz, fig = T, dependson="load-mod", fig.width=6, fig.align='center', fig.height=9----
par(mfrow = c(2,1), 
    font.main = 1)
plot(mod_cov,
     plot_options = list(
       loghazard = T,
       main = "Log-hazard (u,s)",
       xlab = "Time at recurrence",
       ylab = "Time since recurrence"
     ))
plot(mod_cov,
     plot_options = list(
       original = T,
       loghazard = T,
       main = "Log-hazard (t,s)",
       xlab = "Time since randomization",
       ylab = "Time since recurrence"
     ))
par(mfrow = c(1,1))

## ----log10haz, fig = T, dependson="load-mod", fig.width=6, fig.align='center', fig.height=9----
par(mfrow = c(2,1), 
    font.main = 1)
plot(mod_cov,
     plot_options = list(
       log10hazard = T,
       main = "Log10-hazard (u,s)",
       xlab = "Time at recurrence",
       ylab = "Time since recurrence"
     ))
plot(mod_cov,
     plot_options = list(
       original = T,
       log10hazard = T,
       main = "Log10-hazard (t,s)",
       xlab = "Time since randomization",
       ylab = "Time since recurrence"
     ))
par(mfrow = c(1,1))

## ----cut-extra, fig = T, dependson="load-mod", fig.width=6, fig.align='center', fig.height=9----
par(mfrow = c(2,1), 
    font.main = 1)
plot(mod_cov,
     plot_options = list(cut_extrapolated = T,
                         tmax = 3214,
                         main = "Cut extrapolated hazard",
                         xlab = "Time at recurrence",
                         ylab = "Time since recurrence"))

plot(mod_cov,
     plot_options = list(cut_extrapolated = T,
                         tmax = 3214,
                         original = T,
                         main = "Cut extrapolated hazard",
                         xlab = "Time since randomization",
                         ylab = "Time since recurrence"))

## ----colors, fig = T, dependson="load-mod", fig.width=6, fig.align='center', fig.height=4.5----
mycol <- function(nshades){
  colorspace::sequential_hcl(n=nshades, "Blues 3")
}
plot(mod_cov,
     plot_options = list(col_palette = mycol,
                         main = "New colors",
                         xlab = "Time at recurrence",
                         ylab = "Time since recurrence",
                         contour_col = "pink",
                         contour_nlev = 20))

## ----SE, fig = T, chache = T, dependson="load-mod", fig.width=6, fig.align='center', fig.height=4.5----
plot(mod_cov,
     which_plot = "SE",
     plot_options = list(main = "Standard Errors of the hazard",
                         xlab = "Time at recurrence",
                         ylab = "Time since recurrence"))

plot(mod_cov,
     which_plot = "SE",
     plot_options = list(
       loghazard = TRUE,
       main = "Standard Errors of the log-hazard",
       xlab = "Time at recurrence",
       ylab = "Time since recurrence"))

plot(mod_cov,
     which_plot = "SE",
     plot_options = list(
       log10hazard = TRUE,
       main = "Standard Errors of the log10-hazard",
       xlab = "Time at recurrence",
       ylab = "Time since recurrence"))

## ----slices-u, fig = T, chache = T, dependson="load-mod", fig.width=6, fig.align='center', fig.height=4.5----
plot(mod_cov,
     which_plot = "slices",
     where_slices = c(30, 60, 90, 180, 365, 1000, 2000),
     direction = "u",
     plot_options = list(main = "Slices of the hazard",
                         xlab = "Time since recurrence",
                         ylab = "Hazard"))
legend("topright",
       legend = c(30, 60, 90, 180, 365, 1000, 2000), 
       lty = 1, 
       col = grDevices::gray.colors(7))

## ----slices-s, fig = T, chache = T, dependson="load-mod", fig.width=6, fig.align='center', fig.height=4.5----
plot(mod_cov,
     which_plot = "slices",
     where_slices = c(30, 60, 90, 180, 365, 1000, 2000),
     direction = "s",
     plot_options = list(main = "Slices of the hazard",
                         xlab = "Time since randomization",
                         ylab = "Hazard"))
legend("topright",
       legend = c(30, 60, 90, 180, 365, 1000, 2000), 
       lty=1, 
       col= grDevices::gray.colors(7))

## ----cov, fig = T, chache = T, dependson="load-mod", fig.width=6, fig.align='center', fig.height=4.5----
plot(mod_cov,
     which_plot = "covariates",
     plot_options = list(main = "",
                         ylab = "betas"))

