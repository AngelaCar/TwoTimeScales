#' Make a grid of points to evaluate B-splines
#'
#' @inheritParams get_hazard_2d
#' @param model_specifications A list with `'umin'`, `'umax'`, `'smin'`, and `'smax'` used for fitting.
#' @param class_fitmodel The class of the fitted model. Can be `'haz2ts'`, `'haz2tsPGAM'`, or `'haz2tsVCM'`.
#'
#' @return A list of specification for the grid
#' @export
#'
#' @examples
#' make_grid(plot_grid = list(
#'     c(umin = 3, umax = 8.5, du = .1),
#'     c(smin = 0, smax = 7.1, ds = .1)
#'   ),
#'   model_specification = list(umin = 2, umax = 8.5, smin = 0, smax = 7.5),
#'   class_fitmodel = "haz2ts")
#'
make_grid <- function(plot_grid,
                      model_specifications,
                      class_fitmodel,
                      where_slices = NULL,
                      direction = c(NULL, "u", "s"),
                      tmax = NULL,
                      midpoints = FALSE){

  # parameters of desired new grid
  umin <- plot_grid[[1]][1]
  umax <- plot_grid[[1]][2]
  du <- plot_grid[[1]][3]
  smin <- plot_grid[[2]][1]
  smax <- plot_grid[[2]][2]
  ds <- plot_grid[[2]][3]
  if (du <= 0) stop("`du` should be a positive number!")
  if (ds <= 0) stop("`ds` should be a positive number!")

  if(class_fitmodel != "haz2tsVCM"){
    if (umin < model_specifications$umin) {
      umin <- model_specifications$umin
      warning("`umin` is smaller than the lower limit of the domain of Bu. Left boundary adjusted to  =  ", model_specifications$umin)
    }
    if (umax > model_specifications$umax) {
      umax <- model_specifications$umax
      warning("`umax` is larger than the upper limit of the domain of Bu. Right boundary adjusted to  =  ", model_specifications$umax)
    }
    K <- ceiling((umax - umin) / du)

    intu <- seq(umin, umin + K * du, by = du)

    if (smin < model_specifications$smin) {
      smin <- model_specifications$smin
      warning("`smin` is smaller than the lower limit of the domain of Bs. Left boundary adjusted to  =  ", model_specifications$smin)
    }
    if (smax > model_specifications$smax) {
      smax <- model_specifications$smax
      warning("`smax` is larger than the upper limit of the domain of Bs. Right boundary adjusted to  =  ", model_specifications$smax)
    }
    K <- ceiling((smax - smin) / ds)
    ints <- seq(smin, smin + K * ds, by = ds)

    if (midpoints) {
      midu <- intu[-1] - du / 2
      mids <- ints[-1] - ds / 2
      new_grid <- list(
        "intu" = midu,
        "umin" = umin,
        "umax" = umax,
        "du" = du,
        "ints" = mids,
        "smin" = smin,
        "smax" = smax,
        "ds" = ds
      )
    } else {
      new_grid <- list(
        "intu" = intu,
        "umin" = umin,
        "umax" = umax,
        "du" = du,
        "ints" = ints,
        "smin" = smin,
        "smax" = smax,
        "ds" = ds
      )
    }

    # Adjust grid if slices are required
    if (!is.null(where_slices)) {
      if (is.null(direction)) stop("Direction for cutting slices missing.")
      if (direction == "u") {
        if (min(where_slices) < model_specifications$umin | max(where_slices) > model_specifications$umax) {
          stop("Desired cutting points outside of range of `B_u`.")
        } else {
          newu <- unique(sort(c(new_grid$intu, where_slices)))
          new_grid$intu <- newu
        }
      }
      if (direction == "s") {
        if (min(where_slices) < model_specifications$smin | max(where_slices) > model_specifications$smax) {
          stop("Desired cutting points outside of range of `B_s`.")
        } else {
          news <- unique(sort(c(new_grid$ints, where_slices)))
          new_grid$ints <- news
          grid_us <- expand.grid(u = new_grid$intu, s = new_grid$ints)
          grid_us$t <- with(grid_us, u + s)
          t <- unique(grid_us$t)
          if (!is.null(tmax)) {
            new_grid$intt <- t[t <= tmax]
          } else {
            new_grid$intt <- t
          }
          grid_ts <- expand.grid(t = new_grid$intt, s = new_grid$ints)
          grid_ts$u <- grid_ts$t - grid_ts$s
          newu <- unique(grid_ts$u)
          newu <- newu[newu >= new_grid$umin & newu <= new_grid$umax]
          new_grid$intu <- sort(newu)
        }
      }
    }
  }
  return(new_grid)
}


