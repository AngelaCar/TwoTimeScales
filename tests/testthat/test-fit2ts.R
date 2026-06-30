# Tests for function fit2ts

dt <- prepare_data(
  data = reccolon2ts,
  u      = "timer",
  s_out = "timesr",
  events = "status",
  ds = 30
)
m  <- fit2ts(data2ts = dt,
             Bbases_spec = list(bdeg = 3, nseg_s = 10, nseg_u = 8,
                                min_s = 0, max_s = 2730,
                                min_u = 0, max_u = 2300))
mLMM <- fit2ts(data2ts = dt,
               Bbases_spec = list(bdeg = 3, nseg_s = 10, nseg_u = 8,
                                  min_s = 0, max_s = 2730,
                                  min_u = 0, max_u = 2300),
               optim_method = "LMMsolver")

test_that("fit2ts returns correct model class", {
  expect_s3_class(m, "haz2ts")
  expect_s3_class(mLMM, "haz2tsLMM")
})


# ----
test_that("fit2ts Alpha matrix dimensions match B-spline spec", {
  cs <- 10 + 3  # nseg + bdeg for s
  cu <-  8 + 3  # nseg + bdeg for u
  expect_equal(dim(m$optimal_model$Alpha), c(cu, cs))
})

# ----
test_that("smoothing parameters for standard example are the expected ones", {
  expect_equal(round(m$optimal_logrho, 2), c(1.64, -0.33))
})
