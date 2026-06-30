# Tests for function fit1ts

dt <- prepare_data(
  data = reccolon2ts,
  s_out = "timesr",
  events = "status",
  ds = 30
)
m  <- fit1ts(data1ts = dt,
             Bbases_spec = list(bdeg = 3, nseg_s = 10,
                                min_s = 0, max_s = 2730))
mLMM <- fit1ts(data1ts = dt,
               Bbases_spec = list(bdeg = 3, nseg_s = 10,
                                  min_s = 0, max_s = 2730),
               optim_method = "LMMsolver")

# ----
test_that("fit1ts returns haz1ts object with required components", {
  expect_s3_class(m, "haz1ts")
  expect_s3_class(mLMM, "haz1tsLMM")
  expect_true(!is.null(m$optimal_model$alpha))
  expect_true(!is.null(m$optimal_model$SE_alpha))
  expect_true(!is.null(m$optimal_model$aic))
  expect_true(is.finite(m$optimal_model$aic))
})

# ----
test_that("fit1ts Alpha vector dimensions match B-spline spec", {
  cs <- 10 + 3  # nseg + bdeg for s
  expect_equal(length(m$optimal_model$alpha), cs)
})

# ----
test_that("fit1ts ED is within valid range", {
  n_coef <- length(m$optimal_model$alpha)
  expect_true(m$optimal_model$ed > 0)
  expect_true(m$optimal_model$ed < n_coef)
})

