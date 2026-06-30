# Tests for function prepare_data

# --- Conservation of events: 1ts case ---

test_that("prepare_data conserves total events (1ts), and calculates positive exposure", {
  dt <- prepare_data(
    data   = reccolon2ts,
    s_out  = "timesr",
    events = "status",
    ds     = 30
  )
  expect_equal(sum(dt$bindata$y), sum(reccolon2ts$status))
  # Exposure sum should be positive and finite
  expect_true(is.finite(sum(dt$bindata$r)) && sum(dt$bindata$r) > 0)
})

# --- Conservation of events: 2ts case ---

test_that("prepare_data conserves total events (2ts), and calculates positive exposure", {
  dt <- prepare_data(
    data   = reccolon2ts,
    u      = "timer",
    s_out  = "timesr",
    events = "status",
    ds     = 30
  )
  expect_equal(sum(dt$bindata$Y), sum(reccolon2ts$status))
  # Exposure sum should be positive and finite
  expect_true(is.finite(sum(dt$bindata$R)) && sum(dt$bindata$R) > 0)
})

# --- Output structure and class ---
test_that("prepare_data returns correct structure for 1ts case", {
  dt <- prepare_data(
    data = reccolon2ts,
    s_out = "timesr",
    events = "status",
    ds = 30
  )
  expect_s3_class(dt, "data2ts")
  expect_named(dt, c("bins", "bindata"))
  expect_true(all(c("bins_s", "mids", "ns") %in% names(dt$bins)))
  expect_true(all(c("r", "y") %in% names(dt$bindata)))
})

# ----
test_that("prepare_data returns correct structure for 2ts case", {
  dt <- prepare_data(
    data = reccolon2ts,
    u      = "timer",
    s_out = "timesr",
    events = "status",
    ds = 30
  )
  expect_s3_class(dt, "data2ts")
  expect_named(dt, c("bins", "bindata"))
  expect_true(all(c("bins_u", "midu", "nu", "bins_s", "mids", "ns") %in% names(dt$bins)))
  expect_true(all(c("R", "Y") %in% names(dt$bindata)))
})

# ---- Left truncation ----
test_that("left truncation reduces total exposure in 1d", {

  # Setup
  dt_no_trunc <- prepare_data(
    data   = reccolon2ts,
    s_out  = "timesr",
    events = "status",
    ds     = 30
  )

  dt_trunc <- prepare_data(
    data   = reccolon2ts,
    s_in   = "entrys",      # <-- entry times included
    s_out  = "timesr",
    events = "status",
    ds     = 30
  )

  # Expectation: truncation must reduce total exposure
  expect_lt(
    sum(dt_trunc$bindata$r),      # actual (with truncation)
    sum(dt_no_trunc$bindata$r)    # expected to be larger (without)
  )
})


test_that("left truncation reduces total exposure in 2d", {

  # Setup
  dt_no_trunc <- prepare_data(
    data   = reccolon2ts,
    u = "timer",
    s_out  = "timesr",
    events = "status",
    ds     = 30
  )

  dt_trunc <- prepare_data(
    data   = reccolon2ts,
    u = "timer",
    s_in   = "entrys",      # <-- entry times included
    s_out  = "timesr",
    events = "status",
    ds     = 30
  )

  # Expectation: truncation must reduce total exposure
  expect_lt(
    sum(dt_trunc$bindata$R),      # actual (with truncation)
    sum(dt_no_trunc$bindata$R)    # expected to be larger (without)
  )
})
