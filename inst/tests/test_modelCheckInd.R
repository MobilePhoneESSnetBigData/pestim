context("Model Check")

test_that("modelCheckInd", {
  set.seed(1)
  actual_result<-modelCheckInd(nSimPar = 10, nMNO = 29, nReg = 123, fu = list('triang', xMin = 0.2, xMax = 0.25, xMode = 0.22), fv = list('triang', xMin = 115, xMax = 130, xMode = 120), flambda = list('gamma', shape = 21, scale = 123 / 20))
  expect_equal(actual_result[[1]], 29, tolerance = 1e-15, scale = 1)
  expect_equal(actual_result[[2]], 123, tolerance = 1e-15, scale = 1)
  expect_equal(actual_result[[3]], 4.35, tolerance = 1e-15, scale = 1)
  expect_equal(actual_result[[4]], 0.15, tolerance = 1e-15, scale = 1)
  expect_equal(actual_result[[5]], 104.9275, tolerance = 1e-15, scale = 1)
  expect_equal(actual_result[[6]], 0.124765160523187, tolerance = 1e-15, scale = 1)
  expect_equal(actual_result[[7]], 123.85, tolerance = 1e-15, scale = 1)
  expect_equal(actual_result[[8]], 0.147265160523187, tolerance = 1e-15, scale = 1)

  expect_is(actual_result, "data.table")
  expect_type(actual_result, "list")
  expect_type(actual_result[[1]], "double")
  expect_type(actual_result[[2]], "double")
  expect_type(actual_result[[3]], "double")
  expect_type(actual_result[[4]], "double")
  expect_type(actual_result[[5]], "double")
  expect_type(actual_result[[6]], "double")
  expect_type(actual_result[[7]], "double")
  expect_type(actual_result[[8]], "double")
  expect_length(actual_result, 8)
})



