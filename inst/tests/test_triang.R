context("Triangular distribution")


test_that("qtriang", {
  set.seed(1)
  q<-seq(from=0, to=1, by=0.05)
  mode = rnorm(1,15)
  xmin = rnorm(1,10)
  xmax = rnorm(1,25)
  expected_result = c(10.183643324222082, 11.895044688133956, 12.603930343730067, 13.147877438660121, 13.606446052045831, 14.010453110714888,
                      14.375704204071177, 14.731775668541843, 15.101826559104298, 15.487645303540774, 15.891437669317654, 16.315977343186088,
                      16.764834520627172, 17.242738443345601, 17.756184484560670, 18.314523855092769, 18.932108691320710, 19.633098973347128,
                      20.464602954108560, 21.548240039455333, 24.164371387589952)

  actual_result<-qtriang(q, xMin = xmin, xMax = xmax, xMode = mode)


  expect_is(actual_result, "numeric")
  expect_type(actual_result, "double")
  expect_length(actual_result, length(q))


  expect_error(qtriang(q, xMin = xmax, xMax = xmin, xMode = mode), 'xMax must be greater than xMin.')
  expect_error(qtriang(q, xMin = xmin, xMax = xmax, xMode = 2*xmax), 'xMode must be between xMin and xMax.')
  expect_error(qtriang(q, xMin = xmin, xMax = xmin, xMode = mode), 'xMode must be between xMin and xMax.')
  expect_error(qtriang(q, xMin = xmin, xMax = xmin, xMode = xmin), 'NAs are not allowed in subscripted assignments')
  expect_equal(actual_result, expected_result, tolerance = 1e-15, scale = 1)
})


test_that("ptriang", {
  set.seed(1)
  q<-seq(from=0, to=10, by=0.5)
  mode = rnorm(1,7.5)
  xmin = rnorm(1,1)
  xmax = rnorm(1,9)
  expected_result = c(0.000000000000000000, 0.000000000000000000, 0.000000000000000000, 0.002519697668585630, 0.016778546059435958,
                      0.043625617561866219, 0.083060912175876397, 0.135084429901466546, 0.199696170738636603, 0.276896134687386553,
                      0.366684321747716480, 0.469060731919626328, 0.584025365203116098, 0.711578221598185734, 0.849542135042401836,
                      0.951016063213991902, 0.997001636462835905, 1.000000000000000000, 1.000000000000000000, 1.000000000000000000,
                      1.000000000000000000)

  actual_result<-ptriang(q, xMin = xmin, xMax = xmax, xMode = mode)


  expect_is(actual_result, "numeric")
  expect_type(actual_result, "double")
  expect_length(actual_result, length(q))


  expect_error(ptriang(q, xMin = xmax, xMax = xmin, xMode = mode), 'xMax must be greater than xMin.')
  expect_error(ptriang(q, xMin = xmin, xMax = xmax, xMode = 2*xmax), 'xMode must be between xMin and xMax.')
  expect_error(ptriang(q, xMin = xmin, xMax = xmin, xMode = mode), 'xMode must be between xMin and xMax.')
  expect_error(ptriang(q, xMin = xmin, xMax = xmin, xMode = xmin), 'NAs are not allowed in subscripted assignments')
  expect_equal(actual_result, expected_result, tolerance = 1e-15, scale = 1)
})

