context("Triangular distribution")


test_that("qtriang", {
  set.seed(1)
  q<-seq(from=0, to=1, by=0.05)
  mode = rnorm(1,15)
  xmin = rnorm(1,10)
  xmax = rnorm(1,25)
  expected_result = c(10.18364332422208, 11.89504468813396, 12.60393034373007, 13.14787743866012, 13.60644605204583,
                      14.01045311071489, 14.37570420407118, 14.73177566854184, 15.10182655910430, 15.48764530354077,
                      15.89143766931765, 16.31597734318609, 16.76483452062717, 17.24273844334560, 17.75618448456067,
                      18.31452385509277, 18.93210869132071, 19.63309897334713, 20.46460295410856, 21.54824003945533,
                      24.16437138758995)
  actual_result<-qtriang(q, xMin = xmin, xMax = xmax, xMode = mode)


  expect_is(actual_result, "numeric")
  expect_type(actual_result, "double")
  expect_length(actual_result, length(q))


  expect_error(qtriang(q, xMin = xmax, xMax = xmin, xMode = mode), 'xMax must be greater than xMin.')
  expect_error(qtriang(q, xMin = xmin, xMax = xmax, xMode = 2*xmax), 'xMode must be between xMin and xMax.')
  expect_error(qtriang(q, xMin = xmin, xMax = xmin, xMode = mode), 'xMode must be between xMin and xMax.')
  expect_error(qtriang(q, xMin = xmin, xMax = xmin, xMode = xmin), 'NAs are not allowed in subscripted assignments')
  expect_equal(actual_result, expected_result, tolerance = 1e-14, scale = 1)
})


test_that("ptriang", {
  set.seed(1)
  q<-seq(from=0, to=10, by=0.5)
  mode = rnorm(1,7.5)
  xmin = rnorm(1,1)
  xmax = rnorm(1,9)
  expected_result = c(0.00000000000000000, 0.00000000000000000, 0.00000000000000000, 0.00251969766858563,
                      0.01677854605943596, 0.04362561756186622, 0.08306091217587640, 0.13508442990146655,
                      0.19969617073863660, 0.27689613468738655, 0.36668432174771648, 0.46906073191962633,
                      0.58402536520311610, 0.71157822159818573, 0.84954213504240184, 0.95101606321399190,
                      0.99700163646283591, 1.00000000000000000, 1.00000000000000000, 1.00000000000000000,
                      1.00000000000000000)
  actual_result<-ptriang(q, xMin = xmin, xMax = xmax, xMode = mode)


  expect_is(actual_result, "numeric")
  expect_type(actual_result, "double")
  expect_length(actual_result, length(q))


  expect_error(ptriang(q, xMin = xmax, xMax = xmin, xMode = mode), 'xMax must be greater than xMin.')
  expect_error(ptriang(q, xMin = xmin, xMax = xmax, xMode = 2*xmax), 'xMode must be between xMin and xMax.')
  expect_error(ptriang(q, xMin = xmin, xMax = xmin, xMode = mode), 'xMode must be between xMin and xMax.')
  expect_error(ptriang(q, xMin = xmin, xMax = xmin, xMode = xmin), 'NAs are not allowed in subscripted assignments')
  expect_equal(actual_result, expected_result, tolerance = 1e-14, scale = 1)
})
