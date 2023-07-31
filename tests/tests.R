library(testthat)
library(IPTWtrim)

# Test dataGenFun function
# Test 1: Test that dataGenFun returns an error with negative sample size
test_that("dataGenFun throws error with negative n", {
  expect_error(dataGenFun(n = -5,  mu = c(10, 1, 11), sigma.mat = matrix(c(11^2, 4, 4, 4, 3^2, 3, 4, 3, 8^2), 3),
                          eta = c(-1.5,0.1, 0.3, 0.2), theta =  c(110, -12, -0.3, -0.8, -0.2)), "sample size must be a positive integer")
})

# Test 2: Test that dataGenFun returns an error with invalid mu
test_that("dataGenFun throws error with invalid mu", {
  expect_error(dataGenFun(n = 500,  mu = c(10, 1, 11,20), sigma.mat = matrix(c(11^2, 4, 4, 4, 3^2, 3, 4, 3, 8^2), 3),
                          eta = c(-1.5,0.1, 0.3, 0.2), theta =  c(110, -12, -0.3, -0.8, -0.2)), "the mean vector must be a numeric vector of length 3")
})

# Test 3: Test that dataGenFun returns an error with invalid covariance matrix
test_that("dataGenFun throws error with invalid covariance matrix", {
  expect_error(dataGenFun(n = 500,  mu = c(10, 1, 11), sigma.mat = matrix(c(11^2, 4, 4, 4, 3^2, 3, -4, -3, 8^2), 3),
                          eta = c(-1.5,0.1, 0.3, 0.2), theta =  c(110, -12, -0.3, -0.8, -0.2)), "the covariance matrix sigma must be a numeric symmetric 3 by 3 matrix")
})

# Test 4: Test that dataGenFun returns an error with invalid eta
test_that("dataGenFun throws error with invalid eta", {
  expect_error(dataGenFun(n = 500,  mu = c(10, 1, 11), sigma.mat = matrix(c(1, 4, 4, 4, 1, 3, 4, 3, 1), 3),
                          eta = 1.5, theta =  c(110, -12, -0.3, -0.8, -0.2)),"eta must be a numeric vector of length 4")
})

# Test 5: Test that dataGenFun returns an error with invalid theta
test_that("dataGenFun throws error with invalid eta", {
  expect_error(dataGenFun(n = 500,  mu = c(10, 1, 11), sigma.mat = matrix(c(1, 4, 4, 4, 1, 3, 4, 3, 1), 3),
                          eta = c(0,0.1, 0.3, 1), theta =  c(110, -12, 200)),"theta must be a numeric vector of length 5")
})

# Test 6: Test that dataGenFun returns a data frame
test_that("dataGenFun returns a data frame", {
  data <- dataGenFun()
  expect_true(is.data.frame(data))
})

# Test 7: Test that dataGenFun generates data with correct column names
test_that("dataGenFun generates data with correct column names", {
  data <- dataGenFun()
  expect_equal(names(data), c("x1", "x2", "x3", "A", "Y"))
})


# Test measureATEfun function 
# Test 1: Test that measureATEfun returns an error with negative nsim
test_that("measureATEfun throws error with negative nsim", {
  expect_error(measureATEfun(nsim = -5), "the number of iterations must be a positive integer")
})

# Test 2: Test that measureATEfun returns an error with invalid trim.p
test_that("measureATEfun throws error with invalid trim.p", {
  expect_error(measureATEfun(trim.p = 1.1), "the trimming percentile must be a number between 0.5 and 1")
})

# Test 3: Test that measureATEfun returns an error with negative n
test_that("measureATEfun throws error with negative n", {
  expect_error(measureATEfun(n = -500), "the sample size must be a positive integer")
})

# Test 4: Test that measureATEfun returns an error with invalid mu
test_that("measureATEfun throws error with invalid mu", {
  expect_error(measureATEfun(mu = c(10, 1, 11, 20)), "the mean vector must be a numeric vector of length 3")
})

# Test 5: Test that measureATEfun returns an error with invalid covariance matrix
test_that("measureATEfun throws error with invalid covariance matrix", {
  expect_error(measureATEfun(sig.mat = matrix(c(11^2, 4, 4, 4, 3^2, 3, -4, -3, 8^2), 3)),
               "the covariance matrix sigma must be a numeric symmetric 3 by 3 matrix")
})

# Test 6: Test that measureATEfun returns an error with invalid eta
test_that("measureATEfun throws error with invalid eta", {
  expect_error(measureATEfun(eta = 1.5), "eta must be a numeric vector of length 4")
})

# Test 7: Test that measureATEfun returns an error with invalid theta
test_that("measureATEfun throws error with invalid theta", {
  expect_error(measureATEfun(theta = c(110, -12, -0.3)), "theta must be a numeric vector of length 5")
})

# Test 8: Test that measureATEfun returns a numeric matrix
test_that("measureATEfun returns a numeric matrix", {
  result <- measureATEfun()
  expect_true(is.matrix(result))
})

# Test 9: Test that measureATEfun returns a matrix with 1 row and 7 columns
test_that("measureATEfun returns a matrix with 1 row and 7 columns", {
  result <- measureATEfun()
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 7)
})

# Test 10: Test that measureATEfun calculates the bias and mse correctly
test_that("measureATEfun calculates the bias and mse correctly", {
  result <- measureATEfun(nsim = 100, n = 500)
  expect_true(!is.na(result[1, "bias"]))
  expect_true(!is.na(result[1, "mse"]))
})

# Test trimSummary function
# Test 1: Test that trimSummary returns an error with negative nsim
test_that("trimSummary throws error with negative nsim", {
  expect_error(trimSummary(nsim = -5), "the number of iterations must be a positive integer")
})

# Test 2: Test that trimSummary returns an error with negative n
test_that("trimSummary throws error with negative n", {
  expect_error(trimSummary(n = -500), "the sample size must be a positive integer")
})

# Test 3: Test that trimSummary returns an error with invalid mu
test_that("trimSummary throws error with invalid mu", {
  expect_error(trimSummary(mu = c(10, 1, 11, 20)), "the mean vector must be a numeric vector of length 3")
})

# Test 4: Test that trimSummary returns an error with invalid covariance matrix
test_that("trimSummary throws error with invalid covariance matrix", {
  expect_error(trimSummary(sig.mat = matrix(c(11^2, 4, 4, 4, 3^2, 3, -4, -3, 8^2), 3)),
               "the covariance matrix sigma must be a numeric symmetric 3 by 3 matrix")
})

# Test 5: Test that trimSummary returns an error with invalid eta
test_that("trimSummary throws error with invalid eta", {
  expect_error(trimSummary(eta = 'happy'), "eta must be a numeric vector of length 4")
})

# Test 6: Test that trimSummary returns an error with invalid theta
test_that("trimSummary throws error with invalid theta", {
  expect_error(trimSummary(theta = "c(110, -12, -0.3)"), "theta must be a numeric vector of length 5")
})

# Test 7: Test that trimSummary returns an error with invalid trimperc.try
test_that("trimSummary throws error with invalid trimperc.try", {
  expect_error(trimSummary(trimperc.try = c(1.1, 0.5, -0.9)), 
               "trimperc.try must be a numeric list, each number in the list must between 0.5 to 1")
})

# Test 8: Test that trimSummary returns a data frame
test_that("trimSummary returns a data frame", {
  result <- trimSummary()
  expect_true(is.data.frame(result))
})

# Test 9: Test that trimSummary returns a data frame with expected number of rows and columns
test_that("trimSummary returns a data frame with expected number of rows and columns", {
  result <- trimSummary()
  expect_equal(nrow(result), length(c(1, 0.99, 0.95, 0.9, 0.75, 0.5)))
  expect_equal(ncol(result), 10)
})


