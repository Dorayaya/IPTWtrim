#' Summary of ATE Estimation Using Weight Trimming
#'
#' `trimSummary()` performs a simulation study to evaluate the performance of
#'  different weight trimming percentiles. It calculates bias, empirical standard error (EmpSE),
#'  and mean squared error (MSE) for each weight trimming percentile specified in `trimperc.try`.
#'  Additionally, the function computes summary statistics for weights, including mean of trimmed weights,
#'  standard error of mean of trimmed weights, average standard deviation of trimmed weights, average
#'  maximum trimmed weight, and average minimum trimmed weight. The function then generates three
#' plots using `ggplot2` to visualize the relationship between weight trimming percentile
#' and bias, EmpSE, and MSE.
#'
#' @param nsim The number of iterations to perform for each weight trimming percentile.
#' An integer. Default to 1000.
#' @param n The number of observations in each simulated dataset. An integer. Default to 500.
#' @param mu The mean vector of the multivariate normal distribution for confounders x1, x2, x3.
#' A numeric vector of length 3. Default to c(10, 1, 11).
#' @param sig.mat The covariance matrix of the multivariate normal distribution for confounders
#' x1, x2, x3. A 3x3 numeric matrix. Default to a specific matrix.
#' @param eta The coefficients of the logistic regression model for the binary treatment A.
#' A numeric vector of length 4, representing intercept and coefficients for x1, x2, x3.
#' Default to a specific numeric vector.
#' @param theta The coefficients of the linear regression model for the continuous outcome Y.
#' A numeric vector of length 5, representing intercept, coefficients for A, x1, x2, x3.
#' Default to a specific numeric vector.
#' @param trimperc.try A numeric vector specifying the weight trimming percentiles to evaluate.
#' Each value should be between 0.5 and 1 (both inclusive) to represent the desired percentile
#' for weight trimming.
#' @param seed An optional seed value for reproducibility of simulations. Default to 20739377.
#'
#' @return A data frame containing statistics (bias, EmpSE, and MSE) for each weight trimming
#' percentile in `trimperc.try`. The data frame also includes summary statistics for weights:
#' trim.perc (weight trimming percentile), trim.prop (weight trimming proportion),
#' weight.mean (mean of trimmed weights), weight.ese (standard error of mean of trimmed weights),
#' weight.asd (average standard deviation of trimmed weights),
#' weight.max (average maximum trimmed weight), weight.min (average minimum
#' trimmed weight), bias (difference between estimated and true ATE), ese (empirical standard
#' error of the estimator), and mse (mean squared error of the estimator).
#'
#' @importFrom ggplot2 ggplot
#' @importFrom gridExtra grid.arrange
#' @export

trimSummary <- function(nsim = 1000, n = 500, mu = c(10, 1, 11),
                        sig.mat = matrix(c(11^2, 4, 4,
                                           4, 3^2, 3,
                                           4, 3, 8^2), 3),
                        eta = c(-1.5, 0.1, 0.3, 0.2),
                        theta =  c(110, -12, -0.3, -0.8, -0.2),
                        trimperc.try = c(1, 0.99, 0.95, 0.9, 0.75, 0.5),
                        seed = 20739377) {
  
  # Assertions to ensure that inputs are in the correct format and within
  # the valid range.
  
  # nsim must be a positive integer 
  if (!is.numeric(nsim) || length(nsim) != 1 || nsim <= 0 || round(nsim) != nsim) {
    stop("the number of iterations must be a positive integer")
  }
  
  # n must be a positive integer
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || round(n) != n) {
    stop("the sample size must be a positive integer")
  }
  
  # mu must be a numeric vector of length 3
  if (!is.numeric(mu) || length(mu) != 3 ) {
    stop("the mean vector must be a numeric vector of length 3")
  }
  
  ## sigma.mat must be a numeric symmetric 3*3 matrix
  if (!is.numeric(sigma.mat) || !isSymmetric(sigma.mat) || !is.matrix(sigma.mat) ||
      ncol(sigma.mat) != 3 || nrow(sigma.mat) != 3) {
    stop("the covariance matrix sigma must be a numeric symmetric 3 by 3 matrix")
  }
  
  # eta must be a numeric vector of length 4
  if (!is.numeric(eta) || length(eta) != 4 ) {
    stop("eta must be a numeric vector of length 4")
  }
  
  # theta must be a numeric vector of length 5
  if (!is.numeric(theta) || length(theta) != 5 ) {
    stop("theta must be a numeric vector of length 5")
  }
  
  # trimperc.try must be a numeric list, each number in the list must between 0.5 to 1
  if (!is.numeric(trimperc.try) || !all(trimperc.try >= 0.5 & trimperc.try <= 1)){
    stop("trimperc.try must be a numeric list, each number in the list must between 0.5 to 1")
  }
  
  set.seed(seed)
  bias <- c()
  mse <- c()
  weight.mean <- c()
  weight.ese <- c()
  weight.asd <- c()
  weight.max <- c()
  weight.min <- c()

  for (i in 1:length(trimperc.try)) {
    measure.psi1 <- measureATEfun(nsim, trimperc.try[i], n, mu, sig.mat, eta, theta)
    bias[i] <- measure.psi1[1]
    mse[i] <- measure.psi1[2]
    weight.mean[i] <- measure.psi1[3]
    weight.ese[i] <- measure.psi1[4]
    weight.asd[i] <- measure.psi1[5]
    weight.max[i] <- measure.psi1[6]
    weight.min[i] <- measure.psi1[7]
  }

  ese <- sqrt(mse - bias^2)
  trim.proportion <- (1 - trimperc.try) * 2

  data <- data.frame(trim.perc = trimperc.try, trim.prop = trim.proportion, weight.mean = weight.mean,
                     weight.ese = weight.ese, weight.asd = weight.asd,
                     weight.max = weight.max, weight.min = weight.min,
                     bias = bias, ese = ese, mse = mse)

  # Create three plots using ggplot2
  plot_bias <- ggplot2::ggplot(data, ggplot2::aes(x = trim.prop, y = bias)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Weight Trimming Proportion", y = "Bias")

  plot_ese <- ggplot2::ggplot(data, ggplot2::aes(x = trim.prop, y = ese)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Weight Trimming Proportion", y = "EmpSE")

  plot_mse <- ggplot2::ggplot(data, ggplot2::aes(x = trim.prop, y = mse)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Weight Trimming Proportion", y = "MSE")

  # Arrange the three plots side by side using grid.arrange from gridExtra
  combined_plots <- gridExtra::grid.arrange(plot_bias, plot_ese, plot_mse, ncol = 3)

  return(data)
}
