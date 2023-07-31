#' Measure Average Treatment Effect (ATE) Using Trimming
#'
#' `measureATEfun()` estimates the Average Treatment Effect (ATE) using trimming in the
#' weighting step to reduce the influence of extreme weights. It calculates the bias, MSE for estimated ATE,
#' and it gives summary statistics for weights in each iteration. The function generates `nsim`
#' datasets and estimates the ATE for each dataset using a Generalized Linear Model (GLM)
#' with binary treatment (A) and continuous confounders (x1, x2, x3). Trimming is applied
#' to the weights calculated based on the estimated propensity scores to reduce the impact
#' of extreme weights.
#'
#' @param nsim The number of simulations to perform. A positive integer.
#' @param trim.p The percentile of extreme weights to be trimmed on both ends.
#' A numeric value between 0.5 and 1. For example, 0.9 represents trimming the top 10%
#' and bottom 10% of the weights.
#' @param n The number of observations in each simulated dataset. A positive integer.
#' @param mu The mean vector of the multivariate normal distribution for confounders x1, x2, x3.
#' A numeric vector of length 3.
#' @param sig.mat The covariance matrix of the multivariate normal distribution for confounders
#' x1, x2, x3. A 3x3 numeric matrix.
#' @param eta The coefficients of the logistic regression model for the binary treatment A.
#' A numeric vector of length 4, representing intercept and coefficients for x1, x2, x3.
#' @param theta The coefficients of the linear regression model for the continuous outcome Y.
#' A numeric vector of length 5, representing intercept, coefficients for A, x1, x2, x3.
#'
#' @return A matrix containing measures related to ATE estimation using trimming. The
#' columns include bias (difference between estimated and true ATE), mean squared error (MSE)
#' of the estimator, mean of trimmed weights, empirical standard error of mean of trimmed weights,
#' average standard deviation of trimmed weights, average maximum trimmed weight, and average minimum
#' trimmed weight.
#'
#' @importFrom MASS mvrnorm
#' @export

measureATEfun <- function(nsim = 1000, trim.p = 0.95, n = 500,
                          mu = c(10, 1, 11), sig.mat = matrix(c(11^2, 4, 4,
                                                                4, 3^2, 3,
                                                                4, 3, 8^2), 3),
                          eta = c(-1.5,0.1, 0.3, 0.2), theta =  c(110, -12, -0.3, -0.8, -0.2)){
  psi1.trim <- c()
  mean.w.trim <- c()
  min.trim <- c()
  max.trim <- c()
  sd.w.trim <- c()
  for (r in 1:nsim){
    data <- dataGenFun(n, mu, sig.mat, eta, theta)  # Assuming dataGenFun is a function that generates data
    x1 <- data$x1
    x2 <- data$x2
    x3 <- data$x3
    A <- data$A
    Y <- data$Y
    expo.model <- glm(A ~ x1 + x2 + x3, data = data, family = binomial(link = "logit"))
    eta.hat <- expo.model$coefficients
    pi.hat <- (1 + exp(-(eta.hat[1] + x1 * eta.hat[2] + x2 * eta.hat[3] + x3 * eta.hat[4])))^(-1)
    w.un <- data$A / pi.hat + (1 - data$A) / (1 - pi.hat)
    upper_bound <- quantile(w.un, trim.p)
    lower_bound <- quantile(w.un, 1 - trim.p)
    w.trim <- pmax(pmin(w.un, upper_bound), lower_bound)
    mean.w.trim[r] <- mean(w.trim)
    sd.w.trim[r] <- sd(w.trim)
    max.trim[r] <- max(w.trim)
    min.trim[r] <- min(w.trim)
    ### trimmed weighting
    msm.trim <- lm(Y ~ A, data = data, weights = w.trim)
    psi1.trim[r] <- msm.trim$coefficients[2]
  }
  bias <- mean(psi1.trim) - theta[2]
  ese <- sd(psi1.trim)
  mse <- mean((psi1.trim + 12)^2)
  weight.mean <- mean(mean.w.trim)
  weight.ese <- sd(mean.w.trim)
  weight.asd <- mean(sd.w.trim)
  weight.max <- mean(max.trim)
  weight.min <- mean(min.trim)
  measure <- cbind(bias, mse, weight.mean, weight.ese, weight.asd, weight.max, weight.min)
  return(measure)
}
