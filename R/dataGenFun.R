#' Generate Data from a Generalized Linear Model with Binary Outcome
#'
#' `dataGenFun()` generates data from a Generalized Linear Model (GLM) with a binary treatment (A)
#' and confounding variables (x1, x2, x3). The confounding variables are assumed to be generated
#' from a multivariate normal distribution with a specified mean and covariance matrix. The binary
#' treatment A is modeled using a logistic link function with specified coefficients. The continuous
#' outcome Y is generated based on the treatment and confounders using a linear model with normally
#' distributed errors.
#'
#' @param n The number of observations to generate. A positive integer.
#' @param mu The mean vector of the multivariate normal distribution for predictors x1, x2, x3.
#' A numeric vector of length 3.
#' @param sigma.mat The covariance matrix of the multivariate normal distribution for predictors
#' x1, x2, x3. A 3x3 numeric matrix.
#' @param eta The coefficients of the logistic regression model for the binary outcome A.
#' A numeric vector of length 4, representing intercept and coefficients for x1, x2, x3.
#' @param theta The coefficients of the linear regression model for the continuous outcome Y.
#' A numeric vector of length 5, representing intercept, coefficients for A, x1, x2, x3.
#'
#' @return A data frame containing the generated data with columns x1, x2, x3, A, Y.
#'
#' @importFrom MASS mvrnorm
#' @export 

dataGenFun <- function(n = 500,  mu = c(10, 1, 11), sigma.mat = matrix(c(11^2, 4, 4,
                                                                       4, 3^2, 3,
                                                                       4, 3, 8^2), 3),
                       eta = c(-1.5,0.1, 0.3, 0.2), theta =  c(110, -12, -0.3, -0.8, -0.2)){
  mvn <- MASS::mvrnorm(n, mu = mu, Sigma = sigma.mat)
  x1 <- mvn[, 1]
  x2 <- mvn[, 2]
  x3 <- mvn[, 3]
  pi <- 1 / (1 + exp(-(eta[1] + eta[2] * x1 + eta[3] * x2 + eta[4] * x3)))
  A <- rbinom(n, 1, pi)
  Y <- theta[1] + theta[2] * A + theta[3] * x1 + theta[4] * x2 + theta[5] * x3 + rnorm(n, 0, 1)
  data <- as.data.frame(cbind(x1, x2, x3, A, Y))
  return(data)
}
