\name{dataGenFun}
\alias{dataGenFun}
\title{Generate Data from a Specific Model}
\description{
  Generate data based on a specified statistical model. The function simulates a dataset with 'n' observations where the response variable ('Y') is influenced by confounding variables ('x1', 'x2', and 'x3'), and a binary treatment indicator ('A'). The function uses a logistic regression model to determine the probability of treatment assignment ('A') and a linear regression model to generate the outcome variable ('Y'). The predictor variables 'x1', 'x2', and 'x3' are multivariate normally distributed.
}
\usage{
  dataGenFun(n, mu, sig.mat, eta, theta)
}
\arguments{
  \item{n}{
    The number of observations to generate in the dataset. It should be a positive integer.
  }

  \item{mu}{
    A numeric vector of length 3, representing the means of the multivariate normal distribution for the predictor variables 'x1', 'x2', and 'x3', respectively.
  }

  \item{sig.mat}{
    A 3x3 positive definite covariance matrix representing the covariances and variances of the multivariate normal distribution for the predictor variables 'x1', 'x2', and 'x3'.
  }

  \item{eta}{
    A numeric vector of length 4, representing the coefficients of the logistic regression model for calculating the probability of treatment assignment ('A'). The elements correspond to the intercept and the coefficients for 'x1', 'x2', and 'x3', respectively.
  }

  \item{theta}{
    A numeric vector of length 5, representing the coefficients of the linear regression model for generating the outcome variable ('Y'). The elements correspond to the intercept, the coefficient for the treatment indicator 'A', and the coefficients for 'x1', 'x2', and 'x3', respectively.
  }
}
\details{
  The \code{dataGenFun} function generates synthetic data by following these steps:

  - It generates 'n' random samples from a multivariate normal distribution using \code{\link[MASS:mvrnorm]{MASS::mvrnorm}}, with means specified by \code{mu} and covariance matrix specified by \code{sig.mat}.
  - It computes the logistic transformation of a linear combination of the predictor variables 'x1', 'x2', and 'x3', along with corresponding coefficients \code{eta}, to determine the probability of treatment assignment ('A'). The logistic function is used to constrain the probability between 0 and 1.
  - It generates binary treatment indicators ('A') for each observation by sampling from a binomial distribution with probability \code{pi}.
  - It computes the outcome variable ('Y') for each observation using a linear combination of the treatment indicator 'A', confounding variables 'x1', 'x2', and 'x3', and corresponding coefficients \code{theta}. Additionally, a random error term is added to 'Y' by sampling from a normal distribution with mean 0 and standard deviation 1.
  - The generated variables 'x1', 'x2', 'x3', 'A', and 'Y' are combined into a data frame, and the resulting dataset is returned.
}
\value{
  The function returns a data frame containing the generated dataset with 'n' observations and the following columns:

  - \item{x1}{The values of the first confounding variable.}

  - \item{x2}{The values of the second confounding variable.}

  - \item{x3}{The values of the third confounding variable.}

  - \item{A}{The binary treatment indicator (0 or 1) representing the treatment assignment.}

  - \item{Y}{The outcome variable generated based on the specified statistical model.}
}
\author{
  Xiaoya Wang
}
\note{
  - The function relies on the \code{MASS} package to generate multivariate normal samples (\code{\link[MASS:mvrnorm]{MASS::mvrnorm}}) and the \code{\link{rbinom}} and \code{\link{rnorm}} functions for random binomial and normal variates, respectively. Make sure to install and load the \code{MASS} package before using this function.

  - The function assumes that the input arguments \code{mu}, \code{sig.mat}, \code{eta}, and \code{theta} are appropriately specified according to the desired statistical model. Incorrectly specified parameters may lead to unexpected or nonsensical results.

  - It is important to note that the generated data is purely synthetic and is not based on any real-world data source. The user should exercise caution and thoroughly understand the underlying statistical model before using the generated data for any analyses or applications.
}
