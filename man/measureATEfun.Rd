\name{measureATEfun}
\alias{measureATEfun}
\title{Measure Average Treatment Effect (ATE) Using Trimming}
\description{
  `measureATEfun()` estimates the Average Treatment Effect (ATE) using trimming in the weighting step to reduce the influence of extreme weights. It calculates the bias, MSE for estimated ATE, and it gives summary statistics for weights in each iteration. The function generates `nsim` datasets and estimates the ATE for each dataset using a Generalized Linear Model (GLM) with binary treatment (A) and continuous confounders (x1, x2, x3). Trimming is applied to the weights calculated based on the estimated propensity scores to reduce the impact of extreme weights.
}
\usage{
  measureATEfun(nsim = 1000, trim.p = 0.95, n = 500, mu = c(10, 1, 11), sig.mat = matrix(c(11^2, 4, 4, 4, 3^2, 3, 4, 3, 8^2), 3), eta = c(-1.5, 0.1, 0.3, 0.2), theta = c(110, -12, -0.3, -0.8, -0.2))
}
\arguments{
  \item{nsim}{
    The number of simulations to perform. A positive integer.
  }

  \item{trim.p}{
    The percentile of extreme weights to be trimmed on both ends. A numeric value between 0.5 and 1. For example, 0.9 represents trimming the top 10\% and bottom 10\% of the weights.
  }

  \item{n}{
    The number of observations in each simulated dataset. A positive integer.
  }

  \item{mu}{
    The mean vector of the multivariate normal distribution for confounders x1, x2, x3. A numeric vector of length 3.
  }

  \item{sig.mat}{
    The covariance matrix of the multivariate normal distribution for confounders x1, x2, x3. A 3x3 numeric matrix.
  }

  \item{eta}{
    The coefficients of the logistic regression model for the binary treatment A. A numeric vector of length 4, representing intercept and coefficients for x1, x2, x3.
  }

  \item{theta}{
    The coefficients of the linear regression model for the continuous outcome Y. A numeric vector of length 5, representing intercept, coefficients for A, x1, x2, x3.
  }
}
\details{
  The \code{measureATEfun} function generates synthetic data using \code{\link{dataGenFun}} and then estimates the ATE using trimming in the weighting step. The steps include:

  \enumerate{
    \item Generate \code{nsim} random datasets using \code{\link{dataGenFun}} function with specified parameters.
    \item For each dataset, estimate the propensity score for treatment assignment (A) using logistic regression.
    \item Calculate the weights for each observation based on the estimated propensity scores.
    \item Trim the extreme weights by setting them to the corresponding quantiles of the weight distribution.
    \item Estimate the ATE for each dataset using weighted linear regression.
    \item Record the bias, MSE of the estimator, and summary statistics of the trimmed weights for each simulation.
  }
}
\value{
  The function returns a matrix containing measures related to ATE estimation using trimming. The columns include:

  \item{bias}{
    The difference between the estimated and true ATE.
  }

  \item{mse}{
    The mean squared error (MSE) of the ATE estimator.
  }

  \item{weight.mean}{
    The mean of trimmed weights across simulations.
  }

  \item{weight.ese}{
    The empirical standard error of the mean of trimmed weights across simulations.
  }

  \item{weight.asd}{
    The average standard deviation of trimmed weights across simulations.
  }

  \item{weight.max}{
    The average maximum trimmed weight across simulations.
  }

  \item{weight.min}{
    The average minimum trimmed weight across simulations.
  }
}
\examples{
# Generate 1000 simulations with trimming at 95th percentile
result <- measureATEfun(nsim = 1000, trim.p = 0.95, n = 500,
                        mu = c(10, 1, 11),
                        sig.mat = matrix(c(11^2, 4, 4, 4, 3^2, 3, 4, 3, 8^2), 3),
                        eta = c(-1.5, 0.1, 0.3, 0.2),
                        theta = c(110, -12, -0.3, -0.8, -0.2))

# View the results
print(result)
}
