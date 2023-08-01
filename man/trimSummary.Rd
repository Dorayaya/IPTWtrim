\name{trimSummary}
\alias{trimSummary}
\title{Summary of ATE Estimation Using Weight Trimming}
\description{
  `trimSummary()` performs a simulation study to evaluate the performance of different weight trimming percentiles. It calculates bias, empirical standard error (EmpSE), and mean squared error (MSE) for each weight trimming percentile specified in `trimperc.try`. Additionally, the function computes summary statistics for weights, including mean of trimmed weights, standard error of mean of trimmed weights, average standard deviation of trimmed weights, average maximum trimmed weight, and average minimum trimmed weight. The function then generates three plots using `ggplot2` to visualize the relationship between weight trimming percentile and bias, EmpSE, and MSE.
}
\usage{
  trimSummary(nsim = 1000, n = 500, mu = c(10, 1, 11),
              sig.mat = matrix(c(11^2, 4, 4, 4, 3^2, 3, 4, 3, 8^2), 3),
              eta = c(-1.5, 0.1, 0.3, 0.2),
              theta = c(110, -12, -0.3, -0.8, -0.2),
              trimperc.try = c(1, 0.99, 0.95, 0.9, 0.75, 0.5))
}
\arguments{
  \item{nsim}{
    The number of iterations to perform for each weight trimming percentile. An integer. Default to 1000.
  }

  \item{n}{
    The number of observations in each simulated dataset. An integer. Default to 500.
  }

  \item{mu}{
    The mean vector of the multivariate normal distribution for confounders x1, x2, x3. A numeric vector of length 3. Default to c(10, 1, 11).
  }

  \item{sig.mat}{
    The covariance matrix of the multivariate normal distribution for confounders x1, x2, x3. A 3x3 numeric matrix. Default to a specific matrix.
  }

  \item{eta}{
    The coefficients of the logistic regression model for the binary treatment A. A numeric vector of length 4, representing intercept and coefficients for x1, x2, x3. Default to a specific numeric vector.
  }

  \item{theta}{
    The coefficients of the linear regression model for the continuous outcome Y. A numeric vector of length 5, representing intercept, coefficients for A, x1, x2, x3. Default to a specific numeric vector.
  }

  \item{trimperc.try}{
    A numeric vector specifying the weight trimming percentiles to evaluate. Each value should be between 0.5 and 1 (both inclusive) to represent the desired percentile for weight trimming.
  }
  \item{seed}{An optional seed value for reproducibility of simulations. Default to 20739377.}
}
\details{
  The \code{trimSummary} function performs a simulation study to evaluate the performance of different weight trimming percentiles. It generates \code{nsim} datasets for each weight trimming percentile in \code{trimperc.try} and calculates bias, EmpSE, and MSE for each percentile. Additionally, it computes summary statistics for weights, including mean of trimmed weights, standard error of mean of trimmed weights, average standard deviation of trimmed weights, average maximum trimmed weight, and average minimum trimmed weight. The function then generates three plots using \code{ggplot2} to visualize the relationship between weight trimming percentile and bias, EmpSE, and MSE.
}
\value{
  The function returns a data frame containing statistics (bias, EmpSE, and MSE) for each weight trimming percentile in \code{trimperc.try}. The data frame also includes summary statistics for weights:

  \item{trim.perc}{
    Weight trimming percentile.
  }

  \item{trim.prop}{
    Weight trimming proportion (1 - trim.perc) * 2.
  }

  \item{weight.mean}{
    Mean of trimmed weights.
  }

  \item{weight.ese}{
    Standard error of mean of trimmed weights.
  }

  \item{weight.asd}{
    Average standard deviation of trimmed weights.
  }

  \item{weight.max}{
    Average maximum trimmed weight.
  }

  \item{weight.min}{
    Average minimum trimmed weight.
  }

  \item{bias}{
    Difference between estimated and true ATE.
  }

  \item{ese}{
    Empirical standard error of the estimator.
  }

  \item{mse}{
    Mean squared error of the estimator.
  }
}
\examples{
# Evaluate weight trimming at different percentiles
result <- trimSummary(nsim = 1000, n = 500, mu = c(10, 1, 11),
                      sig.mat = matrix(c(11^2, 4, 4, 4, 3^2, 3, 4, 3, 8^2), 3),
                      eta = c(-1.5, 0.1, 0.3, 0.2),
                      theta = c(110, -12, -0.3, -0.8, -0.2),
                      trimperc.try = c(1, 0.99, 0.95, 0.9, 0.75, 0.5))

# View the results
print(result)
}
