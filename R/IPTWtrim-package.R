#' IPTWtrim: Performance Analysis on Inverse Probability of Treatment Weighting
#' (IPTW) Trimming for Causal Inference
#'
#' The IPTWtrim package provides functions for conducting causal inference with inverse probability
#' of treatment weighting (IPTW) and weight trimming. It allows you to estimate the Average Treatment
#' Effect (ATE) using weight trimming based on simulation studies.
#'
#' @docType package
#' @name IPTWtrim-package
#' @aliases IPTWtrim
#' @description
#' The IPTWtrim package contains 3 main functions:
#' \itemize{
#'   \item \code{dataGenFun()}: Generates a dataset with a treatment, confounders and outcomes for simulation studies.
#'     See \link{dataGenFun} for details.
#'   \item \code{measureATEfun()}: Estimates the Average Treatment Effect (ATE) using weight trimming based
#'   on a simulation study, and gives performance measures. See \link{measureATEfun} for details.
#'   \item \code{trimSummary()}: Performs a simulation study to evaluate the performance of weight trimming
#'   with different weight trimming percentiles. It calculates bias, empirical standard error (EmpSE),
#'   and mean squared error (MSE) for each weight trimming percentile.Additionally, the function computes
#'   summary statistics for weights. See \link{trimSummary} for details.
#' }
#'
#' For examples and usage, see \code{vignette("IPTWtrim_handbook")}.
#'
#' @author Xiaoya Wang (\email{x932wang@uwaterloo.ca})
#'
NULL
