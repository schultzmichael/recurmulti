#' recurmulti: A package for estimating Recurrent Multinomial Models
#'
#' The recurmulti package provides functions to fit and evaluated recurrent
#' multinomial models
#'
#' Note that the parameter labels in this package differ from those in the
#' paper.
#' The package uses eta, beta, and zeta for coefficients, corresponding to
#' X, Z, and Q covariate matrices. Eta contains the baseline
#' probability coefficients, analogous to traditional multinomial
#' regression coefficients. Beta contains recurrence coefficients for sequence
#' position-specific covariates, while zeta contains recurrence coefficients
#' for position-independent covariates. This provides memory and performance
#' improvements.
#'
#' Future versions will harmonize the naming of coefficients with the paper.
#' @docType package
#' @name recurmulti
NULL