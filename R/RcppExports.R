# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Gauss-Laguerre quadrature
#'
#' Computes the node and the weights for the Gauss-Laguerre quadrature (integral on the whole real line)
#'
#' @name gauss_laguerre
#' @param N the number of evaluations
#' @return a list containing two numeric vectors of length N, the first one containing the nodes and the second one the weights
#' @importFrom Rcpp evalCpp
#' @export
gauss_laguerre <- function(N) {
    .Call(`_micsr_gauss_laguerre`, N)
}

#' Gauss-Hermitte quadrature
#'
#' Computes the node and the weights for the Gauss-Hermite quadrature (integral on the whole real line)
#'
#' @name gaussian_quad
#' @param N the number of evaluations
#' @return a list containing two numeric vectors of length N, the first one containing the nodes and the second one the weights
#' @export
gauss_hermite <- function(N) {
    .Call(`_micsr_gauss_hermite`, N)
}

punorm0 <- function(x) {
    .Call(`_micsr_punorm0`, x)
}

pbnorm0 <- function(h1, h2, rho) {
    .Call(`_micsr_pbnorm0`, h1, h2, rho)
}

ptnorm0 <- function(z1, z2, z3, rho12, rho13, rho23) {
    .Call(`_micsr_ptnorm0`, z1, z2, z3, rho12, rho13, rho23)
}

#' Compute the probability for the univariate normal function
#'
#' @param z a numeric vector
#' @return a numeric vector
#' @export
punorm <- function(z) {
    .Call(`_micsr_punorm`, z)
}

#' Compute the probability for the bivariate normal function
#'
#' @param z1,z2 two numeric vectors
#' @param rho a numeric vector
#' @return a numeric vector
#' @export
pbnorm <- function(z1, z2, rho) {
    .Call(`_micsr_pbnorm`, z1, z2, rho)
}

ptnormv <- function(z1, z2, z3, rho12, rho13, rho23) {
    .Call(`_micsr_ptnormv`, z1, z2, z3, rho12, rho13, rho23)
}

#' Compute the probability for the trivariate normal function
#'
#' @param z a matrix with three columns
#' @param rho a matrix with three columns
#' @return a numeric vector
#' @export
ptnorm <- function(z, rho) {
    .Call(`_micsr_ptnorm`, z, rho)
}

