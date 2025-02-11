#' Probabilistic Structured Grammatical Evolution in R
#'
#' @description
#' Implementation of Probabilistic Structured Grammatical Evolution (PSGE) in R.
#' PSGE combines Structured Grammatical Evolution with probabilistic grammar-guided
#' genetic programming to evolve solutions to various problems.
#'
#' @details
#' This package provides:
#' * A complete implementation of the PSGE algorithm
#' * Grammar parsing and validation
#' * Probabilistic mapping mechanisms
#' * Evolutionary operators
#' * High-performance C++ core through Rcpp
#'
#' @name Rpsge-package
#' @aliases Rpsge
#' @useDynLib Rpsge, .registration = TRUE
#' @importFrom Rcpp loadModule evalCpp
#' @importFrom methods new
"_PACKAGE"
