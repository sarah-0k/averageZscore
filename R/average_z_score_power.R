#' @include generate_datasets.R zscore_function.R simulation_results.R
NULL

#' Generate simulated data and calculate power of average z-score approach.
#'
#' \code{average_z_score_power} calculates simulated power for two-arm average z-score approach.
#'
#' @param num_sims A numeric value representing the number of simulations to be run.
#'
#' @param num_observations A numeric value representing the number of observations
#' per arm. Currently only two-arms with equal allocation is available.
#'
#' @param num_endpoints A numeric value between 1 and 5 that represents the number
#' of endpoints being used. It is not recommended to have more than 5 endpoints.
#'
#' @param type_vector A character vector of length \code{num_endpoints} that specifies
#' the type of endpoints. Currently normal (\code{n}), log-normal (\code{l}),
#' and binary (\code{b}) are supported.
#'
#' @param mean_vector A numeric vector of length \code{num_endpoints} that specifies
#' the expected mean values in the control arm. Log-normal means represent the
#' mean of follow/baseline values; binary means represent the proportion of
#' events or successes.
#'
#' @param sd_vector A numeric vector of length \code{num_endpoints} that specifies
#' the expected standard deviation (sd) values in the control arm.
#' Log-normal sds represent the sd of log(follow/baseline) values.
#' Binary sds are calculated within the function and should be listed as \code{NA}.
#'
#' @param effect_vector A numeric vector of length \code{num_endpoints} that specifies
#' the expected relative treatment effect for each endpoint (i.e., 1.2 represents
#' a 20% increase in the treatment arm).
#'
#' @param cor_mat A numeric matrix of dimensions \code{num_endpoints} x \code{num_endpoints}
#' that specifies the correlations between all endpoints. The diagnonal of the
#' correlation matrix should equal all 1s. If not specified, no correlations
#' between any endpoints will be assumed. \bold{(optional)}
#'
#' @param direction_vector A numeric vector of length \code{num_endpoints} that specifies
#' the direction of positive treatment effect. All elements should be either a
#' 1 (larger values are better), or a -1 (smaller values are better).
#'
#' @param threshold A numeric value that specifies the p-value cutoff for
#' statistical significance, based on a two-sided test. If not specified, a
#' threshold of 0.05 will be assumed. \bold{(optional)}
#'
#' @param start.seed A numeric value that specifies the starting seed for
#' simulation of the multivariate normal responses. \bold{(optional)}
#'
#' @examples
#' test_scores <- average_z_score_power(num_sims = 100,
#' num_observations = 100,
#' num_endpoints = 5,
#' type_vector = c("n", "n", "l", "n", "b"),
#' mean_vector = c(100, -5, 0.8, 30, 0.3),
#' sd_vector = c(75, 5, 0.8, 20, NA),
#' effect_vector = c(1.2, 1.2, 0.8, 1.2, 0.8),
#' cor_mat = matrix(c(1, -0.3, -0.3, 0, -0.5,
#'                    -0.3, 1, 0, 0, 0,
#'                    -0.3, 0, 1, 0, 0.4,
#'                    0, 0, 0, 1, 0,
#'                    -0.5, 0, 0.4, 0, 1),
#'                  nrow = 5, byrow = T),
#' direction_vector = c(1, -1, -1, 1, -1),
#' threshold = 0.05,
#' start.seed = 123456)
#'
#' test_scores$power
#'
#'
#' @import dplyr
#' @importFrom stats wilcox.test qnorm qlogis
#' @importFrom MASS mvrnorm
#' @export


average_z_score_power <- function(num_sims,
                                  num_observations,
                                  num_endpoints,
                                  type_vector,
                                  mean_vector,
                                  sd_vector,
                                  effect_vector,
                                  cor_mat = NA,
                                  direction_vector,
                                  threshold = 0.05,
                                  start.seed = 123456){

  if(missing(num_sims)){
    return(noquote("Please specify number of simulations (minimum 1)."))
  }
  if(is.na(num_sims) | num_sims < 1){
    return(noquote("Please specify number of simulations (minimum 1)."))
  }

  if(missing(num_observations)){
    return(noquote("Please specify number of observations per arm (minimum 1)."))
  }
  if(is.na(num_observations) | num_observations < 1){
    return(noquote("Please specify number of observations per arm (minimum 1)."))
  }

  if(missing(num_observations)){
    return(noquote("Please specify number of endpoints (minimum 2, maximum 5)."))
  }
  if(is.na(num_endpoints)==T | is.null(num_endpoints)==T){
    return(noquote("Please specify number of endpoints (minimum 2, maximum 5)."))
  }
  if(num_endpoints > 5){
    return(noquote("Please specify no more than 5 endpoints."))
  }
  if(num_endpoints < 2){
    return(noquote("Please specify at least 2 endpoints."))
  }

  if(missing(num_sims)){
    return(noquote("Please specify vector of endpoint types (b, l, or n)."))
  }
  if(length(type_vector) != num_endpoints){
    return(noquote("The length of type vector must equal the number of endpoints."))
  }
  if(sum(type_vector %in% c("b", "l", "n")) != num_endpoints){
    return(noquote("Ensure all elements in type_vector are one of: b, l, or n."))
  }

  if(missing(mean_vector)){
    return(noquote("Please specify vector of means."))
  }
  if(length(mean_vector) != num_endpoints){
    return(noquote("The length of the mean vector must equal the number of endpoints."))
  }
  if(sum(is.na(mean_vector)) > 0){
    return(noquote("Please supply non-NA means."))
  }
  if(sum(mean_vector[which(type_vector == "l")] <= 0) > 0){
    return(noquote("Lognormal means must be positive and non-zero."))
  }
  if(sum(mean_vector[which(type_vector == "b")] <= 0) > 0 |
     sum(mean_vector[which(type_vector == "b")] >= 1) > 0){
    return(noquote("Binary means must be a decimal bounded by (0, 1)."))
  }

  if(missing(sd_vector)){
    return(noquote("Please specify vector of standard deviations (use NA for binary)."))
  }
  if(length(sd_vector) != num_endpoints){
    return(noquote("The length of the standard deviation vector must equal the number of endpoints."))
  }
  if(sum(type_vector %in% c("l", "n")) > 0){
    if(sum(is.na(sd_vector[which(type_vector %in% c("l", "n"))])) > 0 |
       sum(sd_vector[which(type_vector %in% c("l", "n"))] < 0) > 0){
      return(noquote("Please ensure non-negative standard deviations are supplied for non-binary endpoints."))
    }
    if(sum(is.na(sd_vector[which(type_vector %in% c("b"))])) !=
       sum(type_vector %in% c("b"))){
      return(noquote("Please ensure NA standard deviations are supplied for binary endpoints."))
    }
  }

  if(missing(effect_vector)){
    return(noquote("Please specify a vector of treatment effects."))
  }
  if(length(effect_vector) != num_endpoints){
    return(noquote("The length of the treatment effect vector must equal the number of endpoints."))
  }
  if(sum(is.na(effect_vector)) > 0){
    return(noquote("Please supply non-NA effect sizes."))
  }


  if(sum(is.na(cor_mat)) == 0){
    if(dim(cor_mat)[1] != num_endpoints | dim(cor_mat)[2] != num_endpoints){
      return(noquote("The dimensions of the correlation matrix must be equal to the number of endpoints."))
    }
    if(sum(diag(cor_mat) == 1) != num_endpoints){
      return(noquote("The diagonal of the correlation matrix must be all 1s."))
    }
    if(sum(!(min(cor_mat) >= -1 & max(cor_mat) <= 1)) > 0){
      return(noquote("Please ensure all correlations are between -1 and 1."))
    }
    if(isSymmetric(cor_mat) == F){
      return(noquote("Please ensure correlation matrix is symmetric."))
    }
  }


  if(missing(direction_vector)){
    return(noquote("Please specify vector of treatment effect direction (-1, 1)."))
  }
  if(length(direction_vector) != num_endpoints){
    return(noquote("The length of treatment effect direction vector must equal the number of endpoints."))
  }
  if(sum(direction_vector %in% c(-1, 1)) != num_endpoints){
    return(noquote("Ensure all elements in type_vector are one of: -1, 1"))
  }


  if(is.na(threshold)){
    return(noquote("Please specify a p-value threshold for the two-sided test."))
  }
  if(threshold >= 1 | threshold <= 0){
    return(noquote("P-value threshold must be (0, 1)."))
  }


  datalist.sim <- generate_datasets(num_sims = num_sims,
                                    num_observations = num_observations,
                                    num_endpoints = num_endpoints,
                                    type_vector = type_vector,
                                    mean_vector = mean_vector,
                                    sd_vector = sd_vector,
                                    effect_vector = effect_vector,
                                    cor_mat = cor_mat,
                                    start.seed = start.seed)

  datalist.simz <- zscore_function(datalist = datalist.sim,
                                   direction_vector = direction_vector)

  final.results <- simulation_results(datalist = datalist.simz,
                                      type_vector = type_vector,
                                      threshold = threshold)


  return(final.results)
}
