
# Simulate multivariate normal data for correlated normal, lognormal, and
# binary data

generate_datasets <- function(num_sims,
                              num_observations,
                              num_endpoints,
                              type_vector,
                              mean_vector,
                              sd_vector,
                              effect_vector,
                              cor_mat,
                              start.seed){

  datasets <- list()

  bin.vars <- which(type_vector == "b")
  log.vars <- which(type_vector == "l")

  # Mean vectors - update for variable type and include treatment effects
  mean_vector.t <- mean_vector*effect_vector
  mean_vector.c <- mean_vector

  mean_vector.c[bin.vars] <- qlogis(mean_vector[bin.vars])
  mean_vector.t[bin.vars] <- qlogis(mean_vector[bin.vars]*
                                      effect_vector[bin.vars])

  mean_vector.c[log.vars] <- log(mean_vector[log.vars])
  mean_vector.t[log.vars] <- log(mean_vector[log.vars]*
                                   effect_vector[log.vars])

  # Standard deviation vectors - only difference is in binary var
  sd_vector.c <- sd_vector
  sd_vector.t <- sd_vector

  sd_vector.c[bin.vars] <- sqrt(1/(num_observations*mean_vector[bin.vars]*
                                     (1-mean_vector[bin.vars])))
  sd_vector.t[bin.vars] <- sqrt(1/(num_observations*mean_vector[bin.vars]*
                                     effect_vector[bin.vars]*
                                     (1-mean_vector[bin.vars]*
                                        effect_vector[bin.vars])))

  # Sub in 0 cor mat if not given
  if(sum(is.na(cor_mat)) > 0){
    cor_mat <- matrix(0, nrow = num_endpoints, ncol = num_endpoints)
    diag(cor_mat) = 1
  }

  # Variance/covariance matrix given sds and correlations
  vcov_matrix.c <- diag(sd_vector.c)%*%cor_mat%*%diag(sd_vector.c)
  vcov_matrix.t <- diag(sd_vector.t)%*%cor_mat%*%diag(sd_vector.t)


  for (i in 1:num_sims) {
    # Generate correlated normally distributed variables

    # simulate multivariate normal
    set.seed(start.seed + i)
    endpoints.c <- MASS::mvrnorm(n = num_observations,
                                 mu = mean_vector.c,
                                 Sigma = vcov_matrix.c)
    set.seed(start.seed + i + num_observations)
    endpoints.t <- MASS::mvrnorm(n = num_observations,
                                 mu = mean_vector.t,
                                 Sigma = vcov_matrix.t)

    #transform binary back to response var
    if(length(bin.vars) == 1){
      bern.trans.c <- (endpoints.c[,bin.vars] - mean_vector.c[bin.vars]) /
        sd_vector.c[bin.vars]
      bern.tquan.c <- qnorm(mean_vector[bin.vars], lower.tail = F)
      endpoints.c[,bin.vars] <- (bern.trans.c > bern.tquan.c)*1

      bern.trans.t <- (endpoints.t[,bin.vars] - mean_vector.t[bin.vars]) /
        sd_vector.t[bin.vars]
      bern.tquan.t <- qnorm(mean_vector[bin.vars]*effect_vector[bin.vars],
                            lower.tail = F)
      endpoints.t[,bin.vars] <- (bern.trans.t > bern.tquan.t)*1
    }

    if(length(bin.vars) > 1){
      bern.trans.c <- sweep(endpoints.c[,bin.vars], 2, mean_vector.c[bin.vars]) /
        sd_vector.c[bin.vars]
      bern.tquan.c <- qnorm(mean_vector[bin.vars], lower.tail = F)
      endpoints.c[,bin.vars] <- sweep(bern.trans.c, 2, FUN = ">", bern.tquan.c)*1


      bern.trans.t <- sweep(endpoints.t[,bin.vars],2, mean_vector.t[bin.vars]) /
        sd_vector.t[bin.vars]
      bern.tquan.t <- qnorm(mean_vector[bin.vars]*effect_vector[bin.vars],
                            lower.tail = F)
      endpoints.t[,bin.vars] <- sweep(bern.trans.t, 2, FUN = ">", bern.tquan.t)*1
    }


    # Combine all endpoints
    dataset.c <- as.data.frame(endpoints.c) %>%
      mutate(treated = 0,
             id = row_number())

    dataset.t <- as.data.frame(endpoints.t) %>%
      mutate(treated = 1,
             id = row_number() + num_observations)

    dataset <- rbind.data.frame(dataset.c, dataset.t)

    # Store dataset
    datasets[[i]] <- dataset
  }
  return(datasets)
}


