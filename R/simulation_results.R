
# Calculate summarize simulated data results, calculate per-simulation
# Wilcoxon Mann-Whitney p-value, calculate simulated power based on p-val
# threshold (two-sided test)


simulation_results <- function(datalist,
                               type_vector,
                               threshold = 0.05){

  log.vars <- which(type_vector == "l")

  results <- list()

  results$n.sims <- length(datalist)
  results$n <- length(datalist[[1]]$id)
  results$threshold <- threshold

  #report means for data simulations
  results$combined.means <- sapply(datalist, function(df){
    colMeans(df[,-(seq(length(type_vector)+1,length = 2))])
  })

  results$control.means <- sapply(datalist, function(df){
    colMeans(df[df$treated == 0,-(seq(length(type_vector)+1,length = 2))])
  })

  results$treated.means <- sapply(datalist, function(df){
    colMeans(df[df$treated == 1,-(seq(length(type_vector)+1,length = 2))])
  })

  #transform log-normal variables to GMR for reporting
  results$combined.means[log.vars,] <- exp(results$combined.means[log.vars,])
  results$control.means[log.vars,] <- exp(results$control.means[log.vars,])
  results$treated.means[log.vars,] <- exp(results$treated.means[log.vars,])

  #report standard deviations for data simulations
  results$combined.sd <- sapply(datalist, function(df){
    apply(df[,-(seq(length(type_vector)+1,length = 2))], 2, sd)
  })

  results$control.sd <- sapply(datalist, function(df){
    apply(df[df$treated == 0,-(seq(length(type_vector)+1,length = 2))], 2, sd)
  })

  results$treated.sd <- sapply(datalist, function(df){
    apply(df[df$treated == 1,-(seq(length(type_vector)+1,length = 2))], 2, sd)
  })

  #report mean correlations of data simulations
  corlist <- lapply(datalist, function(df){
    cor(df[,c(1:length(type_vector))])
  })

  results$mean.cor <- Reduce("+", corlist)/length(corlist)

  #report p-values and report power
  results$p.vals <- sapply(datalist, function(df){
    wilcox.test(z.average ~ treated, data = df)$p.value})

  results$power <- mean(results$p.vals < threshold)

  return(results)
}
