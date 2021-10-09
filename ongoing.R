
eff2_adasel <- function(res,preferential=F){
  
  result <- res$result
  l_time1 <- nrow(res$summary)
  l_time2 <- nrow(res$selIntsummary)
  if (!preferential){
    samp <- INLA::inla.posterior.sample(500, res$result,
                                           selection = list(time=seq((l_time1-l_time2)+1,l_time1),
                                                            time2b=seq(1,l_time2)))
    posterior <- sapply(samp,function(x) apply(exp(matrix(c(x$latent),ncol=2,byrow=F))^(-1),1,prod))
    median <- apply(posterior,1,quantile, probs=0.5)
    quant0.025 <- apply(posterior,1,quantile, probs=0.025)
    quant0.975 <- apply(posterior,1,quantile, probs=0.975)
    mean <- apply(posterior,1,mean)
    result$effpop2 <- median
    result$effpop2mean <- mean
    result$effpop2_975 <- quant0.975
    result$effpop2_025 <- quant0.025
    result$effpop2summary <- data.frame(time = res$selIntsummary$time, mean = mean,  
                                        quant0.025 = quant0.025, quant0.5 = median, 
                                        quant0.975 = quant0.975)
    # result$effpop2 <- exp(-median)
    # result$effpop2mean <- exp(-mean)
    # result$effpop2_975 <- exp(-quant0.975)
    # result$effpop2_025 <- exp(-quant0.025)
    # result$effpop2summary <- data.frame(time = res$selIntsummary$time, mean = exp(-mean),  
    #                                   quant0.025 = exp(-quant0.025), quant0.5 = exp(-median), 
    #                                   quant0.975 = exp(-quant0.975))

  }
  return(result)
}

selInt_independent <- function(res1,res2){

  #Compute a ratio of two independent BNPR runs
  l_time1 <- nrow(res1$summary)
  l_time2 <- nrow(res2$summary)
  samp1 <- INLA::inla.posterior.sample(500, res1$result,selection = list(time=seq(1,l_time1)))
  posterior1 <- sapply(samp1,function(x) exp(matrix(c(x$latent),ncol=1,byrow=F))^(-1))
  samp2 <- INLA::inla.posterior.sample(500, res2$result,selection = list(time=seq(1,l_time2)))
  posterior2 <- sapply(samp2,function(x) exp(matrix(c(x$latent),ncol=1,byrow=F))^(-1))
  posterior <- posterior2/posterior1
  median <- apply(posterior,1,quantile, probs=0.5)
  quant0.025 <- apply(posterior,1,quantile, probs=0.025)
  quant0.975 <- apply(posterior,1,quantile, probs=0.975)
  mean <- apply(posterior,1,mean)
  result <- list()
  result$selInt_ind <- median
  result$selInt_indmean <- mean
  result$selInt_ind_975 <- quant0.975
  result$selInt_ind_025 <- quant0.025
  result$selInt_indsummary <- data.frame(time = res1$summary$time, mean = mean,  
                                      quant0.025 = quant0.025, quant0.5 = median, 
                                      quant0.975 = quant0.975)
  return(result)
}
 




