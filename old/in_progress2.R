BNPR_sel2 <- function (tree1,tree2, samp_times1,samp_times2, lengthout = 100, prec_alpha = 0.01, 
          prec_beta = 0.01, beta1_prec = 0.001,
          simplify = TRUE, derivative = FALSE, forward = TRUE) 
{

  phy1 <- summarize_phylo(tree1)
  phy1$samp_times <- phy1$samp_times + min(samp_times1)
  phy1$coal_times <- phy1$coal_times + min(samp_times1)
  phy2 <- summarize_phylo(tree2)
  phy2$samp_times <- phy2$samp_times + min(samp_times2)
  phy2$coal_times <- phy2$coal_times + min(samp_times2)
  
    
  result <- infer_coal_samp_selection2(phy1,phy2,lengthout = lengthout, 
                            prec_alpha = prec_alpha, prec_beta = prec_beta, beta1_prec = beta1_prec, 
                            simplify)
  
  #result$samp_times <- phy$samp_times
  #result$n_sampled <- phy$n_sampled
  #result$coal_times <- phy$coal_times
  result$effpop <- exp(-result$result$summary.random$time$`0.5quant`)
  result$effpopmean <- exp(-result$result$summary.random$time$mean)
  result$effpop975 <- exp(-result$result$summary.random$time$`0.025quant`)
  result$effpop025 <- exp(-result$result$summary.random$time$`0.975quant`)
  result$summary <- with(result$result$summary.random$time, 
                         data.frame(time = ID, mean = exp(-mean), sd = sd * exp(-mean), 
                                    quant0.025 = exp(-`0.975quant`), quant0.5 = exp(-`0.5quant`), 
                                    quant0.975 = exp(-`0.025quant`)))

  
  result$beta1 <- result$result$summary.hyperpar[2, "0.5quant"]
  result$beta1summ <- result$result$summary.hyperpar[2,]
  rownames(result$beta1summ) <- "Beta1"
  result$beta1post <- result$result$marginals.hyperpar$`Beta for time2`
  
  return(result)
}





infer_coal_samp_selection2 <- function(phy1,phy2, lengthout=100, prec_alpha=0.01, prec_beta=0.01, 
                            beta1_prec=0.001, simplify = TRUE, events_only = FALSE)
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop('INLA needed for this function to work. Use install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE).',
         call. = FALSE)
  }
  

  coal_times1 <- phy1$coal_times
  coal_times2 <- phy2$coal_times
  n_sampled1 <- phy1$n_sampled
  n_sampled2 <- phy2$n_sampled
  samp_times1 <- phy1$samp_times
  samp_times2 <- phy2$samp_times
  
  grid <- seq(min(c(samp_times1,samp_times2)), max(c(coal_times1,coal_times2)), length.out = lengthout+1)
  
  #if (is.null(n_sampled))
  #  n_sampled <- rep(1, length(samp_times))
  
  coal_data1 <- coal_stats(grid = grid, samp_times = samp_times1, n_sampled = n_sampled1,
                          coal_times = coal_times1)
  if (simplify)
    coal_data1 <- with(coal_data1, condense_stats(time=time, event=event, E=E))
  coal_data1 <- identify_off_grid(coal_data1,samp_times1,coal_times1)
  coal_data2 <- coal_stats(grid = grid, samp_times = samp_times2, n_sampled = n_sampled2,
                           coal_times = coal_times2)
  if (simplify)
    coal_data2 <- with(coal_data2, condense_stats(time=time, event=event, E=E))
  coal_data2 <- identify_off_grid(coal_data2,samp_times2,coal_times2)
  
  
  
  #simplify is used to avoid duplicate in time 

  
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))

  
  
   # data <- with(coal_data, data.frame(y = event, time = time, E_log = E_log))
    
  #  formula <- y ~ -1 + f(time, model="rw1", hyper = hyper, constr = FALSE)
  #  family <- "poisson"

   
    
  data <- joint_coal_stats(coal_data1 = coal_data1, coal_data2 = coal_data2)

  formula <- Y ~ -1 + 
      f(time, model="rw1", hyper = hyper, constr = FALSE) +
      f(time2,copy="time", fixed=FALSE, param=c(0, beta1_prec))
 
  family <- c("poisson", "poisson")

  lc_many <- NULL
  
  mod <- INLA::inla(formula, family = family, data = data,
                    lincomb = lc_many, offset = data$E_log,
                    control.predictor = list(compute=TRUE),
                    control.inla = list(lincomb.derived.only=FALSE))
  
  return(list(result = mod, data = data, grid = grid, x = coal_data1$time))
}



identify_off_grid <- function(coal_data,samp_times,coal_times)
{ 
  if (samp_times[1]>0){
    id <- which.min(abs(coal_data$time-samp_times[1]))
    coal_data$event[1:(id-1)] <- NA
  }
  if (coal_times[length(coal_times)]<coal_data$time[length(coal_data$time)]){
    id <- which.min(abs(coal_data$time-coal_times[length(coal_times)])) #Note that there are double points in $time, that's why id+1
    coal_data$event[(id+1):length(coal_data$event)] <- NA
  }
  return(coal_data)
}


joint_coal_stats <- function(coal_data1, coal_data2)
{
  n1 <- length(coal_data1$time) #Here are the midpts 
  n2 <- length(coal_data2$time) #
  beta0 <- c(rep(0, n1), rep(1, n2))
  E_log <- c(coal_data1$E_log, coal_data2$E_log) #samp_data$E_log does not include NA in the pref_samp, so I am not modifying it here either
  Y <- matrix(c(coal_data1$event, rep(NA, n2), rep(NA, n1), coal_data2$event),
              nrow = n1 + n2, byrow = FALSE)
  #Need to correct Y for the parts of the grid that need to be eccluded
  #Note, there are some parts of Y where they are both NAs
  w <- c(rep(1, n1), rep(1, n2))  #what does the w do? I think that takes into account for the fact that one is at the numerator,
  #the other one at the numerator. 
  
  #we are not correcting the time, i.e. there are times also to part that are not defined
  time  <- c(coal_data1$time, rep(NA, n2))
  time2 <- c(rep(NA, n1), coal_data2$time)
  
  return(list(Y = Y, beta0 = beta0, time = time, time2 = time2, w = w, E_log = E_log))
}



###############################
##### OLD (need for tryout) ###
###############################

coal_stats <- function(grid, samp_times, coal_times, n_sampled = NULL,
                       log_zero = -100)
{
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid)/2 #It is practically the mid point of the grid
  
  #if (is.null(n_sampled))
  #  n_sampled <- rep(1, length(samp_times))
  args <- gen_INLA_args(samp_times = samp_times, n_sampled = n_sampled,
                        coal_times = coal_times) #Looks identical between, just shifted
  #args refer only to the period in which I have coal and samp events. It does not matter the grid. It simply sort them relatively
  
  
  coal_factor <- args$coal_factor
  s <- args$s
  event <- args$event
  
  grid_trimmed <- setdiff(x = grid, y = s) #elements in the grid that are not in the other vector
  sorting <- sort(c(grid_trimmed, s), index.return=TRUE)
  sgrid <- sorting$x #grid + samp_events + coal_events - duplicates (e.g. time=0)
  ordering <- sorting$ix # the numbers before the lenght of the grid are grid points, the rest no
  
  time_index <- cut(x = sgrid[-1], breaks = grid, labels = FALSE) # It tells me which point of sgrid belong to a given interval in the grid
  #it is of the same length of sgrid but without point 0
  time <- field[time_index] #repeats all the midpoint in the grid
  
  event_out <- c(rep(0, length(grid_trimmed)), event)[ordering] #spread out the location of the coal_event (1s) in the sgrid, 
  #grid point and sampling event are 0s
  
  Cfun <- stats::stepfun(x = s, y = c(0, coal_factor, 0), right = TRUE)
  Cvec <- Cfun(sgrid[-1])
  E <- diff(sgrid)*Cvec
  
  E_log = log(E)
  E_log[E == 0] = log_zero
  
  return(data.frame(time = time, event = event_out[-1], E = E, E_log = E_log))
}

gen_INLA_args <- function(samp_times, n_sampled, coal_times)
{
  if (sum(n_sampled) != length(coal_times) + 1)
    stop("Number sampled not equal to number of coalescent events + 1.")
  
  if (length(intersect(coal_times, samp_times)) > 0)
    warning("Coincident sampling event and coalescent event: results may be unpredictable.")
  
  l <- length(samp_times)
  m <- length(coal_times)
  sorting <- sort(c(samp_times, coal_times), index.return=TRUE)
  
  lineage_change <- c(n_sampled, rep(-1, m))[sorting$ix]
  lineages <- utils::head(cumsum(lineage_change), -1) # remove entry for the post-final-coalescent-event open interval
  coal_factor <- lineages*(lineages-1)/2
  
  event <- c(rep(0, l), rep(1, m))[sorting$ix]
  
  return(list(coal_factor=coal_factor, s=sorting$x, event=event, lineages=lineages))
}



condense_stats <- function(time, event, E, log_zero = -100)
{
  result <- stats::aggregate(event ~ time, FUN = sum)
  result$E <- stats::aggregate(E ~ time, FUN = sum)$E
  
  E_log = log(result$E)
  E_log[result$E == 0] = log_zero
  result$E_log <- E_log
  
  return(result)
}

