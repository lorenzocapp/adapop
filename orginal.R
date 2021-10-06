BNPR <- function (data, lengthout = 100, pref = FALSE, prec_alpha = 0.01, 
          prec_beta = 0.01, beta1_prec = 0.001, fns = NULL, log_fns = TRUE, 
          simplify = TRUE, derivative = FALSE, forward = TRUE) 
{
  if (class(data) == "phylo") {
    phy <- summarize_phylo(data)
  }
  else if (all(c("coal_times", "samp_times", "n_sampled") %in% 
               names(data))) {
    phy <- with(data, list(samp_times = samp_times, coal_times = coal_times, 
                           n_sampled = n_sampled))
  }
  result <- infer_coal_samp(samp_times = phy$samp_times, coal_times = phy$coal_times, 
                            n_sampled = phy$n_sampled, fns = fns, lengthout = lengthout, 
                            prec_alpha = prec_alpha, prec_beta = prec_beta, beta1_prec = beta1_prec, 
                            use_samp = pref, log_fns = log_fns, simplify = simplify, 
                            derivative = derivative)
  result$samp_times <- phy$samp_times
  result$n_sampled <- phy$n_sampled
  result$coal_times <- phy$coal_times
  result$effpop <- exp(-result$result$summary.random$time$`0.5quant`)
  result$effpopmean <- exp(-result$result$summary.random$time$mean)
  result$effpop975 <- exp(-result$result$summary.random$time$`0.025quant`)
  result$effpop025 <- exp(-result$result$summary.random$time$`0.975quant`)
  result$summary <- with(result$result$summary.random$time, 
                         data.frame(time = ID, mean = exp(-mean), sd = sd * exp(-mean), 
                                    quant0.025 = exp(-`0.975quant`), quant0.5 = exp(-`0.5quant`), 
                                    quant0.975 = exp(-`0.025quant`)))
  if (derivative) {
    if (forward) 
      ind <- c(1:(lengthout - 1), (lengthout - 1))
    else ind <- c(1, 1:(lengthout - 1))
    result$derivative <- with(result$result$summary.lincomb, 
                              data.frame(time = result$x, mean = -mean[ind], sd = sd[ind], 
                                         quant0.025 = -`0.975quant`[ind], quant0.5 = -`0.5quant`[ind], 
                                         quant0.975 = -`0.025quant`[ind]))
  }
  if (pref) {
    result$beta0 <- result$result$summary.fixed["beta0", 
                                                "0.5quant"]
    result$beta0summ <- result$result$summary.fixed["beta0", 
                                                    ]
    rownames(result$beta0summ) <- "Beta0"
    result$beta1 <- result$result$summary.hyperpar[2, "0.5quant"]
    result$beta1summ <- result$result$summary.hyperpar[2, 
                                                       ]
    rownames(result$beta1summ) <- "Beta1"
  }
  return(result)
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



infer_coal_samp <- function(samp_times, coal_times, n_sampled=NULL, fns = NULL,
                            lengthout=100, prec_alpha=0.01, prec_beta=0.01,
                            beta1_prec=0.001, use_samp = FALSE, log_fns = TRUE,
                            simplify = FALSE, events_only = FALSE,
                            derivative = FALSE)
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop('INLA needed for this function to work. Use install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE).',
         call. = FALSE)
  }
  
  if (min(coal_times) < min(samp_times))
    stop("First coalescent time occurs before first sampling time")
  
  if (max(samp_times) > max(coal_times))
    stop("Last sampling time occurs after last coalescent time")
  
  grid <- seq(min(samp_times), max(coal_times), length.out = lengthout+1)
  
  if (is.null(n_sampled))
    n_sampled <- rep(1, length(samp_times))
  
  coal_data <- coal_stats(grid = grid, samp_times = samp_times, n_sampled = n_sampled,
                          coal_times = coal_times)
  
  if (simplify){ coal_data <- with(coal_data, condense_stats(time=time, event=event, E=E))}
  
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  
  if (!use_samp)
  {
    data <- with(coal_data, data.frame(y = event, time = time, E_log = E_log))
    
    formula <- y ~ -1 + f(time, model="rw1", hyper = hyper, constr = FALSE)
    family <- "poisson"
  }
  else if (use_samp)
  {
    if (events_only)
      samp_data <- samp_stats(grid = grid, samp_times = samp_times)
    else
      samp_data <- samp_stats(grid = grid, samp_times = samp_times,
                              n_sampled = n_sampled)
    
    data <- joint_stats(coal_data = coal_data, samp_data = samp_data)
    
    if (is.null(fns))
    {
      formula <- Y ~ -1 + beta0 +
        f(time, model="rw1", hyper = hyper, constr = FALSE) +
        f(time2, w, copy="time", fixed=FALSE, param=c(0, beta1_prec))
    }
    else
    {
      vals <- NULL
      bins <- sum(data$beta0 == 0)
      for (fni in fns)
      {
        if (log_fns)
          vals <- cbind(vals, c(rep(0, bins), log(fni(samp_data$time))))
        else
          vals <- cbind(vals, c(rep(0, bins), fni(samp_data$time)))
      }
      data$fn <- vals
      
      formula <- Y ~ -1 + beta0 + fn +
        f(time, model="rw1", hyper = hyper, constr = FALSE) +
        f(time2, w, copy="time", fixed=FALSE, param=c(0, beta1_prec))
    }
    
    family <- c("poisson", "poisson")
  }
  else
    stop("Invalid use_samp value, should be boolean.")
  
  if (derivative)
  {
    Imat <- diag(lengthout)
    A <- utils::head(Imat, -1) - utils::tail(Imat, -1)
    field <- grid[-1] - diff(grid)/2
    A <- diag(1/diff(field)) %*% A
    A[A == 0] <- NA
    
    lc_many <- INLA::inla.make.lincombs(time = A)
  }
  else
  {
    lc_many <- NULL
  }
  
  mod <- INLA::inla(formula, family = family, data = data,
                    lincomb = lc_many, offset = data$E_log,
                    control.predictor = list(compute=TRUE),
                    control.inla = list(lincomb.derived.only=FALSE))
  
  return(list(result = mod, data = data, grid = grid, x = coal_data$time))
}



coal_stats <- function(grid, samp_times, coal_times, n_sampled = NULL,
                       log_zero = -100)
{
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid)/2 #It is practically the mid point of the grid
  
  if (is.null(n_sampled))
    n_sampled <- rep(1, length(samp_times))
  args <- gen_INLA_args(samp_times = samp_times, n_sampled = n_sampled,
                        coal_times = coal_times)
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



samp_stats <- function(grid, samp_times, n_sampled = NULL, trim_end = FALSE)
{
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid)/2 #midpoints
  E <- diff(grid) #Delta of the grid
  
  bins <- cut(x = samp_times, breaks = grid, include.lowest = TRUE) #defines bins
  
  if (is.null(n_sampled))
    count <- as.vector(table(bins))
  else
  {
    tab <- stats::aggregate(n_sampled ~ bins, FUN = sum, labels = FALSE)
    count <- rep(0, lengthout)
    count[as.numeric(tab$bins)] <- tab$n_sampled
  }
  
  count[utils::head(grid, -1) >= max(samp_times)] <- NA
  result <- data.frame(time = field, count = count, E = E, E_log = log(E))
  
  if (trim_end)
    result <- result[stats::complete.cases(result),]
  
  return(result)
}

joint_stats <- function(coal_data, samp_data)
{
  n1 <- length(coal_data$time) #Here are the midpts 
  n2 <- length(samp_data$time) #
  beta0 <- c(rep(0, n1), rep(1, n2))
  E_log <- c(coal_data$E_log, samp_data$E_log) #samp_data$E_log does not include NA
  Y <- matrix(c(coal_data$event, rep(NA, n2), rep(NA, n1), samp_data$count),
              nrow = n1 + n2, byrow = FALSE)
  #Note, there are some parts of Y where they are both NAs
  w <- c(rep(1, n1), rep(-1, n2))  #what does the w do? I think that takes into account for the fact that one is at the numerator,
  #the other one at the numerator. 
  time  <- c(coal_data$time, rep(NA, n2))
  time2 <- c(rep(NA, n1), samp_data$time)
  
  return(list(Y = Y, beta0 = beta0, time = time, time2 = time2, w = w, E_log = E_log))
}
