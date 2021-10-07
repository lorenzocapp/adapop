BNPR_sel <- function (tree1,tree2, samp_times1,samp_times2, lengthout = 100,
                      prec_alpha = 0.01, 
                      prec_beta = 0.01, 
                      beta1_prec = 0.001, 
                      parSel=TRUE, 
                      preferential=FALSE,
                      u.truncation=NA,
                      l.truncation=NA,
                      beta0_remove=FALSE,
                      simplify = TRUE, 
                      derivative = FALSE, 
                      forward = TRUE) 
{
  
  if (min(samp_times1)>0){
    warning("Warning: group 1 was not sampled at 0, model may not be identifiable + 
            the current implementation is suboptimal in that case",immediate. = T)
  }
  
  #If parSel is FALSE, I am using adaSel
  phy1 <- summarize_phylo(tree1)
  phy1$samp_times <- phy1$samp_times + min(samp_times1)
  phy1$coal_times <- phy1$coal_times + min(samp_times1)
  phy2 <- summarize_phylo(tree2)
  phy2$samp_times <- phy2$samp_times + min(samp_times2)
  phy2$coal_times <- phy2$coal_times + min(samp_times2)
  
  
  result <- infer_coal_samp_selection(phy1,phy2,lengthout = lengthout, 
                                      prec_alpha = prec_alpha, prec_beta = prec_beta,
                                      beta1_prec = beta1_prec, 
                                      parSel, preferential,
                                      u.truncation,l.truncation,
                                      simplify,beta0_remove)


  
  #result$samp_times <- phy$samp_times
  #result$n_sampled <- phy$n_sampled
  #result$coal_times <- phy$coal_times
  if (!preferential){
      result$effpop <- exp(-result$result$summary.random$time$`0.5quant`)
      result$effpopmean <- exp(-result$result$summary.random$time$mean)
      result$effpop975 <- exp(-result$result$summary.random$time$`0.025quant`)
      result$effpop025 <- exp(-result$result$summary.random$time$`0.975quant`)
      result$summary <- with(result$result$summary.random$time, 
                             data.frame(time = ID, mean = exp(-mean), sd = sd * exp(-mean), 
                                        quant0.025 = exp(-`0.975quant`), quant0.5 = exp(-`0.5quant`), 
                                        quant0.975 = exp(-`0.025quant`)))
    } else {
    result$effpop <- exp(-result$result$summary.random$time1$`0.5quant`)
    result$effpopmean <- exp(-result$result$summary.random$time1$mean)
    result$effpop975 <- exp(-result$result$summary.random$time1$`0.025quant`)
    result$effpop025 <- exp(-result$result$summary.random$time1$`0.975quant`)
    result$summary <- with(result$result$summary.random$time1, 
                           data.frame(time = ID, mean = exp(-mean), sd = sd * exp(-mean), 
                                      quant0.025 = exp(-`0.975quant`), quant0.5 = exp(-`0.5quant`), 
                                      quant0.975 = exp(-`0.025quant`)))
    result$sampInt <- exp(-result$result$summary.random$time3b$`0.5quant`)
    result$sampIntmean <- exp(-result$result$summary.random$time3b$mean)
    result$sampInt975 <- exp(-result$result$summary.random$time3b$`0.975quant`)
    result$sampInt025 <- exp(-result$result$summary.random$time3b$`0.025quant`)
    result$sampIntsummary <- with(result$result$summary.random$time3b,
                                 data.frame(time = ID, mean = exp(mean), sd = sd * exp(mean),
                                            quant0.025 = exp(`0.025quant`), quant0.5 = exp(`0.5quant`),
                                            quant0.975 = exp(`0.975quant`)))
    }
  
    if (parSel){
      if (beta0_remove==FALSE){
        result$beta0 <- result$result$summary.fixed["beta0", "0.5quant"]
        result$beta0summ <- result$result$summary.fixed["beta0",]
        rownames(result$beta0summ) <- "Beta0"
        result$beta0post <- result$result$marginals.fixed$beta0
      }
      result$beta1 <- result$result$summary.hyperpar[2, "0.5quant"]
      result$beta1summ <- result$result$summary.hyperpar[2,]
      rownames(result$beta1summ) <- "Beta1"
      result$beta1post <- result$result$marginals.hyperpar$`Beta for time2`
    } else {
      result$selInt <- exp(-result$result$summary.random$time2b$`0.5quant`)
      result$selIntmean <- exp(-result$result$summary.random$time2b$mean)
      result$selInt975 <- exp(-result$result$summary.random$time2b$`0.975quant`)
      result$selInt025 <- exp(-result$result$summary.random$time2b$`0.025quant`)
      result$selIntsummary <- with(result$result$summary.random$time2b,
                                   data.frame(time = ID, mean = exp(mean), sd = sd * exp(mean),
                                              quant0.025 = exp(`0.025quant`), quant0.5 = exp(`0.5quant`),
                                              quant0.975 = exp(`0.975quant`)))
    }  
  
  
  
  return(result)
}





infer_coal_samp_selection <- function(phy1,phy2, lengthout=100, prec_alpha=0.01, prec_beta=0.01, 
                                      beta1_prec=0.001, parSel=TRUE,preferential=FALSE,
                                      u.truncation=NA,l.truncation=NA, simplify = TRUE,
                                      beta0_remove=FALSE)
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
  if (simplify){coal_data1 <- with(coal_data1, condense_stats(time=time, event=event, E=E))}
  coal_data1 <- identify_off_grid(coal_data1,samp_times1,coal_times1)
  
  coal_data2 <- coal_stats(grid = grid, samp_times = samp_times2, n_sampled = n_sampled2,
                           coal_times = coal_times2)
  if (simplify){coal_data2 <- with(coal_data2, condense_stats(time=time, event=event, E=E))}
  coal_data2 <- identify_off_grid(coal_data2,samp_times2,coal_times2)
  
  
  if (preferential){
    samp_data1 <- samp_stats_adasel(grid = grid, samp_times = samp_times1,
                                n_sampled = n_sampled1)
    samp_data2 <- samp_stats_adasel(grid = grid, samp_times = samp_times2,
                                 n_sampled = n_sampled2)
  }
  
  
  
  
  #If I am using a truncation, this is the place that adjust for it
  if (!is.na(u.truncation)){
    id <- min(which(coal_data1$time>=u.truncation))
    coal_data1 <- coal_data1[1:id,]
    coal_data2 <- coal_data2[1:id,]
    grid <- grid[1:(id+1)]
  }
  if (!is.na(l.truncation)){
    id <- max(which(coal_data1$time<=l.truncation))
    coal_data1 <- coal_data1[id:length(coal_data1$time),]
    coal_data2 <- coal_data2[id:length(coal_data1$time),]
    grid <- grid[(id-1):length(grid)]
  }
  

  
  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
  
  
  if (!preferential){
    data <- joint_coal_stats(coal_data1 = coal_data1, coal_data2 = coal_data2,
                             samp_times2=samp_times2,grid=grid)
  } else {
    data <- joint_coal_stats_adasel(coal_data1 = coal_data1, coal_data2 = coal_data2,
                                    samp_data1 = samp_data1, samp_data2 = samp_data2,
                                    grid=grid,samp_times2=samp_times2)
  
  }
  if (parSel){
    if (!preferential){
        if (beta0_remove){
          formula <- Y ~ -1 + 
            f(time, model="rw1", hyper = hyper, constr = FALSE) +
            f(time2,copy="time", fixed=FALSE, param=c(0, beta1_prec))
        } else {
          formula <- Y ~ -1 + beta0 +
            f(time, model="rw1", hyper = hyper, constr = FALSE) +
            f(time2,copy="time", fixed=FALSE, param=c(0, beta1_prec))
        }
    } else {
      data$time3b=data$time3
      data$time4b=data$time4
      formula <- Y ~ -1 + beta0 +
        f(time1, model="rw1", hyper = hyper, constr = FALSE) +
        f(time2,copy="time1", fixed=FALSE, param=c(0, beta1_prec)) +
        f(time3,w,copy="time1") + f(time3b, model="rw1", hyper = hyper, constr=FALSE) +
        f(time4,w,copy="time2") + f(time4b, w, copy="time3b")
    }
  } else {
    if (!preferential){
      data$beta0<-NULL
      data$time2b=data$time2
      formula <- Y ~ -1 +
        f(time, model="rw1", hyper = hyper, constr = FALSE) +
        f(time2, w, copy="time") + f(time2b, model="rw1",hyper = hyper, constr=FALSE)
    } else {
      data$beta0<-NULL
      data$time2b=data$time2
      data$time3b=data$time3
      data$time4b=data$time4
      data$time4c=data$time4
      formula <- Y ~ -1 +
        f(time1, model="rw1", hyper = hyper, constr = FALSE) +
        f(time2, w, copy="time1") + f(time2b, model="rw1", hyper = hyper, constr=FALSE) +
        f(time3, w, copy="time1") + f(time3b, model="rw1", hyper = hyper, constr=FALSE) +
        f(time4, w, copy="time1") + f(time4b, w, copy="time2b") + f(time4c,copy="time3b")
      
    }
 
  }
  
  if(!preferential){
    family <- c("poisson", "poisson")
  } else {
    family <- c("poisson", "poisson","poisson","poisson")
  }

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


joint_coal_stats <- function(coal_data1, coal_data2,samp_times2,grid)
{
  
  ##Here, I ensure that the field on the second population starts from the minimum sampling points
  id.field <- max(which(utils::tail(grid,-1) <=  min(samp_times2)),0) + 1
  coal_data2 <- coal_data2[(id.field+1):nrow(coal_data2),]

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


joint_coal_stats_adasel <- function(coal_data1, coal_data2,samp_data1,samp_data2,grid,samp_times2)
{

  id.field2 <- max(which(utils::tail(grid,-1) <=  min(samp_times2)),0) + 1
  coal_data2 <- coal_data2[(id.field2+1):nrow(coal_data2),]
  samp_data2 <- samp_data2[(id.field2+1):nrow(samp_data2),]
  

  n1 <- length(coal_data1$time) #Here are the midpts 
  n2 <- length(coal_data2$time) #
  n3 <- length(samp_data1$time)
  n4 <- length(samp_data2$time)
  beta0 <- c(rep(0, n1), rep(1, n2),rep(0,n3),rep(1,n4))
  E_log <- c(coal_data1$E_log, coal_data2$E_log,samp_data1$E_log,samp_data2$E_log) #samp_data$E_log does not include NA in the pref_samp, so I am not modifying it here either
  c1 <- c(coal_data1$event, rep(NA, n2+n3+n4))
  c2 <- c(rep(NA, n1), coal_data2$event, rep(NA,n3+n4))
  c3 <- c(rep(NA,n1+n2),samp_data1$count, rep(NA,n4))
  c4 <- c(rep(NA,n1+n2+n3), samp_data2$count)
  Y <- cbind(c1,c2,c3,c4)
  #Need to correct Y for the parts of the grid that need to be eccluded
  #Note, there are some parts of Y where they are both NAs
  w <- c(rep(1, n1), rep(1, n2),rep(-1,n3),rep(-1,n4))  #what does the w do? I think that takes into account for the fact that one is at the numerator,
  #the other one at the numerator. 
  
  #we are not correcting the time, i.e. there are times also to part that are not defined
  time1 <- c(coal_data1$time, rep(NA, n2+n3+n4))
  time2 <- c(rep(NA, n1), coal_data2$time, rep(NA,n3+n4))
  time3 <- c(rep(NA,n1+n2),samp_data1$time, rep(NA,n4))
  time4 <- c(rep(NA,n1+n2+n3), samp_data2$time)
  
  
  return(list(Y = Y, beta0 = beta0, time1 = time1, time2 = time2, time3=time3,time4=time4, w = w, E_log = E_log))
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


plot_adaSel<-function (BNPR_out, traj = NULL, xlim = NULL, ylim = NULL, nbreaks = 40,
                     lty = 1, lwd = 2, col = "black", main = "", log = "y", ylab = "Sampling Intensity",
                     xlab = "Time", xmarline = 3, axlabs = NULL, traj_lty = 2,
                     traj_lwd = 2, traj_col = col, newplot = TRUE, credible_region = TRUE,
                     heatmaps = TRUE, heatmap_labels = TRUE, heatmap_labels_side = "right",
                     heatmap_labels_cex = 0.7, heatmap_width = 7, yscale = 1, max_samp,
                     ...)
{
  grid = BNPR_out$grid
  if (is.null(xlim)) {
    xlim = c(max(grid), min(grid))
  }
  mask = BNPR_out$x >= min(xlim) & BNPR_out$x <= max(xlim)
  t = BNPR_out$x[mask]
  #t=decimal_date(max_samp)-t
  xlim_dec=xlim
  #xlim=(c(decimal_date(max_samp)-xlim[1],decimal_date(max_samp)-xlim[2]))
  y = BNPR_out$selInt[mask] * yscale
  yhi = BNPR_out$selInt975[mask] * yscale
  ylo = BNPR_out$selInt025[mask] * yscale
  if (newplot) {
    if (is.null(ylim)) {
      ymax = max(yhi)
      ymin = min(ylo)
    }
    else {
      ymin = min(ylim)
      ymax = max(ylim)
    }
    if (heatmaps) {
      yspan = ymax/ymin
      yextra = yspan^(1/10)
      ylim = c(ymin/(yextra^1.35), ymax)
    }
    else {
      ylim = c(ymin, ymax)
    }
    if (is.null(axlabs)) {
      graphics::plot(1, 1, type = "n", log = log, xlab = xlab,
                     ylab = ylab, main = main, xlim = xlim, ylim = ylim,
                     ...)
    }
    else {
      graphics::plot(1, 1, type = "n", log = log, xlab = "",
                     ylab = ylab, main = main, xlim = xlim, ylim = ylim,
                     xaxt = "n", ...)
      graphics::axis(1, at = axlabs$x, labels = axlabs$labs, cex.axis=1,
                     las = 1)
      graphics::mtext(text = xlab, side = 1, line = xmarline)
    }
  }
  if (credible_region) {
    shade_band(x = t, ylo = ylo, yhi = yhi, col = "lightgray")
  }
  if (!is.null(traj)) {
    graphics::lines(t, traj(t), lwd = traj_lwd, lty = traj_lty,
                    col = traj_col)
  }
  if (newplot) {
    if (heatmaps) {
      samps = rep(BNPR_out$samp_times, BNPR_out$n_sampled)
      #samps = decimal_date(max_samp)-samps[samps <= max(xlim_dec) & samps >= min(xlim_dec)]
      samps = samps[samps <= max(xlim_dec) & samps >= min(xlim_dec)]
      coals = BNPR_out$coal_times
      coals = coals[coals <= max(xlim) & coals >= min(xlim)]
      breaks = seq(min(xlim_dec), max(xlim_dec), length.out = nbreaks)
      h_samp = graphics::hist(samps, breaks = breaks, plot = FALSE)
      #h_coal = graphics::hist(coals, breaks = breaks, plot = FALSE)
      hist2heat(h_samp, y = ymin/yextra^0.5, wd = heatmap_width)
      #hist2heat(h_coal, y = ymin/yextra, wd = heatmap_width)
      if (heatmap_labels) {
        if (heatmap_labels_side == "left") {
          lab_x = max(xlim)
          lab_adj = 0
        }
        else if (heatmap_labels_side == "right") {
          lab_x = min(xlim)
          lab_adj = 1
        }
        else {
          warning("heatmap_labels_side not \"left\" or \"right\", defaulting to right")
          lab_x = min(xlim)
          lab_adj = 1
        }
        graphics::text(x = lab_x, y = ymin/(yextra^0.2)+0.0075,
                       labels = "Sampling events", adj = c(lab_adj,
                                                           0), cex = heatmap_labels_cex)
        #graphics::text(x = lab_x, y = ymin/(yextra^1.25),
        #   labels = "Coalescent events", adj = c(lab_adj,
        # 1), cex = heatmap_labels_cex)
      }
    }
  }
  graphics::lines(t, y, lwd = lwd, col = col, lty = lty)
}



plot_BNPR_plus <-function (BNPR_out, traj = NULL, xlim = NULL, ylim = NULL, nbreaks = 40,
                       lty = 1, lwd = 2, col = "black", main = "", log = "y", ylab = "Effective pop size group 1",
                       xlab = "Time", xmarline = 3, axlabs = NULL, traj_lty = 2, parameter="eff1",
                       traj_lwd = 2, traj_col = col, newplot = TRUE, credible_region = TRUE,
                       heatmaps = TRUE, heatmap_labels = TRUE, heatmap_labels_side = "right",
                       heatmap_labels_cex = 0.7, heatmap_width = 7, yscale = 1, max_samp,
                       ...)
{
  grid = BNPR_out$grid
  if (is.null(xlim)) {
    xlim = c(max(grid), min(grid))
  }

  xlim_dec=xlim
  if (parameter=="eff1"){
    mask = BNPR_out$summary$time >= min(xlim) & BNPR_out$summary$time <= max(xlim)
    t = BNPR_out$summary$time[mask]
    y = BNPR_out$effpop[mask] * yscale
    yhi = BNPR_out$effpop975[mask] * yscale
    ylo = BNPR_out$effpop025[mask] * yscale
  } else if (parameter=="eff2"){
   stop("Not available yet in this version of the package! Just flip the order of the populations!")
  } else if (parameter=="sampInt"){
    mask = BNPR_out$sampIntsummary$time >= min(xlim) & BNPR_out$sampIntsummary$time <= max(xlim)
    t = BNPR_out$sampIntsummary$time[mask]
    y = BNPR_out$sampInt[mask] * yscale
    yhi = BNPR_out$sampInt975[mask] * yscale
    ylo = BNPR_out$sampInt025[mask] * yscale
    ylab= "Sampling Intesity"
  } else if (parameter=="selInt"){
    mask = BNPR_out$selIntsummary$time >= min(xlim) & BNPR_out$selIntsummary$time <= max(xlim)
    t = BNPR_out$selIntsummary$time[mask]
    y = BNPR_out$selInt[mask] * yscale
    yhi = BNPR_out$selInt975[mask] * yscale
    ylo = BNPR_out$selInt025[mask] * yscale
    ylab = "Selective intensity"
  }
  
  if (newplot) {
    if (is.null(ylim)) {
      ymax = max(yhi)
      ymin = min(ylo)
    }
    else {
      ymin = min(ylim)
      ymax = max(ylim)
    }
    if (heatmaps) {
      yspan = ymax/ymin
      yextra = yspan^(1/10)
      ylim = c(ymin/(yextra^1.35), ymax)
    }
    else {
      ylim = c(ymin, ymax)
    }
    if (is.null(axlabs)) {
      graphics::plot(1, 1, type = "n", log = log, xlab = xlab,
                     ylab = ylab, main = main, xlim = xlim, ylim = ylim,
                     ...)
    }
    else {
      graphics::plot(1, 1, type = "n", log = log, xlab = "",
                     ylab = ylab, main = main, xlim = xlim, ylim = ylim,
                     xaxt = "n", ...)
      graphics::axis(1, at = axlabs$x, labels = axlabs$labs, cex.axis=1,
                     las = 1)
      graphics::mtext(text = xlab, side = 1, line = xmarline)
    }
  }
  if (credible_region) {
    shade_band(x = t, ylo = ylo, yhi = yhi, col = "lightgray")
  }
  if (!is.null(traj)) {
    graphics::lines(t, traj(t), lwd = traj_lwd, lty = traj_lty,
                    col = traj_col)
  }
  if (newplot) {
    if (heatmaps) {
      samps = rep(BNPR_out$samp_times, BNPR_out$n_sampled)
      #samps = decimal_date(max_samp)-samps[samps <= max(xlim_dec) & samps >= min(xlim_dec)]
      samps = samps[samps <= max(xlim_dec) & samps >= min(xlim_dec)]
      coals = BNPR_out$coal_times
      coals = coals[coals <= max(xlim) & coals >= min(xlim)]
      breaks = seq(min(xlim_dec), max(xlim_dec), length.out = nbreaks)
      h_samp = graphics::hist(samps, breaks = breaks, plot = FALSE)
      #h_coal = graphics::hist(coals, breaks = breaks, plot = FALSE)
      hist2heat(h_samp, y = ymin/yextra^0.5, wd = heatmap_width)
      #hist2heat(h_coal, y = ymin/yextra, wd = heatmap_width)
      if (heatmap_labels) {
        if (heatmap_labels_side == "left") {
          lab_x = max(xlim)
          lab_adj = 0
        }
        else if (heatmap_labels_side == "right") {
          lab_x = min(xlim)
          lab_adj = 1
        }
        else {
          warning("heatmap_labels_side not \"left\" or \"right\", defaulting to right")
          lab_x = min(xlim)
          lab_adj = 1
        }
        graphics::text(x = lab_x, y = ymin/(yextra^0.2)+0.0075,
                       labels = "Sampling events", adj = c(lab_adj,
                                                           0), cex = heatmap_labels_cex)
        #graphics::text(x = lab_x, y = ymin/(yextra^1.25),
        #   labels = "Coalescent events", adj = c(lab_adj,
        # 1), cex = heatmap_labels_cex)
      }
    }
  }
  graphics::lines(t, y, lwd = lwd, col = col, lty = lty)
}

###########From other packages

shade_band = function(x, ylo, yhi, xlim=NULL, col="gray")
{
  if (is.null(xlim))
    xlim = c(0, Inf)
  mask = x >= min(xlim) & x <= max(xlim)
  
  x = x[mask]
  ylo = ylo[mask]
  yhi = yhi[mask]
  
  graphics::polygon(c(x, rev(x)), c(yhi, rev(ylo)), col=col, border=NA)
}

float2gray = function(f)
{
  return(sprintf("#%02x%02x%02x", floor((1-f) * 255), floor((1-f) * 255), floor((1-f) * 255)))
}
hist2heat = function(hist, y, wd)
{
  breaks = hist$breaks
  counts = hist$counts
  upper  = max(counts)
  n = length(counts)
  cols = float2gray(counts / upper)
  graphics::segments(x0 = breaks[1:n], y0=y, x1 = breaks[2:(n+1)], y1=y, lwd=wd, col=cols, lend=1)
}



samp_stats_adasel <- function(grid, samp_times, n_sampled = NULL, trim_end = FALSE)
{
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid)/2
  E <- diff(grid)
  
  bins <- cut(x = samp_times, breaks = grid, include.lowest = TRUE)
  
  if (is.null(n_sampled))
    count <- as.vector(table(bins))
  else
  {
    tab <- stats::aggregate(n_sampled ~ bins, FUN = sum, labels = FALSE)
    count <- rep(0, lengthout)
    count[as.numeric(tab$bins)] <- tab$n_sampled
  }
  
  count[utils::head(grid, -1) >= max(samp_times)] <- NA
  count[utils::tail(grid,-1) <=  min(samp_times)] <- NA #Double Check
  result <- data.frame(time = field, count = count, E = E, E_log = log(E))
  
  if (trim_end) #I am not sure what's this is supposed to do
    result <- result[stats::complete.cases(result),]
  

  return(result)
}

