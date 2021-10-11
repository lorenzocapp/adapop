
#dir <- '/Users/lorenzocappello/Google Drive/Statistics/postdoc Palacios/delta_covid_attempt/selection_model_delta_simulations/'


#setwd(dir)
#setwd("/Users/cappello/Google Drive/postdoc Palacios/stanford_covid/selection_model_trial")

source("functions_inla_selection.R") #where the new functions should be stored 
library("phylodyn")
library("ape")
INLA:::inla.dynload.workaround()

#Parallel variable.
arg<-commandArgs(TRUE)
idx<-as.integer(arg[1])


scez <-seq(1,5)
seedz <- seq(1,100)
precz <- c(0.000001,10)
comb <- expand.grid(scez,seedz,precz)

###########################
###Load the two trees ####
##########################

scenario <- comb[idx,1]
seed <- comb[idx,2]
prec <- comb[idx,3]

filename <- paste("data/sce_",scenario,"_ite_",seed,".rda",sep="")
load(filename)

file.exists(filename)

n <- data$n
scenario <- data$scenario
tree1 <- data$tree1 
tree2 <- data$tree2 
n_sampled1 <- data$n_sampled1
n_sampled2 <- data$n_sampled2
samp_times1 <- data$samp_times1
samp_times2 <- data$samp_times2 
samp_times1 <- samp_times1 - min(samp_times1)
samp_times2 <- samp_times2 - min(samp_times2)
coal_times1last <- data$coal_times1last
coal_times2last <- data$coal_times2last


###Start saving the output 

seed <- idx
n <- 200
scenario <- data$scenario
samp <- "unif"

toadd <- c(seed,n,scenario,samp,prec)

#####Load scenario

source("scenarios.R")



df <- rep(NA,15) 
max2<-0.6*coal_times2last
max1 <- 0.6*coal_times1last
if (max2>1) max2 <-1
if (max1>1) max1 <-1
#############################
####     adaSel noPref  #####
#############################

model <- "adasel_noPref"

res<-BNPR_sel(tree1,tree2,samp_times1,samp_times2,lengthout=100,parSel = F,
                  preferential =F,prec_alpha = prec,prec_beta = prec)

#Selection coeff
table_out<-cbind(res$selIntsummary$time,res$selInt,res$selInt975,res$selInt025)
fun<-function(x){return(covid_exp2(x)/covid_exp1(x))}
env_sel<-envelope(table_out,fun,max2,0,Inf,0)$avg
dev_sel<-dev(table_out,fun,max2,0,Inf,0)$avg
rwd_sel<-relwid(table_out,fun,max2,0,Inf,0)$avg

#Eff1 
table_out<-cbind(res$summary$time,res$effpop,res$effpop975,res$effpop025)
fun<-function(x){covid_exp1(x)}
env_eff1<-envelope(table_out,fun,max1,0,Inf,0)$avg
dev_eff1<-dev(table_out,fun,max1,0,Inf,0)$avg
rwd_eff1<-relwid(table_out,fun,max1,0,Inf,0)$avg


#Eff2 
result <- eff2_adasel(res,preferential=F)
table_out<-cbind(result$effpop2summary$time,result$effpop2,result$effpop2_975,result$effpop2_025)
fun<-function(x){covid_exp2(x)}
env_eff2<-envelope(table_out,fun,max2,0,Inf,0)$avg
dev_eff2<-dev(table_out,fun,max2,0,Inf,0)$avg
rwd_eff2<-relwid(table_out,fun,max2,0,Inf,0)$avg

new <- c(toadd,model,env_sel,dev_sel,rwd_sel,env_eff1,dev_eff1,rwd_eff1,env_eff2,dev_eff2,rwd_eff2)
df <- rbind(df,new)



#############################
####     adaSel adaPref  #####
#############################

model <- "adasel_adaPref"

res<-BNPR_sel(tree1,tree2,samp_times1,samp_times2,lengthout=100,parSel = F,
              preferential =T,prec_alpha = prec,prec_beta = prec)

#Selection coeff
table_out<-cbind(res$selIntsummary$time,res$selInt,res$selInt975,res$selInt025)
fun<-function(x){return(covid_exp2(x)/covid_exp1(x))}
env_sel<-envelope(table_out,fun,max2,0,Inf,0)$avg
dev_sel<-dev(table_out,fun,max2,0,Inf,0)$avg
rwd_sel<-relwid(table_out,fun,max2,0,Inf,0)$avg

#Eff1 
table_out<-cbind(res$summary$time,res$effpop,res$effpop975,res$effpop025)
fun<-function(x){covid_exp1(x)}
env_eff1<-envelope(table_out,fun,max1,0,Inf,0)$avg
dev_eff1<-dev(table_out,fun,max1,0,Inf,0)$avg
rwd_eff1<-relwid(table_out,fun,max1,0,Inf,0)$avg


#Eff2 
result <- eff2_adasel(res,preferential=T)
table_out<-cbind(result$effpop2summary$time,result$effpop2,result$effpop2_975,result$effpop2_025)
fun<-function(x){covid_exp2(x)}
env_eff2<-envelope(table_out,fun,max2,0,Inf,0)$avg
dev_eff2<-dev(table_out,fun,max2,0,Inf,0)$avg
rwd_eff2<-relwid(table_out,fun,max2,0,Inf,0)$avg

new <- c(toadd,model,env_sel,dev_sel,rwd_sel,env_eff1,dev_eff1,rwd_eff1,env_eff2,dev_eff2,rwd_eff2)
df <- rbind(df,new)



#############################
####     parSel noPref  #####
#############################

model <- "parsel_noPref"

res<-BNPR_sel(tree1,tree2,samp_times1,samp_times2,lengthout=100,parSel = T,
              preferential =F,prec_alpha = prec,prec_beta = prec)

#Selection coeff
result <- selInt_parametric(res)
table_out<-cbind(result$selInt_parsummary$time,result$selInt_par,result$selInt_par_975,result$selInt_par_025)
fun<-function(x){return(covid_exp2(x)/covid_exp1(x))}
env_sel<-envelope(table_out,fun,max2,0,Inf,0)$avg
dev_sel<-dev(table_out,fun,max2,0,Inf,0)$avg
rwd_sel<-relwid(table_out,fun,max2,0,Inf,0)$avg

#Eff1 
table_out<-cbind(res$summary$time,res$effpop,res$effpop975,res$effpop025)
fun<-function(x){covid_exp1(x)}
env_eff1<-envelope(table_out,fun,max1,0,Inf,0)$avg
dev_eff1<-dev(table_out,fun,max1,0,Inf,0)$avg
rwd_eff1<-relwid(table_out,fun,max1,0,Inf,0)$avg


#Eff2 
result <- eff2_parametric(res)
table_out<-cbind(result$effpop2summary$time,result$effpop2,result$effpop2_975,result$effpop2_025)
fun<-function(x){covid_exp2(x)}
env_eff2<-envelope(table_out,fun,max2,0,Inf,0)$avg
dev_eff2<-dev(table_out,fun,max2,0,Inf,0)$avg
rwd_eff2<-relwid(table_out,fun,max2,0,Inf,0)$avg

new <- c(toadd,model,env_sel,dev_sel,rwd_sel,env_eff1,dev_eff1,rwd_eff1,env_eff2,dev_eff2,rwd_eff2)
df <- rbind(df,new)




#############################
####     noSel noPref  #####
#############################

model <- "nosel_noPref"


res1<-BNPR(tree1,lengthout=100, grid=res$grid,prec_alpha = prec,prec_beta = prec)
res2<-BNPR(tree2,lengthout=100, grid=res$grid,prec_alpha = prec,prec_beta = prec)



#Selection coeff
result <- selInt_independent(res1,res2)
table_out<-cbind(result$selInt_indsummary$time,result$selInt_ind,result$selInt_ind_975,result$selInt_ind_025)
fun<-function(x){return(covid_exp2(x)/covid_exp1(x))}
env_sel<-envelope(table_out,fun,max2,0,Inf,0)$avg
dev_sel<-dev(table_out,fun,max2,0,Inf,0)$avg
rwd_sel<-relwid(table_out,fun,max2,0,Inf,0)$avg

#Eff1 
table_out<-cbind(res1$summary$time,res1$effpop,res1$effpop975,res1$effpop025)
fun<-function(x){covid_exp1(x)}
env_eff1<-envelope(table_out,fun,max1,0,Inf,0)$avg
dev_eff1<-dev(table_out,fun,max1,0,Inf,0)$avg
rwd_eff1<-relwid(table_out,fun,max1,0,Inf,0)$avg


#Eff2 
table_out<-cbind(res2$summary$time,res2$effpop,res2$effpop975,res2$effpop025)
fun<-function(x){covid_exp2(x)}
env_eff2<-envelope(table_out,fun,max2,0,Inf,0)$avg
dev_eff2<-dev(table_out,fun,max2,0,Inf,0)$avg
rwd_eff2<-relwid(table_out,fun,max2,0,Inf,0)$avg

new <- c(toadd,model,env_sel,dev_sel,rwd_sel,env_eff1,dev_eff1,rwd_eff1,env_eff2,dev_eff2,rwd_eff2)
df <- rbind(df,new)




###############################
####### Save results  #########
###############################


filename <- paste("output/idx_",idx,".rda",sep="")
save(file=filename,output)








