

library(adapref)
n <- 200
scez <-seq(1,5)
seedz <- seq(1,100)
comb <- expand.grid(scez,seedz)

#idx <- 1
save1<-c()
save2 <- c()
for (i in 1:500){
print(i)
scenario <- comb[i,1]
seed <- comb[i,2]
set.seed(seed)

source("scenarios.R")

######Sampling ######
coeff=n/0.4
sampling_traj<-function(t){
  result=rep(0,length(t))
  result<-1
  return(result*coeff)}

samp_times1<-c(0,sampsim_thin(max.n=n-2,traj=sampling_traj,xlim=c(0,3))[-n+1])
n_sampled1 <- c(2,rep(1,n-2)) # simulation like this it will be like this.
samp_times2<-c(0,sampsim_thin(max.n=n-2,traj=sampling_traj,xlim=c(0,3))[-n+1])
n_sampled2 <- c(2,rep(1,n-2)) # simulation like this it will be like this.





#########
#### Sampling
###



simulation1<-coalsim(samp_times=samp_times1-min(samp_times1),n=n_sampled1,traj=covid_exp1)
tree1<-sample_genealogy(simulation1)
tree1<-read.tree(text=tree1$Newick)
#plot(tree1,show.tip.label = FALSE)

simulation2<-coalsim(samp_times2-min(samp_times2),n_sampled2, traj = covid_exp2)
tree2<-sample_genealogy(simulation2)
tree2<-read.tree(text=tree2$Newick)
#plot(tree2,show.tip.label = FALSE)

#################
#### SAVE ######
#################
data <- list()


data$n=n
data$scenario=scenario
data$tree1 <- tree1
data$tree2 <- tree2
data$n_sampled1 <- n_sampled1
data$n_sampled2 <- n_sampled2
data$samp_times1 <- samp_times1
data$samp_times2 <- samp_times2
data$coal_times1last <- max(simulation1$coal_times)
data$coal_times2last <- max(simulation2$coal_times)


save1<-c(save1,data$coal_times1last)
save2<-c(save2,data$coal_times2last)

filename <- paste("data/sce_",scenario,"_ite_",seed,".rda",sep="")
save(file=filename,data)
covid_exp1 <- NULL
covid_exp2 <- NULL
}
