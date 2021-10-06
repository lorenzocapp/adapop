setwd("/Users/lorenzocappello/Google Drive/Statistics/postdoc Palacios/delta_covid_attempt/selection_model_delta_truncated/")

source("functions_inla_selection_truncated.R") #where the new functions should be stored 
source("tajima_inference.R")
library("phylodyn")
library("ape")
library("lubridate")


##################################
###Extract the s and n vector ####
##################################


dev.off()
treeDelta<-read.nexus("../beast2_ver_1&2/tree_singapore_delta")
#treeDelta<-read.nexus("../beast2_ver_3/tree_uk_2_delta")

plot(treeDelta,show.tip.label=FALSE)
axisPhylo()


treeNotdelta<-read.nexus("../beast2_ver_1&2/tree_singapore_notdelta")
#treeNotdelta<-read.nexus("../beast2_ver_3/tree_uk_2_notdelta")

#treeNotdelta<-read.nexus("../beast2/tree_italy_notdelta")

#tree not delta refinement
#treeD <-drop.tip(tree,treeG$tip.label) # This is the D clade
treeNotdelta <- extract.clade(treeNotdelta,153)

plot(treeNotdelta,show.tip.label = FALSE)
axisPhylo()
#plot(treeNotdelta)


data1<-treeDelta
data2 <- treeNotdelta

samp_times1 <- c()
for (label in data1$tip.label){samp_times1 <- c(samp_times1,strsplit(label,"[|]")[[1]][3])}
samp_times1 <- decimal_date(date(samp_times1))
samp_times2 <- c()
for (label in data2$tip.label){samp_times2 <- c(samp_times2,strsplit(label,"[|]")[[1]][3])}
samp_times2 <- decimal_date(date(samp_times2))
lastdate <- max(c(samp_times1,samp_times2))
samp_times1 <- lastdate-samp_times1
samp_times2 <- lastdate-samp_times2






#################################################
### Simulate data from the same pop scenario ####
#################################################


covid_exp_traj1<-function(t){
  result=rep(0,length(t))
  result[t<=0.35]<-10
  result[t>0.35 & t<0.4]<-9.999987e+14*exp(-92.1034*t[t>0.35 & t<0.4]) #50 times smaller
  result[t>=0.4]<-0.1
  return(result)
}


covid_exp_traj2<-function(t){
  result=rep(0,length(t))
  result[t<=0.1696721]<-10
  result[t>0.1696721 & t<0.2196721]<-61218647*exp(-92.1034*t[t>0.1696721 & t<0.2196721]) #50 times smaller
  result[t>=0.2196721]<-0.1
  return(result)
}


covid_exp2<-function(t){
  t=t+0.03825137
  result=rep(0,length(t))
  result[t<=0.35]<-10
  result[t>0.35 & t<0.4]<-9.999987e+14*exp(-92.1034*t[t>0.35 & t<0.4]) #50 times smaller
  result[t>=0.4]<-0.1
  result <- sqrt(result)
  return(result)
}


#why the trajectory
# x<-seq(0,0.5,0.01)
# res1a<-BNPR(data,lengthout=100)
dev.off()
plot_BNPR(res1a,heatmaps=TRUE,heatmap_labels = FALSE,ylab="Ne(t)",main="Effective Pop Size",ylim=c(0.005,50),xlim=c(0.54,0))
lines(x,covid_exp_traj1(x),col="red",lwd=2,type="l")
lines(x,covid_exp2(x),col="blue",lwd=2)

covid_exp1 <- function(t){
  return(covid_exp_traj3,min(samp_times1))
}



simulationG<-coalsim(samp_times=samp_times1,n=n_sampled1,traj=covid_exp2)
treeg<-sample_genealogy(simulationG)
treeGG<-read.tree(text=treeg$Newick)
plot(treeGG,show.tip.label = FALSE)


simulationD<-coalsim(samp_times2-min(samp_times2),n_sampled2, traj = covid_exp_traj1)
treed<-sample_genealogy(simulationD)
treeDD<-read.tree(text=treed$Newick)
plot(treeDD,show.tip.label = FALSE)


# res1a<-BNPR(treeGG,lengthout=100)
# plot_BNPR(res1a,heatmaps=TRUE,heatmap_labels = FALSE,ylab="Ne(t)",main="Effective Pop Size",ylim=c(0.005,50),xlim=c(0.54,0))
# 
# res1a<-BNPR(treeDD,lengthout=100)
# plot_BNPR(res1a,heatmaps=TRUE,heatmap_labels = FALSE,ylab="Ne(t)",main="Effective Pop Size",ylim=c(0.005,50),xlim=c(0.54,0))


res1s <- BNPR_sel(treeGG,treeDD,samp_times1,samp_times2,lengthout = 100,beta0_remove = FALSE)
plot_BNPR(res1s,heatmaps=FALSE,heatmap_labels = FALSE,ylab="Ne(t)",main="Effective Pop Size",ylim=c(0.005,50),xlim=c(0.5,0))

par(mfrow=c(1,2))
plot(res1s$beta0post,type="l",lwd=3,main="Posterior beta0")
abline(v=0,col="red",lwd=2)
plot(res1s$beta1post,type="l",lwd=3,,main="Posterior beta1")
abline(v=1,col="red",lwd=2)

dev.off()
x<-rev(seq(0,0.5,0.001))
plot(x,covid_exp_traj1(x),type="l",lwd=2,xlim=c(0.5,0),log="y")
lines(x,sqrt(covid_exp_traj1(x)),lwd=2,col="red")

#I understand that the priors on the two betas is Gaussian (I assumed centered at zero). I see that we fix the hyperparameters only of one of the two. 

