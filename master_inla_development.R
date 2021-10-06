setwd("/Users/lorenzocappello/Google Drive/Statistics/postdoc Palacios/delta_covid_attempt/selection_model/")

source("functions_inla_selection.R") #where the new functions should be stored 
source("tajima_inference.R")
library("phylodyn")
library("ape")
library("lubridate")


##################################
###Extract the s and n vector ####
##################################


# dev.off()
# treeDelta<-read.nexus("../beast2_ver_3/tree_uk_2_delta")
# 
# plot(treeDelta,show.tip.label=FALSE)
# axisPhylo()
# 
# 
# treeNotdelta<-read.nexus("../beast2_ver_3/tree_uk_2_notdelta")
# 
# #treeNotdelta<-read.nexus("../beast2/tree_italy_notdelta")
# 
# #tree not delta refinement
# #treeD <-drop.tip(tree,treeG$tip.label) # This is the D clade
# treeNotdelta <- extract.clade(treeNotdelta,153)
# 
# plot(treeNotdelta,show.tip.label = FALSE)
# axisPhylo()
# #plot(treeNotdelta)
# 
# 
# data1<-treeDelta
# data2 <- treeNotdelta
# 
# samp_times1 <- c()
# for (label in data1$tip.label){samp_times1 <- c(samp_times1,strsplit(label,"[|]")[[1]][3])}
# samp_times1 <- decimal_date(date(samp_times1))
# samp_times2 <- c()
# for (label in data2$tip.label){samp_times2 <- c(samp_times2,strsplit(label,"[|]")[[1]][3])}
# samp_times2 <- decimal_date(date(samp_times2))
# lastdate <- max(c(samp_times1,samp_times2))
# samp_times1 <- lastdate-samp_times1
# samp_times2 <- lastdate-samp_times2


n=100

coeff=n/0.4
sampling_traj<-function(t){
  result=rep(0,length(t))
  result<-1
  return(result*coeff)}

samp_times1<-c(0,sampsim_thin(max.n=n-2,traj=sampling_traj,xlim=c(0,3))[-n+1])
n_sampled1 <- c(2,rep(1,n-2)) # simulation like this it will be like this.
samp_times2<-c(0,sampsim_thin(max.n=n-2,traj=sampling_traj,xlim=c(0,3))[-n+1])
n_sampled2 <- c(2,rep(1,n-2)) # simulation like this it will be like this.



covid_exp1<-function(t){
  result=rep(0,length(t))
  result[t<=0.35]<-10
  result[t>0.35 & t<0.4]<-9.999987e+14*exp(-92.1034*t[t>0.35 & t<0.4]) #50 times smaller
  result[t>=0.4]<-0.1
  return(result)
}


covid_exp2<-function(t){
  t <- t+0.32
  result=rep(0,length(t))
  result[t<=0.35]<-10
  result[t>0.35 & t<0.4]<-9.999987e+14*exp(-92.1034*t[t>0.35 & t<0.4]) #50 times smaller
  result[t>=0.4]<-0.1
  return(result)
}


simulation1<-coalsim(samp_times=samp_times1-min(samp_times1),n=n_sampled1,traj=covid_exp1)
tree1<-sample_genealogy(simulation1)
tree1<-read.tree(text=tree1$Newick)
#plot(tree1,show.tip.label = FALSE)

simulation2<-coalsim(samp_times2-min(samp_times2),n_sampled2, traj = covid_exp2)
tree2<-sample_genealogy(simulation2)
tree2<-read.tree(text=tree2$Newick)


samp_times2 <- samp_times2



dev.off()
res1s<-BNPR_sel(tree1,tree2,samp_times1,samp_times2,lengthout=100,parSel = F,preferential =F)
plot_BNPR(res1s,heatmaps=FALSE,heatmap_labels = FALSE,ylab="Ne(t)",main="Effective Pop Size",ylim=c(0.005,5000),xlim=c(0.5,0))
plot_BNPR_plus(res1s,heatmaps=FALSE,heatmap_labels = FALSE,ylim=c(0.005,5000),xlim=c(0.5,0),parameter = "selInt")
x <- seq(0,0.5,0.01)
lines(x,covid_exp1(x)/covid_exp2(x),col="red",lwd="2")


res1<-BNPR(tree2,lengthout=100)
plot_BNPR(res1,heatmaps=FALSE,heatmap_labels = FALSE,ylab="Ne(t)",main="Effective Pop Size",ylim=c(0.005,5000),xlim=c(0.5,0))
lines(x,covid_exp2(x),col="red")


plot_adaSel(res1s,heatmaps=FALSE,heatmap_labels = FALSE,ylab="",main="",ylim=c(0.05,50000),xlim=c(0.7,0))
lines(x,covid_exp1(x)/covid_exp2(x),col="red")

#res1s <- BNPR_sel(treeGG,treeDD,samp_times1,samp_times2,lengthout = 100,beta0_remove = FALSE)
#plot_BNPR(res1s,heatmaps=FALSE,heatmap_labels = FALSE,ylab="Ne(t)",main="Effective Pop Size",ylim=c(0.005,50),xlim=c(0.5,0))

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

