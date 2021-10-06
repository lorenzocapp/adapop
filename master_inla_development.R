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
treeDelta<-read.nexus("../beast2_ver_3/tree_uk_2_delta")

plot(treeDelta,show.tip.label=FALSE)
axisPhylo()


treeNotdelta<-read.nexus("../beast2_ver_3/tree_uk_2_notdelta")

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



dev.off()
res1s<-BNPR_sel(data2,data1,samp_times2,samp_times1,lengthout=100,parSel = F,preferential = T)
plot_BNPR(res1s,heatmaps=FALSE,heatmap_labels = FALSE,ylab="Ne(t)",main="Effective Pop Size",ylim=c(0.005,5000),xlim=c(0.2,0))

res1<-BNPR(data1,lengthout=100)
plot_BNPR(res1,heatmaps=FALSE,heatmap_labels = FALSE,ylab="Ne(t)",main="Effective Pop Size",ylim=c(0.005,5000),xlim=c(0.2,0))

res1p<-BNPR_PS_ada(data1,lengthout=100)
plot_BNPR(res1p,heatmaps=FALSE,heatmap_labels = FALSE,ylab="Ne(t)",main="Effective Pop Size",ylim=c(0.005,5000),xlim=c(0.2,0))



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

