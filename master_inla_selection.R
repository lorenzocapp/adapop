
dir <- '/Users/lorenzocappello/Google Drive/Statistics/postdoc Palacios/delta_covid_attempt/selection_model_delta_truncated/'


setwd(dir)
#setwd("/Users/cappello/Google Drive/postdoc Palacios/stanford_covid/selection_model_trial")

source("functions_inla_selection_truncated.R") #where the new functions should be stored 
source("tajima_inference.R")
library("phylodyn")
library("ape")
library("lubridate")


###########################
###Load the two trees ####
##########################



dev.off()
treeDelta<-read.nexus("../beast2_ver_1&2/tree_singapore_delta")
treeDelta<-read.nexus("../beast2_ver_3/tree_uk_2_delta")

plot(treeDelta,show.tip.label=FALSE)
axisPhylo()

treeNotdelta<-read.nexus("../beast2_ver_1&2/tree_singapore_notdelta")
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





res1s<-BNPR_sel(data1,data2,samp_times1,samp_times2,lengthout=100)
plot_BNPR(res1s,heatmaps=FALSE,heatmap_labels = FALSE,ylab="Ne(t)",main="Effective Pop Size",ylim=c(0.005,50),xlim=c(0.54,0))
res1b<-BNPR(data2,lengthout=100)
plot_BNPR(res1b,heatmaps=TRUE,heatmap_labels = FALSE,ylab="Ne(t)",main="Effective Pop Size",ylim=c(0.5,500),xlim=c(0.25,0))

res2s<-BNPR_sel(data2,data1,samp_times2,samp_times1,lengthout=100,parSel = FALSE)
plot_BNPR(res2s,heatmaps=FALSE,heatmap_labels = FALSE,ylab="Ne(t)",main="Effective Pop Size",ylim=c(0.5,500),xlim=c(0.25,0))



par(mfrow=c(2,2))
plot_BNPR(res1s,heatmaps=FALSE,heatmap_labels = FALSE,ylab="Ne(t)",main="",ylim=c(0.05,500),xlim=c(0.5,0))
plot_BNPR(res2s,heatmaps=FALSE,heatmap_labels = FALSE,ylab="",main="",ylim=c(0.05,50000),xlim=c(0.5,0))
lines(res2s$selInt)
plot_adaSel(res2s,heatmaps=FALSE,heatmap_labels = FALSE,ylab="",main="",ylim=c(0.05,50000),xlim=c(0.5,0))


par(mfrow=c(1,2))
plot(res2s$beta0post,type="l",lwd=3,main="Posterior beta0")
abline(v=0,col="red",lwd=2)
plot(res2s$beta1post,type="l",lwd=3,,main="Posterior beta1")
abline(v=1,col="red",lwd=2)


res3s<-BNPR_sel(data1,data2,samp_times1,samp_times2,lengthout=50,u.truncation = 0.25,
                prec_alpha = 0.01)
res4s<-BNPR_sel(data2,data1,samp_times2,samp_times1,lengthout=50,u.truncation = 0.25,
                prec_alpha = 0.000001)

par(mfrow=c(1,2))
plot_BNPR(res3s,heatmaps=FALSE,heatmap_labels = FALSE,ylab="Ne(t)",main="",ylim= c(0.5,50),xlim=c(0.25,0))
plot_BNPR(res4s,heatmaps=FALSE,heatmap_labels = FALSE,ylab="",main="",ylim=c(0.5,500),xlim=c(0.25,0))


par(mfrow=c(1,2))
plot(res4s$beta0post,type="l",lwd=3,main="Posterior beta0")
abline(v=0,col="red",lwd=2)
plot(res4s$beta1post,type="l",lwd=3,,main="Posterior beta1")
abline(v=1,col="red",lwd=2)


#I understand that the priors on the two betas is Gaussian (I assumed centered at zero). I see that we fix the hyperparameters only of one of the two. 

