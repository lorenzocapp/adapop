setwd("/Users/lorenzocappello/Google Drive/Statistics/postdoc Palacios/stanford_covid/selection_model_trial")
#setwd("/Users/cappello/Google Drive/postdoc Palacios/stanford_covid/selection_model_trial")

source("functions_inla_selection.R") #where the new functions should be stored 
source("tajima_inference.R")
library("phylodyn")
library("ape")
library("lubridate")
library(ggplot2)


###########################
###Load the two trees ####
##########################


tree<-read.nexus("../beast2/tree_beast2")
plot(tree,show.tip.label = FALSE)

# ====== x-axis limit and ticks =====
samp.date <- ymd(sapply(strsplit(tree$tip.label, split='|', fixed=TRUE), '[', 3))
lastsd <- max(samp.date) # last sampling date
tot.tr.h <- max(node.depth.edgelength(tree))
plt.xmax <- lastsd
plt.xmin <- as.Date(date_decimal(decimal_date(plt.xmax) - tot.tr.h - 0.06)) # give 3 weeks buffer

plt.breaks.major <- sort(seq(plt.xmax, plt.xmin, by='-4 weeks'))
plt.breaks.minor <- sort(seq(plt.xmax, plt.xmin, by='-2 weeks'))

# ===== plot parameters =====
edge.thk <- 0.8 # tree edge thickness
xlab.size <- 18
xlab.mar <- 15
xtck.lab.size <- 15
xtck.lab.mar <- 5
xtck.lab.vjust <- 0.5


# ===== extract the two trees  =====
treeG <-extract.clade(tree,307) #This is is the G clade
plot(treeG,show.tip.label=FALSE)
treeD <-drop.tip(tree,treeG$tip.label) # This is the D clade
plot(treeD,show.tip.label = FALSE)



samp_times1 <- c()
for (label in treeG$tip.label){samp_times1 <- c(samp_times1,strsplit(label,"[|]")[[1]][3])}
presentG <- max(samp_times1)
samp_times2 <- c()
for (label in treeD$tip.label){samp_times2 <- c(samp_times2,strsplit(label,"[|]")[[1]][3])}
presentD <- max(samp_times2)


g1 <- ggtree(treeD, mrsd=presentD, as.Date=TRUE, size=edge.thk) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.background = element_rect(),
    panel.grid.major   = element_line(color=alpha('grey', alpha=0.3), size=1.5),
    panel.grid.minor   = element_line(color=alpha('grey', alpha=0.2), size=1.2),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank())+
  scale_y_continuous(breaks = NULL, expand=expansion(mult=c(0.03, 0.02))) +
  scale_x_date(breaks=plt.breaks.major, 
               minor_breaks=plt.breaks.minor,
               limits=c(plt.xmin, plt.xmax),
               labels=format(plt.breaks.major, format="%b %d"),
               expand=expansion(mult=c(0, 0.02))) 

g2 <- ggtree(treeG, mrsd=presentG, as.Date=TRUE, size=edge.thk) + 
  theme(axis.text.x=element_text(angle=90, size=xtck.lab.size,
                                 margin=margin(t=xtck.lab.mar),
                                 vjust=xtck.lab.vjust),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.background = element_rect(),
        axis.line.x = element_line(color=plt.xaxis.col, size=1),
        axis.ticks.x = element_line(color=plt.xaxis.col, size=1.5),
        axis.ticks.length=unit(.25, "cm"),
        axis.title.x=element_text(size=xlab.size, face="bold", margin=margin(t=xlab.mar)),
        panel.grid.major   = element_line(color=alpha('grey', alpha=0.3), size=1.5),
        panel.grid.minor   = element_line(color=alpha('grey', alpha=0.2), size=1.2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_y_continuous(breaks = NULL, expand=expansion(mult=c(0.03, 0.02))) +
  scale_x_date(breaks=plt.breaks.major, 
               minor_breaks=plt.breaks.minor,
               limits=c(plt.xmin, plt.xmax),
               labels=format(plt.breaks.major, format="%b %d"),
               expand=expansion(mult=c(0, 0.02))) +
  xlab("Time")


cowplot::plot_grid(plotlist = list(g1, g2), nrow = 2,rel_heights=c(1,6))

#I understand that the priors on the two betas is Gaussian (I assumed centered at zero). I see that we fix the hyperparameters only of one of the two. 
