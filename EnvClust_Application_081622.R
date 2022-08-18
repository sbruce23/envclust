##################################################################################
##################################################################################
## R Code for "Adaptive Clustering and Feature Selection for Categorical Time 
## Series Using Interpretable Frequency-Domain Features" 
## Author: Scott A. Bruce
## Last Updated: August 16, 2022
## R version: 4.2.1
## Tested on Windows machine (Intel i9 10 core 3.3 Ghz processor, 128 GB RAM)
## Description: Code to reproduce application analysis and figures from manuscript.  Run the
## script from start to finish to see all figures.  WARNING: Some of the plots require
## multiple runs of the algorithms and may take some time.  
##
## Data: Data are available for download at https://physionet.org/content/capslpdb/1.0.0/
## Convert .txt files for rbd, ins, plm, n, and nfle to .csv and then you are
## ready to run the script below.
##
## DON'T FORGET TO SET YOUR WORKING DIRECTORY!
##################################################################################
##################################################################################

#######################
## If you want to see the figures immediately, load the .RData file and run the
## code here:
# load('EnvClust_ApplicationResults_081622.RData')
# fig1
# grid::grid.draw(fig7)
# grid::grid.draw(fig8)
# grid::grid.draw(fig9)
########################


rm(list=ls())

#set working directory to the folder containing all files
#setwd("path/to/file")

library(astsa)
library(parallel)
library(foreach)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
library(gridExtra)
library(cowplot)
library(fields)
library(viridis)
library(RColorBrewer)
library(gridBase)
library(abind)
library(fossil)
library(foreach)
library(doRNG)
library(future.apply)
options(future.globals.maxSize= 891289600)

#source R functions needed for running this script
source("EnvClust_Rfunctions_081622.R")

##ggplot settings
hw <- theme_gray()+ theme(
  plot.title=element_text(hjust=0.5),
  plot.subtitle=element_text(hjust=0.5),
  plot.caption=element_text(hjust=-.5),
  
  #  strip.text.y = element_blank(),
  strip.background=element_rect(fill=rgb(.9,.95,1),
                                colour=gray(.5), size=.2),
  
  panel.border=element_rect(fill=FALSE,colour=gray(.70)),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.spacing.x = unit(0.10,"cm"),
  panel.spacing.y = unit(0.05,"cm"),
  
  # axis.ticks.y= element_blank()
  axis.ticks=element_blank(),
  axis.text=element_text(colour="black"),
  axis.text.y=element_text(margin=margin(0,3,0,3)),
  axis.text.x=element_text(margin=margin(-1,0,3,0))
)

######################################################
## Read in data
######################################################
######################################################
rbd <- list()
for(k in 1:22){
  filename = paste("rbd",k,".csv",sep="")
  rbd[[k]] = read.csv(filename, header=TRUE)
}
rbd=lapply(rbd,function(x) x$Sleep.Stage[x$Duration.s.==30])
rbd <- lapply(rbd,function(x) factor(ifelse(x %in% c('W','MT'),'W/MT',x),
                                     levels=c('W/MT','R','S1','S2','S3','S4')));
naidx <- which(unlist(lapply(rbd,function(x) !anyNA(x))))
miscatidx <- which(unlist(lapply(rbd,function(x) sum(c('W/MT','S1','S2','S3','S4','R') %in% unique(x))==6)))
tmp <- table(c(naidx,miscatidx))
idx <- as.numeric(names(tmp)[tmp==2])
rbd_group = (1:22)[idx]; #remove those with NAs or missing a category
######################################################
ins <- list()
for(k in 1:9){
  filename = paste("ins",k,".csv",sep="")
  ins[[k]] = read.csv(filename, header=TRUE)
}
ins=lapply(ins,function(x) x$Sleep.Stage[x$Duration.s.==30])
ins <- lapply(ins,function(x) factor(ifelse(x %in% c('W','MT'),'W/MT',x),
                                     levels=c('W/MT','R','S1','S2','S3','S4')));
naidx <- which(unlist(lapply(ins,function(x) !anyNA(x))))
miscatidx <- which(unlist(lapply(ins,function(x) sum(c('W/MT','S1','S2','S3','S4','R') %in% unique(x))==6)))
tmp <- table(c(naidx,miscatidx))
idx <- as.numeric(names(tmp)[tmp==2])
ins_group = (1:9)[idx]; #remove those with NAs or missing a category
######################################################
plm <- list()
for(k in 1:10){
  filename = paste("plm",k,".csv",sep="")
  plm[[k]] = read.csv(filename, header=TRUE)
}
plm=lapply(plm,function(x) x$Sleep.Stage[x$Duration.s.==30])
plm <- lapply(plm,function(x) factor(ifelse(x %in% c('W','MT'),'W/MT',x),
                                     levels=c('W/MT','R','S1','S2','S3','S4')));
naidx <- which(unlist(lapply(plm,function(x) !anyNA(x))))
miscatidx <- which(unlist(lapply(plm,function(x) sum(c('W/MT','S1','S2','S3','S4','R') %in% unique(x))==6)))
tmp <- table(c(naidx,miscatidx))
idx <- as.numeric(names(tmp)[tmp==2])
plm_group = (1:10)[idx]; #remove those with NAs or missing a category
######################################################
control <- list()
for(k in 1:16){
  filename = paste("n",k,".csv",sep="")
  control[[k]] = read.csv(filename, header=TRUE)
}
control=lapply(control,function(x) x$Sleep.Stage[x$Duration.s.==30])
control[[12]][which(control[[12]]=='1S3')] <- 'S3'
control <- lapply(control,function(x) factor(ifelse(x %in% c('W','MT'),'W/MT',x),
                                     levels=c('W/MT','R','S1','S2','S3','S4')));
naidx <- which(unlist(lapply(control,function(x) !anyNA(x))))
miscatidx <- which(unlist(lapply(control,function(x) sum(c('W/MT','S1','S2','S3','S4','R') %in% unique(x))==6)))
tmp <- table(c(naidx,miscatidx))
idx <- as.numeric(names(tmp)[tmp==2])
control_group = (1:16)[idx]; #remove those with NAs or missing a category
######################################################
nfle <- list()
for(k in 1:40){
  filename = paste("nfle",k,".csv",sep="")
  nfle[[k]] = read.csv(filename, header=TRUE)
}
nfle=lapply(nfle,function(x) x$Sleep.Stage[x$Duration.s.==30])
nfle <- lapply(nfle,function(x) factor(ifelse(x %in% c('W','MT'),'W/MT',x),
                                     levels=c('W/MT','R','S1','S2','S3','S4')));
naidx <- which(unlist(lapply(nfle,function(x) !anyNA(x))))
miscatidx <- which(unlist(lapply(nfle,function(x) sum(c('W/MT','S1','S2','S3','S4','R') %in% unique(x))==6)))
tmp <- table(c(naidx,miscatidx))
idx <- as.numeric(names(tmp)[tmp==2])
nfle_group = (1:40)[idx]; #remove those with NAs or missing a category
######################################################

######################################################
## Data processing and checking
######################################################
######################################################

#subset data for only participants that 
#have all categories represented and no NA
rbd1 <- rbd[rbd_group] #don't lose any
ins1 <- ins[ins_group] #lost 3 since don't observe all categories
plm1 <- plm[plm_group] #lost 1 since don't observe all categories
control1 <- control[control_group] #don't lose any
nfle1 <- nfle[nfle_group] #don't lose any

dat <- c(rbd1, ins1, plm1, control1, nfle1)

group = c(rep('rbd',length(rbd_group)),
          rep('ins',length(ins_group)),
          rep('plm',length(plm_group)),
          rep('control',length(control_group)),
          rep('nfle',length(nfle_group))
          );

yt = list();
for (i in 1:length(dat)){
  yt[[i]] = cat_convert(dat[[i]],'W/MT')
  colnames(yt[[i]]) = c('R','S1','S2','S3','S4')
}
#colnames are R, S1, S2, S3, S4
#test to ensure encoding worked
# test <- lapply(yt,function(x) apply(x,1,function(y) which(y==1)))
# test <- lapply(test,function(x) sapply(x,function(y) ifelse(length(y)==0,0,y)))
# for(i in 1:length(dat)){
#   # print(c(i,sum(diag(table(test[[i]],dat[[i]])))))
#   print(table(test[[i]],dat[[i]]))
# }

######################################################
######################################################
######################################################

#20% to 90% of total sleep time to avoid falling asleep and waking periods
yt=lapply(yt,function(x) x[(floor(nrow(x)*2/10)+1):(floor(nrow(x)*9/10)),])
dat=lapply(dat,function(x) x[(floor(length(x)*2/10)+1):(floor(length(x)*9/10))])

#make sure all categories still represented
idx1=unlist(lapply(yt,function(x) sum(colSums(x)==0)==0))
idx2=unlist(lapply(yt,function(x) sum(rowSums(x)==0)>0))
idx=idx1&idx2;
dat=dat[idx];
yt=yt[idx];
group=group[idx];
table(group) #80 remaining


######################################################
## Figure 1
######################################################
######################################################

#obtain ordering of elements of yt from shortest to longest
ord = order(sapply(yt,length));

#obtain reference frequencies from shortest series
tmp=env.get(yt[[ord[1]]],50);
rf=tmp$freq;

#estimate spectral envelope and scalings
L=floor(sqrt(nrow(yt[[ord[1]]])));std=TRUE;
tmp=envsca.get(yt,L,std);

#plot sample time series from each group
lev = c("W/MT", "R", "S1", "S2", "S3", "S4")
r1 = rbd[[8]]
r1 = as.character(r1)
r1_ind = match(r1, lev)
###############################################################################
i1 = ins[[3]]
i1 = as.character(i1)
i1_ind = match(i1, lev)
###############################################################################
p1 = plm[[2]]
p1 = as.character(p1)
p1_ind = match(p1, lev)
###############################################################################
n1 = nfle[[2]]
n1 = as.character(n1)
n1_ind = match(n1, lev)
###############################################################################
c1 = control[[5]]
c1 = as.character(c1)
c1_ind = match(c1, lev)
###############################################################################
df=data.frame(t=(1:length(r1))*30/(3600),cat=r1_ind)
df$slpcat=rep('RBD',nrow(df))
dflong=df;

r1ts = ggplot()+
  scale_y_continuous(breaks=1:6, minor_breaks = 1:6, labels=lev,limits=c(0.97,6.03),expand = c(0, 0))+
  scale_x_continuous(breaks=1:10,minor_breaks=1:10,expand = c(0, 0))+
  geom_line(data=df,aes(x=t,y=cat))+ggtitle("RBD")+hw+
  theme(axis.title.y = element_text(vjust=-0.5),axis.title.x=element_blank())+coord_fixed(1)+
  # xlab("Time (Hours)")+
  ylab("\n Sleep Stage \n")#+theme(plot.margin=unit(c(.2,.2,.2,-.4),"cm"))
df=data.frame(t=(1:length(p1))*30/(3600),cat=p1_ind)
df$slpcat=rep('PLM',nrow(df))
dflong=rbind(dflong,df);

p1ts = ggplot()+
  scale_y_continuous(breaks=1:6, minor_breaks = 1:6, labels=lev,limits=c(0.97,6.03),expand = c(0, 0))+
  scale_x_continuous(breaks=1:10,minor_breaks=1:10,expand = c(0, 0))+
  geom_line(data=df,aes(x=t,y=cat))+ggtitle("PLM")+hw+
  theme(axis.title.y = element_text(vjust=-0.5))+
  xlab("Time (Hours)")+ylab("\n Sleep Stage \n")+coord_fixed(1)#+theme(plot.margin=unit(c(.2,.2,.2,-.4),"cm"))
df=data.frame(t=(1:length(i1))*30/(3600),cat=i1_ind)
df$slpcat=rep('INS',nrow(df))
dflong=rbind(dflong,df);

i1ts = ggplot()+
  scale_y_continuous(breaks=1:6, minor_breaks = 1:6, labels=lev,limits=c(0.97,6.03),expand = c(0, 0))+
  scale_x_continuous(breaks=1:12,minor_breaks=1:12,expand = c(0, 0))+
  geom_line(data=df,aes(x=t,y=cat))+ggtitle("INS")+hw+coord_fixed(1)+
  theme(axis.title.y = element_blank(),axis.title.x=element_blank())
  #xlab("Time (Hours)")+#ylab("\n Sleep Stage \n")#+theme(plot.margin=unit(c(.2,.2,.2,-.4),"cm"))
df=data.frame(t=(1:length(n1))*30/(3600),cat=n1_ind)
df$slpcat=rep('NFLE',nrow(df))
dflong=rbind(dflong,df);

n1ts = ggplot()+
  scale_y_continuous(breaks=1:6, minor_breaks = 1:6, labels=lev,limits=c(0.97,6.03),expand = c(0, 0))+
  scale_x_continuous(breaks=1:10,minor_breaks=1:10,expand = c(0, 0))+
  geom_line(data=df,aes(x=t,y=cat))+ggtitle("NFLE")+hw+coord_fixed(1)+
  theme(axis.title.y = element_blank())+
  xlab("Time (Hours)")#+ylab("\n Sleep Stage \n")#+theme(plot.margin=unit(c(.2,.2,.2,-.4),"cm"))



fig1 = ggplot()+
  scale_y_continuous(breaks=1:6, minor_breaks = 1:6, labels=lev,limits=c(0.9,6.1),expand = c(0, 0))+
  scale_x_continuous(breaks=1:10,minor_breaks=1:10,expand = c(0, 0))+
  geom_line(data=dflong,aes(x=t,y=cat))+ggtitle("Control")+hw+
  theme(axis.title.y = element_text(vjust=-0.5))+facet_wrap(~slpcat)+
  xlab("Time (Hours)")+ylab("\n Sleep Stage \n")+theme(plot.title=element_blank())+theme(plot.margin=unit(c(.01,.01,.01,.01),'cm'))
fig1

######################################################
## Figure 7 and Table 1
######################################################
######################################################
#cut back to only frequencies below 0.1 since covers most power
mean(sapply(tmp$env_ind,function(x) sum(x[tmp$rf<=0.1])/sum(x))) #89.5% of power in low freqs

tmp$env_ind=lapply(tmp$env_ind,function(x) x[tmp$rf<=0.1])
tmp$sca_ind=lapply(tmp$sca_ind,function(x) x[tmp$rf<=0.1,])
# tmp$env_ind = lapply(tmp$env_ind,function(x) x/sqrt(sum(diff(tmp$rf)[1]*x^2))); #standardize by l2 norm
# tmp$sca_ind = lapply(tmp$sca_ind,function(x) x*sqrt(0.5/0.1)); #standardize by l2 norm
tmp$rf=tmp$rf[tmp$rf<=0.1]

#use k means to select k
#scree plot
k=1:10
seed=12;
plan(multisession);s=500;
test=future_sapply(k,function(x) envsca_kmeans(tmp$rf,tmp$env_ind, tmp$sca_ind, x, s,spar=FALSE, seed)$wcss,future.seed=TRUE);
df=data.frame(objfn=test,k=k)

screep=ggplot(data=df,aes(x=k,y=objfn))+geom_line()+geom_point()+hw+
  scale_x_continuous(breaks=1:max(df$k),minor_breaks=1:max(df$k),expand = c(0, 0),limits=c(0.8,max(k)+.2))+
  labs(title="Scree Plot",x="Number of Clusters (k)",y="Within-Cluster SS")

##use sparse clustering for feature selection

##############################################################
## WARNING: Gap statistic calculations may take a long time (>30 minutes)
## depending on your computing resources 
##############################################################
k=4;s=c(1.1,seq(from=2,to=60,by=2));seed=12;L=floor(sqrt(nrow(yt[[ord[1]]])));std=TRUE;
gpstat=gapstat_s(yt, k, L, s, B=10, seed, std)
save.image('appresults0816.RData') 

gsplot=ggplot(data=data.frame(s=s,gapst=gpstat$gapst),
             aes(x=s,y=gapst))+geom_point()+geom_line()+hw+
  scale_x_continuous(breaks=seq(0,max(s),by=10),
                     minor_breaks=seq(0,max(s),by=10),limits=c(0,max(s)+1),expand = c(0, 0))+
  labs(title="Sparsity Gap Statistic",x="Sparsity Parameter (r)",y="Gap Statistic")

gpstat$s.star #10
s[which(gpstat$gapst==max(gpstat$gapst))[1]] #26
s[which(gpstat$gapst>max(gpstat$gapst)-sd(gpstat$gapst))[1]] #10

##sparse clustering
k=4;
s=10;
spar=TRUE;seed=52;
kmeansfit=envsca_kmeans_wrap(tmp$rf,tmp$env_ind, tmp$sca_ind, k, s, spar, nstart=100,seed)

round(sum(kmeansfit$w.env)/(sum(kmeansfit$w.env)+sum(kmeansfit$w.sca))*100,1) #2.3% envelope
round(sum(kmeansfit$w.sca)/(sum(kmeansfit$w.env)+sum(kmeansfit$w.sca))*100,1) #97.7% scalings

#table 1
apply(table(group,kmeansfit$group_est),2,function(x) x)
table(group)
table(kmeansfit$group_est)
round(apply(table(group,kmeansfit$group_est),2,function(x) x/sum(x)),2)
round(table(group)/80,2)
round(table(kmeansfit$group_est)/sum(kmeansfit$group_est),2)



round(colSums(kmeansfit$w.sca)/sum(kmeansfit$w.sca),2) #70% S2, 21% S3

envwgt=ggplot(data=data.frame(f=kmeansfit$rf,w=kmeansfit$w.env),aes(x=f,y=w))+geom_line()+
  scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0),limits=c(0,0.015))+
  hw+ylab("Weight")+xlab("Frequency")+ggtitle("Spectral Envelope Weights")

dfscas=melt(kmeansfit$w.sca);
names(dfscas)=c("Frequency","Category","Scaling")
dfscas$Frequency=tmp$rf[dfscas$Frequency]
dfscas$Category=c("R", "S1", "S2", "S3", "S4")[dfscas$Category]

scawgt=ggplot(data=dfscas,aes(x=Frequency,y=Scaling,color=Category))+geom_line()+
  scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+
  hw+xlab("Frequency")+ylab("Weight")+ggtitle("Optimal Scaling Weights")+labs(color='Sleep\nStage')


fig7=grid.arrange(grobs=list(screep,gsplot,envwgt,scawgt),
             layout_matrix = rbind(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,4)))
grid::grid.draw(fig7)

######################################################
## Figure 8
######################################################
######################################################

dfenvs=melt(kmeansfit$env_group)
dfenvs$f=rep(tmp$rf,k)
names(dfenvs)=c('Envelope','Cluster','Frequency')
dfenvs$Cluster = factor(dfenvs$Cluster)
envs.plot=ggplot(data=dfenvs[dfenvs$Frequency<=0.1,],aes(x=Frequency,y=Envelope,col=Cluster))+
  scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0),breaks=c(0,2,4,6,8),limits=c(0,8.5))+
  geom_line()+hw+labs(title='Standardized Spectral Envelope',y='')+theme(legend.position= 'bottom');




scadat1=as.data.frame(cbind(melt(kmeansfit$sca_group[[1]]),rep(tmp$rf,5),rep(1:5,each=length(tmp$rf)),rep(1,length(tmp$rf)*5)))[,c(-1,-2)];
names(scadat1)=c('Scaling','f','cat','Cluster');
scadat2=as.data.frame(cbind(melt(kmeansfit$sca_group[[2]]),rep(tmp$rf,5),rep(1:5,each=length(tmp$rf)),rep(2,length(tmp$rf)*5)))[,c(-1,-2)];
names(scadat2)=c('Scaling','f','cat','Cluster');
scadat3=as.data.frame(cbind(melt(kmeansfit$sca_group[[3]]),rep(tmp$rf,5),rep(1:5,each=length(tmp$rf)),rep(3,length(tmp$rf)*5)))[,c(-1,-2)];
names(scadat3)=c('Scaling','f','cat','Cluster');
scadat4=as.data.frame(cbind(melt(kmeansfit$sca_group[[4]]),rep(tmp$rf,5),rep(1:5,each=length(tmp$rf)),rep(4,length(tmp$rf)*5)))[,c(-1,-2)];
names(scadat4)=c('Scaling','f','cat','Cluster');
scadat=rbind(scadat1,scadat2,scadat3,scadat4)
scadat$cat=factor(scadat$cat,levels=1:5,labels=c("R", "S1", "S2", "S3", "S4"))

min(sapply(kmeansfit$sca_group,min)) #-1.22
max(sapply(kmeansfit$sca_group,max)) #0.99
zr=c(-1.24,1.24)

scaplot=ggplot(data=scadat,aes(x=f,y=cat,fill=Scaling))+geom_tile()+hw+
  scale_y_discrete(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  scale_fill_distiller(palette='RdBu',limits=zr)+facet_wrap(~Cluster,labeller=label_both)+
  labs(x="Frequency",y='Sleep Stage',fill='',title='Standardized Optimal Scalings')+theme(legend.position='bottom');


fig8=grid.arrange(grobs=list(envs.plot,scaplot),
             layout_matrix = rbind(c(1,1,1,1),
                                   c(1,1,1,1),
                                   c(1,1,1,1), 
                                   c(1,1,1,1),
                                   c(2,2,2,2),
                                   c(2,2,2,2),
                                   c(2,2,2,2),
                                   c(2,2,2,2),
                                   c(2,2,2,2)))
grid::grid.draw(fig8)

######################################################
## Figure 9 and other metrics
######################################################
######################################################

#proportion of night in various sleep stages overall and by cluster
c(rowMeans(sapply(yt,function(x) colMeans(x))),
  1-sum(rowMeans(sapply(yt,function(x) colMeans(x)))))

c(rowMeans(sapply(yt[kmeansfit$group_est==1],function(x) colMeans(x))),
  1-sum(rowMeans(sapply(yt[kmeansfit$group_est==1],function(x) colMeans(x)))))
c(rowMeans(sapply(yt[kmeansfit$group_est==2],function(x) colMeans(x))),
  1-sum(rowMeans(sapply(yt[kmeansfit$group_est==2],function(x) colMeans(x)))))
c(rowMeans(sapply(yt[kmeansfit$group_est==3],function(x) colMeans(x))),
  1-sum(rowMeans(sapply(yt[kmeansfit$group_est==3],function(x) colMeans(x)))))
c(rowMeans(sapply(yt[kmeansfit$group_est==4],function(x) colMeans(x))),
  1-sum(rowMeans(sapply(yt[kmeansfit$group_est==4],function(x) colMeans(x)))))



###############fuzzy clustering#################
##############################################################
## WARNING: Gap statistic calculations may take a long time (>30 minutes)
## depending on your computing resources 
##############################################################
k=4;
m=seq(from=1.02,10,length.out=25);B=10;seed=6;
gsfuzz2=gapstat_m(tmp$rf, tmp$env_ind, tmp$sca_ind, k, m, B, seed)

gappfuzzk4=ggplot(data=data.frame(m=m,gapst=gsfuzz2$gapst),aes(x=m,y=gapst))+geom_line()+geom_point()+hw+
  scale_x_continuous(breaks=1:max(m),minor_breaks=1:max(m),expand = c(0, 0))+
  labs(title="Fuzziness Parameter Gap Statistic",x="Fuzziness Parameter (m)",y="Gap Statistic")
gappfuzzk4

m[which(gsfuzz2$gapst==max(gsfuzz2$gapst))] #2.52
m[which(gsfuzz2$gapst>max(gsfuzz2$gapst)-sd(gsfuzz2$gapst))[1]] #1.39


k=4;m=1.39;
tmp1=envsca_fuzzyclust_wrap(tmp$rf, tmp$env_ind, tmp$sca_ind, k, m, nstart=100, seed=23)

#explore results
dfbp=melt(tmp1$umat)
dfbp=cbind(dfbp,rep(group,k))
names(dfbp)=c('i','cluster','u','group')
dfbp$cluster= factor(dfbp$cluster)
dfbp$group= factor(dfbp$group,levels=c('control','ins','nfle','plm','rbd'),
                   labels=c('Control','INS','NFLE','PLM','RBD'))

uplot=ggplot(data=dfbp, aes(y=u , x=cluster, fill=group))+ 
  geom_violin()+geom_jitter()+ylim(c(0,1))+facet_wrap(~group,ncol=2)+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))+
  xlab("Cluster") + ylab("Membership")+hw+theme(legend.position = "none")+
  ggtitle("Cluster Membership Degree")
uplot

#arrange into figure
fig9=grid.arrange(gappfuzzk4,uplot,layout_matrix=rbind(c(1,1),
                                                  c(2,2),
                                                  c(2,2)))
grid::grid.draw(fig9)




dfenv=melt(tmp1$env)
names(dfenv)=c('Cluster','Frequency','Envelope')
dfenv$Frequency=tmp$rf[dfenv$Frequency]
dfenv$Cluster = factor(dfenv$Cluster)
env.plot=ggplot(data=dfenv,aes(x=Frequency,y=Envelope,col=Cluster))+
  scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0),breaks=c(0,2,4,6,8),limits=c(0,8.5))+
  geom_line()+hw+labs(title='Standardized Spectral Envelope',y='')+theme(legend.position= 'bottom');
env.plot

scadat1=as.data.frame(cbind(melt(tmp1$sca[1,,]),rep(tmp$rf,5),rep(1:5,each=length(tmp$rf)),rep(1,length(tmp$rf)*5)))[,c(-1,-2)];
names(scadat1)=c('Scaling','f','cat','Cluster');
scadat2=as.data.frame(cbind(melt(tmp1$sca[2,,]),rep(tmp$rf,5),rep(1:5,each=length(tmp$rf)),rep(2,length(tmp$rf)*5)))[,c(-1,-2)];
names(scadat2)=c('Scaling','f','cat','Cluster');
scadat3=as.data.frame(cbind(melt(tmp1$sca[3,,]),rep(tmp$rf,5),rep(1:5,each=length(tmp$rf)),rep(3,length(tmp$rf)*5)))[,c(-1,-2)];
names(scadat3)=c('Scaling','f','cat','Cluster');
scadat4=as.data.frame(cbind(melt(tmp1$sca[4,,]),rep(tmp$rf,5),rep(1:5,each=length(tmp$rf)),rep(4,length(tmp$rf)*5)))[,c(-1,-2)];
names(scadat4)=c('Scaling','f','cat','Cluster');
scadat=rbind(scadat1,scadat2,scadat3,scadat4)
scadat$cat=factor(scadat$cat,levels=1:5,labels=c("R", "S1", "S2", "S3", "S4"))

min(scadat$Scaling)
max(scadat$Scaling)
zr=c(-1.19,1.19)

scaplot=ggplot(data=scadat,aes(x=f,y=cat,fill=Scaling))+geom_tile()+hw+
  scale_y_discrete(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  scale_fill_distiller(palette='RdBu',limits=zr)+facet_wrap(~Cluster,labeller=label_both)+
  labs(x="Frequency",y='Sleep Stage',fill='',title='Standardized Optimal Scalings')+theme(legend.position='bottom');


grid.arrange(grobs=list(env.plot,scaplot),
             layout_matrix = rbind(c(1,1,1,1),
                                   c(1,1,1,1),
                                   c(1,1,1,1), 
                                   c(1,1,1,1),
                                   c(2,2,2,2),
                                   c(2,2,2,2),
                                   c(2,2,2,2),
                                   c(2,2,2,2),
                                   c(2,2,2,2)))


