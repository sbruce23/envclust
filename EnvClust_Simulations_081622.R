##################################################################################
##################################################################################
## R Code for "Adaptive Clustering and Feature Selection for Categorical Time 
## Series Using Interpretable Frequency-Domain Features" 
## Author: Scott A. Bruce
## Last Updated: August 16, 2022
## R version: 4.2.1
## Tested on Windows machine (Intel i9 10 core 3.3 Ghz processor, 128 GB RAM)
## Description: Code to reproduce simulation figures from manuscript.  Run the
## script from start to finish to see all figures.  WARNING: Some of the plots require
## multiple runs of the algorithms and may take some time.  
##
## DON'T FORGET TO SET YOUR WORKING DIRECTORY!
##################################################################################
##################################################################################

#set working directory to the folder containing all files
#setwd("path/to/file")

#######################
## If you want to see the figures immediately, load the .RData file and run the code here:
#######################
# load('EnvClust_SimulationResults_081622.RData')
# grid::grid.draw(fig2)
# grid::grid.draw(fig3)
# fig4
# grid::grid.draw(fig5)
# grid::grid.draw(fig6)
#######################

rm(list=ls())

library(astsa)
library(parallel)
library(foreach)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(reshape2)
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
library(fclust)
options(future.globals.maxSize= 891289600)

#source R functions needed for running this script
source("EnvClust_Rfunctions_081622.R")

#modify ggplot theme settings
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

####################################
## Figure 2
####################################

## Example 1 (4 clusters)
n = 200; #length of time series
nsub = 25 # number of time series per group

set.seed(94)
nobs = n
m = 4   # number of categories

beta1 = cbind(c(0, 3, 1, 1),
              c(0, 1, 3, 1),
              c(0, 1, 1, 3))

beta2 = cbind(c(0, -1, 1, 1),
              c(0, 1, -1, 1),
              c(0, 1, 1, -1))

beta3 = cbind(c(0, -1, 1, 1),
              c(0, 1, 1, 1),
              c(0, 1, 1, 2))

beta4 = cbind(c(0, 2, 1, 1),
              c(0, 1, 1, 1),
              c(0, 1, 1, -1))

# Group 1
yt1 = vector("list",nsub);
index1 = c()
index1[1] = 1
for(sub in 1:nsub){
  yt1[[sub]]=matrix(0,nrow=nobs,ncol=m-1);
  yt1[[sub]][1,1] = 1                            
  zt = matrix(0,nobs,m)
  for (i in 2:nobs){    
    zt[i-1,] = cbind(1, t(yt1[[sub]][i-1,]))
    deno = 1 + (exp(zt[i-1,]%*% beta1[,1]) + exp(zt[i-1,]%*%beta1[,2]) +
                  exp(zt[i-1,]%*% beta1[,3]))
    pi1 = exp(zt[i-1,]%*% beta1[,1])/deno
    pi2 = exp(zt[i-1,]%*% beta1[,2])/deno
    pi3 = exp(zt[i-1,]%*% beta1[,3])/deno
    pi4 = 1/deno
    
    index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
    if (index1[i]==4){
      yt1[[sub]][i,1:(m-1)] = rep(0,m-1)
    } else {
      yt1[[sub]][i,index1[i]]=1
    }
  }
}

# Group 2
yt2 = vector("list",nsub);
index1 = c()
index1[1] = 1
for(sub in 1:nsub){
  yt2[[sub]]=matrix(0,nrow=nobs,ncol=m-1);
  yt2[[sub]][1,1] = 1                            
  zt = matrix(0,nobs,m)
  for (i in 2:nobs){    
    zt[i-1,] = cbind(1, t(yt2[[sub]][i-1,]))
    deno = 1 + (exp(zt[i-1,]%*% beta2[,1]) + exp(zt[i-1,]%*%beta2[,2]) +
                  exp(zt[i-1,]%*% beta2[,3]))
    pi1 = exp(zt[i-1,]%*% beta2[,1])/deno
    pi2 = exp(zt[i-1,]%*% beta2[,2])/deno
    pi3 = exp(zt[i-1,]%*% beta2[,3])/deno
    pi4 = 1/deno
    
    index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
    if (index1[i]==4){
      yt2[[sub]][i,1:(m-1)] = rep(0,m-1)
    } else {
      yt2[[sub]][i,index1[i]]=1
    }
  }
}

# Group 3
yt3 = vector("list",nsub);
index1 = c()
index1[1] = 1
for(sub in 1:nsub){
  yt3[[sub]]=matrix(0,nrow=nobs,ncol=m-1);
  yt3[[sub]][1,1] = 1                            
  zt = matrix(0,nobs,m)
  for (i in 2:nobs){    
    zt[i-1,] = cbind(1, t(yt3[[sub]][i-1,]))
    deno = 1 + (exp(zt[i-1,]%*% beta3[,1]) + exp(zt[i-1,]%*%beta3[,2]) +
                  exp(zt[i-1,]%*% beta3[,3]))
    pi1 = exp(zt[i-1,]%*% beta3[,1])/deno
    pi2 = exp(zt[i-1,]%*% beta3[,2])/deno
    pi3 = exp(zt[i-1,]%*% beta3[,3])/deno
    pi4 = 1/deno
    
    index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
    if (index1[i]==4){
      yt3[[sub]][i,1:(m-1)] = rep(0,m-1)
    } else {
      yt3[[sub]][i,index1[i]]=1
    }
  }
}

# Group 4
yt4 = vector("list",nsub);
index1 = c()
index1[1] = 1
for(sub in 1:nsub){
  yt4[[sub]]=matrix(0,nrow=nobs,ncol=m-1);
  yt4[[sub]][1,1] = 1                            
  zt = matrix(0,nobs,m)
  for (i in 2:nobs){    
    zt[i-1,] = cbind(1, t(yt4[[sub]][i-1,]))
    deno = 1 + (exp(zt[i-1,]%*% beta4[,1]) + exp(zt[i-1,]%*%beta4[,2]) +
                  exp(zt[i-1,]%*% beta4[,3]))
    pi1 = exp(zt[i-1,]%*% beta4[,1])/deno
    pi2 = exp(zt[i-1,]%*% beta4[,2])/deno
    pi3 = exp(zt[i-1,]%*% beta4[,3])/deno
    pi4 = 1/deno
    
    index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
    if (index1[i]==4){
      yt4[[sub]][i,1:(m-1)] = rep(0,m-1)
    } else {
      yt4[[sub]][i,index1[i]]=1
    }
  }
}

yt = c(yt1,yt2,yt3,yt4);
group_true = c(rep(1,nsub), rep(2,nsub),rep(3,nsub),rep(4,nsub));

#plot time series, envelopes, and scalings
cts = matrix(sapply(yt,function(x) apply(x,1,function(y) ifelse(length(which(y==1))==0,4,which(y==1)))),ncol=1);
cts=cbind(cts,rep(1:n,nsub*4),rep(1:4,each=n*nsub),rep(1:(nsub*4),each=n))
colnames(cts)=c('Category','Time','Cluster','Series')
cts=as.data.frame(cts)
cts$Cluster = factor(cts$Cluster,levels=1:4);
cts$Series = factor(cts$Series,levels=1:(nsub*4));

tsplot=ggplot(data=cts[cts$Series %in% c(1,nsub+1,2*nsub+1,3*nsub+1),],aes(x=Time,y=Category,color=Cluster))+geom_line()+
  scale_y_continuous(expand=c(0,0),limits=c(0.5,4.5),breaks=1:4)+
  scale_x_continuous(expand = c(0, 0),limits=c(1,n),breaks=c(1,seq(0,n-25,by=25)[-1]))+hw+
  xlab("Time")+ylab("Category")+labs(color = "Cluster")+facet_wrap(~Cluster,ncol=2)+theme(legend.position="top")

#estimate spectral envelope and optimal scalings
L=floor(sqrt(n));std=TRUE;
tmp=envsca.get(yt,L,std);

envmat=as.data.frame(unlist(tmp$env_ind));
envmat=cbind(envmat,rep(tmp$rf,nsub*4),rep(1:4,each=nsub*length(tmp$rf)),rep(1:(nsub*4),each=length(tmp$rf)))
names(envmat)=c('env','f','k','i');
envmat$k = factor(envmat$k,levels=1:4)

scamat=as.data.frame(unlist(tmp$sca_ind));
scamat=cbind(scamat,rep(rep(1:3,each=length(tmp$rf)),nsub*4),
             rep(tmp$rf,3*nsub*4),rep(1:4,each=3*nsub*length(tmp$rf)),
             rep(1:(nsub*4),each=3*length(tmp$rf)))
names(scamat)=c('env','s','f','k','i');
scamat$k = factor(scamat$k,levels=1:4)
scamat$s = factor(scamat$s,levels=1:3)

#cluster mean envelope and scalings
envmeans=aggregate(envmat$env,by=envmat[,c('f','k')],FUN=mean)
envplot=ggplot(data=envmeans,aes(x=f,y=x,color=k))+
  geom_line()+hw+scale_y_continuous(expand=c(0,0),limits=c(0,4))+
  scale_x_continuous(expand = c(0, 0),limits=c(0,0.5),breaks=c(0,.1,.2,.3,.4))+
  xlab("Frequency (Hz)")+ylab("Standardized Spectral Envelope")+labs(color = "Cluster")+theme(legend.position="none")

scameans=aggregate(scamat$env,by=scamat[,c('s','f','k')],FUN=mean)
scaplot=ggplot(data=scameans,aes(x=f,y=x,color=k))+hw+scale_x_continuous(expand=c(0,0),limits=c(0,.5),breaks=c(0,.1,.2,.3,.4))+
  scale_y_continuous(expand=c(0,0),limits=c(-1.25,1.25))+geom_line()+facet_wrap(~s,nrow=1)+
  xlab("Frequency (Hz)")+ylab("Standardized Scalings")+labs(color = "Cluster")+theme(legend.position="none")


#arrange into figure
fig2=grid.arrange(tsplot,envplot,scaplot,
             layout_matrix=rbind(c(1,1),c(1,1),c(1,1),c(2,2),c(2,2),c(3,3),c(3,3)))
fig2

####################################
## Figure 3
####################################
####################################
## WARNING: Gap statistic calculations may take a long time (>30 minutes)
## depending on computing resources available
####################################

#scree plot
k=1:10
seed=12;

##regular k means
plan(multisession);
test=future_sapply(k,function(x) envsca_kmeans(tmp$rf,tmp$env_ind, tmp$sca_ind, x, s=25,spar=FALSE, seed)$wcss,future.seed=TRUE);
df=data.frame(objfn=test,k=k)

screep=ggplot(data=df,aes(x=k,y=objfn))+geom_line()+geom_point()+hw+
  scale_x_continuous(breaks=1:max(df$k),minor_breaks=1:max(df$k),expand = c(0, 0),limits=c(0.8,max(k)+.2))+
  labs(title="Scree Plot",x="Number of Clusters (k)",y="Within-Cluster SS")

##gap stat to select sparsity param
k=4;s=c(1.1,seq(from=2,to=24,by=1));seed=49;std=TRUE;
set.seed(seed)
gsdat=gapstat_s(yt, k, L, s, B=10, seed, std)

gsdat$s.star
s[which(gsdat$gapst==max(gsdat$gapst))[1]] 
s[which(gsdat$gapst>(max(gsdat$gapst)-sd(gsdat$gapst)))[1]] 

gapstatplot=ggplot(data=data.frame(s=s,gapst=gsdat$gapst),aes(x=s,y=gapst))+geom_point()+geom_line()+hw+
  scale_x_continuous(breaks=seq(0,max(s),by=10),minor_breaks=seq(0,max(s),by=10),limits=c(0,max(s)+1),expand = c(0, 0))+
  labs(title="Sparsity Gap Statistic",x="Sparsity Parameter (r)",y="Gap Statistic")

##sparse clustering
k=4;s=7;
spar=TRUE;seed=8;
sparclus_results=envsca_kmeans(tmp$rf,tmp$env_ind, tmp$sca_ind, k, s, spar, seed)

sum(sparclus_results$w.env)/(sum(sparclus_results$w.env)+sum(sparclus_results$w.sca)) #89.5%
sum(sparclus_results$w.sca)/(sum(sparclus_results$w.env)+sum(sparclus_results$w.sca)) #10.5%


envwgt=ggplot(data=data.frame(f=sparclus_results[[2]],w=sparclus_results[[3]]),aes(x=f,y=w))+geom_line()+
  scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0),limits=c(-0.01,0.35))+
  hw+ylab("Weight")+xlab("Frequency")+ggtitle("Spectral Envelope Weights")

dfscas=melt(sparclus_results[[4]]);
names(dfscas)=c("Frequency","Category","Scaling")
dfscas$Frequency=tmp$rf[dfscas$Frequency]
dfscas$Category=c(1,2,3)[dfscas$Category]
dfscas$Category=factor(dfscas$Category,levels=c(1,2,3))

scawgt=ggplot(data=dfscas,aes(x=Frequency,y=Scaling,color=Category))+geom_line()+
  scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0),limits=c(-0.001,0.05))+
  hw+xlab("Frequency")+ylab("Weight")+ggtitle("Optimal Scaling Weights")+theme(legend.position = c(.8,.8))


##replicated simulations for sparse setting
############################################
R=100 #replications
n = 200; #length of time series

set.seed(94)

#start cluster
ncore <- detectCores();
c1 <- makeCluster(max(1,ncore-2));
registerDoParallel(c1);
options(scipen=999)

for (nsub in c(10,25,50)){# number of time series per group
  print(nsub);
  out1=foreach(r=1:R, .combine='rbind',.packages=c('astsa','fossil')) %dorng% {
    print(r);
    nobs = n
    m = 4   # number of categories
    
    beta1 = cbind(c(0, 3, 1, 1),
                  c(0, 1, 3, 1),
                  c(0, 1, 1, 3))
    
    beta2 = cbind(c(0, -1, 1, 1),
                  c(0, 1, -1, 1),
                  c(0, 1, 1, -1))
    
    beta3 = cbind(c(0, -1, 1, 1),
                  c(0, 1, 1, 1),
                  c(0, 1, 1, 2))
    
    beta4 = cbind(c(0, 2, 1, 1),
                  c(0, 1, 1, 1),
                  c(0, 1, 1, -1))
    
    # Group 1
    yt1 = vector("list",nsub);
    index1 = c()
    index1[1] = 1
    for(sub in 1:nsub){
      yt1[[sub]]=matrix(0,nrow=nobs,ncol=m-1);
      yt1[[sub]][1,1] = 1                            
      zt = matrix(0,nobs,m)
      for (i in 2:nobs){    
        zt[i-1,] = cbind(1, t(yt1[[sub]][i-1,]))
        deno = 1 + (exp(zt[i-1,]%*% beta1[,1]) + exp(zt[i-1,]%*%beta1[,2]) +
                      exp(zt[i-1,]%*% beta1[,3]))
        pi1 = exp(zt[i-1,]%*% beta1[,1])/deno
        pi2 = exp(zt[i-1,]%*% beta1[,2])/deno
        pi3 = exp(zt[i-1,]%*% beta1[,3])/deno
        pi4 = 1/deno
        
        index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
        if (index1[i]==4){
          yt1[[sub]][i,1:(m-1)] = rep(0,m-1)
        } else {
          yt1[[sub]][i,index1[i]]=1
        }
      }
    }
    
    # Group 2
    yt2 = vector("list",nsub);
    index1 = c()
    index1[1] = 1
    for(sub in 1:nsub){
      yt2[[sub]]=matrix(0,nrow=nobs,ncol=m-1);
      yt2[[sub]][1,1] = 1                            
      zt = matrix(0,nobs,m)
      for (i in 2:nobs){    
        zt[i-1,] = cbind(1, t(yt2[[sub]][i-1,]))
        deno = 1 + (exp(zt[i-1,]%*% beta2[,1]) + exp(zt[i-1,]%*%beta2[,2]) +
                      exp(zt[i-1,]%*% beta2[,3]))
        pi1 = exp(zt[i-1,]%*% beta2[,1])/deno
        pi2 = exp(zt[i-1,]%*% beta2[,2])/deno
        pi3 = exp(zt[i-1,]%*% beta2[,3])/deno
        pi4 = 1/deno
        
        index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
        if (index1[i]==4){
          yt2[[sub]][i,1:(m-1)] = rep(0,m-1)
        } else {
          yt2[[sub]][i,index1[i]]=1
        }
      }
    }
    
    # Group 3
    yt3 = vector("list",nsub);
    index1 = c()
    index1[1] = 1
    for(sub in 1:nsub){
      yt3[[sub]]=matrix(0,nrow=nobs,ncol=m-1);
      yt3[[sub]][1,1] = 1                            
      zt = matrix(0,nobs,m)
      for (i in 2:nobs){    
        zt[i-1,] = cbind(1, t(yt3[[sub]][i-1,]))
        deno = 1 + (exp(zt[i-1,]%*% beta3[,1]) + exp(zt[i-1,]%*%beta3[,2]) +
                      exp(zt[i-1,]%*% beta3[,3]))
        pi1 = exp(zt[i-1,]%*% beta3[,1])/deno
        pi2 = exp(zt[i-1,]%*% beta3[,2])/deno
        pi3 = exp(zt[i-1,]%*% beta3[,3])/deno
        pi4 = 1/deno
        
        index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
        if (index1[i]==4){
          yt3[[sub]][i,1:(m-1)] = rep(0,m-1)
        } else {
          yt3[[sub]][i,index1[i]]=1
        }
      }
    }
    
    # Group 4
    yt4 = vector("list",nsub);
    index1 = c()
    index1[1] = 1
    for(sub in 1:nsub){
      yt4[[sub]]=matrix(0,nrow=nobs,ncol=m-1);
      yt4[[sub]][1,1] = 1                            
      zt = matrix(0,nobs,m)
      for (i in 2:nobs){    
        zt[i-1,] = cbind(1, t(yt4[[sub]][i-1,]))
        deno = 1 + (exp(zt[i-1,]%*% beta4[,1]) + exp(zt[i-1,]%*%beta4[,2]) +
                      exp(zt[i-1,]%*% beta4[,3]))
        pi1 = exp(zt[i-1,]%*% beta4[,1])/deno
        pi2 = exp(zt[i-1,]%*% beta4[,2])/deno
        pi3 = exp(zt[i-1,]%*% beta4[,3])/deno
        pi4 = 1/deno
        
        index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
        if (index1[i]==4){
          yt4[[sub]][i,1:(m-1)] = rep(0,m-1)
        } else {
          yt4[[sub]][i,index1[i]]=1
        }
      }
    }
    
    yt = c(yt1,yt2,yt3,yt4);
    group_true = c(rep(1,nsub), rep(2,nsub),rep(3,nsub),rep(4,nsub));
    
    #estimate spectral envelope and scalings
    L=floor(sqrt(n));std=TRUE;
    tmp=envsca.get(yt,L,std);
    
    
    #regular kmeans
    k=4;s=25;
    spar=FALSE;seed=sample(1:1000,1);
    tmp1=envsca_kmeans(tmp$rf,tmp$env_ind, tmp$sca_ind, k, s, spar, seed)
    
    ##sparse clustering
    k=4;s=7;
    spar=TRUE;seed=sample(1:1000,1);
    tmp2=envsca_kmeans(tmp$rf,tmp$env_ind, tmp$sca_ind, k, s, spar, seed)
    
    #return rand indices
    ri.reg=adj.rand.index(tmp1$group_est,group_true);
    ri.sp=adj.rand.index(tmp2$group_est,group_true);
    
    return(c(ri.reg,ri.sp))
    
  }
  save.image(paste('EnvClust_ex1_R100n200nsub',nsub,'.RData',sep=''));
}
stopCluster(c1);

load('EnvClust_ex1_R100n200nsub10.RData')
outlong=out1[1:R,];
load('EnvClust_ex1_R100n200nsub25.RData')
outlong=rbind(outlong,out1[1:R,])
load('EnvClust_ex1_R100n200nsub50.RData')
outlong=rbind(outlong,out1[1:R,])

colnames(outlong) = c("K Means","Sparse K Means")
rownames(outlong)=rep("",nrow(outlong))
outlong=cbind(melt(outlong),rep(rep(c(10,25,50),each=R),2),rep(1:R,6))[,-1]
colnames(outlong)=c("ClusteringAlgorithm","AdjustedRandIndex","N_k","Replicate")

randplot=ggplot(data=as.data.frame(outlong),aes(x=ClusteringAlgorithm,y=AdjustedRandIndex,fill=ClusteringAlgorithm))+geom_boxplot()+
  hw+facet_wrap(~N_k,labeller = label_both)+
  labs(x="Clustering Algorithm",y="Adjusted Rand Index",title="Adjusted Rand Index by Cluster Size")+theme(legend.position="none")

fig3=grid.arrange(grobs=list(screep,gapstatplot,envwgt,scawgt,randplot),
             layout_matrix = rbind(c(1,2),
                                   c(3,4),
                                   c(5,5))
)
fig3

##get mean ARIs
round(mean(outlong$AdjustedRandIndex[outlong$ClusteringAlgorithm=="K Means" & outlong$N_k == 10]),2)
round(mean(outlong$AdjustedRandIndex[outlong$ClusteringAlgorithm=="K Means" & outlong$N_k == 25]),2)
round(mean(outlong$AdjustedRandIndex[outlong$ClusteringAlgorithm=="K Means" & outlong$N_k == 50]),2)

round(mean(outlong$AdjustedRandIndex[outlong$ClusteringAlgorithm=="Sparse K Means" & outlong$N_k == 10]),2)
round(mean(outlong$AdjustedRandIndex[outlong$ClusteringAlgorithm=="Sparse K Means" & outlong$N_k == 25]),2)
round(mean(outlong$AdjustedRandIndex[outlong$ClusteringAlgorithm=="Sparse K Means" & outlong$N_k == 50]),2)

#test for location shifts in ARIs
wilcox.test(outlong$AdjustedRandIndex[outlong$ClusteringAlgorithm=="K Means" & outlong$N_k == 10],
            outlong$AdjustedRandIndex[outlong$ClusteringAlgorithm=="Sparse K Means" & outlong$N_k == 10],
            alternative='two.sided')

wilcox.test(outlong$AdjustedRandIndex[outlong$ClusteringAlgorithm=="K Means" & outlong$N_k == 25],
            outlong$AdjustedRandIndex[outlong$ClusteringAlgorithm=="Sparse K Means" & outlong$N_k == 25],
            alternative='two.sided')

wilcox.test(outlong$AdjustedRandIndex[outlong$ClusteringAlgorithm=="K Means" & outlong$N_k == 50],
            outlong$AdjustedRandIndex[outlong$ClusteringAlgorithm=="Sparse K Means" & outlong$N_k == 50],
            alternative='two.sided')


####################################
## Figure 4
####################################

R=100 #replications
nsub=10 #cluster size
set.seed(675)

#start cluster
ncore <- detectCores();
c1 <- makeCluster(max(1,ncore-2));
registerDoParallel(c1);
options(scipen=999)

for (n in c(200,400,1000)){# length of time series
  
  print(n);
  out1=foreach(r=1:R, .combine='rbind',.packages=c('astsa','fossil')) %dorng% {
    print(r);
    nobs = n
    m = 4   # number of categories

    beta1 = cbind(c(0, 3, 1, 1),
                  c(0, 1, 3, 1),
                  c(0, 1, 1, 3))
    
    beta2 = cbind(c(0, -1, 1, 1),
                  c(0, 1, -1, 1),
                  c(0, 1, 1, -1))
    
    beta3 = cbind(c(0, -1, 1, 1),
                  c(0, 1, 1, 1),
                  c(0, 1, 1, 2))
    
    beta4 = cbind(c(0, 2, 1, 1),
                  c(0, 1, 1, 1),
                  c(0, 1, 1, -1))
    
    # Group 1
    yt1 = vector("list",nsub);
    index1 = c()
    index1[1] = 1
    for(sub in 1:nsub){
      if(sub>(nsub-2)){
        eps=rbind(rep(0,3),matrix(runif(9,0,1),ncol=3))
        beta=beta1+eps
      }else{beta=beta1}
      yt1[[sub]]=matrix(0,nrow=nobs,ncol=m-1);
      yt1[[sub]][1,1] = 1 
      
      zt = matrix(0,nobs,m)
      for (i in 2:nobs){    
        zt[i-1,] = cbind(1, t(yt1[[sub]][i-1,]))
        deno = 1 + (exp(zt[i-1,]%*% beta[,1]) + exp(zt[i-1,]%*%beta[,2]) +
                      exp(zt[i-1,]%*% beta[,3]))
        pi1 = exp(zt[i-1,]%*% beta[,1])/deno
        pi2 = exp(zt[i-1,]%*% beta[,2])/deno
        pi3 = exp(zt[i-1,]%*% beta[,3])/deno
        pi4 = 1/deno
        
        index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
        if (index1[i]==4){
          yt1[[sub]][i,1:(m-1)] = rep(0,m-1)
        } else {
          yt1[[sub]][i,index1[i]]=1
        }
      }
    }
    
    # Group 2
    yt2 = vector("list",nsub);
    index1 = c()
    index1[1] = 1
    for(sub in 1:nsub){
      if(sub>(nsub-2)){
        eps=rbind(rep(0,3),matrix(runif(9,0,1),ncol=3))
        beta=beta2+eps
      }else{beta=beta2}
      yt2[[sub]]=matrix(0,nrow=nobs,ncol=m-1);
      yt2[[sub]][1,1] = 1 
      zt = matrix(0,nobs,m)
      for (i in 2:nobs){    
        zt[i-1,] = cbind(1, t(yt2[[sub]][i-1,]))
        deno = 1 + (exp(zt[i-1,]%*% beta[,1]) + exp(zt[i-1,]%*%beta[,2]) +
                      exp(zt[i-1,]%*% beta[,3]))
        pi1 = exp(zt[i-1,]%*% beta[,1])/deno
        pi2 = exp(zt[i-1,]%*% beta[,2])/deno
        pi3 = exp(zt[i-1,]%*% beta[,3])/deno
        pi4 = 1/deno
        
        index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
        if (index1[i]==4){
          yt2[[sub]][i,1:(m-1)] = rep(0,m-1)
        } else {
          yt2[[sub]][i,index1[i]]=1
        }
      }
    }
    
    # Group 3
    yt3 = vector("list",nsub);
    index1 = c()
    index1[1] = 1
    for(sub in 1:nsub){
      if(sub>(nsub-2)){
        eps=rbind(rep(0,3),matrix(runif(9,0,1),ncol=3))
        beta=beta3+eps
      }else{beta=beta3}
      yt3[[sub]]=matrix(0,nrow=nobs,ncol=m-1);
      yt3[[sub]][1,1] = 1  
      zt = matrix(0,nobs,m)
      for (i in 2:nobs){    
        zt[i-1,] = cbind(1, t(yt3[[sub]][i-1,]))
        deno = 1 + (exp(zt[i-1,]%*% beta[,1]) + exp(zt[i-1,]%*%beta[,2]) +
                      exp(zt[i-1,]%*% beta[,3]))
        pi1 = exp(zt[i-1,]%*% beta[,1])/deno
        pi2 = exp(zt[i-1,]%*% beta[,2])/deno
        pi3 = exp(zt[i-1,]%*% beta[,3])/deno
        pi4 = 1/deno
        
        index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
        if (index1[i]==4){
          yt3[[sub]][i,1:(m-1)] = rep(0,m-1)
        } else {
          yt3[[sub]][i,index1[i]]=1
        }
      }
    }
    
    # Group 4
    yt4 = vector("list",nsub);
    index1 = c()
    index1[1] = 1
    for(sub in 1:nsub){
      if(sub>(nsub-2)){
        eps=rbind(rep(0,3),matrix(runif(9,0,1),ncol=3))
        beta=beta4+eps
      }else{beta=beta4}
      yt4[[sub]]=matrix(0,nrow=nobs,ncol=m-1);
      yt4[[sub]][1,1] = 1   
      zt = matrix(0,nobs,m)
      for (i in 2:nobs){    
        zt[i-1,] = cbind(1, t(yt4[[sub]][i-1,]))
        deno = 1 + (exp(zt[i-1,]%*% beta[,1]) + exp(zt[i-1,]%*%beta[,2]) +
                      exp(zt[i-1,]%*% beta[,3]))
        pi1 = exp(zt[i-1,]%*% beta[,1])/deno
        pi2 = exp(zt[i-1,]%*% beta[,2])/deno
        pi3 = exp(zt[i-1,]%*% beta[,3])/deno
        pi4 = 1/deno
        
        index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
        if (index1[i]==4){
          yt4[[sub]][i,1:(m-1)] = rep(0,m-1)
        } else {
          yt4[[sub]][i,index1[i]]=1
        }
      }
    }
    
    yt = c(yt1,yt2,yt3,yt4);
    group_true = c(rep(1,nsub), rep(2,nsub),rep(3,nsub),rep(4,nsub));
  
    #estimate spectral envelope and scalings
    L=floor(sqrt(n));std=TRUE;
    tmp=envsca.get(yt,L,std);
    
    #regular kmedoids
    k=4;s=25;
    spar=FALSE;seed=sample(1:1000,1);
    tmp1=envsca_kmedoids(tmp$rf,tmp$env_ind, tmp$sca_ind, k, s, spar, seed)
    
    #regular kmeans
    k=4;s=25;
    spar=FALSE;seed=sample(1:1000,1);
    tmp2=envsca_kmeans(tmp$rf,tmp$env_ind, tmp$sca_ind, k, s, spar, seed)
    
    #return rand indices
    ri.medoid=adj.rand.index(tmp1$group_est,group_true);
    ri.mean=adj.rand.index(tmp2$group_est,group_true);
    
    print(c(ri.medoid,ri.mean))
    
    return(c(ri.medoid,ri.mean))
    
  }
  save.image(paste('EnvClust_ex2_R100nsub10n',n,'.RData',sep=''));
}
stopCluster(c1);

load('EnvClust_ex2_R100nsub10n200.RData')
outlong2=out1[1:R,];
load('EnvClust_ex2_R100nsub10n400.RData')
outlong2=rbind(outlong2,out1[1:R,])
load('EnvClust_ex2_R100nsub10n1000.RData')
outlong2=rbind(outlong2,out1[1:R,])

colnames(outlong2) = c("K Medoids","K Means")
rownames(outlong2)=rep("",nrow(outlong2))
outlong2=cbind(melt(outlong2),rep(rep(c(200,400,1000),each=R),2),rep(1:R,6))[,-1]
colnames(outlong2)=c("ClusteringAlgorithm","AdjustedRandIndex","T","Replicate")

fig4=ggplot(data=as.data.frame(outlong2),aes(x=ClusteringAlgorithm,y=AdjustedRandIndex,fill=ClusteringAlgorithm))+geom_boxplot()+
  hw+facet_wrap(~T,labeller = label_both)+
  labs(x="Clustering Algorithm",y="Adjusted Rand Index",title="Adjusted Rand Index by Series Length")+theme(legend.position="none")
fig4

##get mean ARIs
round(mean(outlong2$AdjustedRandIndex[outlong2$ClusteringAlgorithm=="K Means" & outlong2$T == 200]),2)
round(mean(outlong2$AdjustedRandIndex[outlong2$ClusteringAlgorithm=="K Means" & outlong2$T == 400]),2)
round(mean(outlong2$AdjustedRandIndex[outlong2$ClusteringAlgorithm=="K Means" & outlong2$T == 1000]),2)

round(mean(outlong2$AdjustedRandIndex[outlong2$ClusteringAlgorithm=="K Medoids" & outlong2$T == 200]),2)
round(mean(outlong2$AdjustedRandIndex[outlong2$ClusteringAlgorithm=="K Medoids" & outlong2$T == 400]),2)
round(mean(outlong2$AdjustedRandIndex[outlong2$ClusteringAlgorithm=="K Medoids" & outlong2$T == 1000]),2)

#test for location shifts in ARIs
wilcox.test(outlong2$AdjustedRandIndex[outlong2$ClusteringAlgorithm=="K Means" & outlong2$T == 200],
            outlong2$AdjustedRandIndex[outlong2$ClusteringAlgorithm=="K Medoids" & outlong2$T == 200],
            alternative='two.sided')

wilcox.test(outlong2$AdjustedRandIndex[outlong2$ClusteringAlgorithm=="K Means" & outlong2$T == 400],
            outlong2$AdjustedRandIndex[outlong2$ClusteringAlgorithm=="K Medoids" & outlong2$T == 400],
            alternative='two.sided')

wilcox.test(outlong2$AdjustedRandIndex[outlong2$ClusteringAlgorithm=="K Means" & outlong2$T == 1000],
            outlong2$AdjustedRandIndex[outlong2$ClusteringAlgorithm=="K Medoids" & outlong2$T == 1000],
            alternative='two.sided')


####################################
## Figure 5
####################################

##case 2, 3 clusters with switching series
## 

n = 400; #length of time series
nsub = 25 # number of time series per group

set.seed(95) 
nobs = n
m = 4   # number of categories

beta1 = cbind(c(0, 2.5, 1, 1),
              c(0, 1, 1, 1),
              c(0, 1, 1, 2.5))
beta2 = cbind(c(0, 1, 1, 2),
              c(0, 1, 1, 1),
              c(0, 2, 1, 1))

# Group 1
yt1 = vector("list",nsub);
index1 = c()
index1[1] = 1
for(sub in 1:nsub){
  yt1[[sub]]=matrix(0,nrow=nobs,ncol=m-1);
  yt1[[sub]][1,1] = 1                            
  zt = matrix(0,nobs,m)
  for (i in 2:nobs){    
    zt[i-1,] = cbind(1, t(yt1[[sub]][i-1,]))
    deno = 1 + (exp(zt[i-1,]%*% beta1[,1]) + exp(zt[i-1,]%*%beta1[,2]) +
                  exp(zt[i-1,]%*% beta1[,3]))
    pi1 = exp(zt[i-1,]%*% beta1[,1])/deno
    pi2 = exp(zt[i-1,]%*% beta1[,2])/deno
    pi3 = exp(zt[i-1,]%*% beta1[,3])/deno
    pi4 = 1/deno
    
    index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
    if (index1[i]==4){
      yt1[[sub]][i,1:(m-1)] = rep(0,m-1)
    } else {
      yt1[[sub]][i,index1[i]]=1
    }
  }
}

# Group 2
yt2 = vector("list",nsub);
index1 = c()
index1[1] = 1
for(sub in 1:nsub){
  yt2[[sub]]=matrix(0,nrow=nobs,ncol=m-1);
  yt2[[sub]][1,1] = 1                            
  zt = matrix(0,nobs,m)
  for (i in 2:nobs){    
    zt[i-1,] = cbind(1, t(yt2[[sub]][i-1,]))
    deno = 1 + (exp(zt[i-1,]%*% beta2[,1]) + exp(zt[i-1,]%*%beta2[,2]) +
                  exp(zt[i-1,]%*% beta2[,3]))
    pi1 = exp(zt[i-1,]%*% beta2[,1])/deno
    pi2 = exp(zt[i-1,]%*% beta2[,2])/deno
    pi3 = exp(zt[i-1,]%*% beta2[,3])/deno
    pi4 = 1/deno
    
    index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
    if (index1[i]==4){
      yt2[[sub]][i,1:(m-1)] = rep(0,m-1)
    } else {
      yt2[[sub]][i,index1[i]]=1
    }
  }
}

# Group 3 (switching)
yt3 = vector("list",floor(5));
index1 = c()
index1[1] = 1
for(sub in 1:floor(5)){
  rsam=runif(n=1,min=0.75,max=0.75);
  yt3[[sub]]=matrix(0,nrow=nobs,ncol=m-1);
  yt3[[sub]][1,1] = 1
  zt = matrix(0,nobs,m)
  for (i in 2:floor(nobs*rsam)){
    zt[i-1,] = cbind(1, t(yt3[[sub]][i-1,]))
    deno = 1 + (exp(zt[i-1,]%*% beta1[,1]) 
                + exp(zt[i-1,]%*%beta1[,2]) +
                  exp(zt[i-1,]%*% beta1[,3]))
    pi1 = exp(zt[i-1,]%*% beta1[,1])/deno
    pi2 = exp(zt[i-1,]%*% beta1[,2])/deno
    pi3 = exp(zt[i-1,]%*% beta1[,3])/deno
    pi4 = 1/deno

    index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
    if (index1[i]==4){
      yt3[[sub]][i,1:(m-1)] = rep(0,m-1)
    } else {
      yt3[[sub]][i,index1[i]]=1
    }
  }
  for (i in floor(nobs*rsam):nobs){
    zt[i-1,] = cbind(1, t(yt3[[sub]][i-1,]))
    deno = 1 + (exp(zt[i-1,]%*% beta2[,1])
                + exp(zt[i-1,]%*%beta2[,2]) +
                  exp(zt[i-1,]%*% beta2[,3]))
    pi1 = exp(zt[i-1,]%*% beta2[,1])/deno
    pi2 = exp(zt[i-1,]%*% beta2[,2])/deno
    pi3 = exp(zt[i-1,]%*% beta2[,3])/deno
    pi4 = 1/deno
    
    index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
    if (index1[i]==4){
      yt3[[sub]][i,1:(m-1)] = rep(0,m-1)
    } else {
      yt3[[sub]][i,index1[i]]=1
    }
  }
}

yt = c(yt1,yt2,yt3);
group_true = c(rep(1,nsub), rep(2,nsub),rep("Switching",floor(5)));

#plot time series, envelopes, and scalings
cts = matrix(sapply(yt,function(x) apply(x,1,function(y) ifelse(length(which(y==1))==0,4,which(y==1)))),ncol=1);
cts=data.frame(cts,rep(1:n,length(group_true)),rep(c('1','2','Switching'),times=c(nsub*n,nsub*n,floor(5)*n)),rep(1:length(group_true),each=n))
colnames(cts)=c('Category','Time','Cluster','Series')
cts=as.data.frame(cts)
cts$Cluster = factor(cts$Cluster,levels=c('1','2','Switching'));
cts$Series = factor(cts$Series,levels=1:(nsub*3));

tsplot=ggplot(data=cts[cts$Series %in% c(1,nsub+1,2*nsub+1,3*nsub+1),],aes(x=Time,y=Category,color=Cluster))+geom_line()+
  scale_y_continuous(expand=c(0,0),limits=c(0.5,4.5),breaks=1:4)+
  scale_x_continuous(expand = c(0, 0),limits=c(1,n),breaks=c(1,seq(0,n-50,by=50)[-1]))+hw+
  xlab("Time")+ylab("Category")+labs(color = "Cluster")+facet_wrap(~Cluster,ncol=1)+theme(legend.position="top")

L=floor(sqrt(n));std=TRUE;
tmp=envsca.get(yt,L,std);

envmat=as.data.frame(unlist(tmp$env_ind));
envmat=cbind(envmat,rep(tmp$rf,length(group_true)),
             rep(1:3,times=c(nsub*length(tmp$rf),nsub*length(tmp$rf),floor(5)*length(tmp$rf))),
                 rep(1:length(group_true),each=length(tmp$rf)))
names(envmat)=c('env','f','k','i');
envmat$k = factor(envmat$k,levels=1:3)


scamat=as.data.frame(unlist(tmp$sca_ind));
scamat=cbind(scamat,rep(rep(1:3,each=length(tmp$rf)),length(group_true)),
             rep(tmp$rf,3*length(group_true)),rep(1:3,times=c(3*nsub*length(tmp$rf),3*nsub*length(tmp$rf),3*floor(5)*length(tmp$rf))),
             rep(1:length(group_true),each=3*length(tmp$rf)))
names(scamat)=c('sca','s','f','k','i');
scamat$k = factor(scamat$k,levels=1:3)
scamat$s = factor(scamat$s,levels=1:3)


#average by cluster
envmeans=aggregate(envmat$env,by=envmat[,c('f','k')],FUN=mean)
envplot=ggplot(data=envmeans,aes(x=f,y=x,color=k))+
  geom_line()+hw+scale_y_continuous(expand=c(0,0),limits=c(0,3.5))+
  scale_x_continuous(expand = c(0, 0),limits=c(0,0.5),breaks=c(0,.1,.2,.3,.4))+
  xlab("Frequency (Hz)")+ylab("Standardized Spectral Envelope")+labs(color = "Cluster")+theme(legend.position="none")


scameans=aggregate(scamat$sca,by=scamat[,c('s','f','k')],FUN=mean)
scaplot=ggplot(data=scameans,aes(x=f,y=x,color=k))+hw+scale_x_continuous(expand=c(0,0),limits=c(0,.5),breaks=c(0,.1,.2,.3,.4))+
  scale_y_continuous(expand=c(0,0),limits=c(-1.25,1.25))+geom_line()+facet_wrap(~s,nrow=1)+
  xlab("Frequency (Hz)")+ylab("Standardized Scalings")+labs(color = "Cluster")+theme(legend.position="none")

#arrange into figure
fig5=grid.arrange(tsplot,envplot,scaplot,layout_matrix=rbind(c(1,1),c(1,1),c(1,1),c(1,1),
                                                        c(2,2),c(2,2),
                                                        c(3,3),c(3,3)))
grid::grid.draw(fig5)


####################################
## Figure 6
####################################
####################################
## WARNING: Gap statistic calculations may take a long time (>30 minutes)
## depending on computing resources available
####################################

#scree plot
k=1:10
seed=54;

#regular k means
plan(multisession);
test=future_sapply(k,function(x) envsca_kmeans(tmp$rf,tmp$env_ind, tmp$sca_ind, x, s=25,spar=FALSE, seed)$wcss,future.seed=TRUE);
df=data.frame(objfn=test,k=k)

screep=ggplot(data=df,aes(x=k,y=objfn))+geom_line()+geom_point()+hw+
  scale_x_continuous(breaks=1:max(df$k),minor_breaks=1:max(df$k),expand = c(0, 0),limits=c(0.8,max(k)+.2))+
  labs(title="Scree Plot",x="Number of Clusters (k)",y="Within-Cluster SS")

#select k=2  and then use gapstat to find m
k=2;
m=seq(from=1.02,7,length.out=25);B=10;seed=6;
gpst_m=gapstat_m(tmp$rf, tmp$env_ind, tmp$sca_ind, k, m, B, seed)

gapfuzzdf=data.frame(m=m,gapst=gpst_m$gapst);
gappfuzz2=ggplot(data=gapfuzzdf,aes(x=m,y=gapst))+geom_line()+geom_point()+hw+
  scale_x_continuous(breaks=seq(1,max(m),by=1),minor_breaks=seq(1,max(m),by=0.5),expand = c(0, 0))+
  labs(title="Fuzziness Gap Statistic",x="Fuzziness Parameter (m)",y="Gap Statistic")

gapfuzzdf$m[which(gapfuzzdf$gapst==max(gapfuzzdf$gapst))[1]] #3.511667
gapfuzzdf$m[which(gapfuzzdf$gapst>(max(gapfuzzdf$gapst)-sd(gapfuzzdf$gapst)))[1]] #1.7675

#fuzzy clustering with specific m
k=2;m=1.7675;
tmp2=envsca_fuzzyclust_wrap(tmp$rf, tmp$env_ind, tmp$sca_ind, k, m, nstart=50, seed=23)


dfbpsim=melt(tmp2$umat)
dfbpsim=cbind(dfbpsim,rep(group_true,k))
names(dfbpsim)=c('i','cluster','u','group')
dfbpsim$cluster= factor(dfbpsim$cluster)
dfbpsim$group= factor(dfbpsim$group)

uplot=ggplot(data=dfbpsim[dfbpsim$cluster==1,], aes(y=u , x=cluster, fill=group))+ 
  geom_violin()+geom_jitter()+facet_wrap(~group)+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))+
  xlab("") + ylab("Membership")+hw+theme(legend.position = "none")+theme(axis.text.x=element_blank())+
  ggtitle("First Cluster Membership Degree")


##replicated simulations for fuzzy clustering
############################################
R=100 #replications
n = 400; #length of time series

set.seed(94)
for (mfuzz in c(1.7675,1.1)){
  for (nsub in c(10,25,50,200)){# number of time series per group
    print(nsub);
    out1=foreach(r=1:R, .combine='rbind',.packages=c('astsa','fossil','fclust','future.apply')) %do% {
      print(r);
      nobs = n
      m = 4   # number of categories
      
      beta1 = cbind(c(0, 2.5, 1, 1),
                    c(0, 1, 1, 1),
                    c(0, 1, 1, 2.5))
      beta2 = cbind(c(0, 1, 1, 2),
                    c(0, 1, 1, 1),
                    c(0, 2, 1, 1))
      
      # Group 1
      yt1 = vector("list",nsub);
      index1 = c()
      index1[1] = 1
      for(sub in 1:nsub){
        yt1[[sub]]=matrix(0,nrow=nobs,ncol=m-1);
        yt1[[sub]][1,1] = 1                            
        zt = matrix(0,nobs,m)
        for (i in 2:nobs){    
          zt[i-1,] = cbind(1, t(yt1[[sub]][i-1,]))
          deno = 1 + (exp(zt[i-1,]%*% beta1[,1]) + exp(zt[i-1,]%*%beta1[,2]) +
                        exp(zt[i-1,]%*% beta1[,3]))
          pi1 = exp(zt[i-1,]%*% beta1[,1])/deno
          pi2 = exp(zt[i-1,]%*% beta1[,2])/deno
          pi3 = exp(zt[i-1,]%*% beta1[,3])/deno
          pi4 = 1/deno
          
          index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
          if (index1[i]==4){
            yt1[[sub]][i,1:(m-1)] = rep(0,m-1)
          } else {
            yt1[[sub]][i,index1[i]]=1
          }
        }
      }
      
      # Group 2
      yt2 = vector("list",nsub);
      index1 = c()
      index1[1] = 1
      for(sub in 1:nsub){
        yt2[[sub]]=matrix(0,nrow=nobs,ncol=m-1);
        yt2[[sub]][1,1] = 1                            
        zt = matrix(0,nobs,m)
        for (i in 2:nobs){    
          zt[i-1,] = cbind(1, t(yt2[[sub]][i-1,]))
          deno = 1 + (exp(zt[i-1,]%*% beta2[,1]) + exp(zt[i-1,]%*%beta2[,2]) +
                        exp(zt[i-1,]%*% beta2[,3]))
          pi1 = exp(zt[i-1,]%*% beta2[,1])/deno
          pi2 = exp(zt[i-1,]%*% beta2[,2])/deno
          pi3 = exp(zt[i-1,]%*% beta2[,3])/deno
          pi4 = 1/deno
          
          index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
          if (index1[i]==4){
            yt2[[sub]][i,1:(m-1)] = rep(0,m-1)
          } else {
            yt2[[sub]][i,index1[i]]=1
          }
        }
      }
      
      # Group 3 (switching)
      yt3 = vector("list",floor(5));
      index1 = c()
      index1[1] = 1
      for(sub in 1:floor(5)){
        rsam=runif(n=1,min=0.75,max=0.75);
        yt3[[sub]]=matrix(0,nrow=nobs,ncol=m-1);
        yt3[[sub]][1,1] = 1
        zt = matrix(0,nobs,m)
        for (i in 2:floor(nobs*rsam)){
          zt[i-1,] = cbind(1, t(yt3[[sub]][i-1,]))
          deno = 1 + (exp(zt[i-1,]%*% beta1[,1]) 
                      + exp(zt[i-1,]%*%beta1[,2]) +
                        exp(zt[i-1,]%*% beta1[,3]))
          pi1 = exp(zt[i-1,]%*% beta1[,1])/deno
          pi2 = exp(zt[i-1,]%*% beta1[,2])/deno
          pi3 = exp(zt[i-1,]%*% beta1[,3])/deno
          pi4 = 1/deno
          
          index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
          if (index1[i]==4){
            yt3[[sub]][i,1:(m-1)] = rep(0,m-1)
          } else {
            yt3[[sub]][i,index1[i]]=1
          }
        }
        for (i in floor(nobs*rsam):nobs){
          zt[i-1,] = cbind(1, t(yt3[[sub]][i-1,]))
          deno = 1 + (exp(zt[i-1,]%*% beta2[,1])
                      + exp(zt[i-1,]%*%beta2[,2]) +
                        exp(zt[i-1,]%*% beta2[,3]))
          pi1 = exp(zt[i-1,]%*% beta2[,1])/deno
          pi2 = exp(zt[i-1,]%*% beta2[,2])/deno
          pi3 = exp(zt[i-1,]%*% beta2[,3])/deno
          pi4 = 1/deno
          
          index1[i] = which(rmultinom(1,1, c(pi1,pi2,pi3, pi4))==1)
          if (index1[i]==4){
            yt3[[sub]][i,1:(m-1)] = rep(0,m-1)
          } else {
            yt3[[sub]][i,index1[i]]=1
          }
        }
      }
      
      yt = c(yt1,yt2,yt3);
      group_true = c(rep(1,nsub), rep(2,nsub),rep(3,floor(5)));
      
      #estimate spectral envelope and scalings
      L=floor(sqrt(n));std=TRUE;
      tmp=envsca.get(yt,L,std);
      
      #fuzzy clustering - fit
      k=2;seed=sample(1:1000,1);
      tmp1=envsca_fuzzyclust_wrap(tmp$rf, tmp$env_ind, tmp$sca_ind, k, mfuzz, nstart=1, seed=seed)

      #return rand indices
      ri.reg=ARI.F(group_true,tmp1$umat)
      print(ri.reg)
      return(c(ri.reg))
      
    }
    save.image(paste('EnvClust_ex3_R100n400nsub',nsub,'mfuzz',mfuzz,'.RData',sep=''));
  }
}

load('EnvClust_ex3_R100n400nsub10mfuzz1.7675.RData')
outlong3=out1[1:R];
# load('EnvClust_ex3_R100n400nsub25mfuzz1.7675.RData')
# outlong3=c(outlong3,out1[1:R])
load('EnvClust_ex3_R100n400nsub50mfuzz1.7675.RData')
outlong3=c(outlong3,out1[1:R])
load('EnvClust_ex3_R100n400nsub200mfuzz1.7675.RData')
outlong3=c(outlong3,out1[1:R])
load('EnvClust_ex3_R100n400nsub10mfuzz1.1.RData')
outlong3=c(outlong3,out1[1:R])
# load('EnvClust_ex3_R100n400nsub25mfuzz1.1.RData')
# outlong3=c(outlong3,out1[1:R])
load('EnvClust_ex3_R100n400nsub50mfuzz1.1.RData')
outlong3=c(outlong3,out1[1:R])
load('EnvClust_ex3_R100n400nsub200mfuzz1.1.RData')
outlong3=c(outlong3,out1[1:R])

outlong3=cbind(c(rep(1.7675,length(outlong3)/2),rep(1.1,length(outlong3)/2)),
              outlong3,
              rep(rep(c(10,50,200),each=R),2),rep(1:R,8))
colnames(outlong3)=c("m","AdjustedRandIndex","N_k","Replicate")
outlong3=as.data.frame(outlong3)
outlong3$N_k = factor(outlong3$N_k,levels=c(10,50,200))
outlong3$m = factor(outlong3$m,levels=c(1.7675,1.1))
outlong3$AdjustedRandIndex=as.numeric(outlong3$AdjustedRandIndex)
randplot2=ggplot(data=outlong3,aes(x=m,y=AdjustedRandIndex,fill=m))+geom_boxplot()+
  hw+ylim(c(0.3,1))+facet_wrap(~N_k,labeller = label_both,nrow=1)+theme(legend.position="none")+
  labs(x="Fuzziness Parameter (m)",y="Adjusted Rand Index",title="Adjusted Rand Index by Fuzziness Parameter and Cluster Size")

##get mean ARIs
round(mean(outlong3$AdjustedRandIndex[outlong3$m==1.7675 & outlong3$N_k == 10]),2)
round(mean(outlong3$AdjustedRandIndex[outlong3$m==1.7675 & outlong3$N_k == 50]),2)
round(mean(outlong3$AdjustedRandIndex[outlong3$m==1.7675 & outlong3$N_k == 200]),2)

round(mean(outlong3$AdjustedRandIndex[outlong3$m==1.1 & outlong3$N_k == 10]),2)
round(mean(outlong3$AdjustedRandIndex[outlong3$m==1.1 & outlong3$N_k == 50]),2)
round(mean(outlong3$AdjustedRandIndex[outlong3$m==1.1 & outlong3$N_k == 200]),2)

fig6=grid.arrange(screep,gappfuzz2,uplot,randplot2,
             layout_matrix=rbind(c(1,1,2,2),c(3,3,3,3),c(4,4,4,4)))
grid::grid.draw(fig6)
