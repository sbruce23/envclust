##################################################################################
##################################################################################
## R Functions for "Adaptive Clustering and Feature Selection for Categorical Time 
## Series Using Interpretable Frequency-Domain Features" 
## Author: Scott A. Bruce
## Last Updated: August 16, 2022
## R version: 4.2.1
## Tested on Windows machine (Intel i9 10 core 3.3 Ghz processor, 128 GB RAM)
##################################################################################
##################################################################################

##################
#convert a categorical ts into multivariate indicator process
##################
cat_convert <- function(xt,refcat){
  stage = sort(unique(xt))
  nobs = length(xt)
  idx = which(stage==refcat)
  stage = as.factor(levels(stage)[c(stage[-idx],stage[idx])])
  yt = matrix(0,nobs,length(stage)-1)
  for (j in 1:length(stage)-1){
    yt[,j] = (xt==stage[j])*1
  }
  return(yt)
}

##################
#interpolation function
##################
intp <- function(xl,xm,xu,yl,yu){
  if(xu==xl){
    ym=yu;
  } else{
    tmp = (xm - xl)/(xu-xl);
    ym = yl*(1-tmp) + yu*tmp;
  }
  return(ym);
}

##################
#wrapper for interpolation function
##################
intp.wrap <- function(x,freq,fref,fyy_re){
  fa = max(which(freq<=fref[x]));
  fb = min(length(freq),fa+1);
  return(intp(xl=freq[fa],xm=fref[x],xu=freq[fb],
              yl=fyy_re[,,fa],yu=fyy_re[,,fb]));
} 

##################
#estimation function for envelope and scaling with interpolation
##################
env.get <- function(yt,L,rf=NULL,std=FALSE){
  
  #smoothed estimate of the periodogram matrix
  fyy = mvspec(yt, spans = c(L,L), plot=FALSE,kernel="fejer");
  fyy_re = Re(fyy$fxx);
  freq = fyy$freq;
  nb = dim(yt)[2];
  
  if(!is.null(rf)){
    #interpolate fit at reference frequencies
    fyy_re = array(sapply(1:length(rf),function(x) intp.wrap(x,freq,rf,fyy_re)),
                   dim=c(nb,nb,length(rf)));
    freq=rf;
  }
  
  #estimates of spectral envelope and scalings
  Q = diag(nb);
  num = fyy$n.used;	
  ev = apply(fyy_re,3,function(x) eigen(2*Q%*%x%*%Q/num,symmetric=TRUE));
  specenv = sapply(ev,function(x) x$values[1]);
  
  
  b = sapply(ev,function(x) Q%*%x$vectors[,1]);
  beta = apply(b,2,function(x) x/sqrt(sum(x^2)));

  #standardize envelope and scalings
  # if (std){specenv = specenv/sum(diff(rf)[1]*specenv)}
  ncat=nrow(beta)+1;
  if (std){
    specenv = specenv/sqrt(sum(diff(rf)[1]*specenv^2)); #standardize by l2 norm
    beta = beta/sqrt(diff(rf)[1]*length(rf)); #standardize by l2 norm
    
  } 
  
  #output results
  output = list(freq = freq, 
                envelope = matrix(specenv,ncol=1), 
                scale = t(beta));
  return(output)
}   

##################
#estimation function for group envelope and scaling by averaging individual estimates
##################
group_env <- function(env_ind,sca_ind){
  
  nsub = length(env_ind);
  ncat=dim(sca_ind[[1]])[2]+1;
  nf=length(env_ind[[1]]);
  specenv=rowMeans(matrix(unlist(env_ind),nrow=nf));
  beta=array(unlist(sca_ind),dim=c(nf,ncat-1,nsub));
  beta=apply(beta,c(1,2),mean);
  
  output = list(envelope = specenv, scale = beta)
  return(output)
}

##################
#function to estimate individual envelopes and scalings
##################
envsca.get <- function(yt,L,std){
  #obtain ordering of elements of yt from shortest to longest
  ord = order(sapply(yt,length));
  
  #obtain reference frequencies from shortest series
  tmp=env.get(yt[[ord[1]]],L);
  rf=tmp$freq;
  
  #obtain individual envelopes and scalings
  tmp = lapply(yt,function(x) env.get(x,L,rf,std));
  env_ind = lapply(tmp,function(x) x$envelope);
  sca_ind = lapply(tmp,function(x) x$scale);
  
  return(list(rf=rf,env_ind=env_ind,sca_ind=sca_ind));
}

##################
#function for k medoids clustering
##################
envsca_kmedoids <- function(rf, env_ind, sca_ind, k, s, spar, seed=123){
  
  set.seed(seed)
  
  #initialize weights (sum of weights squared should equal 1)
  nf=length(rf);
  ncat=dim(sca_ind[[1]])[2]+1;
  if (spar){
    w.env=matrix(1/sqrt(nf*ncat),nrow=nf,ncol=1);
    w.sca=matrix(1/sqrt(nf*ncat),nrow=nf,ncol=ncat-1);
  } else{
    w.env=matrix(1,nrow=nf,ncol=1);
    w.sca=matrix(1,nrow=nf,ncol=ncat-1);
  }
  
  #initialize random starting point for clustering
  nobs=length(env_ind);
  classes = 1:k;
  meds_tmp=sample(nobs,k,replace=FALSE);
  env=env_ind[meds_tmp];
  sca=sca_ind[meds_tmp];
  print(meds_tmp)
  
  #compute distances and assign to clusters accordingly
  idx = expand.grid(1:nobs,classes);
  g = matrix(mapply(function(x,y) sum(w.env*abs(env_ind[[x]]-env[[y]])^1) +
                      sum(w.sca*abs(sca_ind[[x]]-sca[[y]])^1/(1)),
                    x=idx[,1],y=idx[,2]),nrow=nobs,ncol=k);
  group_tmp = apply(g,1,function(x) which(x==min(x)));
  
  #compute cost
  cost_tmp = sum(apply(g,1, min))
  print(cost_tmp)
  
  #while loop to continue until assignments don't change
  stp=0;
  while(stp==0){
    
    #step 1) with fixed weights, run k-medoids
    ########################################
    stp.km=0;
    while(stp.km==0){
      
      #compute distances for all swaps of each obs and each center
      costall=matrix(NA,ncol=1,nrow=nrow(idx))
      for (i in 1:nrow(idx)){
        
        #propose each obs as new kth medoid
        meds_prop=meds_tmp;
        meds_prop[idx[i,2]] = idx[i,1];
        
        #get new proposed medoids
        env_prop=env_ind[meds_prop];
        sca_prop=sca_ind[meds_prop];
        
        #compute distances
        g_prop = matrix(mapply(function(x,y) sum(w.env*abs(env_ind[[x]]-env_prop[[y]])^1) +
                                 sum(w.sca*abs(sca_ind[[x]]-sca_prop[[y]])^1/(1)),
                               x=idx[,1],y=idx[,2]),nrow=nobs,ncol=k);
        
        #compute cost
        cost_prop = sum(apply(g_prop,1, min))
        
        costall[i]=cost_prop
        
      }
      
      #if swap with lower cost exists, update medoids and group memberships accordingly
      if(min(costall)<cost_tmp){
        minidx=which(costall==min(costall))[1]
        meds_tmp[idx[minidx,2]] =  idx[minidx,1];
        print(meds_tmp)
        env=env_ind[meds_tmp];
        sca=sca_ind[meds_tmp];
        
        g = matrix(mapply(function(x,y) sum(w.env*abs(env_ind[[x]]-env[[y]])^1) +
                            sum(w.sca*abs(sca_ind[[x]]-sca[[y]])^1/(1)),
                          x=idx[,1],y=idx[,2]),nrow=nobs,ncol=k);
        
        cost_tmp = sum(apply(g,1, min))
        print(cost_tmp)
        group_tmp =  apply(g,1,function(x) which(x==min(x)));
        
        
        #handle empty clusters
        if (length(unique(group_tmp))<k){
          #identify empty clusters
          miss = (1:k)[!(1:k %in% unique(group_tmp))]; 
          #randomly assign point(s) to empty cluster(s)
          group_tmp[sample(1:nobs,length(miss))] = miss;
        }
        
      } else{
        stp.km=1;
      }
      
    }
    
    #step 2) with fixed clusters, update weights
    ########################################
    
    #get pairwise distances (probably faster way to do this...)
    idx = expand.grid(1:nobs,1:nobs,1:nf);
    darr.env = array(mapply(function(x,y,z) sum((env_ind[[x]][z]-env_ind[[y]][z])^2),
                            x=idx[,1],y=idx[,2],z=idx[,3]),dim=c(nobs,nobs,nf));
    
    idx = expand.grid(1:nobs,1:nobs,1:nf,1:(ncat-1));
    darr.sca = array(mapply(function(w,x,y,z) sum((sca_ind[[w]][y,z]-sca_ind[[x]][y,z])^2/(1)),
                            w=idx[,1],x=idx[,2],y=idx[,3],z=idx[,4]),dim=c(nobs,nobs,nf,ncat-1));
    
    #get array of cluster membership
    idx = expand.grid(group_tmp,group_tmp,classes);
    clarr=array(apply(idx,1,function(x) as.numeric(all(x==x[1]))),dim=c(nobs,nobs,k));
    
    #total cluster ss
    a.env.1 = apply(darr.env,3,function(x) (1/nobs)*sum(x)); #sum of all distances for each frequency
    a.sca.1 = apply(darr.sca,c(3,4),function(x) (1/nobs)*sum(x)); #sum of all distances for each frequency and categories (-1)
    
    #within cluster ss
    nclus=table(sort(group_tmp));
    idx=expand.grid(1:nf,classes);
    a.env.2 = matrix(mapply(function(x,y) sum(darr.env[,,x]*clarr[,,y])*(1/nclus[y]),
                            x=idx[,1],y=idx[,2]),nrow=nf,ncol=k);
    idx=expand.grid(1:nf,1:(ncat-1),classes);
    a.sca.2 = array(mapply(function(x,y,z) sum(darr.sca[,,x,y]*clarr[,,z]*(1/nclus[z])),
                           x=idx[,1],y=idx[,2],z=idx[,3]),dim=c(nf,ncat-1,k));
    
    #between cluster ss
    a.env=a.env.1-rowSums(a.env.2);
    a.sca=a.sca.1-apply(a.sca.2,c(1,2),sum);
    
    #soft thresholding on a
    avecpos=pmax(c(a.env,a.sca),0);
    delta=0;
    sthres=sign(avecpos)*pmax(abs(avecpos)-delta,0);
    wnew=sthres/sqrt(sum(sthres^2));
    while (sum(abs(wnew))>s){
      # print(c(sum(abs(wnew)),delta))
      
      delta=delta+0.0001;
      sthres=sign(avecpos)*pmax(abs(avecpos)-delta,0);
      wnew=sthres/sqrt(sum(sthres^2));
    }
    
    #update weights
    w.env.new=wnew[1:nf];
    w.sca.new=matrix(wnew[(nf+1):length(wnew)],nrow=nf,ncol=ncat-1);
    
    # plot(x=rf,y=w.env,type='l',ylim=c(0,0.5));
    # lines(x=rf,y=w.env.new,type='l',col='red');
    # 
    # plot(x=rf,y=w.sca[,1],type='l',ylim=c(0,0.5));
    # lines(x=rf,y=w.sca.new[,1],type='l',col='red');
    
    #check stopping criterion
    chg=(sum(abs(w.env.new-w.env))+sum(abs(w.sca.new-w.sca)))/(sum(abs(w.env))+sum(abs(w.sca)));
    if(chg<.0001 | spar==FALSE){
      objcrit= sum(w.env*a.env)+sum(w.sca*a.sca); #get objective criterion with current weights
      wcss= sum(w.env*rowSums(a.env.2))+sum(w.sca*apply(a.sca.2,c(1,2),sum)); #weighted within cluster ss
      stp=1; #stop loop
      
    } else{
      w.env=w.env.new;
      w.sca=w.sca.new;
    }
  }
  
  #get final group envelope and scalings
  env=env_ind[meds_tmp];
  sca=sca_ind[meds_tmp];
  
  #return cluster membership, weights, and estimated cluster envelope and scalings
  return(list(group_est=group_tmp,rf=rf,w.env=w.env,w.sca=w.sca,
              env_group=env,sca_group=sca,env_ind=env_ind,sca_ind=sca_ind,
              objcrit=objcrit,wcss=wcss))
  
} 

##################
# function for k means clustering
##################
envsca_kmeans <- function(rf, env_ind, sca_ind, k, s, spar, seed=123){
  
  set.seed(seed)
  
  #initialize weights (sum of weights squared should equal 1)
  nf=length(rf);
  ncat=dim(sca_ind[[1]])[2]+1;
  if (spar){
    w.env=matrix(1/sqrt(nf*ncat),nrow=nf,ncol=1);
    w.sca=matrix(1/sqrt(nf*ncat),nrow=nf,ncol=ncat-1);
  } else{
    w.env=matrix(1,nrow=nf,ncol=1);
    w.sca=matrix(1,nrow=nf,ncol=ncat-1);
  }
  
  #initialize random starting point for clustering
  nobs=length(env_ind);
  classes = 1:k;
  group_tmp=sample(classes,nobs,replace=TRUE);
  
  #handle empty clusters
  if (length(unique(group_tmp))<k){
    #identify empty clusters
    miss = (1:k)[!(1:k %in% unique(group_tmp))]; 
    #randomly assign point(s) to empty cluster(s)
    group_tmp[sample(1:nobs,length(miss))] = miss;
  }
  
  #while loop to continue until assignments don't change
  stp=0;
  while(stp==0){
  
    #step 1) with fixed weights, run kmeans
    ########################################
    stp.km=0;
    while(stp.km==0){
      
      #handle empty clusters
      if (length(unique(group_tmp))<k){
        #identify empty clusters
        miss = (1:k)[!(1:k %in% unique(group_tmp))]; 
        #randomly assign point(s) to empty cluster(s)
        group_tmp[sample(1:nobs,length(miss))] = miss;
      }
      
      #get group envelope and scalings
      tmp = lapply(classes,function(x) group_env(env_ind[group_tmp==x],sca_ind[group_tmp==x]));
      env = lapply(tmp,function(x) x$envelope);
      sca = lapply(tmp,function(x) x$scale);
      
      #compute distances
      idx = expand.grid(1:nobs,classes);
      g = matrix(mapply(function(x,y) sum(w.env*(env_ind[[x]]-env[[y]])^2) +
                          sum(w.sca*(sca_ind[[x]]-sca[[y]])^2/(1)),
                        x=idx[,1],y=idx[,2]),nrow=nobs,ncol=k);
      gmin = sum(apply(g,1,function(x) min(x)));

      #print(gmin);
      
      #reassign each observation to group with most similar measure
      group_new =  apply(g,1,function(x) which(x==min(x)));
      
      #handle empty clusters
      if (length(unique(group_new))<k){
        #identify empty clusters
        miss = (1:k)[!(1:k %in% unique(group_new))]; 
        #randomly assign point(s) to empty cluster(s)
        group_new[sample(1:nobs,length(miss))] = miss;
      }
      
      #stop if group assignment doesn't change
      #otherwise update group assignments and keep going
      if (identical(group_tmp,group_new)){
        stp.km=1;
      } else{
        group_tmp=group_new;
      } 
      
    }
  
    #step 2) with fixed clusters, update weights
    ########################################
    
    #get pairwise distances (probably faster way to do this...)
    idx = expand.grid(1:nobs,1:nobs,1:nf);
    darr.env = array(mapply(function(x,y,z) sum((env_ind[[x]][z]-env_ind[[y]][z])^2),
                            x=idx[,1],y=idx[,2],z=idx[,3]),dim=c(nobs,nobs,nf));

    idx = expand.grid(1:nobs,1:nobs,1:nf,1:(ncat-1));
    darr.sca = array(mapply(function(w,x,y,z) sum((sca_ind[[w]][y,z]-sca_ind[[x]][y,z])^2/(1)),
                            w=idx[,1],x=idx[,2],y=idx[,3],z=idx[,4]),dim=c(nobs,nobs,nf,ncat-1));
    
    #get array of cluster membership
    idx = expand.grid(group_tmp,group_tmp,classes);
    clarr=array(apply(idx,1,function(x) as.numeric(all(x==x[1]))),dim=c(nobs,nobs,k));
    
    #total cluster ss
    a.env.1 = apply(darr.env,3,function(x) (1/nobs)*sum(x)); #sum of all distances for each frequency
    a.sca.1 = apply(darr.sca,c(3,4),function(x) (1/nobs)*sum(x)); #sum of all distances for each frequency and categories (-1)
    
    #within cluster ss
    nclus=table(sort(group_tmp));
    idx=expand.grid(1:nf,classes);
    a.env.2 = matrix(mapply(function(x,y) sum(darr.env[,,x]*clarr[,,y])*(1/nclus[y]),
                            x=idx[,1],y=idx[,2]),nrow=nf,ncol=k);
    idx=expand.grid(1:nf,1:(ncat-1),classes);
    a.sca.2 = array(mapply(function(x,y,z) sum(darr.sca[,,x,y]*clarr[,,z]*(1/nclus[z])),
                            x=idx[,1],y=idx[,2],z=idx[,3]),dim=c(nf,ncat-1,k));
    
    #between cluster ss
    a.env=a.env.1-rowSums(a.env.2);
    a.sca=a.sca.1-apply(a.sca.2,c(1,2),sum);
    
    #soft thresholding on a
    avecpos=pmax(c(a.env,a.sca),0);
    delta=0;
    sthres=sign(avecpos)*pmax(abs(avecpos)-delta,0);
    wnew=sthres/sqrt(sum(sthres^2));
    while (sum(abs(wnew))>s){
      # print(c(sum(abs(wnew)),delta))
      
      delta=delta+0.0001;
      sthres=sign(avecpos)*pmax(abs(avecpos)-delta,0);
      wnew=sthres/sqrt(sum(sthres^2));
    }
    
    #update weights
    w.env.new=wnew[1:nf];
    w.sca.new=matrix(wnew[(nf+1):length(wnew)],nrow=nf,ncol=ncat-1);
    
    # plot(x=rf,y=w.env,type='l',ylim=c(0,0.5));
    # lines(x=rf,y=w.env.new,type='l',col='red');
    # 
    # plot(x=rf,y=w.sca[,1],type='l',ylim=c(0,0.5));
    # lines(x=rf,y=w.sca.new[,1],type='l',col='red');
    
    #check stopping criterion
    chg=(sum(abs(w.env.new-w.env))+sum(abs(w.sca.new-w.sca)))/(sum(abs(w.env))+sum(abs(w.sca)));
    if(chg<.0001 | spar==FALSE){
      objcrit= sum(w.env*a.env)+sum(w.sca*a.sca); #get objective criterion with current weights
      wcss= sum(w.env*rowSums(a.env.2))+sum(w.sca*apply(a.sca.2,c(1,2),sum)); #weighted within cluster ss
      stp=1; #stop loop
      
    } else{
      w.env=w.env.new;
      w.sca=w.sca.new;
    }
  }
  
  #get final group envelope and scalings
  tmp = lapply(classes,function(x) group_env(env_ind[group_tmp==x],sca_ind[group_tmp==x]));
  env = lapply(tmp,function(x) x$envelope);
  sca = lapply(tmp,function(x) x$scale);
  
  #return cluster membership, weights, and estimated cluster envelope and scalings
  return(list(group_est=group_tmp,rf=rf,w.env=w.env,w.sca=w.sca,
              env_group=env,sca_group=sca,env_ind=env_ind,sca_ind=sca_ind,
              objcrit=objcrit,wcss=wcss))
      
} 

##################
#gapstat to choose sparsity parameter s
##################
gapstat_s <- function(yt, k, L, s, B=10, seed=123, std=TRUE){
  
  spar=TRUE;
  
  #get individual estimates of envelope and scalings
  tmp=envsca.get(yt,L,std);
  
  #permute observations independently within each feature B times
  #set seed for replication
  set.seed(seed);
  nf=length(tmp$rf);
  ncat=dim(tmp$sca_ind[[1]])[2]+1;
  nobs=length(tmp$env_ind);
  ns=length(s);
  
  # envp=array(unlist(tmp$env_ind),dim=c(nf,nobs,B,ns));
  # scap=array(unlist(tmp$sca_ind),dim=c(nf,ncat-1,nobs,B,ns));
  envp=array(unlist(tmp$env_ind),dim=c(nf,nobs,B));
  scap=array(unlist(tmp$sca_ind),dim=c(nf,ncat-1,nobs,B));
  
  # envp2=apply(envp,c(1,3,4),function(x) x[sample(nobs)]);
  # scap2=apply(scap,c(1,2,4,5),function(x) x[sample(nobs)]);
  envp2=apply(envp,c(1,3),function(x) x[sample(nobs)]);
  scap2=apply(scap,c(1,2,4),function(x) x[sample(nobs)]);
  
  # envp3=lapply(seq(B), function(x) lapply(seq(ns), function(y) lapply(seq(nobs), function(z) envp2[z,,x,y])));
  # scap3=lapply(seq(B), function(x) lapply(seq(ns), function(y) lapply(seq(nobs), function(z) scap2[z,,,x,y])));
  envp3=lapply(seq(B), function(x) lapply(seq(ns), function(y) lapply(seq(nobs), function(z) envp2[z,,x])));
  scap3=lapply(seq(B), function(x) lapply(seq(ns), function(y) lapply(seq(nobs), function(z) scap2[z,,,x])));
  
  #run the cluster algo on each permuted set for each value of s and produce objective fn
  idx=expand.grid(1:B,s);
  seed1=sample(1:1000,size=nrow(idx));
  plan(multisession(workers=max(1,availableCores()-4),gc=TRUE));
  # plan(multisession(workers=2,gc=TRUE));
  # print('CPU go now!');
  out1=future_mapply(function(x,y,z) envsca_kmeans(tmp$rf,envp3[[x]][[which(s==y)]],scap3[[x]][[which(s==y)]],
                                                      k,y,spar,z)$objcrit ,x=idx[,1],y=idx[,2],z=seed1,future.seed=TRUE);

  #run cluster algo for original set for each value of s
  seed2=sample(1:1000,size=length(s));
  out2=future_mapply(function(y,z) envsca_kmeans(tmp$rf,tmp$env_ind,tmp$sca_ind,k,y,spar,z)$objcrit ,y=s,z=seed2,future.seed=TRUE); 

  #produce gap statistics
  gapst = numeric(length(s));
  for (i in 1:length(s)){
    gapst[i] = log(out2[i]) - mean(log(out1[idx[,2]==s[i]]));
  }

  #choose s = argmax(s) gapst 
  #or s as the smallest value within 1sd of log(objcrit.out1) at s of the largest gap stat
  # s.star = s[which(gapst==max(gapst))[1]];
  s.star=s[which(gapst>(max(gapst)-sd(gapst)))[1]]; 
  
  return(list(s.star=s.star,gapst.star=max(gapst),gapst=gapst));
}

##################
#function for fuzzy clustering
##################
envsca_fuzzyclust <- function(rf, env_ind, sca_ind, k, m=1.2, seed=123){
  
  nobs=length(env_ind);
  nf=length(rf);
  ncat=dim(sca_ind[[1]])[2]+1;
  classes = 1:k;
  
  #initialize membership matrix (sum of membership degrees should equal 1)
  umat=matrix(sample(1:100,size=nobs*k,replace=TRUE),nrow=nobs,ncol=k);
  umat=apply(umat,2,function(x) x/rowSums(umat));
  
  stp=0;
  while(stp==0){
    
    #get group envelope and scalings
    tmp=matrix(unlist(env_ind),nrow=nf,ncol=nobs);
    env=matrix(apply(tmp,1,function(x) (x%*%umat^m)/colSums(umat^m)),ncol=nf);
    tmp=array(unlist(sca_ind),dim=c(nf,ncat-1,nobs));
    sca=array(apply(tmp,c(1,2),function(x) (x%*%umat^m)/colSums(umat^m)),dim=c(length(classes),nf,ncat-1));
    
    # plot(x=rf,y=env[1,],type='l')
    # lines(x=rf,y=env[2,],col='red')
    
    #update membership mat
    idx = expand.grid(1:nobs,classes);
    g = matrix(mapply(function(x,y) sum(env_ind[[x]]-env[y,])^2 +
                        sum((sca_ind[[x]]-sca[y,,])^2),
                      x=idx[,1],y=idx[,2]),nrow=nobs,ncol=k);
    
    umatnew=t(matrix(apply(g,1,function(x) 1/rowSums((x%*%t(1/x))^(1/(m-1)))),ncol=nobs));
    
    #get new objective function
    objfn=sum((umatnew^m)*g);
    
    #stop if membership mat doesn't change much
    #otherwise update membership mat and keep going
    chg=sum(abs(umatnew-umat))/sum(abs(umat));
    # print(c(objfn,chg));
    if(chg<.0001){
      stp=1; #stop loop
      
    } else{
      umat=umatnew;
    }
  }
  return(list(umat=umat,env=env,sca=sca,objfn=objfn))
}

##################
#second version of fuzzy clustering that produces between cluster SS
##################
envsca_fuzzyclust2 <- function(rf, env_ind, sca_ind, k, m=1.2, seed=123){
  
  nobs=length(env_ind);
  nf=length(rf);
  ncat=dim(sca_ind[[1]])[2]+1;
  classes = 1:k;
  
  #initialize membership matrix (sum of membership degrees should equal 1)
  umat=matrix(sample(1:100,size=nobs*k,replace=TRUE),nrow=nobs,ncol=k);
  umat=apply(umat,2,function(x) x/rowSums(umat));
  
  stp=0;
  while(stp==0){
    
    #get group envelope and scalings
    tmp=matrix(unlist(env_ind),nrow=nf,ncol=nobs);
    env=matrix(apply(tmp,1,function(x) (x%*%umat^m)/colSums(umat^m)),ncol=nf);
    tmp=array(unlist(sca_ind),dim=c(nf,ncat-1,nobs));
    sca=array(apply(tmp,c(1,2),function(x) (x%*%umat^m)/colSums(umat^m)),dim=c(length(classes),nf,ncat-1));
    
    # plot(x=rf,y=env[1,],type='l')
    # lines(x=rf,y=env[2,],col='red')
    
    #update membership mat
    idx = expand.grid(1:nobs,classes);
    g = matrix(mapply(function(x,y) sum(env_ind[[x]]-env[y,])^2 +
                        sum((sca_ind[[x]]-sca[y,,])^2),
                      x=idx[,1],y=idx[,2]),nrow=nobs,ncol=k);
    
    umatnew=t(matrix(apply(g,1,function(x) 1/rowSums((x%*%t(1/x))^(1/(m-1)))),ncol=nobs));
    
    #get new objective function
    objfn=sum((umatnew^m)*g);
    

    
    #stop if membership mat doesn't change much
    #otherwise update membership mat and keep going
    chg=sum(abs(umatnew-umat))/sum(abs(umat));
    # print(c(objfn,chg));
    if(chg<.0001){
      stp=1; #stop loop
      
      #get pairwise distances (probably faster way to do this...)
      idx = expand.grid(1:nobs,1:nobs,1:nf);
      darr.env = array(mapply(function(x,y,z) sum((env_ind[[x]][z]-env_ind[[y]][z])^2),
                              x=idx[,1],y=idx[,2],z=idx[,3]),dim=c(nobs,nobs,nf));
      
      idx = expand.grid(1:nobs,1:nobs,1:nf,1:(ncat-1));
      darr.sca = array(mapply(function(w,x,y,z) sum((sca_ind[[w]][y,z]-sca_ind[[x]][y,z])^2),
                              w=idx[,1],x=idx[,2],y=idx[,3],z=idx[,4]),dim=c(nobs,nobs,nf,ncat-1));
      
      #get array of cluster membership
      clarr=array(sapply(classes,function(x) umatnew[,x]%*%t(umatnew[,x])),dim=c(nobs,nobs,k))
      
      #total cluster ss
      a.env.1 = apply(darr.env,3,function(x) (1/nobs)*sum(x)); #sum of all distances for each frequency
      a.sca.1 = apply(darr.sca,c(3,4),function(x) (1/nobs)*sum(x)); #sum of all distances for each frequency and categories (-1)
      
      #within cluster ss
      nclus=colSums(umatnew);
      idx=expand.grid(1:nf,classes);
      a.env.2 = matrix(mapply(function(x,y) sum(darr.env[,,x]*clarr[,,y])*(1/nclus[y]),
                              x=idx[,1],y=idx[,2]),nrow=nf,ncol=k);
      idx=expand.grid(1:nf,1:(ncat-1),classes);
      a.sca.2 = array(mapply(function(x,y,z) sum(darr.sca[,,x,y]*clarr[,,z]*(1/nclus[z])),
                             x=idx[,1],y=idx[,2],z=idx[,3]),dim=c(nf,ncat-1,k));
      
      #between cluster ss
      a.env=a.env.1-rowSums(a.env.2);
      a.sca=a.sca.1-apply(a.sca.2,c(1,2),sum);
      objcrit= sum(a.env)+sum(a.sca); #get objective criterion
      
    } else{
      umat=umatnew;
    }
  }
  return(list(umat=umat,env=env,sca=sca,objfn=objfn,objcrit=objcrit))
}

##################
#wrapper for fuzzy clustering algorithm to run multiple starts
##################
envsca_fuzzyclust_wrap <- function(rf, env_ind, sca_ind, k, m=1.2, nstart=100, seed=123){
  
  set.seed(seed)
  seedlist=sample(1:1000,size=nstart);
  
  #run multiple starts
  plan(multisession);
  tmp=future_lapply(1:nstart,function(x) envsca_fuzzyclust(rf,env_ind,sca_ind, k, m, seedlist[x]),future.seed=TRUE)
  
  #return start with min objective function
  objfn=unlist(lapply(tmp,function(x) x$objfn))
  return(tmp[[which(objfn==min(objfn))[1]]]);
}

##################
#gap statistic to choose fuzziness parameter m
##################
gapstat_m <- function(rf, env_ind, sca_ind, k, m=1.2, B=10, seed=123){
  
  set.seed(seed)

  nf=length(rf);
  ncat=dim(sca_ind[[1]])[2]+1;
  nobs=length(env_ind);
  nm=length(m);
  
  envp=array(unlist(tmp$env_ind),dim=c(nf,nobs,B));
  scap=array(unlist(tmp$sca_ind),dim=c(nf,ncat-1,nobs,B));
  
  envp2=apply(envp,c(1,3),function(x) x[sample(nobs)]);
  scap2=apply(scap,c(1,2,4),function(x) x[sample(nobs)]);
  
  envp3=lapply(seq(B), function(x) lapply(seq(nm), function(y) lapply(seq(nobs), function(z) envp2[z,,x])));
  scap3=lapply(seq(B), function(x) lapply(seq(nm), function(y) lapply(seq(nobs), function(z) scap2[z,,,x])));
  
  #run the cluster algo on each permuted set for each value of s and produce objective fn
  idx=expand.grid(1:B,m);
  seed1=sample(1:1000,size=nrow(idx));
  
  plan(multisession(workers=max(1,availableCores()-2),gc=TRUE));
  out1=future_mapply(function(x,y,z) envsca_fuzzyclust2(rf,envp3[[x]][[which(m==y)]],scap3[[x]][[which(m==y)]],
                                                      k,y,z)$objcrit ,x=idx[,1],y=idx[,2],z=seed1,future.seed=TRUE);
  
 
  #run cluster algo for original set for each value of m
  seed2=sample(1:1000,size=length(m));
  out2=future_mapply(function(y,z) envsca_fuzzyclust2(rf,env_ind,sca_ind,k,y,z)$objcrit,y=m,z=seed2,future.seed=TRUE); 
  
  #produce gap statistics
  gapst = numeric(length(m));
  for (i in 1:length(m)){
    gapst[i] = log(out2[i]) - mean(log(out1[idx[,2]==m[i]]));
  }
  
  # plot(x=m,y=gapst,type='b')
  
  #choose m = argmax(m) gapst 
  #or m as the smallest value within 1sd of log(objfn.out1) at m of the largest gap stat
  m.star = m[which(gapst==max(gapst))[1]];
  
  return(list(gapst=gapst,m.star=m.star))
}