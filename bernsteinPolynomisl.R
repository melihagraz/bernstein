
rm(list = ls())
library(huge)
dt<-read.table("datar.txt")
nsim=1000
for (m in 1:nsim){
  
  load(file= paste("adj", m, ".rda", sep=""))
  adj.true
  
  datar<-dt[(20*m-19):(20*m),]
  datar<-as.matrix(datar)
  
  
  size.all<-dim(datar)
  n<-size.all[1]
  p<-size.all[2]
  b<-matrix(0,n,p)
  for (i1 in 1:p){
    y<-datar[,i1]
    a=rep(0,n)
    z=rep(0,n)
    x=rep(0,n)
    for (j in 1:n){
      x[j]=(j-1)/(n-1)
      for (i in 1:n){
        a[i]=choose(n-1, (i-1))*(x[j]^(i-1))*((1-x[j])^(n-i))*y[i]
        b[j,i1]<-sum(a)	
      }
    }	
  }
  sz<-matrix(0,n,p)
  for (i1 in 1:p){
    y<-datar[,i1]
    n=length(y)
    a=rep(0,n)
    for (j in 1:n){
      x[j]=(j-1)/(n-1)
      for (i in 1:n){
        a[i]=exp(-(n-1)*x[j])*((n-1)*x[j])^(i-1)/factorial(i-1)*y[i]
        sz[j,i1]<-sum(a)	
      }
    }
  }
  
  test.out1= huge(datar, method = "mb")
  r.lambda<-huge.select(test.out1,criterion="ric")$opt.lambda
  out1<-huge(datar, method = "mb",lambda=r.lambda)
  r.path<-out1$path
  
  
  test.out2= huge(b, method = "mb")
  b.lambda<-huge.select(test.out2,criterion="ric")$opt.lambda
  out2<-huge(b, method = "mb",lambda=b.lambda)
  b.path<-out2$path
  
  test.out3= huge(sz, method = "mb")
  sz.lambda<-huge.select(test.out3,criterion="ric")$opt.lambda
  out3<-huge(sz, method = "mb",lambda=sz.lambda)
  sz.path<-out3$path
  
  AUCr<-huge.roc(path=r.path,theta=adj.true)$AUC
  AUCb<-huge.roc(path=b.path,theta=adj.true)$AUC
  AUCsz<-huge.roc(path=sz.path,theta=adj.true)$AUC
  
  write(c(m,AUCr),"AUCr-1000gene.txt",ncol=2,append=T)
  write(c(m,AUCb),"AUCb-1000gene.txt",ncol=2,append=T)
  write(c(m,AUCsz),"AUCsz-1000gene.txt",ncol=2,append=T)
  
  F1r<-huge.roc(path=r.path,theta=adj.true)$F1
  F1b<-huge.roc(path=b.path,theta=adj.true)$F1
  F1sz<-huge.roc(path=sz.path,theta=adj.true)$F1
  
  write(c(m,F1r),"F1r-1000gene.txt",ncol=2,append=T)
  write(c(m,F1b),"F1b-1000gene.txt",ncol=2,append=T)
  write(c(m,F1sz),"F1sz-1000gene.txt",ncol=2,append=T)
  
  tpr<-huge.roc(path=r.path,theta=adj.true)$tp
  tpb<-huge.roc(path=b.path,theta=adj.true)$tp
  tpsz<-huge.roc(path=sz.path,theta=adj.true)$tp
  
  fpr<-huge.roc(path=r.path,theta=adj.true)$fp
  fpb<-huge.roc(path=b.path,theta=adj.true)$fp
  fpsz<-huge.roc(path=sz.path,theta=adj.true)$fp 
  
  write(c(m,tpr),"tpr-1000gene.txt",ncol=2,append=T)
  write(c(m,tpb),"tpb-1000gene.txt",ncol=2,append=T)
  write(c(m,tpsz),"tpsz-1000gene.txt",ncol=2,append=T)
  
  write(c(m,fpr),"fpr-1000gene.txt",ncol=2,append=T)
  write(c(m,fpb),"fpb-1000gene.txt",ncol=2,append=T)
  write(c(m,fpsz),"fpsz-1000gene.txt",ncol=2,append=T)
}