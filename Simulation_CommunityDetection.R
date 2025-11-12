rm(list=ls())
library(irlba)
library(aricode)  #NMI
library(mclust)  #ARI
library(RSpectra)
source("HSC_function.R")


#Case 1
N<-400
K<-2
R<-2

p<-0.2
lambda<-0.5
a<-0.001
b<-1

m1<-1.2
m2<-1
sd<-1

Z<-matrix(0,N,K)

for (j in 1:K){
	Z[((j-1)*N/K+1):(j*N/K),j]<-1
}
B<-p*lambda*diag(K)+p*(1-lambda)*rep(1,K)%*%t(rep(1,K))
M<-(m1-m2)*diag(K)+m2*rep(1,K)%*%t(rep(1,K))
EX<-Z%*%M
EA<-Z%*%B%*%t(Z)

true_labels <- rep(1:K, each = N/K)
times<-2
result0<-replicate(times,simul(N,K,R,p,lambda,a,b,m1,m2,sd,Z,B,M,EX,EA,true_labels))
mean_result<-apply(result0,1,mean)
mean_result <- data.frame(t(mean_result))

MR_result<-data.frame(HSCwithCov=mean_result[[12]],
                      CASC=mean_result[[9]],
                      EdgeSC=mean_result[[3]],
                      MotifSC=mean_result[[4]],
                      CovC=mean_result[[1]])

NMI_result<-data.frame(HSCwithCov=mean_result[[27]],
                       CASC=mean_result[[24]],
                       EdgeSC=mean_result[[18]],
                       MotifSC=mean_result[[19]],
                       CovC=mean_result[[16]])

ARI_result<-data.frame(HSCwithCov=mean_result[[40]],
                       CASC=mean_result[[37]],
                       EdgeSC=mean_result[[31]],
                       MotifSC=mean_result[[32]],
                       CovC=mean_result[[29]])

#Case 2
N<-400
K<-2
R<-2

p<-0.18
lambda<-0.18
a<-0.001
b<-1

m1<-1.2
m2<-0.2
sd<-0.4

Z<-matrix(0,N,K)

for (j in 1:K){
	Z[((j-1)*N/K+1):(j*N/K),j]<-1
}
B<-p*lambda*diag(K)+p*(1-lambda)*rep(1,K)%*%t(rep(1,K))
M<-(m1-m2)*diag(K)+m2*rep(1,K)%*%t(rep(1,K))
EX<-Z%*%M
EA<-Z%*%B%*%t(Z)

true_labels <- rep(1:K, each = N/K)
times<-100

result0<-replicate(times,simul(N,K,R,p,lambda,a,b,m1,m2,sd,Z,B,M,EX,EA,true_labels))
mean_result<-apply(result0,1,mean)
mean_result


#Case 3
N<-600
K<-4
R<-4

p<-0.2
lambda<-0.4
a<-0.1
b<-1.1

m1<-1.5
m2<-0.5
sd<-0.6

R.seq <-seq(1,10,1)
mymat <- matrix(data=NA,nrow=length(R.seq),ncol=38)
i=1
for (R in R.seq){
Z<-matrix(0,N,K)

for (j in 1:K){
	Z[((j-1)*N/K+1):(j*N/K),j]<-1
}
B<-p*lambda*diag(K)+p*(1-lambda)*rep(1,K)%*%t(rep(1,K))
M<-(m1-m2)*diag(K)+m2*rep(1,K)%*%t(rep(1,K))
EX<-Z%*%M
EA<-Z%*%B%*%t(Z)

true_labels <- rep(1:K, each = N/K)
times<-100
result0<-replicate(times,simul(N,K,R,p,lambda,a,b,m1,m2,sd,Z,B,M,EX,EA,true_labels))
mean_result<-apply(result0,1,mean)
mymat[i,] <- c(R,mean_result)
print(mymat[i,])
i<-i+1
}