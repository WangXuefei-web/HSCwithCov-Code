rm(list=ls())
library(irlba)
library(aricode)  #NMI
library(mclust)  #ARI
library(ggplot2)
library(gridExtra)
library(grid)
library(showtext)

#--------------------------
specClust <- function(adjMat, ifcov=T, covMat, alpha, nBlocks, method = "regLaplacian",
                      rowNorm = T, nIter = 20, verbose = F) {
  
  # eigs does not support dsCMatrix type at this time
  # Matrix has Namespace problems when using dsCMatrix
  adjMat = as(adjMat, "dgCMatrix")
  
  similarityMat = getSimilarityMat(adjMat, method)
  
  if(ifcov==T){
    similarityMat = as(similarityMat + alpha*covMat%*%t(covMat), "dgCMatrix")
  }
  
  # eigsDecomp = eigs(similarityMat, nBlocks + 3)
  
  #nu??????????????��??��??nv??????????????��??��
  eigsDecomp = irlba(similarityMat, nu = nBlocks + 1, nv = 0, maxit=2000, work = 3*(nBlocks + 1)) 
  eigsDecomp = list(vectors = eigsDecomp$u,  values = eigsDecomp$d)
  
  if(rowNorm == T) {
    eigsDecomp$vectors[,1:nBlocks] = eigsDecomp$vectors[,1:nBlocks] /
      sqrt(rowSums(eigsDecomp$vectors[,1:nBlocks]^2))
    
    # if there were rows of zeros need to handle NaN's
    eigsDecomp$vectors[is.nan(eigsDecomp$vectors)] = 0
  }
  
  # kmeansResult = bigkmeans(eigsDecomp$vectors[,1:nBlocks], nBlocks,nstart = nIter)
  kmeansResult = kmeans(eigsDecomp$vectors[,1:nBlocks], nBlocks,nstart = nIter)
  kmeansResult$tot.withinss
  
  if(verbose == T) {
    return( list(cluster = kmeansResult$cluster,
                 wcss = kmeansResult$tot.withinss,
                 eigenVals = eigsDecomp$values[1:(nBlocks+1)],
                 eigenVecs = eigsDecomp$vectors[,1:(nBlocks+1)]) )
  } else {
    return(kmeansResult$cluster)
  }
}


getTuningRange = function(graphMatrix, covariates, nBlocks,
                          assortative) {
  
  nCov = ncol(covariates)
  
  singValCov = svd(covariates, nu = min(nBlocks, nCov))$d
  
  if(assortative == T) {
    # eigenValGraph = eigs(graphMatrix, nBlocks + 2, which = "LR",
    #     opts = list(retvec = F))$values #eigenvalues only
    
    eigenValGraph = irlba(graphMatrix, nu = nBlocks + 1, nv = 0, maxit=2000, work = 3*(nBlocks + 1))$d
    
    if(nCov > nBlocks) {
      hmax = eigenValGraph[1]/(singValCov[nBlocks]^2 - singValCov[nBlocks+1]^2) 
    } else {
      hmax = eigenValGraph[1]/singValCov[nCov]^2 
    }
    hmin = (eigenValGraph[nBlocks] - eigenValGraph[nBlocks + 1])/singValCov[1]^2
  } else {
    # eigenValGraph = eigs(graphMatrix, nBlocks + 2,
    #     opts = list(retvec = F))$values #eigenvalues only
    # eigenValGraph = sort(eigenValGraph^2, decreasing=T)
    
    eigenValGraph = (irlba(graphMatrix, nu = nBlocks + 1, nv = 0, maxit=2000, work = 3*(nBlocks + 1))$d)^2
    
    if(nCov > nBlocks) {
      hmax = eigenValGraph[1]/(singValCov[nBlocks]^2 - singValCov[nBlocks+1]^2) 
    } else {
      hmax = eigenValGraph[1]/singValCov[nCov]^2 
    }
    hmin = (eigenValGraph[nBlocks] - eigenValGraph[nBlocks + 1])/singValCov[1]^2
  }
  
  # return
  list( hmax = hmax, hmin = hmin )
}

getSimilarityMat <- function(adjMat, method) {
  if(method == "regLaplacian") {
    rSums = Matrix::rowSums(adjMat)
    tau = mean(rSums)
    normMat = Diagonal(length(rSums), 1/sqrt(rSums + tau))
    return(normMat %*% adjMat %*% normMat)
  }
  else if(method == "laplacian") {
    rSums = Matrix::rowSums(adjMat)
    normMat = Diagonal(length(rSums), 1/sqrt(rSums))
    return(normMat %*% adjMat %*% normMat)
  }
  else if (method == "CASCregLaplacian"){
    rSums = Matrix::rowSums(adjMat)
    tau = mean(rSums)
    normMat = Diagonal(length(rSums), 1/sqrt(rSums + tau))
    return(normMat %*% adjMat %*% normMat %*% normMat %*% adjMat %*% normMat)
  }
  else if(method == "adjacency"){
    return(adjMat)
  }
  else {
    stop(paste("Error: method =", method, "Not valid"))
  }
}

misClustRate <- function(clusters, nMembers) {
  
  nNodes = sum(nMembers)
  nBlocks = length(nMembers)
  nMisClustNodes = 0
  uniqueClusters = unique(clusters)
  clusterCounts = matrix(0, ncol = nBlocks, nrow = nBlocks)
  clusterLabels = rep(0, nBlocks)
  
  #get label counts for each cluster
  for(i in 1:nBlocks) {
    
    clustStart = sum(nMembers[1:i]) - sum(nMembers[i]) + 1
    clustEnd = sum(nMembers[1:i])
    
    for(j in uniqueClusters) {
      clusterCounts[j,i] = sum(j == clusters[clustStart:clustEnd]) 
    }    
  }
  
  #determine cluster label based on counts
  clusterCountsTemp = clusterCounts
  for(i in 1:nBlocks) {
    
    maxCoor = t(which(clusterCountsTemp ==
                        max(clusterCountsTemp), arr.ind = T))
    clusterLabels[maxCoor[2]] = maxCoor[1]
    clusterCountsTemp[maxCoor[1], ] = rep(-1, nBlocks)
    clusterCountsTemp[, maxCoor[2]] = rep(-1, nBlocks)
    
  }
  
  for(i in 1:nBlocks) {
    nMisClustNodes = nMisClustNodes + sum(clusterCounts[-clusterLabels[i],i])
  }  
  
  return( nMisClustNodes/nNodes )
}

tuningOnlyalpha = function(graphMat, covMat, alpha_n=20, method = "adjacency") {
  true_labels <- rep(1:K, each = N/K)
  rangehTuning <- getTuningRange(graphMat, covMat, nBlocks=K, assortative=T)
  hTuningSeq = seq(rangehTuning$hmin, rangehTuning$hmax,length.out = alpha_n)
  
  a_wcssVec = vector(length = alpha_n)
  nmi_value<- vector(length = alpha_n)
  for(i in 1:alpha_n) {
    res_cluster<- specClust(graphMat, alpha=hTuningSeq[i], ifcov=T, covMat, nBlocks=K, method,rowNorm = T, nIter = 20, verbose = T)
    a_wcssVec[i] = res_cluster$wcss
    nmi_value[i] <- NMI(true_labels, res_cluster$cluster)
  }
  hOpt = hTuningSeq[which.min(a_wcssVec)]
  return(list(bestalpha=hOpt,
              hTuningSeq=hTuningSeq,
              a_wcssVec=a_wcssVec,
              nmi_value=nmi_value))
}

tuningPara<-function(delta_n=20,alpha_n=20,verbose=F,adjMat, wadjMat, covMat,method = "adjacency"){
  true_labels <- rep(1:K, each = N/K)
  deltaseq<- seq(0,1,length.out = delta_n)
  d_wcssVec = vector(length = delta_n)
  d_gapValK = vector(length = delta_n)
  d_gapValKplus1 = vector(length = delta_n)
  error_rate<- vector(length = delta_n)
  bestalphaseq<- vector(length = delta_n)
  nmi_value<- vector(length = delta_n)
  Valgaps<- vector(length = delta_n)
  for (j in 1:delta_n){
    delta=deltaseq[j]
    madjMat=delta*adjMat+(1-delta)*wadjMat
    
    alpha1=tuningOnlyalpha(graphMat=madjMat, covMat, alpha_n, method)
    bestalphaseq[j] = alpha1$bestalpha
    
    res_cluster<- specClust(madjMat, alpha=bestalphaseq[j], ifcov=T, covMat, nBlocks=K, method,rowNorm = T, nIter = 20, verbose = T)
    
    d_wcssVec[j] = res_cluster$wcss
    d_gapValK[j] = res_cluster$eigenVals[K]
    d_gapValKplus1[j] = res_cluster$eigenVals[K+1]
    Valgaps[j]<-(d_gapValK[j]-d_gapValKplus1[j])/d_gapValK[j]
    error_rate[j]<-misClustRate(res_cluster$cluster, rep(N/K,K))
    nmi_value[j] <- NMI(true_labels, res_cluster$cluster)
  }
  wcssgaps<-abs((d_wcssVec[1:(delta_n-1)]-d_wcssVec[2:delta_n])/d_wcssVec[1:(delta_n-1)])
  #return
  if (verbose==F){
    return(list(bestAlpha=bestalphaseq[which.min(d_wcssVec)],  #wcssgaps
                bestDelta=deltaseq[which.min(d_wcssVec)]))
  }else{
    return(list(bestAlpha=bestalphaseq[which.min(d_wcssVec)],
                bestDelta=deltaseq[which.min(d_wcssVec)],
                deltaseq=deltaseq,
                bestalphaseq=bestalphaseq,
                wcssgaps=wcssgaps,
                d_wcssVec=d_wcssVec,
                d_gapValK=d_gapValK,
                d_gapValKplus1=d_gapValKplus1,
                Valgaps=Valgaps,
                error_rate=error_rate,
                nmi_value=nmi_value))
  }
}

minmax <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

#--------------------------
set.seed(1234)

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
covMat<-matrix(sapply(as.vector(EX), function(mu) rnorm(1, mean = mu, sd=sd)),N,R)
#covMat<-apply(covMat,2,minmax)
EA<-Z%*%B%*%t(Z)
adjMat<-matrix(0,N,N)
adjMat[upper.tri(adjMat)]<-runif(N*(N-1)/2,0.001,1)*rbinom(N*(N-1)/2,1,prob=as.vector(EA[upper.tri(EA)]))
adjMat[lower.tri(adjMat)]<-t(adjMat)[lower.tri(t(adjMat))]

##-------motif matrix--------
AMat<-adjMat
AMat[AMat!=0]<-1
wadjMat<-adjMat%*%AMat*AMat+AMat%*%adjMat*AMat+AMat%*%AMat*adjMat
wadjMat[wadjMat!=0]<-wadjMat[wadjMat!=0]/3

delta_n=15
alpha_n=15
res_Para<-tuningPara(delta_n,alpha_n,verbose=T,adjMat, wadjMat, covMat,method = "adjacency")
delta<-res_Para$bestDelta
delta<-res_Para$deltaseq[which.min(res_Para$d_wcssVec)]
res_alpha<-tuningOnlyalpha(graphMat=delta*adjMat+(1-delta)*wadjMat, covMat, alpha_n, method = "adjacency")

deltaseq=res_Para$deltaseq
bestalphaseq=res_Para$bestalphaseq
d_wcss<-res_Para$d_wcssVec

alpahseq=res_alpha$hTuningSeq
a_wcss<-res_alpha$a_wcssVec

# dataframe
data1 <- data.frame(
  delta = deltaseq,
  normalized_wcss = (d_wcss - min(d_wcss)) / (max(d_wcss) - min(d_wcss)),
  nmi_value = res_Para$nmi_value
)

data2 <- data.frame(
  alpha = alpahseq,
  normalized_wcss = (a_wcss - min(a_wcss)) / (max(a_wcss) - min(a_wcss)),
  nmi_value = res_alpha$nmi_value
)

showtext_auto()

# ͼ1: Delta NMI
plot1<-ggplot(data1, aes(x = delta)) +
	geom_line(aes(y = normalized_wcss, color = "WCSS"), linewidth = 0.6) +
	geom_point(aes(y = normalized_wcss, color = "WCSS"), , size = 1) +
	
	geom_line(aes(y = nmi_value, color = "NMI"), linewidth = 0.6) +
	geom_point(aes(y = nmi_value, color = "NMI"), size = 1) +
	scale_color_manual(values = c("WCSS" = "black", "NMI" = "red")) +
	labs(title=expression(paste("with ", zeta, " selected to minimize WCSS")),
			 y = "",
			 color = "") +
	theme_minimal(base_size = 15)+
  theme(
    plot.title = element_text(hjust = 0.5, size = 10.5, face = "bold"), 
    axis.text = element_text(size = 10.5),      
    axis.title = element_text(size = 12),       
    legend.position = "bottom", 
    plot.margin = margin(0.5, 0.5, 0, 0, "cm"),
    legend.margin = margin(t = -10),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)  # ???Ӻ?ɫ?߿?
  )

delta_wcssmin<-round(res_Para$deltaseq[which.min(res_Para$d_wcssVec)],2)
# ͼ2: Alpha NMI
plot2 <- ggplot(data2, aes(x = alpha)) +
	geom_line(aes(y = normalized_wcss, color = "WCSS"), linewidth = 0.6) +
	geom_point(aes(y = normalized_wcss, color = "WCSS"), , size = 1) +
	
	geom_line(aes(y = nmi_value, color = "NMI"), linewidth = 0.6) +
	geom_point(aes(y = nmi_value, color = "NMI"), size = 1) +
	scale_color_manual(values = c("WCSS" = "black", "NMI" = "red")) +
	labs(title = bquote(delta == .(delta_wcssmin)),
			 x = expression(zeta),
			 y = "",
			 color = "") +
	theme_minimal(base_size = 15)+
  theme(
    plot.title = element_text(hjust = 0.5, size = 10.5, face = "bold"), 
    axis.text = element_text(size = 10.5),      
    axis.title = element_text(size = 12),       
    legend.position = "bottom", 
    plot.margin = margin(0.5, 0.5, 0, 0, "cm"),
    legend.margin = margin(t = -10),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)  
  )

# arrange
grid.arrange(plot1, plot2, ncol = 2)

res_Para$bestAlpha
res_Para$bestDelta

res_Para$deltaseq[which.min(res_Para$d_wcssVec)]
res_Para$bestalphaseq[which.min(res_Para$d_wcssVec)]
