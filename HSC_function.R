#rm(list=ls())
library(irlba)
library(aricode)  #NMI
library(mclust)  #ARI
library(RSpectra)


NMI<-aricode::NMI

misClustRate2 <- function(clusters, nMembers) {
	
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
												max(clusterCountsTemp), arr.ind = T)) #arr.ind = TRUE
		clusterLabels[maxCoor[2]] = maxCoor[1]
		clusterCountsTemp[maxCoor[1], ] = rep(-1, nBlocks)
		clusterCountsTemp[, maxCoor[2]] = rep(-1, nBlocks)
		
	}
	
	for(i in 1:nBlocks) {
		nMisClustNodes = nMisClustNodes + sum(clusterCounts[-clusterLabels[i],i])
	}  
	
	return( list(misRate=nMisClustNodes/nNodes,
							 clusterLabels = clusterLabels))
}


map_labels <- function(est_labels, clusterLabels) {
	correct_est_label <- sapply(est_labels, function(x) which(clusterLabels == x))
  return(correct_est_label)
}

getSimilarityMat <- function(adjMat, method) {
	adjMat = as(adjMat, "CsparseMatrix")
	
	if(method == "regLaplacian") {
		rSums = Matrix::rowSums(adjMat)
		tau = mean(rSums)
		normMat = Diagonal(length(rSums), 1/sqrt(rSums + tau))   
		return(normMat %*% adjMat %*% normMat)
	}
	else if(method == "laplacian") {
		rSums = Matrix::rowSums(adjMat)
		normMat = Diagonal(length(rSums), 1/sqrt(rSums + tau))
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


getTuningRange<-function(adjMat,covMat,K,method){
	graphMatrix = getSimilarityMat(adjMat, method)
	nCov = ncol(covMat)
	
	covMat2<-covMat%*%t(covMat)

	cov_eig = svd(covMat2, nu = min(K, nCov))$d
#	net_eig = svd(graphMatrix, K+1)$d
	#net_eig = RSpectra::svds(graphMatrix, K+1)$d
	net_eig = irlba(graphMatrix, nu = K + 1, nv = K + 1, maxit=2000, work = 3*(K + 1))$d
	
	if(nCov > K) {
		hmax = net_eig[1]/(cov_eig[K] - cov_eig[K+1]) 
	} else {
		hmax = net_eig[1]/cov_eig[nCov]
	}
	hmin = (net_eig[K] - net_eig[K + 1])/cov_eig[1]
	return(list(alpha_min=hmin, alpha_max=hmax))
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


tuningPara<-function(delta_n=20,alpha_n=20,verbose=F,adjMat,wadjMat,covMat,method = "adjacency",true_labels,misclust_labels){
	deltaseq<- seq(0,1,length.out = delta_n)
	hOpt<- vector(length = delta_n)
	
	d_wcssVec = vector(length = delta_n)
	d_gapValK = vector(length = delta_n)
	
	error_rate<- vector(length = delta_n)
	nmi_value<- vector(length = delta_n)
	Valgaps<- vector(length = delta_n)
	
	doubleCov=covMat%*%t(covMat)
	for (j in 1:delta_n){
		delta=deltaseq[j]
		madjMat=delta*adjMat+(1-delta)*wadjMat
		madjMat<-getSimilarityMat(madjMat, method)
		rangehTuning <- getTuningRange(madjMat, covMat, K=K, method=method)
		hTuningSeq = seq(rangehTuning$alpha_min, rangehTuning$alpha_max,length.out = alpha_n)
		
		a_wcssVec = vector(length = alpha_n)
		for(i in 1:alpha_n) {
			alpha=hTuningSeq[i]
			targetMat<-madjMat+alpha*doubleCov
			res_cluster<- specClust(targetMat,K,nIter = 20, verbose = T)
			a_wcssVec[i] = res_cluster$wcss
		}
		hOpt[j] = hTuningSeq[which.min(a_wcssVec)]
		
		targetMat = madjMat+hOpt[j]*doubleCov
		res_cluster<- specClust(targetMat,K,nIter = 20, verbose = T)
		
		d_wcssVec[j] = res_cluster$wcss
		d_gapValK[j] = res_cluster$eigenVals[K]
		Valgaps[j]<-(res_cluster$eigenVals[K]-res_cluster$eigenVals[K+1])/res_cluster$eigenVals[K]
		
		error_rate[j]<-misClustRate(res_cluster$cluster, misclust_labels)
		nmi_value[j] <- NMI(true_labels, res_cluster$cluster)
	}
	
	wcssgaps<-abs((d_wcssVec[1:(delta_n-1)]-d_wcssVec[2:delta_n])/d_wcssVec[1:(delta_n-1)])
	#return
	if (verbose==F){
		return(list(bestAlpha=hOpt[which.min(wcssgaps)],  #wcssgaps
								bestDelta=deltaseq[which.min(wcssgaps)]))
	}else{
		return(list(bestAlpha=hOpt[which.min(wcssgaps)],
								bestDelta=deltaseq[which.min(wcssgaps)],
								wcssgaps=wcssgaps,
								alphaOPt=hOpt,
								d_wcssVec=d_wcssVec,
								d_gapValK=d_gapValK,
								Valgaps=Valgaps,
								error_rate=error_rate,
								nmi_value=nmi_value))
	}
}

tuningOnlyalpha = function(graphMat, covMat, alpha_n=20, method = "adjacency") {
	graphMat=getSimilarityMat(graphMat, method)
	
	rangehTuning <- getTuningRange(graphMat, covMat, K=K, method=method)
	hTuningSeq = seq(rangehTuning$alpha_min, rangehTuning$alpha_max,length.out = alpha_n)
	a_wcssVec = vector(length = alpha_n)
	
	doubleCov=covMat%*%t(covMat)
	for(i in 1:alpha_n) {
		targetMat<-graphMat+hTuningSeq[i]*doubleCov
		res_cluster<- specClust(targetMat,K,nIter = 20, verbose = T)
		a_wcssVec[i] = res_cluster$wcss
	}
	hOpt = hTuningSeq[which.min(a_wcssVec)]
	return(hOpt)
}


tuningOnlydelta = function(adjMat, wadjMat, delta_n=20, method = "adjacency",verbose=F,misclust_labels) {
	deltaseq<- seq(0,1,length.out = delta_n)
	d_wcssVec = vector(length = delta_n)
	eigValK = vector(length = delta_n)
	error_rate<- vector(length = delta_n)
	for(i in 1:delta_n) {
		delta=deltaseq[i]
		madjMat=delta*adjMat+(1-delta)*wadjMat
		targetMat<-getSimilarityMat(madjMat, method)
		res_cluster <- specClust(targetMat,K,nIter = 20, verbose = T)
		d_wcssVec[i] = res_cluster$wcss
		eigValK[i]<-res_cluster$eigenVals[K]
		error_rate[i]<-misClustRate(res_cluster$cluster, misclust_labels)
	}
	wcssgaps<-abs((d_wcssVec[1:(delta_n-1)]-d_wcssVec[2:delta_n])/d_wcssVec[1:(delta_n-1)])
	bestdelta = deltaseq[which.min(wcssgaps)]
	if (verbose==T){
		return(list(bestdelta=bestdelta, d_wcssVec=d_wcssVec, wcssgaps=wcssgaps,
								eigValK=eigValK, error_rate=error_rate))
	}else{
		return(bestdelta)
	}
}

#min_max
minmax <- function(x) {
	return((x - min(x)) / (max(x) - min(x)))
}

specClust<-function(targetMat,K,nIter = 20, verbose = F){
	targetMat = as(targetMat, "CsparseMatrix")
	
	#	eigsDecomp = svd(targetMat, K+1)   
	eigsDecomp = irlba(targetMat, nu = K + 1, nv = K + 1, maxit=2000, work = 3*(K + 1)) 
	
	eig_vector<-eigsDecomp$u[,1:K]   
	eig_value<-eigsDecomp$d[1:(K+1)]
	

	eig_vector = eig_vector / sqrt(rowSums(eig_vector^2))
	# if there were rows of zeros need to handle NaN's
	eig_vector[is.na(eig_vector) | is.nan(eig_vector) | is.infinite(eig_vector)] <- 0
	
	kmeansResult = kmeans(eig_vector, K, nstart = nIter)
	
	if(verbose == T) {
		return( list(cluster = kmeansResult$cluster,
								 wcss = kmeansResult$tot.withinss,
								 eigenVals = eig_value))
	} else {
		return(kmeansResult$cluster)
	}
}


simul<-function(N,K,R,p,lambda,a,b,m1,m2,sd,Z,B,M,EX,EA,true_labels){
	covMat<-matrix(sapply(as.vector(EX), function(mu) rnorm(1, mean = mu, sd=sd)),N,R)
	adjMat<-matrix(0,N,N)
	adjMat[upper.tri(adjMat)]<-runif(N*(N-1)/2,a,b)*rbinom(N*(N-1)/2,1,prob=as.vector(EA[upper.tri(EA)]))
	adjMat[lower.tri(adjMat)]<-t(adjMat)[lower.tri(t(adjMat))]
	
	##motif matrix
	AMat<-adjMat
	AMat[AMat!=0]<-1
	wadjMat<-adjMat%*%AMat*AMat+AMat%*%adjMat*AMat+AMat%*%AMat*adjMat
	#weight<-(AMat%*%AMat*AMat)*3
	wadjMat[wadjMat!=0]<-wadjMat[wadjMat!=0]/3
	
	doubleCov<-covMat%*%t(covMat)
	
	true_labels<-rep(1:K, each = N/K)
	misclust_labels=rep(N/K,K)

	##XX,kmeans
	XX_kmeans1 = scale(covMat, center = F, scale = sqrt(Matrix::colSums(covMat^2))) 
	XX_kmeans <- kmeans(XX_kmeans1, centers = K,nstart = 20)$cluster
	
	##XX 
	XX<-specClust(targetMat=doubleCov,K,nIter = 20, verbose = F)
	
	##edge, A
	We <- specClust(targetMat=adjMat,K,nIter = 20, verbose = F)
	
	##motif, A
	Wm <- specClust(targetMat=wadjMat,K,nIter = 20, verbose = F)
	
	##mixed, A 
	delta_m<-tuningOnlydelta(adjMat, wadjMat, delta_n=20, method = "adjacency",verbose=F,misclust_labels)
	Mixed_W <- specClust(targetMat = delta_m*adjMat+(1-delta_m)*wadjMat,K,nIter = 20, verbose = F)
	
	##edge, A + alpha * X %*% Matrix::t(X) 
	alpha_e<-tuningOnlyalpha(graphMat=adjMat, covMat, alpha_n=20, method = "adjacency")
	We_X <- specClust(targetMat = adjMat+alpha_e*doubleCov, K,nIter = 20, verbose = F)
	
	##motif, A + alpha * X %*% Matrix::t(X)
	alpha_m<-tuningOnlyalpha(graphMat=wadjMat, covMat, alpha_n=20, method = "adjacency")
	Wm_X <- specClust(targetMat = wadjMat+alpha_m*doubleCov, K,nIter = 20, verbose = F)
	
	##mixed, LA + alpha * X %*% Matrix::t(X)
	##L mixed matrix
	LadjMat<-getSimilarityMat(adjMat, method = "regLaplacian")
	LwadjMat<-getSimilarityMat(wadjMat, method = "regLaplacian")
	LbeyondPara<-tuningPara(delta_n=20,alpha_n=20,verbose=F,LadjMat,LwadjMat,covMat,method = "adjacency",true_labels,misclust_labels)
	Ldelta=LbeyondPara$bestDelta
	Lalpha=LbeyondPara$bestAlpha
	
	LmadjMat<-Ldelta*LadjMat+(1-Ldelta)*LwadjMat
	Mixed_LW_X <- specClust(targetMat = LmadjMat+Lalpha*doubleCov, K,nIter = 20, verbose = F)
	
	##casc_We_X
	alpha_casc_We<-tuningOnlyalpha(graphMat=adjMat, covMat, alpha_n=20, method = "CASCregLaplacian")
	casc_adjMat<-getSimilarityMat(adjMat, method = "CASCregLaplacian")
	casc_We_X<- specClust(targetMat = casc_adjMat+alpha_casc_We*doubleCov, K,nIter = 20, verbose = F)
	
	##casc_Wm_X
	alpha_casc_Wm<-tuningOnlyalpha(graphMat=wadjMat, covMat, alpha_n=20, method = "CASCregLaplacian")
	casc_wadjMat<-getSimilarityMat(wadjMat, method = "CASCregLaplacian")
	casc_Wm_X<- specClust(targetMat = casc_wadjMat+alpha_casc_Wm*doubleCov, K,nIter = 20, verbose = F)
	
	##casc_W_X
	cascPara<-tuningPara(delta_n=20,alpha_n=20,verbose=F,adjMat,wadjMat,covMat,method = "CASCregLaplacian",true_labels,misclust_labels)
	casc_d=cascPara$bestDelta
	casc_a=cascPara$bestAlpha
	casc_madjMat<-getSimilarityMat(casc_d*adjMat+(1-casc_d)*wadjMat, method = "CASCregLaplacian")
	
	casc_W_X<- specClust(targetMat = casc_madjMat+casc_a*doubleCov, K,nIter = 20, verbose = F)
	
	##mixed, A + alpha * X %*% Matrix::t(X)
	##mixed matrix
	beyondPara<-tuningPara(delta_n=20,alpha_n=20,verbose=F,adjMat,wadjMat,covMat,method = "adjacency",true_labels,misclust_labels)
	delta=beyondPara$bestDelta
	alpha=beyondPara$bestAlpha
	madjMat<-delta*adjMat+(1-delta)*wadjMat
	
	Mixed_W_X <- specClust(targetMat = madjMat+alpha*doubleCov, K,nIter = 20, verbose = F)
	
	##error_rate
	XX_kmeans_error<-misClustRate(XX_kmeans, misclust_labels)
	XX_error<-misClustRate(XX, misclust_labels)
	We_error<-misClustRate(We, misclust_labels)
	Wm_error<-misClustRate(Wm, misclust_labels)
	Mixed_W_error<-misClustRate(Mixed_W, misclust_labels)
	
	We_X_error<-misClustRate(We_X, misclust_labels)
	Wm_X_error<-misClustRate(Wm_X, misclust_labels)
	Mixed_LW_X_error<-misClustRate(Mixed_LW_X, misclust_labels)
	casc_We_X_error<-misClustRate(casc_We_X, misclust_labels)
	casc_Wm_X_error<-misClustRate(casc_Wm_X, misclust_labels)
	casc_W_X_error<-misClustRate(casc_W_X, misclust_labels)
	Mixed_W_X_error<-misClustRate(Mixed_W_X, misclust_labels)
	
	##nmi
	XX_kmeans_nmi<- NMI(true_labels, XX_kmeans)
	XX_nmi<-NMI(true_labels,XX)
	We_nmi<-NMI(true_labels,We)
	Wm_nmi<-NMI(true_labels,Wm)
	Mixed_W_nmi<-NMI(true_labels,Mixed_W)
	
	We_X_nmi<-NMI(true_labels,We_X)
	Wm_X_nmi<-NMI(true_labels,Wm_X)
	Mixed_LW_X_nmi<-NMI(true_labels,Mixed_LW_X)
	casc_We_X_nmi<-NMI(true_labels,casc_We_X)
	casc_Wm_X_nmi<-NMI(true_labels,casc_Wm_X)
	casc_W_X_nmi<-NMI(true_labels,casc_W_X)
	Mixed_W_X_nmi<-NMI(true_labels,Mixed_W_X)
	
	##ari
	XX_kmeans_ari<- adjustedRandIndex(true_labels, XX_kmeans)
	XX_ari<-adjustedRandIndex(true_labels,XX)
	We_ari<-adjustedRandIndex(true_labels,We)
	Wm_ari<-adjustedRandIndex(true_labels,Wm)
	Mixed_W_ari<-adjustedRandIndex(true_labels,Mixed_W)
	
	We_X_ari<-adjustedRandIndex(true_labels,We_X)
	Wm_X_ari<-adjustedRandIndex(true_labels,Wm_X)
	Mixed_LW_X_ari<-adjustedRandIndex(true_labels,Mixed_LW_X)
	casc_We_X_ari<-adjustedRandIndex(true_labels,casc_We_X)
	casc_Wm_X_ari<-adjustedRandIndex(true_labels,casc_Wm_X)
	casc_W_X_ari<-adjustedRandIndex(true_labels,casc_W_X)
	Mixed_W_X_ari<-adjustedRandIndex(true_labels,Mixed_W_X)
	
	return(c(XX_kmeans_error,XX_error, We_error, Wm_error, Mixed_W_error, 
	         We_X_error, Wm_X_error,Mixed_LW_X_error, casc_We_X_error, casc_Wm_X_error,
	         casc_W_X_error , Mixed_W_X_error,delta,alpha,100,
					 XX_kmeans_nmi,XX_nmi, We_nmi, Wm_nmi, Mixed_W_nmi,
					 We_X_nmi, Wm_X_nmi,Mixed_LW_X_nmi, casc_We_X_nmi,casc_Wm_X_nmi,
					 casc_W_X_nmi,  Mixed_W_X_nmi,100,XX_kmeans_ari,XX_ari, 
					 We_ari, Wm_ari, Mixed_W_ari, We_X_ari, Wm_X_ari,
					 Mixed_LW_X_ari, casc_We_X_ari, casc_Wm_X_ari, casc_W_X_ari,  Mixed_W_X_ari))
}

realnet<-function(adjMat,covMat){
	doubleCov<-covMat%*%t(covMat)

	##XX,kmeans
	XX_kmeans1 = scale(covMat, center = F, scale = sqrt(Matrix::colSums(covMat^2))) 
	XX_kmeans <- kmeans(XX_kmeans1, centers = K,nstart = 20)$cluster
	
	##XX 
	XX<-specClust(targetMat=doubleCov,K,nIter = 20, verbose = F)
	
	##edge, A
	We <- specClust(targetMat=adjMat,K,nIter = 20, verbose = F)
	
	##motif, A
	Wm <- specClust(targetMat=wadjMat,K,nIter = 20, verbose = F)
	
	##mixed, A 
	delta_m<-tuningOnlydelta(adjMat, wadjMat, delta_n=20, method = "adjacency",verbose=F,misclust_labels)
	Mixed_W <- specClust(targetMat = delta_m*adjMat+(1-delta_m)*wadjMat,K,nIter = 20, verbose = F)
	
	##edge, A + alpha * X %*% Matrix::t(X) 
	alpha_e<-tuningOnlyalpha(graphMat=adjMat, covMat, alpha_n=20, method = "adjacency")
	We_X <- specClust(targetMat = adjMat+alpha_e*doubleCov, K,nIter = 20, verbose = F)
	
	##motif, A + alpha * X %*% Matrix::t(X)
	alpha_m<-tuningOnlyalpha(graphMat=wadjMat, covMat, alpha_n=20, method = "adjacency")
	Wm_X <- specClust(targetMat = wadjMat+alpha_m*doubleCov, K,nIter = 20, verbose = F)
	
	##mixed, LA + alpha * X %*% Matrix::t(X)
	##L mixed matrix
	LadjMat<-getSimilarityMat(adjMat, method = "regLaplacian")
	LwadjMat<-getSimilarityMat(wadjMat, method = "regLaplacian")
	LbeyondPara<-tuningPara(delta_n=20,alpha_n=20,verbose=F,LadjMat,LwadjMat,covMat,method = "adjacency",true_labels,misclust_labels)
	Ldelta=LbeyondPara$bestDelta
	Lalpha=LbeyondPara$bestAlpha
	
	LmadjMat<-Ldelta*LadjMat+(1-Ldelta)*LwadjMat
	Mixed_LW_X <- specClust(targetMat = LmadjMat+Lalpha*doubleCov, K,nIter = 20, verbose = F)
	
	##casc_We_X
	alpha_casc_We<-tuningOnlyalpha(graphMat=adjMat, covMat, alpha_n=20, method = "CASCregLaplacian")
	casc_adjMat<-getSimilarityMat(adjMat, method = "CASCregLaplacian")
	casc_We_X<- specClust(targetMat = casc_adjMat+alpha_casc_We*doubleCov, K,nIter = 20, verbose = F)
	
	##casc_Wm_X
	alpha_casc_Wm<-tuningOnlyalpha(graphMat=wadjMat, covMat, alpha_n=20, method = "CASCregLaplacian")
	casc_wadjMat<-getSimilarityMat(wadjMat, method = "CASCregLaplacian")
	casc_Wm_X<- specClust(targetMat = casc_wadjMat+alpha_casc_Wm*doubleCov, K,nIter = 20, verbose = F)
	
	##casc_W_X
	cascPara<-tuningPara(delta_n=20,alpha_n=20,verbose=F,adjMat,wadjMat,covMat,method = "CASCregLaplacian",true_labels,misclust_labels)
	casc_d=cascPara$bestDelta
	casc_a=cascPara$bestAlpha
	casc_madjMat<-getSimilarityMat(casc_d*adjMat+(1-casc_d)*wadjMat, method = "CASCregLaplacian")
	
	casc_W_X<- specClust(targetMat = casc_madjMat+casc_a*doubleCov, K,nIter = 20, verbose = F)
	
	##mixed, A + alpha * X %*% Matrix::t(X)
	##mixed matrix
	beyondPara<-tuningPara(delta_n=20,alpha_n=20,verbose=F,adjMat,wadjMat,covMat,method = "adjacency",true_labels,misclust_labels)
	delta=beyondPara$bestDelta
	alpha=beyondPara$bestAlpha
	madjMat<-delta*adjMat+(1-delta)*wadjMat
	
	Mixed_W_X <- specClust(targetMat = madjMat+alpha*doubleCov, K,nIter = 20, verbose = F)
	
	##error_rate
	XX_kmeans_error<-misClustRate(XX_kmeans, misclust_labels)
	XX_error<-misClustRate(XX, misclust_labels)
	We_error<-misClustRate(We, misclust_labels)
	Wm_error<-misClustRate(Wm, misclust_labels)
	Mixed_W_error<-misClustRate(Mixed_W, misclust_labels)
	
	We_X_error<-misClustRate(We_X, misclust_labels)
	Wm_X_error<-misClustRate(Wm_X, misclust_labels)
	Mixed_LW_X_error<-misClustRate(Mixed_LW_X, misclust_labels)
	casc_We_X_error<-misClustRate(casc_We_X, misclust_labels)
	casc_Wm_X_error<-misClustRate(casc_Wm_X, misclust_labels)
	casc_W_X_error<-misClustRate(casc_W_X, misclust_labels)
	Mixed_W_X_error<-misClustRate(Mixed_W_X, misclust_labels)
	
	##nmi
	XX_kmeans_nmi<- NMI(true_labels, XX_kmeans)
	XX_nmi<-NMI(true_labels,XX)
	We_nmi<-NMI(true_labels,We)
	Wm_nmi<-NMI(true_labels,Wm)
	Mixed_W_nmi<-NMI(true_labels,Mixed_W)
	
	We_X_nmi<-NMI(true_labels,We_X)
	Wm_X_nmi<-NMI(true_labels,Wm_X)
	Mixed_LW_X_nmi<-NMI(true_labels,Mixed_LW_X)
	casc_We_X_nmi<-NMI(true_labels,casc_We_X)
	casc_Wm_X_nmi<-NMI(true_labels,casc_Wm_X)
	casc_W_X_nmi<-NMI(true_labels,casc_W_X)
	Mixed_W_X_nmi<-NMI(true_labels,Mixed_W_X)
	
	##ari
	XX_kmeans_ari<- adjustedRandIndex(true_labels, XX_kmeans)
	XX_ari<-adjustedRandIndex(true_labels,XX)
	We_ari<-adjustedRandIndex(true_labels,We)
	Wm_ari<-adjustedRandIndex(true_labels,Wm)
	Mixed_W_ari<-adjustedRandIndex(true_labels,Mixed_W)
	
	We_X_ari<-adjustedRandIndex(true_labels,We_X)
	Wm_X_ari<-adjustedRandIndex(true_labels,Wm_X)
	Mixed_LW_X_ari<-adjustedRandIndex(true_labels,Mixed_LW_X)
	casc_We_X_ari<-adjustedRandIndex(true_labels,casc_We_X)
	casc_Wm_X_ari<-adjustedRandIndex(true_labels,casc_Wm_X)
	casc_W_X_ari<-adjustedRandIndex(true_labels,casc_W_X)
	Mixed_W_X_ari<-adjustedRandIndex(true_labels,Mixed_W_X)
	
	return(list(value=c(XX_kmeans_error,XX_error, We_error, Wm_error, Mixed_W_error,
	                    We_X_error, Wm_X_error,Mixed_LW_X_error, casc_We_X_error, casc_Wm_X_error,
	                    casc_W_X_error , Mixed_W_X_error,delta,alpha,100,
											XX_kmeans_nmi,XX_nmi, We_nmi, Wm_nmi, Mixed_W_nmi,
											We_X_nmi, Wm_X_nmi,Mixed_LW_X_nmi, casc_We_X_nmi,casc_Wm_X_nmi,
											casc_W_X_nmi,  Mixed_W_X_nmi,100,XX_kmeans_ari,XX_ari,
											We_ari, Wm_ari, Mixed_W_ari, We_X_ari, Wm_X_ari,
											Mixed_LW_X_ari, casc_We_X_ari, casc_Wm_X_ari, casc_W_X_ari,  Mixed_W_X_ari),
							Mixed_W_X_label = Mixed_W_X,
							casc_We_X_label = casc_We_X,
							XX_kmeans_label = XX_kmeans,
							We_label = We,
							Wm_label = Wm))
}


