#rm(list=ls())
library(irlba)
library(aricode)  #NMI
library(mclust)  #ARI
library(RSpectra) #specClust


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
	
	#wcssgaps<-abs((d_wcssVec[1:(delta_n-1)]-d_wcssVec[2:delta_n])/d_wcssVec[1:(delta_n-1)])
	#return
	if (verbose==F){
		return(list(bestAlpha=hOpt[which.min(d_wcssVec)],  #wcssgaps
								bestDelta=deltaseq[which.min(d_wcssVec)]))
	}else{
		return(list(bestAlpha=hOpt[which.min(d_wcssVec)],
								bestDelta=deltaseq[which.min(d_wcssVec)],
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
	#wcssgaps<-abs((d_wcssVec[1:(delta_n-1)]-d_wcssVec[2:delta_n])/d_wcssVec[1:(delta_n-1)])
	bestdelta = deltaseq[which.min(d_wcssVec)]
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
	
	covMat <- scale(covMat, center = TRUE, scale = TRUE)
	doubleCov<-covMat%*%t(covMat)
	
	true_labels<-rep(1:K, each = N/K)
	misclust_labels=rep(N/K,K)

	##XX,kmeans
	XX_kmeans <- kmeans(covMat, centers = K,nstart = 20)$cluster
	
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
	#XX_kmeans1 = scale(covMat, center = F, scale = sqrt(Matrix::colSums(covMat^2))) 
	XX_kmeans <- kmeans(covMat, centers = K,nstart = 20)$cluster
	
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

# generate adjMat
generate_adjMat_complex <- function(N,K,p,B1=0.3,B2=0.9,lambda,a,b,Z){
  B0 <- runif(K,B1,B2)
  B <- p*matrix(B0 %x% B0,nrow=K,ncol=K,byrow=TRUE)
  offdiag.flag <- matrix(TRUE, nrow = K, ncol = K)
  diag(offdiag.flag) <- FALSE
  B[offdiag.flag] <- B[offdiag.flag]*(1-lambda)
  
  EA<-Z%*%B%*%t(Z)
  
  adjMat<-matrix(0,N,N)
  adjMat[upper.tri(adjMat)]<-runif(N*(N-1)/2,a,b)*rbinom(N*(N-1)/2,1,prob=as.vector(EA[upper.tri(EA)]))
  adjMat[lower.tri(adjMat)]<-t(adjMat)[lower.tri(t(adjMat))]
  
  return(adjMat)
}


# generate varies covmat
#' Generate Covariate Matrix with Various Distributions and Shapes
#'
#' @param N Number of nodes
#' @param K Number of communities
#' @param R Dimension of covariates
#' @param type Type of distribution: "spherical", "non-spherical", "non-convex"
#' @param m1,m2 Mean parameters for covariates
#' @param sd Base standard deviation
#' @param var_ratio Variance ratio between communities (for non-spherical)
#' @param shape_params Additional parameters for non-convex shapes
#' @return List containing covariate matrix and true labels
#'
generate_covMat <- function(N, K, R, Z,
                            type = "spherical",
                            m1, m2 , 
                            sd , var_ratio ,
                            shape_params = list()) {
  true_labels <- rep(1:K, each = N/K)
  if (type == "spherical") {
    # Original spherical clusters with equal variance
    M <- (m1 - m2) * diag(K) + m2 * matrix(1, K, K)
    EX <- Z %*% M
    covMat <- matrix(
      sapply(as.vector(EX), function(mu) rnorm(1, mean = mu, sd = sd)),
      N, R
    )
    
  } else if (type == "non-spherical") {
    # Clusters with different variances (elliptical distributions)
    M <- (m1 - m2) * diag(K) + m2 * matrix(1, K, K)
    
    covMat <- matrix(0, N, R)
    for (j in 1:K) {
      community_indices <- which(true_labels == j)
      n_in_community <- length(community_indices)
      community_mean <- M[j, ]
      
      # Create different covariance matrices for each community
      if (j == 1) {
        # Community 1: Spherical with base variance
        cov_matrix <- diag(R) * sd^2
      } else if (j == 2) {
        # Community 2: Elliptical with larger variance in first dimension
        cov_matrix <- diag(R) * (sd * var_ratio)^2
        cov_matrix[1, 1] <- (sd * var_ratio * 1.5)^2
      } else {
        # Community 3: Correlated features with moderate variance
        cov_matrix <- diag(R) * (sd * var_ratio)^2
        if (R >= 2) {
          cov_matrix[1, 2] <- cov_matrix[2, 1] <- 0.7 * (sd * var_ratio)^2
        }
      }
      
      # Generate multivariate normal data
      if (n_in_community > 0) {
        covMat[community_indices, ] <- MASS::mvrnorm(
          n = n_in_community,
          mu = community_mean,
          Sigma = cov_matrix
        )
      }
    }
    
  } else if (type == "non-convex") {
    # Non-convex shapes (moons, circles, etc.)
    shape_type <- ifelse(is.null(shape_params$shape_type), "moons", shape_params$shape_type)
    noise <- ifelse(is.null(shape_params$noise), 0.1, shape_params$noise)
    
    covMat <- matrix(0, N, R)
    
    if (shape_type == "moons") {
      # Generate interleaving moons (2 communities)
      for (j in 1:min(K, 2)) {
        community_indices <- which(true_labels == j)
        n_in_community <- length(community_indices)
        
        if (n_in_community > 0) {
          if (j == 1) {
            # First moon
            angles <- runif(n_in_community, 0, pi)
            radius <- 1 + rnorm(n_in_community, 0, noise)
            x <- radius * cos(angles) + m1
            y <- radius * sin(angles) + m1
          } else {
            # Second moon (shifted and inverted)
            angles <- runif(n_in_community, pi, 2*pi)
            radius <- 1 + rnorm(n_in_community, 0, noise)
            x <- radius * cos(angles) + m2 + 0.5
            y <- radius * sin(angles) + m2 - 0.5
          }
          
          covMat[community_indices, 1] <- x
          covMat[community_indices, 2] <- y
        }
      }
      
      # For K > 2, make remaining communities spherical
      if (K > 2) {
        for (j in 3:K) {
          community_indices <- which(true_labels == j)
          n_in_community <- length(community_indices)
          if (n_in_community > 0) {
            covMat[community_indices, 1] <- rnorm(n_in_community, m1 + 1.5, sd)
            covMat[community_indices, 2] <- rnorm(n_in_community, m2 + 1.5, sd)
          }
        }
      }
      
    } else if (shape_type == "circles") {
      # Concentric circles
      for (j in 1:K) {
        community_indices <- which(true_labels == j)
        n_in_community <- length(community_indices)
        
        if (n_in_community > 0) {
          angles <- runif(n_in_community, 0, 2*pi)
          radius <- j * 0.8 + rnorm(n_in_community, 0, noise)
          
          x <- radius * cos(angles) + m1
          y <- radius * sin(angles) + m2
          
          covMat[community_indices, 1] <- x
          covMat[community_indices, 2] <- y
        }
      }
    }
    
    # Add remaining dimensions as noise if R > 2
    if (R > 2) {
      for (r in 3:R) {
        covMat[, r] <- rnorm(N, mean = 0, sd = noise)
      }
    }
  }
  
  return(list(
    covMat = covMat,
    true_labels = true_labels,
    Z = Z,
    type = type
  ))
}

#' Visualize Covariate Matrix
#'
#' @param data_matrix N x R matrix of covariates
#' @param labels Vector of true cluster labels
#' @param title Plot title
#' @param dims Which dimensions to plot (default: first 2 or 3)
#'
visualize_embeddings <- function(data_matrix, labels, title = "", dims = NULL) {
  if (is.null(dims)) {
    dims <- 1:min(3, ncol(data_matrix))
  }
  
  if (length(dims) == 2) {
    # 2D visualization
    df <- data.frame(
      x = data_matrix[, dims[1]],
      y = data_matrix[, dims[2]],
      cluster = as.factor(labels)
    )
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = cluster)) +
      ggplot2::geom_point(alpha = 0.7, size = 2) +
      ggplot2::labs(title = title, x = paste("Dimension", dims[1]), y = paste("Dimension", dims[2])) +
      ggplot2::theme_minimal() +
      ggplot2::scale_color_viridis_d()
    
    print(p)
    return(p %>% plotly::layout(title = ""))
    
  } else if (length(dims) == 3) {
    # 3D visualization
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("plotly package required for 3D visualization. Installing...")
      install.packages("plotly")
    }
    
    df <- data.frame(
      x = data_matrix[, dims[1]],
      y = data_matrix[, dims[2]],
      z = data_matrix[, dims[3]],
      cluster = as.factor(labels)
    )
    
    
    p <- plotly::plot_ly(df, x = ~x, y = ~y, z = ~z, color = ~cluster, 
                         colors = RColorBrewer::brewer.pal(length(unique(labels)), "Set1"),
                         type = 'scatter3d', mode = 'markers',
                         marker = list(size = 3, opacity = 0.7)) %>%
      plotly::layout(title = title,
                     scene = list(xaxis = list(title = ""),
                                  yaxis = list(title = ""),
                                  zaxis = list(title = "")),
                     showlegend = FALSE)
    
    print(p)
    return(p %>% plotly::layout(title = ""))
  } else {
    stop("Please specify 2 or 3 dimensions for visualization")
  }
}

simul_complexCov <- function(adjMat,covMat,N,K){
  true_labels<-rep(1:K, each = N/K)
  misclust_labels=rep(N/K,K)
  
  # motif matrix
  AMat<-adjMat
  AMat[AMat!=0]<-1
  wadjMat<-adjMat%*%AMat*AMat+AMat%*%adjMat*AMat+AMat%*%AMat*adjMat
  #weight<-(AMat%*%AMat*AMat)*3
  wadjMat[wadjMat!=0]<-wadjMat[wadjMat!=0]/3
  
  # Covariance similarity
  covMat = scale(covMat, center = T, scale = T)
  doubleCov<-covMat%*%t(covMat)
  
  ##X,kmeans
  X_only <- kmeans(covMat, centers = K,nstart = 20)$cluster
  
  ##edge, A
  We <- specClust(targetMat=adjMat,K,nIter = 20, verbose = F)
  
  ##motif, A
  Wm <- specClust(targetMat=wadjMat,K,nIter = 20, verbose = F)
  
  
  ##casc_We_X
  alpha_casc_We<-tuningOnlyalpha(graphMat=adjMat, covMat, alpha_n=20, method = "CASCregLaplacian")
  casc_adjMat<-getSimilarityMat(adjMat, method = "CASCregLaplacian")
  casc_We_X<- specClust(targetMat = casc_adjMat+alpha_casc_We*doubleCov, K,nIter = 20, verbose = F)
  
  
  ##mixed, A + alpha * X %*% Matrix::t(X)
  ##mixed matrix
  beyondPara<-tuningPara(delta_n=20,alpha_n=20,verbose=F,adjMat,wadjMat,covMat,method = "adjacency",true_labels,misclust_labels)
  delta=beyondPara$bestDelta
  alpha=beyondPara$bestAlpha
  madjMat<-delta*adjMat+(1-delta)*wadjMat
  Lmat <- madjMat+alpha*doubleCov
  Mixed_W_X <- specClust(targetMat = Lmat, K,nIter = 20, verbose = F)
  
  ##error_rate
  X_only_error<-misClustRate(X_only, misclust_labels)
  We_error<-misClustRate(We, misclust_labels)
  Wm_error<-misClustRate(Wm, misclust_labels)
  casc_We_X_error<-misClustRate(casc_We_X, misclust_labels)
  Mixed_W_X_error<-misClustRate(Mixed_W_X, misclust_labels)
  
  ##nmi
  X_only_nmi<- NMI(true_labels, X_only)
  We_nmi<-NMI(true_labels,We)
  Wm_nmi<-NMI(true_labels,Wm)
  casc_We_X_nmi<-NMI(true_labels,casc_We_X)
  Mixed_W_X_nmi<-NMI(true_labels,Mixed_W_X)
  
  ##ari
  X_only_ari<- adjustedRandIndex(true_labels, X_only)
  We_ari<-adjustedRandIndex(true_labels,We)
  Wm_ari<-adjustedRandIndex(true_labels,Wm)
  casc_We_X_ari<-adjustedRandIndex(true_labels,casc_We_X)
  Mixed_W_X_ari<-adjustedRandIndex(true_labels,Mixed_W_X)
  
  re <- list(metric = c(delta,alpha,
                        X_only_error, We_error, Wm_error, casc_We_X_error, Mixed_W_X_error,
                        X_only_nmi,We_nmi, Wm_nmi, casc_We_X_nmi, Mixed_W_X_nmi, 
                        X_only_ari, We_ari, Wm_ari,casc_We_X_ari,  Mixed_W_X_ari), Lmat=Lmat)
  
  
  return(re)
  
  
}

#' Enhanced Spectral Clustering with Multiple Hyperparameter Selection and Clustering Methods
#'
#' @param targetMat Similarity matrix for spectral clustering
#' @param K Number of clusters
#' @param nIter Maximum number of iterations for clustering algorithms
#' @param verbose Whether to print progress messages
#' @param clust_method Clustering method: "kmeans", "silhouette", "dbindex"
#' @param n_init Number of random initializations for stability
#' @param criterion Hyperparameter selection criterion: "wcss", "silhouette", "dbindex"
#' @return Cluster assignments and clustering metrics
#'
specClust_enhanced <- function(targetMat, K, nIter = 20, verbose = FALSE,
                               clust_method = "kmeans", n_init = 20, 
                               criterion = "wcss") {
  
  targetMat = as(targetMat, "CsparseMatrix")
  
  # Compute eigenvectors for spectral embedding
  eigsDecomp = irlba(targetMat, nu = K + 1, nv = K + 1, maxit = 2000, work = 3 * (K + 1))
  eig_vector <- eigsDecomp$u[, 1:K]
  eig_value <- eigsDecomp$d[1:(K + 1)]
  
  # Normalize rows to unit length
  eig_vector = eig_vector / sqrt(rowSums(eig_vector^2))
  eig_vector[is.na(eig_vector) | is.nan(eig_vector) | is.infinite(eig_vector)] <- 0
  
  # Apply selected clustering method
  clustering_result <- switch(clust_method,
                              "kmeans" = {
                                # Standard k-means with multiple random initializations (built-in CV)
                                kmeans(eig_vector, centers = K, nstart = n_init, iter.max = nIter)
                              },
                              "silhouette" = {
                                # Clustering optimized for silhouette score
                                cluster_optimize_silhouette(eig_vector, K, n_init, nIter)
                              },
                              "dbindex" = {
                                # Clustering optimized for Davies-Bouldin index
                                cluster_optimize_dbindex(eig_vector, K, n_init, nIter)
                              },
                              {
                                # Default to standard k-means
                                warning("Unknown clustering method. Using standard kmeans.")
                                kmeans(eig_vector, centers = K, nstart = n_init, iter.max = nIter)
                              }
  )
  
  if (verbose) {
    # Calculate additional metrics from kmeans result
    metrics <- calculate_clustering_metrics(eig_vector, clustering_result$cluster, criterion)
    
    return(list(
      cluster = clustering_result$cluster,
      wcss = clustering_result$tot.withinss,
      eigenVals = eig_value,
      eig_vector = eig_vector,
      clustering_method = clust_method,
      metrics = metrics
    ))
  } else {
    return(clustering_result$cluster)
  }
}

#' Calculate Clustering Metrics from kmeans Result
#'
#' @param data Input data matrix (eigenvectors)
#' @param clusters Cluster assignments from kmeans
#' @param criterion Metric to calculate
#' @return Calculated metric value
#'
calculate_clustering_metrics <- function(data, clusters, criterion = "wcss") {
  switch(criterion,
         "wcss" = {
           # Directly use WCSS from kmeans result (no need to recalculate)
           # This would be calculated from the clustering result elsewhere
           NULL  # Placeholder - actual WCSS comes from kmeans object
         },
         "silhouette" = {
           # Silhouette score
           if (requireNamespace("cluster", quietly = TRUE)) {
             sil <- cluster::silhouette(clusters, dist(data))
             mean(sil[, 3])  # Return average silhouette width
           } else {
             warning("cluster package not available for silhouette calculation.")
             NA
           }
         },
         "dbindex" = {
           # Davies-Bouldin index
           if (requireNamespace("clusterSim", quietly = TRUE)) {
             clusterSim::index.DB(data, clusters)$DB
           } else {
             warning("clusterSim package not available for DB index calculation.")
             NA
           }
         }
  )
}

#' Clustering Optimized for Silhouette Score
#'
#' @param data Input data matrix
#' @param K Number of clusters
#' @param n_init Number of random initializations
#' @param nIter Maximum iterations
#' @return Clustering result with best silhouette score
#'
cluster_optimize_silhouette <- function(data, K, n_init = 20, nIter = 20) {
  best_silhouette <- -Inf
  best_result <- NULL
  
  for (i in 1:n_init) {
    current_result <- kmeans(data, centers = K, nstart = 1, iter.max = nIter)
    current_silhouette <- calculate_clustering_metrics(data, current_result$cluster, "silhouette")
    
    if (!is.na(current_silhouette) && current_silhouette > best_silhouette) {
      best_silhouette <- current_silhouette
      best_result <- current_result
    }
  }
  
  # If no valid silhouette scores, fall back to standard kmeans
  if (is.null(best_result)) {
    best_result <- kmeans(data, centers = K, nstart = n_init, iter.max = nIter)
  }
  
  return(best_result)
}

#' Clustering Optimized for Davies-Bouldin Index
#'
#' @param data Input data matrix
#' @param K Number of clusters
#' @param n_init Number of random initializations
#' @param nIter Maximum iterations
#' @return Clustering result with best DB index
#'
cluster_optimize_dbindex <- function(data, K, n_init = 20, nIter = 20) {
  best_dbindex <- Inf
  best_result <- NULL
  
  for (i in 1:n_init) {
    current_result <- kmeans(data, centers = K, nstart = 1, iter.max = nIter)
    current_dbindex <- calculate_clustering_metrics(data, current_result$cluster, "dbindex")
    
    if (!is.na(current_dbindex) && current_dbindex < best_dbindex) {
      best_dbindex <- current_dbindex
      best_result <- current_result
    }
  }
  
  # If no valid DB indices, fall back to standard kmeans
  if (is.null(best_result)) {
    best_result <- kmeans(data, centers = K, nstart = n_init, iter.max = nIter)
  }
  
  return(best_result)
}

#' Enhanced Hyperparameter Tuning with Multiple Criteria
#'
#' @param delta_n Number of delta values to try
#' @param alpha_n Number of alpha values to try
#' @param verbose Whether to print progress
#' @param adjMat Adjacency matrix
#' @param wadjMat Weighted adjacency matrix
#' @param covMat Covariate matrix
#' @param method Similarity matrix method
#' @param true_labels True cluster labels
#' @param misclust_labels Misclustering labels
#' @param criterion Hyperparameter selection criterion
#' @param clust_method Clustering method
#' @param n_init Number of initializations
#' @return Best hyperparameters and evaluation metrics
#'
tuningPara_enhanced <- function(delta_n = 20, alpha_n = 20, verbose = FALSE,
                                adjMat, wadjMat, covMat, method = "adjacency",
                                true_labels, misclust_labels,
                                criterion = "wcss", clust_method = "kmeans", n_init = 20) {
  
  deltaseq <- seq(0, 1, length.out = delta_n)
  hOpt <- vector(length = delta_n)
  
  d_criterionVec <- vector(length = delta_n)  # Store criterion values for each delta
  d_wcssVec <- vector(length = delta_n)
  d_gapValK <- vector(length = delta_n)
  error_rate <- vector(length = delta_n)
  nmi_value <- vector(length = delta_n)
  Valgaps <- vector(length = delta_n)
  
  doubleCov <- covMat %*% t(covMat)
  
  for (j in 1:delta_n) {
    delta <- deltaseq[j]
    madjMat <- delta * adjMat + (1 - delta) * wadjMat
    madjMat <- getSimilarityMat(madjMat, method)
    rangehTuning <- getTuningRange(madjMat, covMat, K = K, method = method)
    hTuningSeq <- seq(rangehTuning$alpha_min, rangehTuning$alpha_max, length.out = alpha_n)
    
    a_criterionVec <- vector(length = alpha_n)
    for (i in 1:alpha_n) {
      alpha <- hTuningSeq[i]
      targetMat <- madjMat + alpha * doubleCov
      res_cluster <- specClust_enhanced(targetMat, K, nIter = 20, verbose = TRUE,
                                        clust_method = clust_method, n_init = n_init,
                                        criterion = criterion)
      
      # Use appropriate criterion for selection
      if (criterion == "wcss") {
        a_criterionVec[i] <- res_cluster$wcss
      } else if (criterion == "silhouette") {
        a_criterionVec[i] <- res_cluster$metrics
      } else if (criterion == "dbindex") {
        a_criterionVec[i] <- res_cluster$metrics
      }
    }
    
    # Select best alpha based on criterion
    if (criterion == "silhouette") {
      best_idx <- which.max(a_criterionVec)
    } else {
      best_idx <- which.min(a_criterionVec)
    }
    hOpt[j] <- hTuningSeq[best_idx]
    d_criterionVec[j] <- a_criterionVec[best_idx]  # Store best criterion value for this delta
    
    # Get final results with selected alpha
    targetMat <- madjMat + hOpt[j] * doubleCov
    res_cluster <- specClust_enhanced(targetMat, K, nIter = 20, verbose = TRUE,
                                      clust_method = clust_method, n_init = n_init,
                                      criterion = criterion)
    
    d_wcssVec[j] <- res_cluster$wcss
    d_gapValK[j] <- res_cluster$eigenVals[K]
    Valgaps[j] <- (res_cluster$eigenVals[K] - res_cluster$eigenVals[K + 1]) / res_cluster$eigenVals[K]
    error_rate[j] <- misClustRate(res_cluster$cluster, misclust_labels)
    nmi_value[j] <- NMI(true_labels, res_cluster$cluster)
  }
  
  # Select best delta using the same criterion
  if (criterion == "silhouette") {
    best_delta_idx <- which.max(d_criterionVec)
  } else {
    best_delta_idx <- which.min(d_criterionVec)
  }
  
  if (verbose == FALSE) {
    return(list(
      bestAlpha = hOpt[best_delta_idx],
      bestDelta = deltaseq[best_delta_idx]
    ))
  } else {
    return(list(
      bestAlpha = hOpt[best_delta_idx],
      bestDelta = deltaseq[best_delta_idx],
      criterion_values = d_criterionVec,
      alphaOpt = hOpt,
      d_wcssVec = d_wcssVec,
      d_gapValK = d_gapValK,
      Valgaps = Valgaps,
      error_rate = error_rate,
      nmi_value = nmi_value,
      criterion_used = criterion
    ))
  }
}

#' Unified Hyperparameter Tuning Function (Replaces tuningPara, tuningOnlydelta, tuningOnlyalpha)
#'
#' @param adjMat Adjacency matrix
#' @param wadjMat Weighted adjacency matrix (optional)
#' @param covMat Covariate matrix (optional)
#' @param delta_n Number of delta values
#' @param alpha_n Number of alpha values
#' @param method Similarity matrix method
#' @param tuning_type Type of tuning: "both", "delta_only", "alpha_only"
#' @param criterion Hyperparameter selection criterion
#' @param clust_method Clustering method
#' @param n_init Number of initializations
#' @param verbose Whether to print progress
#' @param true_labels True labels for evaluation
#' @param misclust_labels Misclustering labels
#' @return Tuning results
#'
tuning_unified <- function(adjMat, wadjMat = NULL, covMat = NULL,
                           delta_n = 20, alpha_n = 20, method = "adjacency",
                           tuning_type = "both", criterion = "wcss", 
                           clust_method = "kmeans", n_init = 20, verbose = FALSE,
                           true_labels = NULL, misclust_labels = NULL) {
  
  if (tuning_type == "delta_only") {
    # Tune only delta (no covariates)
    deltaseq <- seq(0, 1, length.out = delta_n)
    d_criterionVec <- vector(length = delta_n)
    d_wcssVec <- vector(length = delta_n)
    eigValK <- vector(length = delta_n)
    error_rate <- vector(length = delta_n)
    
    for (i in 1:delta_n) {
      delta <- deltaseq[i]
      madjMat <- delta * adjMat + (1 - delta) * wadjMat
      targetMat <- getSimilarityMat(madjMat, method)
      res_cluster <- specClust_enhanced(targetMat, K, nIter = 20, verbose = TRUE,
                                        clust_method = clust_method, n_init = n_init,
                                        criterion = criterion)
      
      # Store criterion value based on the specified criterion
      if (criterion == "wcss") {
        d_criterionVec[i] <- res_cluster$wcss
      } else if (criterion == "silhouette") {
        d_criterionVec[i] <- res_cluster$metrics
      } else if (criterion == "dbindex") {
        d_criterionVec[i] <- res_cluster$metrics
      }
      
      d_wcssVec[i] <- res_cluster$wcss
      eigValK[i] <- res_cluster$eigenVals[K]
      if (!is.null(misclust_labels)) {
        error_rate[i] <- misClustRate(res_cluster$cluster, misclust_labels)
      }
    }
    
    # Select best delta using the specified criterion
    if (criterion == "silhouette") {
      best_delta_idx <- which.max(d_criterionVec)
    } else {
      best_delta_idx <- which.min(d_criterionVec)
    }
    bestdelta <- deltaseq[best_delta_idx]
    
    if (verbose) {
      return(list(bestDelta = bestdelta, d_criterionVec = d_criterionVec, d_wcssVec = d_wcssVec,
                  eigValK = eigValK, error_rate = error_rate))
    } else {
      return(bestdelta)
    }
    
  } else if (tuning_type == "alpha_only") {
    # Tune only alpha (fixed graph structure)
    graphMat <- getSimilarityMat(adjMat, method)
    rangehTuning <- getTuningRange(graphMat, covMat, K = K, method = method)
    hTuningSeq <- seq(rangehTuning$alpha_min, rangehTuning$alpha_max, length.out = alpha_n)
    a_criterionVec <- vector(length = alpha_n)
    
    doubleCov <- covMat %*% t(covMat)
    for (i in 1:alpha_n) {
      targetMat <- graphMat + hTuningSeq[i] * doubleCov
      res_cluster <- specClust_enhanced(targetMat, K, nIter = 20, verbose = TRUE,
                                        clust_method = clust_method, n_init = n_init,
                                        criterion = criterion)
      
      if (criterion == "wcss") {
        a_criterionVec[i] <- res_cluster$wcss
      } else if (criterion == "silhouette") {
        a_criterionVec[i] <- res_cluster$metrics
      } else if (criterion == "dbindex") {
        a_criterionVec[i] <- res_cluster$metrics
      }
    }
    
    if (criterion == "silhouette") {
      hOpt <- hTuningSeq[which.max(a_criterionVec)]
    } else {
      hOpt <- hTuningSeq[which.min(a_criterionVec)]
    }
    
    return(hOpt)
    
  } else {
    # Tune both delta and alpha (default)
    return(tuningPara_enhanced(delta_n, alpha_n, verbose, adjMat, wadjMat, covMat, 
                               method, true_labels, misclust_labels, criterion, 
                               clust_method, n_init))
  }
}

#' Enhanced Simulation with Complex Covariates and Multiple Methods
#'
#' @param adjMat Adjacency matrix
#' @param covMat Covariate matrix
#' @param N Number of nodes
#' @param K Number of clusters
#' @param criterion Hyperparameter selection criterion
#' @param clust_method Clustering method
#' @param n_init Number of initializations
#' @return Simulation results with various metrics
#'
simul_complexCov_enhanced <- function(adjMat, covMat, N, K,
                                      criterion = "wcss", clust_method = "kmeans", n_init = 20) {
  
  true_labels <- rep(1:K, each = N / K)
  misclust_labels <- rep(N / K, K)
  
  # Motif matrix calculation
  AMat <- adjMat
  AMat[AMat != 0] <- 1
  wadjMat <- adjMat %*% AMat * AMat + AMat %*% adjMat * AMat + AMat %*% AMat * adjMat
  wadjMat[wadjMat != 0] <- wadjMat[wadjMat != 0] / 3
  
  # Covariance similarity with scaling
  covMat <- scale(covMat, center = TRUE, scale = TRUE)
  doubleCov <- covMat %*% t(covMat)
  
  # Different clustering approaches
  X_only <- specClust_enhanced(doubleCov, K, nIter = 20, verbose = FALSE,
                               clust_method = clust_method, n_init = n_init)
  
  We <- specClust_enhanced(adjMat, K, nIter = 20, verbose = FALSE,
                           clust_method = clust_method, n_init = n_init)
  
  Wm <- specClust_enhanced(wadjMat, K, nIter = 20, verbose = FALSE,
                           clust_method = clust_method, n_init = n_init)
  
  # CASC with enhanced tuning
  alpha_casc_We <- tuning_unified(
    adjMat = adjMat, covMat = covMat, alpha_n = 20,
    tuning_type = "alpha_only", criterion = criterion,
    clust_method = clust_method, n_init = n_init
  )
  
  casc_adjMat <- getSimilarityMat(adjMat, method = "CASCregLaplacian")
  casc_We_X <- specClust_enhanced(
    casc_adjMat + alpha_casc_We * doubleCov, K, nIter = 20, verbose = FALSE,
    clust_method = clust_method, n_init = n_init
  )
  
  # Mixed method with enhanced tuning
  beyondPara <- tuning_unified(
    adjMat = adjMat, wadjMat = wadjMat, covMat = covMat,
    delta_n = 20, alpha_n = 20, tuning_type = "both",
    criterion = criterion, clust_method = clust_method, n_init = n_init,
    true_labels = true_labels, misclust_labels = misclust_labels
  )
  
  delta <- beyondPara$bestDelta
  alpha <- beyondPara$bestAlpha
  madjMat <- delta * adjMat + (1 - delta) * wadjMat
  Lmat <- madjMat + alpha * doubleCov
  Mixed_W_X <- specClust_enhanced(Lmat, K, nIter = 20, verbose = FALSE,
                                  clust_method = clust_method, n_init = n_init)
  
  # Calculate evaluation metrics
  X_only_error <- misClustRate(X_only, misclust_labels)
  We_error <- misClustRate(We, misclust_labels)
  Wm_error <- misClustRate(Wm, misclust_labels)
  casc_We_X_error <- misClustRate(casc_We_X, misclust_labels)
  Mixed_W_X_error <- misClustRate(Mixed_W_X, misclust_labels)
  
  X_only_nmi <- NMI(true_labels, X_only)
  We_nmi <- NMI(true_labels, We)
  Wm_nmi <- NMI(true_labels, Wm)
  casc_We_X_nmi <- NMI(true_labels, casc_We_X)
  Mixed_W_X_nmi <- NMI(true_labels, Mixed_W_X)
  
  X_only_ari <- adjustedRandIndex(true_labels, X_only)
  We_ari <- adjustedRandIndex(true_labels, We)
  Wm_ari <- adjustedRandIndex(true_labels, Wm)
  casc_We_X_ari <- adjustedRandIndex(true_labels, casc_We_X)
  Mixed_W_X_ari <- adjustedRandIndex(true_labels, Mixed_W_X)
  
  re <- list(
    metric = c(X_only_error, We_error, Wm_error, casc_We_X_error, Mixed_W_X_error,
               X_only_nmi, We_nmi, Wm_nmi, casc_We_X_nmi, Mixed_W_X_nmi,
               X_only_ari, We_ari, Wm_ari, casc_We_X_ari, Mixed_W_X_ari),
    Lmat = Lmat,
    hyperparams = list(delta = delta, alpha = alpha),
    criterion = criterion,
    clust_method = clust_method
  )
  
  return(re)
}
