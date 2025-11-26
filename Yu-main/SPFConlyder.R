###Distance-based clustering of functional data with derivative principal component analysis
###Implement the clustering procedure only using the derivative information
### the parameter setting is same as SPFCder.
kCFConlyder = function(y, yder, t, k = 3, kSeed = 123, maxIter = 125, fvethreshold1=0.9,
                optnsSW = list( methodMuCovEst = 'smooth', FVEthreshold = 0.90, methodBwCov = 'GCV', methodBwMu = 'GCV'), 
                optnsCS = list( methodMuCovEst = 'smooth', FVEthreshold = 0.70, methodBwCov = 'GCV', methodBwMu = 'GCV'),
                optnsderSW=list(method="DPC",kernelType="epan"),optnsderCS=list(method="DPC",kernelType="epan")){ 
  
  if( (k <2) || (floor(length(y)*0.5) < k) ){
    warning("The value of 'k' is outside [2, 0.5*N]; reseting to 3.")
  } 
  if(maxIter <1){
    stop("Please allow at least 1 iteration.")
  }
  
  ## First FPCA
  DPC<-FPCAder(FPCA(y, t, optnsSW),optnsderSW)
  N <- length(y)
  if( DPC$optns$dataType == 'Sparse' ){
    stop(paste0("The data has to be 'Dense' for kCFC to be relevant; the current dataType is : '", fpcaObjY$optns$dataType,"'!") )
  }
  
  ## Initial clustering and cluster-associated FPCAs
  if(!is.null(kSeed)){
    set.seed(kSeed)
  }
  
  K.min<-min(which(cumsum(DPC$lambdaDer)/sum(DPC$lambdaDer)>fvethreshold1))
  initialClustering <- kmeans(as.matrix(DPC$xiDer[,1:K.min]), centers = k, algorithm = "MacQueen", iter.max = maxIter)

  clustConf0 <- as.factor(initialClustering$cluster)
  indClustIds <- lapply(levels(clustConf0), function(u) which(clustConf0 == u) )
  if( any( min( sapply( indClustIds, length)) <= c(3)) ){
    stop(paste0("kCFC stopped during the initial k-means step. The smallest cluster has three (or less) curves. " ,
                "Consider using a smaller number of clusters (k) or a different random seed (kSeed)."))
  }

  listOfFPCAderobjs<- lapply(indClustIds, function(u) FPCAder(FPCA(y[u], t[u], optnsCS),optnsderCS))
  
  ## Iterative clustering
  ymat <- List2Mat(y,t); 
  ymatder<-List2Mat(yder,t)

  convInfo <- "None"
  clustConf <- list() 
  
  for(j in seq_len(maxIter)){ 
    
    # Get new costs and relevant cluster configuration
    iseCosts       <- sapply(listOfFPCAderobjs, function(u) GetISEfromFPCAonlyder(u,yder,t,ymatder))
    clustConf[[j]] <- as.factor(apply(iseCosts, 1, which.min))
    
    # Check that clustering progressed reasonably 
    #ie. Still having k clster AND the minimum cluster size is reasonable 
    if( (length(unique(clustConf[[j]])) < k) || any( min(summary(clustConf[[j]])) <= c(0.01 * N,3))){
      convInfo <- ifelse( length(unique(clustConf[[j]])) < k , "LostCluster", "TinyCluster")
      break;
    }
    # Check if algorithm converged
    if( (j>=1) && any(sapply(clustConf[1:(j-1)], function(u) all(u == clustConf[[j]]))) ){
      convInfo <- "WeMadeIt!"
      break;
    } 
    
    indClustIds       <- lapply(levels(clustConf[[j]]), function(u) which(clustConf[[j]] == u) )
    listOfFPCAderobjs<- lapply(indClustIds, function(u) FPCAder(FPCA(y[u], t[u], optnsCS),optnsderCS) )
    curvesThatChanged <- ifelse(j > 1, sum(!( as.numeric(clustConf[[j]])  == as.numeric(clustConf[[j-1]] ))),
                                sum(!( as.numeric(clustConf[[j]])  == as.numeric(clustConf0))))
  } 
  
  if(convInfo == 'None'){
    warning(paste0( 'FkC did not converge after maxIter = ', maxIter, ' iterations. ', curvesThatChanged, ' curve(s) are undecided.'))
  }
  if(convInfo == 'TinyCluster'){
    warning(paste0("kCFCder did not fully converge. It stopped because the smallest cluster has ",
                   "less than 1% of the samples' curves. Consider using a smaller number of clusters."))
  } 
  if(convInfo == 'LostCluster'){
    warning(paste0("kCFCder did not fully converge. It stopped because it 'lost a cluster'. Consider using a smaller number of  clusters."))
  }
  
  kCFCderobj <-  list(cluster = clustConf[[j]], fpcaderList=listOfFPCAderobjs, iterToConv = j, prevConf = clustConf, clustConf0 =  clustConf0)
  class(kCFCderobj) <- 'kCFCderobj'
  return( kCFCderobj )
}  

GetISEfromFPCAonlyder = function(fpcaObj,yder,t,ymatder){
  # First get the fitted curves for all the sample based on the mu/phi/lambda/sigma2
  # of 'fpcaObj' and then calculate their associated ISE; 'iseCost' is a n-dim vector.

    K.min<-min(which(cumsum(fpcaObj$lambdaDer)/sum(fpcaObj$lambdaDer)>0.9))
 
    numIntResultsder<-mapply(function(yvecder,tvec) 
    GetINScores(yvecder,tvec,optns= fpcaObj$optns,fpcaObj$obsGrid,mu=fpcaObj$muDer,lambda =fpcaObj$lambdaDer[1:K.min],phi =as.matrix(fpcaObj$phiDer[,1:K.min]),sigma2 = fpcaObj$sigma2),yder,t)
  fittedYmatder=List2Mat(numIntResultsder[3,],t)
  iseCostder<-apply((fittedYmatder - ymatder)^2, 1, function(y) {notNA <- !is.na(y);  trapzRcpp(X = fpcaObj$obsGrid[notNA], Y = y   [notNA])})
  iseCost<-iseCostder 
  return( iseCost )
}


GetINScores <- function(yvec, tvec, optns,obsGrid, mu, lambda, phi, sigma2=NULL){
  if(length(lambda) != ncol(phi)){
    stop('No. of eigenvalues is not the same as the no. of eigenfunctions.')
  }
  
  #tau = sort(unique(signif( unlist(t),14 ))) # get observed time grid
  ranget <- diff(range(tvec))
  mu= approx(obsGrid,mu,tvec)$y
  cy = yvec - mu
  phi = apply(phi,2,function(phivec){return(approx(obsGrid,phivec,tvec)$y)})
  
  xiEst = matrix(0,length(lambda)) 
  # Get Scores xiEst
  for(i in 1:length(lambda)){
    temp = cy * phi[,i]
    xiEst[i,1] = trapzRcpp(X = tvec[!is.na(temp)], Y = temp[!is.na(temp)])
    if (optns[['shrink']] && !is.null(sigma2)) {
      xiEst[i,1] <- xiEst[i,1] * lambda[i] / 
        (lambda[i] + ranget * sigma2 / length(tvec))
    }
  }
  
  # Get Fitted Y: n by p matrix on observed time grid
  fittedY = mu + t(phi %*% xiEst)
  
  ret = list('xiEst' = xiEst,xiVar=matrix(NA, length(lambda), length(lambda)), 'fittedY' = fittedY)
  
  return(ret)
  
}



# This function converts dense regular functional input lists
# to a matrix for easier dense case implementation
##########################################################################
# Input:  - y: list of n dense regular observed p-dim functional objects
##########################################################################
# Output: - ymat: n by p matrix containing all functional data
##########################################################################

List2Mat <- function(y,t){
  n = length(y)
  obsGrid = sort(unique(unlist(t)))
  ymat = matrix( rep(NA, n * length(obsGrid)), nrow = n, byrow = TRUE)
  
  for (i in 1:n){
    ymat[i, is.element(obsGrid, t[[i]])] = y[[i]]   
  }
  return(ymat)
}
