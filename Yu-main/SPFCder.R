#' Distance-based clustering of functional data with derivative principal component analysis.
#' 
#' @param y A list of \emph{n} vectors containing the observed values for each individual. Missing values specified by \code{NA}s are  #  supported for dense case (\code{dataType='dense'}).
#' @param yder A list of \emph{n} vectors containing the derivatives values for each individual.  
#' @param t A list of \emph{n} vectors containing the observation time points for each individual corresponding to y.
#' @param k A scalar defining the number of clusters to define; default 3. Values that define very small clusters (eg.cluster size  # #   <=3) will potentially err.
#' @param wt A scalar defining the weight of derivatives information; range over (0,1); default 0.5.
#' @param kSeed A scalar valid seed number to ensure replication; default: 123
#' @param maxIter A scalar defining the maximum number of iterations allowed; default 20, common for both the simple kmeans initially  #   and the subsequent k-centres
#' @param optnsSW A list of options control parameters specified by \code{list(name=value)} to be used for sample-wide FPCA; by  #   # #  default: "list( methodMuCovEst ='smooth', FVEthreshold= 0.90, methodBwCov = 'GCV', methodBwMu = 'GCV' )". See `Details in ?FPCA'.
#' @param optnsCS A list of options control parameters specified by \code{list(name=value)} to be used for cluster-specific FPCA; by  #  default:  "list( methodMuCovEst ='smooth', FVEthreshold= 0.70, methodBwCov = 'GCV', methodBwMu = 'GCV' )". See `Details in ?FPCA'.
#' @param optnsderSW; optnsderCS Alist of options control parameters soecified by \code{list(name=value)} to be used for FPCAder; by # #  default:(method="DPC",kernel="epan")
#' @return A list containing the following fields:
#' \item{cluster}{A vector of levels 1:k, indicating the cluster to which each curve is allocated.} 
#' \item{fpcaList}{A list with the fpcaObj for each separate cluster.} 
#' \item{iterToConv}{A number indicating how many iterations where required until convergence.} 



# (90, 70) -----> (99, 99)
kCFCder = function(y, yder, t, k = 3, wt = 0.5, kSeed = 123, maxIter = 125, fvethreshold1=0.90,
                optnsSW = list( methodMuCovEst = 'smooth', FVEthreshold = 0.99, methodBwCov = 'GCV', methodBwMu = 'GCV'), 
                optnsCS = list( methodMuCovEst = 'smooth', FVEthreshold = 0.99, methodBwCov = 'GCV', methodBwMu = 'GCV'),
                optnsderSW=list(method="DPC",kernelType="epan"),optnsderCS=list(method="DPC",kernelType="epan")){ 
  
  if( (k <2) || (floor(length(y)*0.5) < k) ){
    warning("The value of 'k' is outside [2, 0.5*N]; reseting to 3.")
  } 
  if(maxIter <1){
    stop("Please allow at least 1 iteration.")
  }
  
  ## First FPCA
  DPC<-FPCAder(FPCA(y, t, optnsSW), optnsderSW)
  N <- length(y)
  if( DPC$optns$dataType == 'Sparse' ){
    stop(paste0("The data has to be 'Dense' for kCFC to be relevant; the current dataType is : '", fpcaObjY$optns$dataType,"'!") )
  }
  
  ## Initial clustering and cluster-associated FPCAs
  if(!is.null(kSeed)){
    set.seed(kSeed)
  }
   
  K.min<-min(which(cumsum(DPC$lambdaDer)/sum(DPC$lambdaDer)>fvethreshold1))
  initialClustering <- kmeans( cbind(DPC$xiEst,DPC$xiDer[,1:K.min]), centers = k, algorithm = "MacQueen", iter.max = maxIter)
  dist.wt <- as.matrix(dist(DPC$xiEst)) + as.matrix(dist(DPC$xiDer[,1:K.min]))
  
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
    iseCosts       <- sapply(listOfFPCAderobjs, function(u) GetISEfromFPCAder(u, y, yder, t, ymat, ymatder,wt))
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
  
  kCFCderobj <-  list(cluster = clustConf[[j]], fpcaderList=listOfFPCAderobjs, iterToConv = j, prevConf = clustConf, clustConf0 =  clustConf0, dis=dist.wt)
  class(kCFCderobj) <- 'kCFCderobj'
  return( kCFCderobj )
}  

GetISEfromFPCAder = function(fpcaObj,y,yder,t,ymat,ymatder,wt){
  # First get the fitted curves for all the sample based on the mu/phi/lambda/sigma2
  # of 'fpcaObj' and then calculate their associated ISE; 'iseCost' is a n-dim vector.
    K.min<-min(which(cumsum(fpcaObj$lambdaDer)/sum(fpcaObj$lambdaDer)>0.9))

    numIntResults <- mapply(function(yvec,tvec)
    GetINScores(yvec, tvec,optns= fpcaObj$optns,fpcaObj$obsGrid,mu = fpcaObj$mu,lambda =fpcaObj$lambda ,phi =fpcaObj$phi,sigma2=fpcaObj$sigma2),y,t)
   
    numIntResultsder<-mapply(function(yvecder,tvec) 
    GetINScores(yvecder,tvec,optns= fpcaObj$optns,fpcaObj$obsGrid,mu=fpcaObj$muDer,lambda =fpcaObj$lambdaDer[1:K.min],phi=as.matrix(fpcaObj$phiDer[,1:K.min]),sigma2 = fpcaObj$sigma2),yder,t)
  
  fittedYmat = List2Mat(numIntResults[3,],t)
  fittedYmatder=List2Mat(numIntResultsder[3,],t)
  iseCostfpc <- apply((fittedYmat - ymat)^2, 1, function(y) {notNA <- !is.na(y);  trapzRcpp(X = fpcaObj$obsGrid[notNA], Y = y [notNA])}) 
  iseCostder<-apply((fittedYmatder - ymatder)^2, 1, function(y) {notNA <- !is.na(y);  trapzRcpp(X = fpcaObj$obsGrid[notNA], Y = y   [notNA])})
  iseCost<-(1-wt)*iseCostfpc + wt*iseCostder 
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