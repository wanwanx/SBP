require(tidyverse)
require(data.table)
require(igraph)
require(Matrix)
require(RColorBrewer)


#' Main function to call for SEK
#'
#' Need to compute the average A, laplacian L, eigenvector U, supervised R
#' @param U matrix with columns are leading eigenvector of Laplacian L. 
#' @param R Separation preference matrix, in our case the correlation matrix between outcome and pair (i,j)
#' @param K Number of cluster. 
#' @param lambda A number as penalty factor, default is 0 
#' @param nstart Number of random initialization, default is 1
#' @param max.iter Maximum number of iteration if not converged
#' @param tol Tolerance level pre-set for convergence, default is 1e-4 
#' @param trace Logical constant for intermediate output, default if FALSE. 
#' @return Result list, including estimated center, cluster assignment, objective function values, and corresponding lambda 
SBP <- function(U,R,K,lambda=0,nstart=1,max.iter=20,tol=1e-4,trace=FALSE){
  result <- list()
  # nstart is the number of random initizalizations. 
  for(m in 1:nstart){
    result[[m]] <- SupervisedKmeans(U=U,R=R,K=K,lambda=lambda,
                                            max.iter=max.iter,tol=tol,
                                            trace=trace)
  }
  objs <- unlist(lapply(result,function(x){min(x$obj)}))
  # pick the one with smallest objective function
  m.min <- which.min(objs)
  return(result[[m.min]])
}


#' Core function for performing SEK for one randomization
#'
#' @param U matrix with columns are leading eigenvector of Laplacian L. 
#' @param R Separation preference matrix, in our case the correlation matrix between outcome and pair (i,j)
#' @param K Number of cluster. 
#' @param lambda A number as penalty factor, default is 0 
#' @param max.iter Maximum number of iteration if not converged
#' @param tol Tolerance level pre-set for convergence, default is 1e-4 
#' @param trace Logical constant for intermediate output, default if FALSE. 
#' @return Result list, including estimated center, cluster assignment, objective function values, and corresponding lambda 
SupervisedKmeans <- function(U,R,K,lambda=0,max.iter=20,tol=1e-4,trace=FALSE){
  
  # p is the number of voxels. 
  p <- nrow(U)
  
  #### initialize by Kmeans
  km <- kmeans(U,centers=K,nstart=1,algorithm="Lloyd")
  clusters  <- km$cluster
  centers <- km$centers
  # Z is the matrix of cluster assignment, p x k
  Z <- matrix(0,nrow=p,ncol=K)
  Z[cbind(1:p,clusters)] <- 1
  Z <- Matrix(Z,sparse=TRUE)
  obj <- NULL
  URow2Norm <- rowSums(U^2)
  for(iter in 1:max.iter){
    
    D.z.mat <- outer(URow2Norm, rowSums(centers^2),"+") -2*U%*%t(centers)  
    ## calculate the squared distance between each row of U and each center: 
    ## ||u_i - mu_k||^2 = ||u_i||^2+||mu_k||^2-2<u_i,mu_k>

    co.cluster <- Z%*%t(Z)
    
    current.obj <- sum(D.z.mat[cbind(1:p,clusters)]) + lambda*sum(co.cluster*R)/2  
    ## We do this to print out the objective. It is not needed for the algorithm to proceed.
    
    if(trace) print(paste("The objective before iteration",iter,":",current.obj,"..."))
    
    obj <- c(obj,current.obj)
    
    if((iter>1)&&(current.obj>(1-tol)*obj[iter-1])){
      if(trace)       print("Convergence!")
      break
    }
    
    for(i in 1:p){
      R.z.i <- R[i,]%*%Z 
      ##### R score with the K clusters
      Dist.z.i <- D.z.mat[i,]
      min.k <- which.min(as.numeric(Dist.z.i+lambda*R.z.i))
      Z[i,clusters[i]] <- 0
      clusters[i] <- min.k
      Z[i,min.k] <- 1
    }
    dk <- colSums(Z)
    centers <- t(Z)%*%U/dk
  }
  return(list(centers=as.matrix(centers),cluster=clusters,obj=obj,lambda=lambda))
}



