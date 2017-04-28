poly.kernel <- function(v1, v2=v1, p=2) { 
	
	((as.matrix(v1) %*% t(v2))+1)^p
	
	}
		

inverse <- function(TrainObjects, eps = 1e-28){
	
	X <- poly.kernel(TrainObjects)
	
	eig <- eigen(X, symmetric = T)
	
	P <- eig[[2]]
	
	lambda <- eig[[1]]
	
	ind <- lambda > eps
	
	lambda[ind] <- 1/lambda[ind]
	
	lambda[!ind] <- 0
	
	ans <- P %*% diag(lambda) %*% t(P)
	
	return(ans)
}


KernelRidgeReg <- function(TrainObjects,TrainLabels,TestObjects,lambda){
 	
  X <- TrainObjects
 	 	
  y <- TrainLabels
 		 						
  kernel <- poly.kernel(X)
  
  design.mat <- cbind(1, kernel)
  
  I <- rbind(0, cbind(0, kernel))
  
  M <- crossprod(design.mat) + lambda*I
  #crossprod is just x times  traspose of x, just looks neater in my openion
  
  M.inv <- inverse(M)
  #solve function in R returns singular matrix, this one avoids it.
  
  k <- as.matrix(diag(poly.kernel(cbind(TrainObjects,TrainLabels))))
  #Removing diag still gives the same MSE, but will output a vector of prediction.
  
  Labels <- rbind(0,as.matrix(TrainLabels))
    
  y.hat <- t(Labels) %*% M.inv %*% rbind(0,k)
  
  y.true <- Y.test
  
  MSE <-mean((y.hat - y.true)^2) 
    
  return(list(MSE=MSE,y.hat=y.hat))
	
}

KernelRidgeReg(X.train,Y.train,X.test,0.5)
