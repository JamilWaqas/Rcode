AAL <- function(y, experts, lambda, w0 = NULL, training = NULL) {
  experts <- as.matrix(experts)
  N <- ncol(experts)
  T <- nrow(experts)
  
  # Uniform initial weight vector if unspecified
  if (is.null(w0)) {
    w0 <- matrix(1/N, ncol = N)
  }
  if (sum(is.na(experts)) > 0) {
    warning("There are not allowed NA's in expert advice")
  }
  
  w <- matrix(0, ncol = N, nrow = T)
  g<-(1/N)
  
  if (!is.null(training)) {
    At <- training$At
    bt <- training$bt
  } else {
    At <- diag(1, N)
    bt <- matrix(0, nrow = N)
  }
  
  for (t in 1:T) {
    
    Dt<- diag(sqrt(abs(g)),N)
    At<- At + (experts[t,] %*% t(experts[t,])) + (ginv(Dt)*lambda)
    w[t,]<-  ginv(At)%*%bt  #solve(At,bt) # 
    g<- as.numeric(t(experts[t,])%*%(w[t,]))   #sqrt(abs(w[t,]))
    bt<- bt + y[t]*experts[t,]
  }
  # w[1,] = w0
  
  object <- list(model = "Lasso", loss.type = list(name = "square"), coefficients = solve(At, bt))
  
  object$parameters <- list(lambda = lambda)
  object$weights <- w
  object$prediction <-  apply(experts * w, 1, sum)
  
  object$training <- list(At = At, bt = bt)
  
  return(object)
} 
