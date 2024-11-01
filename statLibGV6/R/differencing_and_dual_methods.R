#' Compute kth order differencing matrix
#' 
#' @param k order of the differencing matrix
#' @param n Number of time points
#' @export
myGetDkn <- function(k, n) {
  
  if(k==1){
    D <- diag(-1, n, n)#rep(-1, n)
    #D <- diag(d)
    for (row in 1:n-1){  D[row,row+1] <-1}  
  }
  else if(k==2){
    D <- diag(-1, n, n)#rep(-1, n)
    #D <- diag(d)
    for (row in 1:n-1){  
      D[row,row+1] <- 2
    }
    for (row in 1:n-2){  
      D[row,row+2] <- -1
      if (row!=1){D[row, 1] <- 0 }
    }
  }
  return(D)
}

#' Compute KKT residual
#' 
#' @param y response
#' @param b primal variable
#' @param theta primal variable
#' @param v dual variable
#' @param D differencing matrix
#' @param lambda regularization parameter
#' @export
kkt_residual <- function(y, b, theta, v, D, lambda) {
  k1 <- theta - D%*%b
  k2 <- b - y + t(D)%*%v
  
  n <-length(theta)
  z <- vector("numeric",n)
  for(j in 1:n){
    if(theta[j]>0.0){
      z[j]<- abs(v[j] -  lambda)
    }
    else if(theta[j]<0.0){
      z[j]<- abs(v[j] +  lambda)
    }
  else{
   z[j] <- sign(v[j])*(abs(v[j]) - lambda) 
  }
  }
  k3 <- z
  KKT<- c(max(abs(k1)),max(abs(k2)),max(abs(k3)))
  return(max(KKT))

}


#' Compute dual objective
#' 
#' @param y response
#' @param D differencing matrix
#' @param v dual
#' @param lambda range variable
#' @export
dual <- function(y, v, D, lambda) {
  #check if v satisfies constraint
  cond <- inf_norm(v,lambda)
  temp <- t(D) %*% v
  diff <- (y - temp)^2
  norm_2 <- sum(diff)
  final <- (0.5)*norm_2
  return(final)
}

#' Compute dual objective gradient
#' 
#' @param y response
#' @param D differencing matrix
#' @param v dual
#' @param lambda range variable
#' @export
dual_grad <- function(y, v, D, lambda) {
  temp <- t(D)%*%v
  diff <- temp - y
  grad_temp <-(t(D)%*% v - y)
  return(D%*%grad_temp)
}


#' Compute duality gap
#' 
#' @param y response
#' @param b primal variable
#' @param v dual variable
#' @param D differencing matrix
#' @param lambda regularization parameter
#' @export
duality_gap <- function(y, b, v, D, lambda) {
  #compute primal objective
  diff1 <- y-b
  err <- sum(diff1^2)
  
  theta = D %*% b 
  reg <- sum(abs(theta))
  
  primal <- (0.5)*err + lambda*reg
  
  #compute dual objective
  y_l2 <- (0.5)* sum(y^2)
  diff2 <- y - (t(D)%*%v)
  diff2_l2 <- (0.5)* sum(diff2^2)
  dual <- y_l2 - diff2_l2
  duality_gap <- primal - dual
  return(duality_gap)
}

#' Compute projection to a closed interval -k,k
#' 
#' @param z variable to be projected
#' @param lambda range variable
#' @export
project <- function(z, lambda) {
  
  if (lambda < z){
    p <- lambda  
  } 
  else if ((0.0-lambda) > z){
    p <- (0.0-lambda)
  } 
  else{
    p <- z
  }
  return(p)
}

#' Compute infinity norm to closed interval
#' 
#' @param v variable to be projected
#' @param lambda range variable
#' @export
inf_norm <- function(v, lambda) {
  v_abs <- abs(v)
  m <- max(v_abs)
  if (m < lambda){return (TRUE)}
  else{ return (FALSE)}
}


#' Compute infinity norm to closed interval
#' 
#' @param v variable to be projected
#' @param lambda range variable
#' @export
nuclear_norm <- function(D){
  ss <- svd(D)
  dd <-ss$d
  nuc <- sum(dd^2)
  return(nuc)
}

#' Compute infinity norm to closed interval
#' 
#' @param v variable to be projected
#' @param lambda range variable
#' @export
compute_rk <- function(y, v, D, k){
  A <- t(D)
  n <- length(v)
  temp <- vector("numeric",n)
  for (i in 1:n){
    if(i != k){
      A_k <- A[,i]*v[i]
      temp <- temp + A_k
    }
  }
  final <- y - temp
  return(final)
}

#' Solve trend-filtering by coordinate descent on dual problem
#' 
#' @param y response
#' @param k order of differencing matrix
#' @param v Initial dual variable
#' @param lambda regularization parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
trend_filter_cd <- function(y, k, v, lambda=10, max_iter=100, tol=1e-3) {
  converged = FALSE
  n = length(y)
  
  #set the step size
  D <- myGetDkn(k, n)
  dd<-nuclear_norm(D)
  t<-1/nuclear_norm(D)
  A <- t(D)
  
  x <- matrix(0,nrow=n, ncol=max_iter+1)#c(0,max_iter+1)
  f_x <- c(0,max_iter+1)
  
  x[,1] <-v
  f_x[1] <- dual(y,v,D,lambda)
  
  #to keep track of relative changes, kkt residual change, and duality gaps
  rel_f <-  c(0,max_iter)
  rel_x <-  c(0,max_iter)
  kkt_res <- c(0,max_iter)
  dual_gap <- c(0,max_iter)
  
  cc<-2
  while(converged == FALSE){
    iterate <- n
    for(k in iterate){
      #compute r_k
      r_k <- compute_rk(y, x[1:n,cc], D, k)
      
      #compute D_k/||D||_k (transposed)
      A_k<-A[1:n,k]
      a <- sqrt(sum(A_k^2))
      v_a <- A_k / a # normalized vector
      
      #compute z_k = <r_k, v_a>
      z_k <- r_k %*% v_a
      
      #update kth v term
      x[k,cc] <- project(z_k,lambda)
    }
    b <- y - (t(D)%*%x[,cc])
    theta<- D%*%b
    
    next_x <- x[1:n,cc]
    f_x[cc] <- dual(y,next_x,D,lambda)
    rel_f[cc] <- sqrt((f_x[cc] - f_x[cc-1])^2)  
    rel_x[cc] <-  sqrt(sum( (x[,cc] - x[,cc-1])^2))
    kkt_res[cc] <- kkt_residual(y, b, theta, v, D, lambda) 
    dual_gap[cc] <- duality_gap(y, b, x[,cc], D, lambda)
    
    #compute criterion & if satisfied
    if (rel_f[cc] < tol){converged=TRUE}
    cc<- cc+1
  }
  x_final <- x[1:n,cc]
  return(list(x_final,f_x, rel_f, rel_x, kkt_res, dual_gap))
}


#' Solve trend-filtering by proximal gradient on dual problem
#' 
#' @param y response
#' @param k order of differencing matrix
#' @param v Initial dual variable
#' @param lambda regularization parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
trend_filter_pg <- function(y, k, v, lambda=10, max_iter=100, tol=1e-3) {
  
  converged = FALSE
  n = length(y)
  
  #set the step size
  D <- myGetDkn(k, n)
  t<-1/nuclear_norm(D)
  A <- t(D)
  
  x <- matrix(0,nrow=n, ncol=max_iter+1)#c(0,max_iter+1)
  f_x <- c(0,max_iter+1)
  
  x[,1] <- v
  f_x[1] <- dual(y,v,D,lambda)
  
  #to keep track of relative changes, kkt residual change, and duality gaps
  rel_f <-  c(0,max_iter)
  rel_x <-  c(0,max_iter)
  kkt_res <- c(0,max_iter)
  dual_gap <- c(0,max_iter)
  
    
  #set the step size
  t<-1/nuclear_norm(D)
  
  cc<-2
  while(converged == FALSE){
    DD_inv <-inv( D%*%t(D))
    x_k <- DD_inv %*% ( x[,cc] - D %*%y)
    temp <- x[1:n,cc] - t*x_k
    
    for(x in temp){
      x[k,cc] <- project(x,lambda)}
    
    b <- y - (t(D)%*%x[,cc])
    theta<- D%*%b
    
    next_x <- x[,cc]
    f_x[cc] <- dual(y,next_x,D,lambda) #dual <- function(y, v, D, lambda)
    rel_f[cc] <- f_x[cc] - f_x[cc-1]  
    rel_x[cc] <-  x[cc] - x[cc-1]
    kkt_res[cc] <- kkt_residual(y, b, theta, v, D, lambda) 
    dual_gap[cc] <- duality_gap(y, b, x[,cc], D, lambda)
    cc<- cc+1
    if(cc==maxiter){converged = TRUE}
  }
  x_final <- x[,cc]
  length(x_final)
  
  return(list(x_final,f_x, rel_f, rel_x, kkt_res, dual_gap))
}




