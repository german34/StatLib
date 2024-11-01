#' Gradient Step
#'
#' @param gradf handle to function that returns gradient of objective function #' @param x current parameter estimate
#' @param t step-size
#' @export
gradient_step <- function(gradf, t, x) {
  x_new <-x - (t* gradf(x))
  return (x_new )
  }


#' Gradient Descent (Fixed Step-Size)
#'
#' @param fx handle to function that returns objective function values
#' @param gradf handle to function that returns gradient of objective function #' @param x0 initial parameter estimate
#' @param t step-size
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
gradient_descent_fixed <- function(fx, gradf, t, x0, max_iter=1e2, tol=1e-3) {
  'Returns:
  The final iterate value
  The objective function values
  The 2-norm of the gradient values
  The relative change in the function values
  The relative change in the iterate values'
  eps=0.03
  x <- matrix(0,nrow=2, ncol=max_iter+1) #c(0,max_iter+1)
  f_x <- vector("numeric", max_iter+1)
  x[,1] <-x0
  f_x[1] <- fx(x0)

  df_l2 <- c(0,max_iter+1)
  df_l2[1]< sum(gradf(x0)^2) #|| * ||_2^2

  rel_f <-  c(0,max_iter)
  rel_x <-  c(0,max_iter)

  #while not converged
  converged=FALSE
  i=2
  while(converged==FALSE){ #for( i in 2:maxiter+1){
    #x_new <- gradient_step(gradf, t, x_old)
    x[,i]<- gradient_step(gradf, t, x[,i-1])
    #print('new step')
    #print(x[,i])

    #f_new <- fx(x_new)
    f_x[i] <- fx(x[,i])

    #df_new <- gradf(x_new)
    df_new <- gradf(x[,i])
    df_l2[i]<- norm(df_new,type="2")^2

    f_diff <- f_x[i] - f_x[i-1]
    x_diff <- x[,i] - x[,i-1]
    rel_f[i-1] <-  norm(f_diff,type="2")^2 /norm(f_x[i],type="2")^2 #c(0,max_iter)
    rel_x[i-1] <-  norm(x_diff,type="2")^2 /norm(x[i],type="2")^2 #c(0,max_iter)
    if( i==max_iter ){converged=TRUE} #other stopping crietrions (df_l2[i]< eps*df_l2[1]) || (df_l2[i]<tol
    i=i+1
  }

  x_final <- x[,max_iter+1]
  return(list(x_final,f_x, df_l2, rel_f, rel_x))
}



#' Backtracking
#'
#' @param fx handle to function that returns objective function values
#' @param x current parameter estimate
#' @param t current step-size
#' @param df the value of the gradient of objective function evaluated at the current x
#' @param alpha the backtracking parameter
#' @param beta the decrementing multiplier
#' @export
#'
backtrack <- function(fx, t, x, df, alpha=0.5, beta=0.9,counter=1) {
  t_k <- t
  #f(x) - gamma*alpha_k * ||df||_2^2
  df_temp <- norm(df , type="2")^2
  fval <- fx(x)
  max_iter = 40
  converged=FALSE

  iter=1

  while (converged==FALSE){
    t_k = beta*t_k
    x_temp <- x - t_k*df
    res = fval - alpha*t_k*df_temp #for checking armijo condition
    new_fx <- fx(x_temp)
    if ((res > new_fx)||(max_iter<iter)){converged=TRUE} #checking...
    iter<-iter+1
  }
  #if(verbose ==TRUE){print(res)}
  #if(counter>=28){print(res)}#verbose=TRUE}
  #print('***final_step_size')
  #print(t_k)
  return(t_k)
}


#' Gradient Descent (Backtracking Step-Size)
#'
#' @param fx handle to function that returns objective function values
#' @param gradf handle to function that returns gradient of objective function #' @param x0 initial parameter estimate
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
gradient_descent_backtrack <- function(fx, gradf, x0, max_iter=1e2, tol=1e-3) {
  #print('Max iteration: ')
  #print(max_iter)
  x <- matrix(0,nrow=2, ncol=max_iter+1)#c(0,max_iter+1)
  f_x <- c(0,max_iter+1)
  converged = FALSE

  x[1:2,1] <- x0

  f_x <- c(0,max_iter+1)
  x[,1] <-x0
  t <- 1
  f_x[1] <- fx(x0)

  df_l2 <- vector("numeric",max_iter+1)
  df_l2[1]<-sum(gradf(x0)^2) #|| * ||_2^2

  rel_f <-  c(0,max_iter)
  rel_x <-  c(0,max_iter)

  #while not converged
  i=2
  while(converged==FALSE){ #for( i in 2:max_iter+1){
    #print('iterate ---- ')
    #print(i-1)

    df_xk <- gradf(x[,i-1])
    t_k <-backtrack(fx, t, x[,i-1], df_xk, alpha=0.5, beta=0.9,counter=i)

    #x_new <- gradient_step(gradf, t, x_old)
    x[,i]<- gradient_step(gradf, t_k, x[,i-1])
    #print(x[,i])

    #f_new <- fx(x_new)
    f_x[i] <- fx(x[,i])

    #df_new <- gradf(x_new)
    df_new <- gradf(x[,i])
    df_l2[i]<- norm(df_new,type="2")^2

    f_diff <- f_x[i] - f_x[i-1]
    x_diff <- x[,i] - x[,i-1]
    rel_f[i-1] <-  norm(f_diff,type="2")^2 /norm(f_x[i],type="2")^2 #c(0,max_iter)
    rel_x[i-1] <-  norm(x_diff,type="2")^2 /norm(x[i],type="2")^2 #c(0,max_iter)
    if( max_iter+1==i ){converged=TRUE} #other stopping crietrions (df_l2[i]< eps*df_l2[1]) || (df_l2[i]<tol
    i=i+1
  }
  #print('Final iterate....')
  x_final <- x[,max_iter+1]
  #print(x_final)

  return(list(x_final,f_x, df_l2, rel_f, rel_x))
}


#' Gradient Descent
#'
#' @param fx handle to function that returns objective function values
#' @param gradf handle to function that returns gradient of objective function #' @param x0 initial parameter estimate
#' @param t step-size
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
gradient_descent <- function(fx, gradf, x0, t=NULL, max_iter=1e2, tol=1e-3) {

  if(is.null(t) == TRUE) #then use backtrack
  {
    gradient_descent_backtrack(fx, gradf, x0, max_iter=max_iter, tol=tol)
  }
  else{
    #otherwise use fixed
    gradient_descent_fixed(fx, gradf, t, x0, max_iter=max_iter, tol=tol)
  }
}


#' Objective Function for Logistic Regression
#'
#' @param y binary response
#' @param X design matrix
#' @param beta regression coefficient vector
#' @param lambda regularization parameter
#' @export
fx_logistic <- function(y, X, beta, lambda=0) {

  lb <- 0
  for (i in 1:length(y)){
    y_i <- -y[i]
    x_i = X[i,]
    temp <- x_i[1]*beta[1] + x_i[2]*beta[2]
    comp1 <- y_i * temp
    comp2 <- log(1 + exp(temp))
    lb <- lb + comp1 + comp2
  }
  if (lambda > 0){
    comp3 <- (lambda/2.0)* norm(beta,type="2")^2
    lb <- lb+t(comp3)
  }
  return(t(lb))
}


#' Gradient for Logistic Regession
#'
#' @param y binary response
#' @param X design matrix
#' @param beta regression coefficient vector
#' @param lambda regularization parameter
#' @export
gradf_logistic <- function(y, X, beta, lambda=0) {

  n=nrow(X)
  p=ncol(X)
  #print(X[1,1:2])

  lb <-vector("numeric",p)#c(0,p)
  for (i in 1:n){
    y_i <- -y[i]
    x_i = X[i,]
    temp <- x_i[1]*beta[1] + x_i[2]*beta[2] #temp <- x_i %*% beta
    #print(temp)
    #print('----')
    mult <- y_i + 1/(1+exp(-temp))
    lb<-lb+(mult*t(x_i))
  }
  #print(lb)

  if (lambda > 0){
    comp3 <- lambda * beta
    #print(lb)
    #print(comp3)
    #print(comp3)
    lb <- lb+t(comp3)
  }
  return(t(lb))
}









