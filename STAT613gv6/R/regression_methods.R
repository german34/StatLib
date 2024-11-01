#' Plot restrictions
#'
#' \code{ridge_regression} returns the ridge regression coefficient estimates
#' for a sequence of regularization parameter values
#'
#' @param y response variable
#' @param X design matrix
#' @param lambda vector of tuning parameters
#' @export
#'

#usethis::use_package(caret)
#usethis::use_tidy_description(caret)

get_ridge_regression <- function(y, X, lambda){
  #First grab num of lambda, cols and rows for X
  n_lambda = length(lambda)
  n=nrow(X)
  p=ncol(X)

  #Take SVD of X
  ss <- svd(X)
  dd <- ss$d
  D <- diag(dd)
  U <- ss$u
  V <- ss$v
  DtD <- (as.matrix(t(D))%*%as.matrix(D))

  U_t <- t(U)
  V_t <- t(V)

  #create empty matrix for each ridge estimate (correspond to columns)
  mat1 <- matrix(0, nrow=p, ncol=n_lambda)
  #temp <- matrix(0,nrow=n,ncol=n)

  for (i in 1:n_lambda){
    #choose correct lambda
    temp_lambda = lambda[i]
    #print(temp_lambda)

    #print(temp_lambda)
    I_n <- temp_lambda*n*diag(n)
    #print(I_n[1:3,1:3])

    #Compute inverse of (D^t D + lambda I)
    M <- (DtD + I_n)#as.matrix(DtD) + as.matrix(I_n)

    M_inv = solve(M)
    UtY <- U_t %*% y#as.matrix(U_t) %*% as.matrix(y)

    V_M_inv  <- V %*% M_inv %*% t(D) #as.matrix(V) %*% as.matrix(M_inv) %*% t(D)



    #for(j in 1:n){
    #  b_temp <- dd[j]^2/(dd[j]^2 + temp_lambda) * UtY[i]*
    #}
    #compute B for this lambda
    #B_temp <- V %*% M_inv %*% t(D) %*% t(U)
    #B_final<- B_temp %*% y
    B_final <- V_M_inv %*% UtY #as.matrix(V_M_inv) %*% as.matrix(UtY)
    #print(length(B_final))

    #print('----')
    #add as a column to B
    mat1[,i] <- B_final

  }
  #return B_hat for each lambda
  return (mat1)
}


#' Generalized Cross Validation
#'
#' \code{gcv} returns the leave-one-out
#' for a sequence of regularization parameter values. #'
#' @param y response variables
#' @param X design matrix
#' @param lambda vector of tuning parameters
#' @export

gcv <- function(y, X, lambda){
  #First check how many lambda
  n_lambda = length(lambda)
  n=nrow(X)
  p=ncol(X)

  #Take SVD of X
  ss <- svd(X)
  dd <- ss$d
  D <- diag(dd)
  U <- ss$u
  V <- ss$v
  U_t <- t(U)
  V_t <- t(V)

  #for each lambda we will compute the GCV and return the minimum of all (or each??)
  gcv_temp = c(0, n_lambda) #vector("numeric",n_lambda)

  for (i in 1:n_lambda){
    #choose correct lambda
    temp_lambda = lambda[i]
    dof = 0

    #First compute S_lambda (using SVD decomp)
    #i.e.-> U D ( D^tD + lambda n I )^{-1} D^t U^t
    M <- (t(D)%*%D + temp_lambda*n*diag(n))
    M_inv = solve(M)
    S_lambda <- U %*% D %*% M_inv %*% t(D) %*% t(U)
    #as.matrix(U) %*% as.matrix(D) %*% as.matrix(M_inv) %*% as.matrix(t(D)) %*% as.matrix(t(U))


    #second compute y hat
    y_hat <- S_lambda %*% y#as.matrix(y)
    #print(y_hat[1])

    #compute degrees of freedom (a.k.a. take trace of S_lambda)
    dof <- sum(diag(S_lambda))#tr(S_lambda)
    #print(dof)
    #print(y)
    #Take the difference between y_true and y_hat
    res <- y-y_hat

    #compute gcv
    total <- 0.0
    for (j in 1:length(res)){
      temp <- (res[j]/(1-(dof/n)))^2
      total <- total + temp
    }
    total <- total / n
    gcv_temp[i] <- total

    #"faste"r way to compute sum
    #total <- (res/ (1-(dof/n)  )^2
    #gcv_temp[i] <- (1/n) * sum(total)

  }
  #return(min(gcv_temp))
  return(gcv_temp)
}



#' Leave One Out
#'
#' \code{leave_one_out} returns the leave-one-out prediction error estimates #' for a sequence of regularization parameter values.
#'
#' @param y response variables
#' @param X design matrix
#' @param lambda vector of tuning parameters
#' @param B ridge regression coefficient estimates for seq of reg/tuning parameters
#' @export
#'
#'
leave_one_out <- function(y, X, lambda) {
  n_lambda = length(lambda)
  print("Length of lambda is ....")
  print(n_lambda)
  n=nrow(X)
  p=ncol(X)

  #Take SVD of X
  ss <- svd(X)
  dd <- ss$d
  D <- diag(dd)
  U <- ss$u
  V <- ss$v

  B = get_ridge_regression(y, X, lambda)

  #create vector to hold loo(nu) for each nu
  LOO <- vector("numeric",n_lambda)
  #print(length(LOO))

  for (i in 1:n_lambda){
    temp_lambda <- lambda[i]

    #***compute h(lambda) diagonals
    #I_n <-temp_lambda *n* diag(n)
    #temp <- (t(D)%*%D + I_n)
    #D_inv = solve(temp)
    #H_lambda <- U %*% D %*% D_inv %*% t(D) %*% t(U)
    #H_lambda <- diag(d_row) - H_lambda

    #compute y_hat = X * B(nu)
    y_hat <- X %*% B[,i] # as.matrix(X) %*% as.matrix(B[,i])

    #compute y_k - y_k^{-k}(nu)
    res <- y - y_hat#[,i]

    total <- 0.0


    for (k in 1:n){
      #k residual value
      num_k <- res[k,1]

      #compute h_k(nu)
      hh_k <- 0.0
      for(j in 1:n){
        hh_k <- hh_k + (dd[j]^2/(dd[j]^2+temp_lambda*n) * (U[k,j]^2))
      }

      #adding each individual (y_k - y^{-k}_k)^2
      total <- total + (num_k/(1 - hh_k))^2
    }
    #print(total/n)
    LOO[i] = total/n
  }
  return(LOO)
}

#' K-fold Cross Validation
#'
#' \code{k_fold_cv} returns K-fold cross-validation prediction error estimates #' for a sequence of regularization parameter values.
#'
#' @param y response variables
#' @param X design matrix
#' @param k number of folds
#' @param lambda vector of tuning parameters
#' @param seed seed for random number generator
#' @export

k_fold_cv <- function(y, X, lambda, k=5, seed=12345) {
  set.seed(seed) #<- Always set a seed for reproducibility
  n_lambda = length(lambda)

  #Split into different sets using create folds
  flds <- createFolds(y, k, list = TRUE, returnTrain = TRUE)


  avg_err <- vector("numeric",n_lambda) #vector("numeric",n_lambda)
  for(j in 1:k){
      #First obtain subset:
      indx <- unlist(flds[j])#,use.names=FALSE)

      #specify training set, i.e. -C_k
      y_k <- y[indx]
      X_k <- X[indx,]

      #specify testing set, i.e., C_k
      y_test<-y[-indx]
      X_test<-X[-indx,]

      #First compute B^{ - C_k } for with respect to lambdas
      B_k = get_ridge_regression(y_k, X_k, lambda)


      for(i in 1:n_lambda) {
        #Compute X{C_k} * B^{- C_k}
        y_pred <- X_test %*% B_k[,i]#as.matrix(X_test) %*% as.matrix(B_k[,i])

        #take difference & check for nans just in case
        res <- y_test - y_pred
        if (any(is.na(res))){break}
        #temp <- sum(res^2)/dim(X_k)[1]
        #print(sum(res^2))
        #print(norm(res,type="2")^2)
        avg_err[i] = avg_err[i]+(norm(res,type="2")^2/dim(X_test)[1])
        #print(avg_err[i])
      }
  }

  errors <- avg_err / k
  #print(errors)
  return(errors)
}




