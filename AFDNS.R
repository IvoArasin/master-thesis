#################################################################################
### Arbitrage-Free Dynamic Nelson Siegel Model (Christensen et al. 2009) in R ###
### Ivo L. Arasin                                                             ###
### ivo.arasin@aol.com ivolovis.arasin@student.kuleuven.be                    ###
### March 2024                                                                ###
### Disclaimer: No warranty of any kind                                       ###
#################################################################################
library(readxl)
require(xts)
library(expm)
# Import Data
data <- read_excel("/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/zero_rates08til23.xlsx")
# Define maturity vector (this one is normalized to years, not months, as is commonly done)
m <- c(1/365, 7/365, 14/365, 1/12, 2/12, 3/12, 6/12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15)
N <- 3
dt <- 1/12
scaler <- 100 # Scale data such that 1=100% and 0.01=1%
data <- as.matrix(subset(data, select=-c(1, 21, 22)))/scaler
train_data <- data[1:152,] # 152/191 is approximately 80%, such that you have an 80/20 train/test-split

# Factor Loadings
factor_loadings <- function(l=0.72,m=c(1/365, 7/365, 14/365, 1/12, 2/12, 3/12, 6/12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15))
{
  column1 <- rep.int(1,length(m))
  column2 <- (1 - exp(-l * m))/(l * m)
  column3 <- (1 - exp(-l * m) - l*m*exp(-l*m))/(l*m)
  
  lambmat <- cbind(column1,column2,column3)
  return(lambmat)
} 

# Yield Adjustment Term (see Christensen 2009 p.20 at the bottom)
yield_adjustment_term <- function(sigma11, sigma22, sigma33, lambda, maturity){
  term1 <- -sigma11^2*((maturity^2)/(6))
  term2 <- -sigma22^2*(1/(2*lambda^2)-(1/lambda^3)*((1-exp(-lambda*maturity))/maturity)+(1)/(4*lambda^3)*((1-exp(-2*lambda*maturity))/(maturity)) )
  term3 <- -sigma33^2*( (1)/(2*lambda^2) + (1)/(lambda^2)*exp(-lambda*maturity)-(1)/(4*lambda)*maturity*exp(-2*lambda*maturity)-
                          (3)/(4*lambda^2)*exp(-2*lambda*maturity)-(2)/(lambda^3)*(1-exp(-lambda*maturity))/(maturity)+(5)/(8*lambda^3)*(1-exp(-2*lambda*maturity))/(maturity))
  y_adj <- term1+term2+term3
  return(y_adj)
}

initialize <- function(K_P, Sigma, delta_t){
  error <- FALSE
  K_P_eigen <- eigen(K_P)
  
  K_P_eigenvalues <- K_P_eigen$values
  K_P_eigenvectors <- K_P_eigen$vectors

  # Ensure K_P remains stationary by negating all values that lie outside the unit circle or are negative
  # For reference, see Christensen 2009, p.14 footnote
  if(min(K_P_eigenvalues)<0){
    print("Real part of Eigenvalue not positive")
    error <- TRUE
  }
  
  # Compute conditional covariance matrix Q and unconditional variance P 
  # For reference, see  Eq.21 & Eq.22 in Caldeira et al. (2016) with title
  # "Forecasting the yield curve with the arbitrage-free dynamic Nelson–Siegel model: Brazilian evidence"
  inv_K_P_eigenvectors <- solve(K_P_eigenvectors)
  Omega <- inv_K_P_eigenvectors%*%Sigma%*%t(Sigma)%*%t(inv_K_P_eigenvectors)
  
  Q_middle_part <- matrix(NA, N,N)
  P_init_middle_part <- matrix(NA, N,N)
  for(i in 1:N){
    for(j in 1:N){
      eigenSum <- (K_P_eigenvalues[i]+K_P_eigenvalues[j])
      P_init_middle_part[i,j] <- (Omega[i,j]/eigenSum)
      Q_middle_part[i,j] <- Omega[i,j]*(1-exp(-eigenSum*delta_t))/eigenSum
    }
  }

  P_init <- K_P_eigenvectors%*%P_init_middle_part%*%t(K_P_eigenvectors)
  Q <- K_P_eigenvectors%*%Q_middle_part%*%t(K_P_eigenvectors)
  
  list("P_init"=P_init, "Q"=Q, "error"=error)
}

total_iterations <<- 0
iter_count <<- 0
KalmanFilter_AFDNS <- function(parameters, dataset=train_data, returnOnlyLikelihood=TRUE, N=3, maturity_vector=m, h=0){
  iter_count <<- iter_count + 1
  if(iter_count>=100){
    print(parameters)
    iter_count <<- 0
  }
  total_iterations <<- total_iterations + 1
  error <- FALSE
  T <- nrow(dataset)
  nrCols <- ncol(dataset)
  y <- dataset
  
  # Define static variables
  theta <- c(NA, NA, NA)
  K_P <- diag(3)
  Sigma <- diag(3)
  Q <- diag(3)
  H <- diag(nrCols)
  B <- matrix(NA, nrCols, N)
  yield_adj <- rep(NA, nrCols)
  delta_t <- 1/12
  logLikelihood <- 0
  
  # Assign static variables
  l <- parameters[1]
  
  # convert Phi1 parameters back to K_P matrix parameters
  K_P1 <- -12*log(parameters[2])
  K_P2 <- -12*log(parameters[3])
  K_P3 <- -12*log(parameters[4])
  K_P <- diag(c(K_P1, K_P2, K_P3))
  
  #K_P <- diag(parameters[2:4]) # thsese lines that are commented out are just an alternative specification. Results are the same
  # but the method with commented out lines seems slower (just based on my own observation)
  theta <- parameters[5:7]
  Sigma <- abs(diag(parameters[8:10]))
  
  H <- diag(parameters[11:length(parameters)]^2)
  yield_adj <- yield_adjustment_term(sigma11=Sigma[1,1], sigma22=Sigma[2,2], sigma33=Sigma[3,3], lambda=l, m=maturity_vector)
  B <- factor_loadings(l, maturity_vector)
  
  Phi1 <- diag(parameters[2:4])
  #Phi1 <- expm(-K_P*delta_t)
  Phi1_eigen <- eigen(Phi1)
  Phi1_eigenvals <- Phi1_eigen$values

  # Ensure autoregressive process is stationary by raising error if transition matrix has unit root (i.e. eigenvalue >= 1)
  # For reference, see Christensen 2009, p.14 footnote. Here "Phi1" corresponds to "A" in the footnote
  if(max(abs(Re(Phi1_eigenvals)))>=1){
    print("Phi1 has Eigenvalue >1")
    return(10000000000000)
    break
  }
  
  Phi0 <- (diag(N) - Phi1)%*%theta
  
  result_object <- initialize(K_P, Sigma, delta_t)
  error <- result_object$error
  if(error){
    print("K_P is not positive definite!")
    return(10000000000000)
  }
  Q <- result_object$Q
  
  # Define recursive KF variables
  x.t <- matrix(NA, T+1, N)
  x.tt <- matrix(NA, T, N)
  P.t <- array(NA, c(N,N,T+1))
  P.tt <- array(NA, c(N,N,T))
  prediction <- matrix(NA, T, nrCols)
  filtered_pred <- matrix(NA, T, nrCols)
  v <- matrix(NA, T, nrCols)
  v1 <- matrix(NA, T+1, nrCols)
  v2 <- matrix(NA, T, nrCols)
  
  P.t[,,1] <- result_object$P_init
  x.t[1,] <- theta
  
  Q_eigen <- eigen(Q)
  Q_eigvals <- Q_eigen$values
  
  H_eigen <- eigen(H)
  H_eigvals <- H_eigen$values
  
  # Ensure matrices are positive (semi-)definite, i.e. have real-valued positive entries on their diagonal
  # Constraint is not strictly necessary since values for H and Q are squared and are forced to be absolute (e.g. abs(x))
  if(max(H - t(H))> 1e-9 | any(Q_eigvals < -1e-9) | any(H_eigvals < -1e-9)){
    return(10000000000)
    break
  }
  
  # Kalman Filter Recursion
  for(t in 1:T){
    
    prediction[t,] <- B%*%x.t[t,] + yield_adj
    v[t,] <- as.numeric(y[t,]) - prediction[t,]
    
    F <- B%*%P.t[,,t]%*%t(B)+H
    
    tryCatch({
      detF <- det(F)
      log_detF <- log(detF)
      invF <- solve(F)
    }, warning = function(w) {
      print("Warning with manipulation of F")
      return(10000000000000)
      break
    }, error = function(e) {
      print("Error with manipulation of F")
      return(10000000000000)
      break
    }, finally = {
      # Update Step
      
      x.tt[t,] <- x.t[t,] + P.t[,,t]%*%t(B)%*%invF%*%v[t,]
      P.tt[,,t] <- P.t[,,t] - P.t[,,t]%*%t(B)%*%invF%*%B%*%P.t[,,t]
      
      filtered_pred[t,] <- B%*%x.tt[t,] + yield_adj
      v2[t,] <- as.numeric(y[t,]) - filtered_pred[t,]
      
      # Prediction Step
      x.t[t+1,] <- Phi0+Phi1%*%x.tt[t,]
      P.t[,,t+1] <- Phi1%*%P.tt[,,t]%*%t(Phi1)+Q
      
      # Evaluate Log Likelihood 
      logLikelihood <- logLikelihood + (-0.5*nrCols*log(2*pi)-0.5*log_detF-0.5*t(v[t,])%*%invF%*%v[t,])
      
    })
    
  }
  
  # Forecasting
  y_t_h <- matrix(NA, h, ncol(data))
  if(h>0){
    recurstion_term <- diag(3)
    #xtt <- x.tt[T,]
    for(i in 1:h){
      #xtt <- Phi0+Phi1%*%xtt
      
      # Forecasting as described in Christensen 2009 p.24
      recurstion_term <- recurstion_term + ifelse(i>1, 1, 0)*Phi1**(i-1)
      x.t_pred <- recurstion_term%*%Phi0 + Phi1**(i)%*%x.tt[T,]

      # measurement-equation for prediction
      y_t_h[i, ] <- B%*%x.t_pred+yield_adj
    }
    #y_t_h[i, ] <- B%*%xtt+yield_adj
  }
  
  if(returnOnlyLikelihood){
    -logLikelihood
  }
  else{
  list("x.t"=x.t, "x.tt"=x.tt, "P.t"=P.t, "P.tt"=P.tt,
       "v"=v, "v1"=v1, "v2"=v2, "NeglogLikelihood"=-logLikelihood,
       "Q"=Q, "H"=H, "K_P"=K_P, "theta"=theta, "B"=B, "Sigma"=Sigma, "Phi0"=Phi0, "Phi1"=Phi1, "yield_adj"=yield_adj, "prediction"=prediction, "filtered_pred"=filtered_pred, "y_t_h"=y_t_h)
  }
}

para_init <- c(
  0.35,
  0.99, 0.98, 0.95,
  0.03, -0.018, -0.0255,
  0.01, 0.01, 0.01,
  
  0.01,  0.01,  0.01,  0.01,  0.01,0.01,  0.01,  0.01,  0.01,  0.01,
  0.01,  0.01,  0.01,  0.01,  0.01,0.01,  0.01,  0.01,  0.01,  0.01
)
################################################################
# Automatic optimization with multiple self-initiated restarts #
################################################################
automatic_optimization_test <- automatic_optimization(para_init, KalmanFilter_AFDNS, train_data, maxiter=30)

automatic_optimization <- function(para_init, model, input_data, maxiter=10){
  maxiter <- maxiter
  parameter_length <- length(para_init)
  model_parameters <- matrix(NA, maxiter, parameter_length)
  optim_statistics <- matrix(NA, maxiter, 3)
  colnames(optim_statistics) <- c("convergence", "nrOfIterations", "NegLogLike")
  filtered_error_RMSE <- matrix(NA, maxiter, ncol(input_data))
  parameters <- para_init
  
  prevLogLike <- 0
  logLike <- 2
  i<-0
  while(abs(prevLogLike-logLike)>0.1 & i<maxiter){
    i<-i+1
    
    # optimize model parameters
    optim_values<-optim(parameters,model,control = list(trace=1, maxit=10000))
    
    #convergence
    optim_statistics[i,1] <- optim_values$convergence
    
    #iterations
    optim_statistics[i,2] <- total_iterations
    
    # NegLogLike
    prevLogLike <- logLike
    logLike <- model(para=parameters, dataset=input_data, return=TRUE, h=0)
    optim_statistics[i,3] <- logLike
    
    # model parameters
    model_parameters[i, ] <- optim_values$par
    parameters <- optim_values$par
    
    # filtered RMSE
    model_output <- model(parameters=parameters, dataset=input_data, returnOnlyLikelihood=FALSE, h=0)
    
    filtered_error_RMSE[i, ] <- sqrt(colMeans((model_output$v2*100)^2)) # multiplied with 100 due to scaling for AFDNS data
    
    total_iterations <<- 0
    intermediate_result <- list("optim_statistics"=optim_statistics, "model_parameters"=model_parameters, "filtered_error_RMSE"=filtered_error_RMSE)
    write.csv(intermediate_result, "intermediate_result.csv", row.names = TRUE)
  }
  
  result <- list("model_parameters"=model_parameters)
  result <- write.csv(result, "optimization_file_AFDNS.csv", row.names=TRUE)
}

# basic optimization
#afns_model_optim <- optim(par=para_init, fn=KalmanFilter_AFDNS, dataset=data[1:152,], returnOnlyLikelihood=TRUE, N=3, maturity_vector=m, control=list(trace=1, maxit=5000))


# sandbox: some tests and (visual) checks as well as forecasting analysis
param_set <- as.matrix(read.csv2("/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/optimized_files/RollingWindowForcasts/AFNDS/AFDNS_152.csv"))
params <- as.numeric(param_set[10,5:34])

# Apply model with optimized parameters
AFNS_object <- KalmanFilter_AFDNS(parameter=reinit, dataset=data[1:172,], returnOnlyLikelihood=FALSE, N=3, maturity_vector=m, h=12)
AFNS_object$NeglogLikelihood

round(colMeans(abs(AFNS_object$v2))*100,6)
plot(m,AFNS_object$yield_adj*100,type="l")
plot(AFNS_object$x.t[,1], type="l")

round(colMeans(abs(AFNS_object$v))*100,6)
plot(m, colMeans(AFNS_object$filtered_pred)*100, type="l", col="green")
lines(m, colMeans(data[1:152,])*100, type="p", col="blue")
AFNS_object$theta
AFNS_object$Phi0

time <- 152
ahead <- 12
preds <- KalmanFilter_AFDNS(parameter=reinit, dataset=data[1:time,], returnOnlyLikelihood=FALSE, N=3, maturity_vector=m, h=12)
maxyield <- max(data[(time+ahead),]*100)
minyield <- min(data[(time+ahead),]*100)
plot(m, data[(time+ahead),]*100, ylim=c(-1,3))
lines(m, data[time,]*100, type="l", col="red")
lines(m,preds$y_t_h[ahead,]*100, type="l", col="blue")

#randomWalk <- as.numeric(sqrt(((data[(time+ahead),]-data[(time),])*100)^2))
pred_RMSE <- as.numeric(sqrt(((data[(time+ahead),]-preds$y_t_h[ahead,])*100)^2))
for(i in 1:length(pred_RMSE)){
  cat(pred_RMSE[i], "\n", sep="")
}
