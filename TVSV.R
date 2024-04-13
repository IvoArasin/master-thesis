#################################################################################
### DNS with time-varying shape and volatility (Koopman et al. 2010) in R     ###
### Ivo L. Arasin                                                             ###
### ivo.arasin@aol.com ivolovis.arasin@student.kuleuven.be                    ###
### March 2024                                                                ###
### Disclaimer: No warranty of any kind                                       ###
#################################################################################
library(readxl)
data <- read_excel("/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/zero_rates08til23.xlsx")
m <- c(1/365, 7/365, 14/365, 1/12, 2/12, 3/12, 6/12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15)
N <- 4
scaler <- 1
data <- as.matrix(subset(data, select=-c(1, 21, 22)))/scaler

# Factor Loadings for observation-matrix B
factor_loadings <- function(l,m)
{
  column1 <- rep.int(1,length(m))
  column2 <- (1 - exp(-l * m))/(l * m)
  column3 <- (1 - exp(-l * m) - l*m*exp(-l*m))/(l*m)
  
  B <- cbind(column1,column2,column3)
  B
}

# parameters gathering meta-information during MLE estimation
print_param_counter <<- 0
total_iterations <<- 0

# Dynamic Nelson Siegel model with time-varying shape and volatility
KalmanFilterTVSV <- function(para,Y,lik=TRUE, h=0) {
  # print parameters every 100th iteration
  if(print_param_counter >=100){
    print(para)
    print_param_counter <<- 0
  }
  print_param_counter <<- print_param_counter +1
  total_iterations <<- total_iterations +1
  
  W <- ncol(Y)
  T <- nrow(Y)
  if(is.null(W) | is.null(T)){
    W <- length(Y)
    T <- 1
  }
  
  # Create vectors and matrices
  mu	<- matrix(NA,N,1) # Mean vector
  phi<- diag(N+1) # Vector Autoregressive coeffient matrix VAR(1)	
  H	<- diag(W) # Variance matrix of residuals
  Q	<- diag(N+1) # Transition covariance matrix of residuals
  
  # Assign measurement-noise matrix
  H_start_ <- 13
  H_end_ <- 13 + W - 1
  H <- diag(para[H_start_:H_end_])^2
  
  # Transition Matrix A (=Phi)
  phi[1,1] <- para[1] # level autoregressive factor
  phi[2,2] <- para[2] # slope autoregressive factor
  phi[3,3] <- para[3] # curvature autoregressive factor
  phi[4,4] <- para[4] # shape parameter autoregressive factor
  phi[5,5] <- 0
  
  # Ensure stationarity of autoregressive factors
  diagA <- eigen(phi[1:4,1:4])
  diagA_eigenvalues <- diagA$values
  if(max(abs(Re(diagA_eigenvalues)))>=1){
    print("Phi has Eigenvalue >1")
    return(10000000000000)
    #break
  }
  
  # Mean vector
  mu[1]<-para[5] # level
  mu[2]<-para[6] # slope
  mu[3]<-para[7] # curvature
  mu[4]<-para[8] # shape parameter
  
  # Process noise
  Q[1,1] <- para[9]  # level
  Q[2,2] <- para[10]  # slope
  Q[3,3] <- para[11] #curvature
  Q[4,4] <- para[12] # shape parameter
  Q[5,5] <- 0
  
  # Take square of process-noise (thereby also ensuring positivity)
  Q <- Q %*% t(Q) 
  
  # Gamma vector
  Gamma <- abs(para[(H_end_+1):(H_end_+W)])
  
  gamma0 <- 0.0001 #fixed for identifiability
  gamma1 <- para[length(para)-1]
  gamma2 <- para[length(para)]
  
  # Set up empty arrays to save recursive parameters
  v1   <- matrix(NA,T,W)			  
  v2   <- matrix(NA,T,W)
  a.tt <- matrix(NA, T, N+1)
  a.t  <- matrix(NA, (T+1), N+1)
  P.tt <- array(NA, c(T, N+1, N+1))
  P.t  <- array(NA, c((T+1), N+1, N+1))
  h_t <- rep(T, 0)
  
  # Start state vector and variance matrix
  a.t[1, 1] <- mu[1] # mean of level
  a.t[1, 2] <- mu[2] # mean of slope
  a.t[1, 3] <- mu[3] # mean of curvature
  a.t[1, 4] <- mu[4] # mean of shape parameter 
  a.t[1, 5] <- 0
  
  # Start variance matrix
  lyapunov<-function(N,phi,Q){
    matrix(solve(diag(N^2) - kronecker(phi,phi)) %*% matrix(Q,(N^2),1),N,N)
  }
  tryCatch(
    {
      cannot_be_done <- 0
      P.t[1,1:4,1:4] <-lyapunov(N=N,phi=phi[1:4,1:4],Q=Q[1:4,1:4]) # Start variance matrix. Pt0
    },
    warning = function(w){
      print(w)
      print("Difficult to invert lyapunov!")
      cannot_be_done <- 1
      logLik<- -1000000000000000
      
    },
    error=function(e){
      print(e)
      print("lyapunov could not be inverted!")
      cannot_be_done <- 1
      logLik<- -1000000000000000
      return(1000000000000000)
    },
    finally={
      if(cannot_be_done==1){
        return(1000000000000000)
      }
    }
  )
  
  # Initialize h_1
  P.t[1, 5, 5] <- gamma0/(1-gamma1-gamma2)
  
  P.t[1,5, 1:4 ] <- rep(0, 4)
  P.t[1,1:4, 5 ] <- rep(0, 4)
  
  # Initial log-likelihood	
  logLik <- 0
  
  for (t in 1:T) 
  {
    # Gamma Restrictions
    if(gamma1<=0 || gamma2<=0 || gamma1+gamma2>=1){
      print(paste("Gamma error: ", (gamma1+gamma2)))
      logLik<- -1000000000000000
      break
      
    }else{
      
      NS_factor_term <- factor_loadings(a.t[t,4], m) %*% a.t[t, 1:3] # Gamma isn't considered here because the
      # predicted value of the common volatility term is in zero in expectation, which would result in
      # Gamma always being multiplied with zero, meaning its omission is just more efficient
      
      # This if-clause allows the function to run on one single observation
      if(T==1){
        Y_data <- as.numeric(Y)
      }else{
        Y_data <- (as.numeric(Y[t, ]))
      }
      
      # Prediction error
      v <- Y_data - NS_factor_term
      
      # Jacobian, i.e. first derivative of the observation model
      jacobian_EKF <- function(beta2, beta3, l,m){
        column1 <- rep.int(1,length(m))
        column2 <- (1 - exp(-l * m))/(l * m)
        column3 <- (1 - exp(-l * m) - l*m*exp(-l*m))/(l*m)
        column4 <- (( beta2*exp(-m*l)*(m*l-exp(m*l)+1) + beta3*exp(-m*l)*(m**2*l**2+m*l-exp(m*l)+1))/(m*l**2))
        jac_dataframe <- cbind(column1,column2,column3, column4)
        jac_dataframe
      }
      
      jac <- jacobian_EKF(a.t[t, 2], a.t[t, 3], a.t[t, 4], m)
      
      jac_w_Gamma <- cbind(jac, Gamma)
      F <- jac_w_Gamma %*% P.t[t, ,] %*% t(jac_w_Gamma) + H
      
      detF <- det(F)
      # Ensure invertibility of F
      if(is.na(detF) || is.nan(detF) || is.infinite(detF) || abs(detF)<1e-150 || detF<0){
        print(paste("Determinant of F can't be determined on iteration ",t))
        logLik<- -1000000000000000
        break
        
      }
      else{
        tryCatch(
          {
            cannot_be_done <<- 0
            
            F.inv  <- solve(F)
          },
          warning = function(w){
            #print(w)
            print(paste("Difficult to invert F on iteration ",t))
            cannot_be_done <<- 1
            logLik<- -1000000000000000
          },
          error=function(e){
            #print(e)
            print(paste("Can't invert F on iteration ",t))
            cannot_be_done <<- 1
            logLik<- -1000000000000000
          },
          finally={
            if(cannot_be_done==1){
              break
            }
          }
        )
        
        # Log-likelihood
        logLik <- logLik-0.5*(length(v))*log(2*pi)-0.5*log(detF)-0.5*t(v)%*%F.inv%*%v # constructed via the prediction error decomposition
        
        # Updating the state vector and its variance matrix
        a.tt[t, ]   <- a.t[t, ] +  P.t[t, , ] %*% t(jac_w_Gamma) %*% F.inv %*% v
        P.tt[t, , ] <- P.t[t, , ] - P.t[t, , ] %*% t(jac_w_Gamma) %*% F.inv %*% jac_w_Gamma %*% P.t[t, , ]
        
        # Filtered errors
        NS_factor_term <- factor_loadings(a.tt[t,4], m) %*% a.tt[t, 1:3]
        NS_factor_term_w_Gamma <- cbind(NS_factor_term, Gamma)
        v1[t, ]	<- NS_factor_term + Gamma * a.tt[t, 5]
        v2[t, ] <- Y_data - v1[t, ]
        
        # Predicting the state vector and its estimation error covariance matrix
        a.t[t + 1, ]  <- phi %*% a.tt[t, ] + (diag(N+1) - phi) %*% c(mu, 0)
        
        h_t[t] <- Q[5,5]
        
        # Update common volatility component h_t via GARCH process
        Q[5,5] <- gamma0 + gamma1*(a.tt[t, 5]^2+P.tt[t, 5, 5]) + gamma2 * h_t[t]
        
        P.t[t + 1, ,] <- phi %*% P.tt[t, ,] %*% t(phi) + diag(5) %*% Q %*% diag(5)
      }
    }
  }
  
  # Forecasting
  y_t_h <- matrix(NA, h, W)
  if(t > T-1 & h>0){
    recurstion_term <- diag(N+1)
    for(i in 1:h){
      recurstion_term <- recurstion_term + ifelse(i>1, 1, 0)*phi**(i-1)
      a_t_pred <- recurstion_term%*%(diag(N+1)-phi)%*%c(mu, 0) + phi**(i)%*%a.tt[T, ]
      y_t_h[i, ] <- factor_loadings(a_t_pred[4], m)%*%a_t_pred[1:3]
    }
  }
  
  if(lik)
  {
    as.numeric(-logLik)
  }else{
    return(list(a.tt=a.tt,a.t=a.t,P.tt=P.tt,P.t=P.t,h_t=h_t,v2=v2,v1=v1, y_t_h=y_t_h, negLogLik=-logLik))
  }
}

########################################################################################
##########     Model Estimation
########################################################################################

TVSV_twoStepInitGuess <- c(
  0.97, 0.98, 0.94, 0.95,
  
  1.9, -1.5, -2.3, 0.35,
  
  0.2, 0.2, 0.2, 0.2,
  
  0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
  0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
  
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  0.4, 0.4)
# Manually optimize model
TVSV_model_optim <- optim(TVSV_twoStepInitGuess,KalmanFilterTVSV,control = list(REPORT=5), Y=data, lik=TRUE, h=0, method="BFGS")
optim_paras <- TVSV_model_optim$par

# Reconstruct fitted yields
TVSV_model <- KalmanFilterTVSV(optim_paras, Y=data, lik=FALSE, h=12)
time_ <- 178
recon <- factor_loadings(TVSV_model$a.tt[time_,4], m)%*%TVSV_model$a.tt[time_,1:3]+optim_paras[33:52]*TVSV_model$a.tt[time_,5]
plot(m, data[(time_),], type="l")
lines(m, recon, type="l", lty=2, col="blue")

# Automatically optimize model in recursive fashion (Iterative Nelder-Mead)
myDirectory = "/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/optimized_files/RollingWindowForcasts"
automated_rollingWindow(para_init=TVSV_twoStepInitGuess, model=KalmanFilterTVSV, maxiter=2, maxFuncEvals=50, data=data, directoryPath=myDirectory, method_="Nelder-Mead")

# This function can be used to recursively optimize expanding windows
automated_rollingWindow <- function(para_init, model, start_=1, end_=0, data=data, maxiter=10, maxFuncEvals=100000, directoryPath="", method_="Nelder-Mead"){
  starting_values <- para_init
  start <- start_
  if(end_==0){
    end_ <- nrow(data)
  }
  for(i in 1:length(end_)){
    end <- end_[i]
    window_data <- data[start:end,]
    starting_values <- automatic_optimization(starting_values, model, window_data, maxiter=maxiter, window_nr=i, directoryPath=directoryPath, maxFuncEvals=maxFuncEvals, method_=method_)
    preds <- model(para=starting_values, Y=window_data, lik=FALSE,h=12)
    preds <- preds$y_t_h
    preds_h12 <- preds[12,]
    preds_h6 <- preds[6,]
    preds_h1 <- preds[1,]
    predictions <- rbind(preds_h12, preds_h6, preds_h1)
    myDirectory = directoryPath
    write.csv(predictions, paste0(myDirectory,"predictions_window_", i, "DNS.csv"), row.names=TRUE)
  }
}

automatic_optimization <- function(para_init, model, input_data, maxiter=10, window_nr=0, directoryPath="", maxFuncEvals=maxFuncEvals, method_="Nelder-Mead"){
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
  if(method_=="BFGS"){
    method_ <- "BFGS"
    maxiter <- 1
  }
  while(abs(prevLogLike-logLike)>0.1 & i<maxiter){
    i<-i+1
    
    # optimize model parameters
    optim_values<-optim(parameters,model,control = list(trace=1, maxit=maxFuncEvals), Y=input_data, method=method_)
    
    #convergence
    optim_statistics[i,1] <- optim_values$convergence
    
    #iterations
    optim_statistics[i,2] <- total_iterations
    
    # model parameters
    model_parameters[i, ] <- optim_values$par
    parameters <- optim_values$par
    
    # NegLogLike
    prevLogLike <- logLike
    logLike <- model(para=parameters, Y=input_data, lik=TRUE)
    optim_statistics[i,3] <- logLike
    
    # filtered RMSE
    model_output <- model(para=parameters, Y=input_data, lik=FALSE)
    
    filtered_error_RMSE[i, ] <- sqrt(colMeans((model_output$v2)^2))
    
    total_iterations <<- 0
    
    myDirectory = directoryPath
    fileName <- paste0(myDirectory,"TVSV_intermediate_estimate.csv")
    intermediate_result <- list("optim_statistics"=optim_statistics, "model_parameters"=model_parameters, "filtered_error_RMSE"=filtered_error_RMSE)
    write.csv(intermediate_result, file=fileName, row.names = TRUE)
  }
  fileName_finalOutput <- paste0(myDirectory,"TVSV_", window_nr, "_final.csv")
  result <- list("optim_statistics"=optim_statistics, "model_parameters"=model_parameters, "filtered_error_RMSE"=filtered_error_RMSE)
  write.csv(result, fileName_finalOutput, row.names=TRUE)
  return(optim_values$par)
}

# Function to convert MLE parameters between differently scaled data
para_scaler <- function(para){
  scaled_para <- para
  scaled_para[5:7] <- para[5:7]/100
  scaled_para[9:11] <- para[9:11]/100
  scaled_para[13:32] <- para[13:32]/100
  scaled_para[33:52] <- para[33:52]/100
  return(scaled_para)
}
