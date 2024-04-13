library(readxl)
# import data
data <- read_excel("/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/zero_rates08til23.xlsx")
# assign maturity vector / vector of tenors
m <- c(1/365, 7/365, 14/365, 1/12, 2/12, 3/12, 6/12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15)
# specify the number of latent variables your state-space model uses
N <- 4
# if you want to scale the data, set the scaler
scaler <- 1
# remove some unwanted columns (e.g. date column and the two interpolated columns)
data <- as.matrix(subset(data, select=-c(1, 21, 22)))/scaler

# Factor Loading Matrix B
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

# Dynamic Nelson Siegel model with time-varying shape parameter (DNS-TVS)
independent_DNS_TVS <- function(para,Y,lik=TRUE, h=0) {
  if(print_param_counter >=100){
    print(para)
    print_param_counter <<- 0
  }
  print_param_counter <<- print_param_counter +1
  total_iterations <<- total_iterations +1
  W <- ncol(Y)
  T <- nrow(Y)
  
  # Create vectors and matrices
  mu	<- matrix(NA,N,1) # Mean vector
  phi<- diag(N) # Create transition-Matrix A (also commonly referred to as Phi)
  H	<- diag(W) # Create measurement-noise matrix
  Q	<- diag(N) # Create process-noise matrix
  
  # Vector autoregressive coeffient matrix: VAR(1)
  phi[1,1] <- para[1]
  phi[2,2] <- para[2]
  phi[3,3] <- para[3]
  phi[4,4] <- para[4]
  
  diagA <- eigen(phi)
  diagA_eigenvalues <- diagA$values
  
  # Ensure autoregressive components of latent factors remains stationary 
  if(max(abs(Re(diagA_eigenvalues)))>=1){
    print("Phi has Eigenvalue >1")
    return(10000000000000)
  }
  
  # Assign mean vector Mu (= \mu)
  mu[1:4]<-para[5:8]
  
  # Assign process-noise matrix
  Q[1,1] <- para[9]
  Q[2,2] <- para[10]
  Q[3,3] <- para[11]
  Q[4,4] <- para[12]
  
  # Square process noise matrix (also ensures non-negativity)
  Q <- Q %*% t(Q) 
  
  # Assign measurement-noise matrix
  H_end_ <- length(para)
  H_start_ <- H_end_ - (W-1)
  H <- diag(para[H_start_:H_end_])^2
  
  #Set up dataframes for parameter collection
  v1   <- matrix(NA,T,W)			  
  v2   <- matrix(NA,T,W)
  a.tt <- matrix(NA, T, N)
  a.t  <- matrix(NA, (T+1), N)
  P.tt <- array(NA, c(T, N, N))
  P.t  <- array(NA, c((T+1), N, N))
  
  # Initialize predicted state vector with state-means Mu (=\mu)
  a.t[1, ]  <- mu
  
  # Initialize estimation-error matrix
  lyapunov<-function(N,phi,Q){
    matrix(solve(diag(N^2) - kronecker(phi,phi)) %*% matrix(Q,(N^2),1),N,N)
  }  
  
  tryCatch(
    {
      cannot_be_done <- 0
      P.t[1, ,] <-lyapunov(N=N,phi=phi,Q=Q)
    },
    warning = function(w){
      print(w)
      print("Difficult to evaluate lyapunov!")
      cannot_be_done <- 1
      #logLik<- -1000000000000000
      return(1000000000000000)
      
    },
    error=function(e){
      print(e)
      print("Lyapunov could not be evaluated!")
      cannot_be_done <- 1
      #logLik<- -1000000000000000
      return(1000000000000000)
    },
    finally={
      if(cannot_be_done==1){
        break
      }
    }
    
  )
  
  # Initial log-likelihood	
  logLik <- 0
  
  # Start actual Extended Kalman filter recursion
  for (t in 1:T) 
  { 
    # Prediction Error
    NS_factor_term <- factor_loadings(exp(a.t[t,4]), m) %*% a.t[t, 1:3]
    v <- (as.numeric(Y[t, ])) - NS_factor_term
    
    # Jacobian Matrix as defined in original paper by Koopman et al. (2010)
    jacobian <- function(beta2, beta3, l,m){
      column1 <- rep.int(1,length(m))
      column2 <- (1 - exp(-l * m))/(l * m)
      column3 <- (1 - exp(-l * m) - l*m*exp(-l*m))/(l*m)
      # column4 is multiplied once more with "l" since we estimate the log of the shape parameter to ensure positivity
      column4 <- (beta2*((exp(-m*l))/(l)-(1-exp(-m*l))/(m*l**2))+beta3*(-(exp(-m*l)*(exp(m*l)-m**2*l**2-m*l-1))/(m*l**2)))*l
      jac_dataframe <- cbind(column1,column2,column3, column4)
      jac_dataframe
    }
    jac_w_params <- jacobian(a.t[t, 2], a.t[t, 3], exp(a.t[t, 4]), m)
    
    F <-  jac_w_params%*% P.t[t, ,] %*% t(jac_w_params) + H
    
    detF <- det(F)
    
    # Some safety checks to ensure F can be inverted
    if(is.na(detF) || is.nan(detF) || is.infinite(detF) || abs(detF)<1e-150 || detF<0){
      print(detF)
      print("detF cannot be determined")
      logLik<- -1000000000000000
      break
    }
    else{
      cannot_be_done <<- 0
      tryCatch(
        {
          F.inv  <- solve(F)
        },
        warning = function(w){
          print(w)
          print("Difficult to invert F!")
          cannot_be_done <<- 1
          logLik<- -1000000000000000
          
        },
        error=function(e){
          print(e)
          print("F could not be inverted!")
          cannot_be_done <<- 1
          logLik<- -1000000000000000
        },
        finally={
          if(cannot_be_done==1){
            return(10000000000000)
          }
        }
        
      )
      
      # Calculate Log-Likelihood
      logLik <- logLik-0.5*(length(v))*log(2*pi)-0.5*log(detF)-0.5*t(v)%*%F.inv%*%v # constructed via the prediction error decomposition
      
      # Kalman Gain
      Kalman_Gain <- P.t[t, , ] %*% t(jac_w_params) %*% F.inv
      
      # Updating the state vector and its estimation error-matrix
      a.tt[t, ]   <- a.t[t, ] +  Kalman_Gain %*% v
      P.tt[t, , ] <- P.t[t, , ] - Kalman_Gain %*% jac_w_params %*% P.t[t, , ]
      
      # Filtered estimates and filtered errors
      # The shape parameter estimate is exponentiated because we estimate its log(.) to ensure positivity
      v1[t, ]	<- factor_loadings(exp(a.tt[t,4]), m) %*% a.tt[t, 1:3]
      v2[t, ] <- (as.numeric(Y[t, ])) - v1[t, ]
      
      # Predicting the state vector and estimation error-matrix
      a.t[t + 1, ]  <- phi %*% a.tt[t, ] + (diag(N) - phi) %*% mu  
      P.t[t + 1, ,] <- phi %*% P.tt[t, ,] %*% t(phi) + Q
    }
  }
  
  # Forecasting
  y_t_h <- matrix(NA, h, W)
  y_t_h_fix <- matrix(NA, h, W)
  if(t > T-1 & h>0){
    recurstion_term <- diag(N)
    for(i in 1:h){
      recurstion_term <- recurstion_term + ifelse(i>1, 1, 0)*phi**(i-1)
      a_t_pred <- recurstion_term%*%(diag(N)-phi)%*%mu + phi**(i)%*%a.tt[T, ]
      
      # Forecasting all latent parameters
      y_t_h[i, ] <- factor_loadings(exp(a_t_pred[4]), m)%*%a_t_pred[1:3]
      
      # Forecasting only level, slope and curvature while keeping the shape parameter
      # fixed at its last known value
      y_t_h_fix[i, ] <- factor_loadings(exp(a.tt[T,4]), m)%*%a_t_pred[1:3]
    }
  }
  
  if(lik)
  {
    as.numeric(-logLik)
  }else{
    return(list(a.tt=a.tt,a.t=a.t,P.tt=P.tt,P.t=P.t,v2=v2,v1=v1, y_t_h=y_t_h, y_t_h_fix=y_t_h_fix, negLogLik=-logLik))
  }
}


# To convert optimized parameters between estimations using differently scaled data
para_scaler <- function(parameters){
  parameters_scaledDown <- parameters
  parameters_scaledDown[5:7] <- parameters_scaledDown[5:7]/100
  parameters_scaledDown[9:11] <- parameters_scaledDown[9:11]/100
  parameters_scaledDown[13:32] <- parameters_scaledDown[13:32]/100
  parameters_scaledDown
}

para_init3 <- c( # This one works very well along with the bounds
  0.98, 0.97, 0.93, 0.95,
  3.2, -1.6, -2.4, -1.05,
  0.26, 0.31, 0.5, 0.2,
  
  0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
  0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
)


lower_ <- c(0.1, 0.1, 0.1, 0.1,
            2, -10, -10, -5,
            0.01, 0.01, 0.01, 0.01,
            
            0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0)

upper_ <- c(0.999999, 0.999999, 0.999999, 0.999999,
            5, 0, 0, 3,
            1, 1, 1, 1,
            
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

TVS_lbfgsb_calmPeriod_window_6 <- optim(TVS_twoStepFullSample_inits,independent_DNS_TVS,Y=data[25:145,],lik=TRUE,control = list(REPORT=2), method="L-BFGS-B", lower=lower_, upper=upper_)
TVS_lbfgsb_calmPeriod_window_6$par


###############################################################
###############################################################
###############################################################

# Automatic Window Estimation
para_init_TVS <- c(
  0.97, 0.98, 0.94, 0.95,
  
  1.9, -1.5, -2.3, log(0.35),
  
  0.2, 0.2, 0.2, 0.4,
  
  0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
  0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)

# Single optimization run
optimized_TVS_model_object <- optim(para_init_TVS,independent_DNS_TVS,control = list(trace=1, maxit=100000), Y=data, method="BFGS")
optim_para <- optimized_TVS_model_object$par

# Apply model 
DNS_model_fitted <- independent_DNS_TVS(para=optim_para, Y=data, lik=FALSE, h=12)

# Reconstruct fitted yields
time_ <- 90
reconstruction <- factor_loadings(exp((DNS_model_fitted$a.tt[time_,4])),m)%*% DNS_model_fitted$a.tt[time_,1:3]
par(mfrow=c(1,1))
plot(m, data[time_,], type="l", main="DNS-TVV")
lines(m, reconstruction, type="l", lty=2, col="blue")

# Automatically optimize model in recursive fashion (Iterative Nelder-Mead)
myDirectory = "/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/optimized_files/RollingWindowForcasts"
automated_rollingWindow(para_init=para_init_TVS, model=independent_DNS_TVS, maxiter=15, maxFuncEvals=2000, data=data, directoryPath=myDirectory, method_="BFGS")

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
    preds <- model(para=starting_values, Y=window_data, lik=FALSE, forecast_horizon=12)
    preds <- preds$y_t_h
    preds_h12 <- preds[12,]
    preds_h6 <- preds[6,]
    preds_h1 <- preds[1,]
    predictions <- rbind(preds_h12, preds_h6, preds_h1)
    myDirectory = directoryPath
    write.csv(predictions, paste0(myDirectory,"predictions_window_", i, "TVS.csv"), row.names=TRUE)
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
    fileName <- paste0(myDirectory,"TVS_calmTime_window_intermediate.csv")
    intermediate_result <- list("optim_statistics"=optim_statistics, "model_parameters"=model_parameters, "filtered_error_RMSE"=filtered_error_RMSE)
    write.csv(intermediate_result, file=fileName, row.names = TRUE)
  }
  fileName_finalOutput <- paste0(myDirectory,"TVS_calmTime_window_", window_nr, "_final.csv")
  #result <- list("model_parameters"=model_parameters)
  result <- list("optim_statistics"=optim_statistics, "model_parameters"=model_parameters, "filtered_error_RMSE"=filtered_error_RMSE)
  write.csv(result, fileName_finalOutput, row.names=TRUE)
  return(optim_values$par)
}
