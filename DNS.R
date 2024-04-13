library(readxl)
# import data
data <- read_excel("/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/zero_rates08til23.xlsx")
# assign maturity vector / vector of tenors
m <- c(1/365, 7/365, 14/365, 1/12, 2/12, 3/12, 6/12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15)
# specify the number of latent variables your state-space model uses
N <- 3
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

# Dynamic Nelson Siegel model
independent_DNS <- function(para,Y,lik=TRUE, forecast_horizon=0) {
  if(print_param_counter >=100){
    print(para)
    print_param_counter <<- 0
  }
  print_param_counter <<- print_param_counter +1
  total_iterations <<- total_iterations <<- 0 +1
  T <- nrow(Y)
  W <- ncol(Y)
  l <- para[1]
  pars<-list()

  # Create vectors and matrices
  mu	<- matrix(NA,N,1) # Mean vector
  phi<- diag(N) # Create transition-Matrix A (also commonly referred to as Phi)
  H	<- diag(W) # Create measurement-noise matrix
  Q	<- diag(N) # Create process-noise matrix
  
  # Create and assign factor loading matrix B
  B	<- factor_loadings(l,m) 
  
  # Assign measurement-noise matrix
  H_end_ <- length(para)
  H_start_ <- H_end_ - (W-1)
  H <- diag(para[H_start_:H_end_])^2
  #H <- diag(para[11:30])^2
  
  # Assign transition-matrix A
  phi[1,1] <- para[2]
  phi[2,2] <- para[3]
  phi[3,3] <- para[4]
  
  # Assign mean vector Mu (= \mu)
  mu[1]<-para[5]
  mu[2]<-para[6]
  mu[3]<-para[7]
  
  # Assign process-noise matrix
  Q[1,1] <- para[8]
  Q[2,2] <- para[9]
  Q[3,3] <- para[10]
  
  # Square process noise matrix (also ensures non-negativity)
  Q <- Q %*% t(Q) 
  
  #Set up dataframes for parameter collection
  v1   <- matrix(NA,T,W)
  v2   <- matrix(NA,T,W)
  a.tt <- matrix(NA, T, N) # Filtered state vector
  a.t  <- matrix(NA, (T+1), N) # Predicted state vector
  P.tt <- array(NA, c(T, N, N)) # Filtered estimation-error matrix
  P.t  <- array(NA, c((T+1), N, N)) # Predicted estimation-error matrix
  
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
  
  # Start actual Kalman filter recursion
  for (t in 1:T) 
  {
    # Prediction Error
    v <- (as.numeric(Y[t, ])) - B %*% a.t[t, ]
    
    F <- B %*% P.t[t, ,] %*% t(B) + H
    
    detF <- det(F)
    
    # Some safety checks to ensure F can be inverted
    if(detF<=1e-150 || is.na(detF) || is.nan(detF) || is.infinite(detF) || detF<0){
      print("break reached in KFilter")
      logLik<- -1000000000000000; break
    }
    else{
      # Try Catch block to prevent model estimation from being interrupted by an error
      # If inversion of F fails, return a high likelihood so optimizer knows it's headed the wrong way
      cannot_be_done <<- 0
      tryCatch(
        {
          F.inv  <- solve(F)
        },
        warning = function(w){
          print(w)
          print("Difficult to invert F!")
          cannot_be_done <<- 1
        },
        error=function(e){
          print(e)
          print("F could not be inverted!")
          cannot_be_done <<- 1
        },
        finally={
          if(cannot_be_done==1){
            return(10000000000000)
            break
          }
        }
        
      )
      
      # Calculate Log-Likelihood
      logLik <- logLik-0.5*(W)*log(2*pi)-0.5*log(detF)-0.5*t(v)%*%F.inv%*%v
      }
    # Updating the state vector and its estimation error-matrix
    a.tt[t, ]   <- a.t[t, ] +  P.t[t, , ] %*% t(B) %*% F.inv %*% v
    P.tt[t, , ] <- P.t[t, , ] - P.t[t, , ] %*% t(B) %*% F.inv %*% B %*% P.t[t, , ]
    
    # Filtered estimates and filtered errors
    v1[t, ]	<- B %*% a.tt[t, ]
    v2[t, ] <- (as.numeric(Y[t, ])) - B %*% a.tt[t, ]
    
    # Predicting the state vector and estimation error-matrix
    a.t[t + 1, ]  <- phi %*% a.tt[t, ] + (diag(N) - phi) %*% mu  
    P.t[t + 1, ,] <- phi %*% P.tt[t, ,] %*% t(phi) + Q
  }

  # Forecasting
  y_t_h <- matrix(NA, forecast_horizon, W)
  if(t > T-1 & forecast_horizon>0){
    recurstion_term <- diag(N)
    for(i in 1:forecast_horizon){
      recurstion_term <- recurstion_term + ifelse(i>1, 1, 0)*phi**(i-1)
      a_t_pred <- recurstion_term%*%(diag(N)-phi)%*%mu + phi**(i)%*%a.tt[T, ]
      y_t_h[i, ] <- B%*%a_t_pred
    }
  }
  
  if(lik)
  {
    as.numeric(-logLik)
  }else{
    return(list(a.tt=a.tt,a.t=a.t,P.tt=P.tt,P.t=P.t,v2=v2,v1=v1, y_t_h=y_t_h, negLokLig=-logLik))
  }
}

################################################################
# Automatic optimization with multiple self-initiated restarts #
################################################################

para_init <- c(
  0.72,
  0.98, 0.97, 0.96,
  2.2, -1.9, -3.3,
  0.2, 0.2, 0.2,
  
  0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
  0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
  
)

# Single optimization run
optimized_DNS_model_object <- optim(para_init,independent_DNS,control = list(trace=1, maxit=100000), Y=data, method="BFGS")
optim_para <- optimized_DNS_model_object$par

# Apply model 
DNS_model_fitted <- independent_DNS(para=optim_para, Y=data, lik=FALSE, forecast_horizon=12)

# Reconstruct fitted yields
time_ <- 130
B <- factor_loadings(optim_para[1],m)
recon <- B%*%DNS_model_fitted$a.tt[time_,1:3]
plot(m,data[time_,], type="l")
lines(m,recon, type="l", col="blue", lty=2)

# Automatically optimize model in recursive fashion (Iterative Nelder-Mead)
myDirectory = "/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/optimized_files/RollingWindowForcasts"
automated_rollingWindow(para_init=para_init, model=independent_DNS, maxiter=15, maxFuncEvals=2000, data=data, directoryPath=myDirectory, method_="BFGS")

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
    preds <- model(para=starting_values, Y=window_data, lik=FALSE,forecast_horizon=12)
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
    fileName <- paste0(myDirectory,"DNS_calmTime_window_intermediate.csv")
    intermediate_result <- list("optim_statistics"=optim_statistics, "model_parameters"=model_parameters, "filtered_error_RMSE"=filtered_error_RMSE)
    write.csv(intermediate_result, file=fileName, row.names = TRUE)
  }
  fileName_finalOutput <- paste0(myDirectory,"DNS_calmTime_window_", window_nr, "_final.csv")
  #result <- list("model_parameters"=model_parameters)
  result <- list("optim_statistics"=optim_statistics, "model_parameters"=model_parameters, "filtered_error_RMSE"=filtered_error_RMSE)
  write.csv(result, fileName_finalOutput, row.names=TRUE)
  return(optim_values$par)
}
