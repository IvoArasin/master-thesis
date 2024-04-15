library(readxl)
# import data
data <- read_excel("/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/zero_rates08til23.xlsx")
# create maturity vector
m <- c(1/365, 7/365, 14/365, 1/12, 2/12, 3/12, 6/12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15)
N <- 3
scaler <- 1
data <- as.matrix(subset(data, select=-c(1, 21, 22)))/scaler

# Factor Loading matrix B
factor_loadings <- function(l,m)
{
  column1 <- rep.int(1,length(m))
  column2 <- (1 - exp(-l * m))/(l * m)
  column3 <- (1 - exp(-l * m) - l*m*exp(-l*m))/(l*m)
  
  B <- cbind(column1,column2,column3)
  B
} 

# Hyper parameters to create some statistics during optimization
print_param_counter <<- 0
total_iterations <<- 0

# DNS with time-varying common volatility component
independent_DNS_TVV <- function(para,Y,lik=TRUE, h=0) {
  if(print_param_counter >=100){
    print(para)
    print_param_counter <<- 0
  }
  print_param_counter <<- print_param_counter +1
  total_iterations <<- total_iterations +1
  W <- ncol(Y)
  T <- nrow(Y)
  l<- para[1]
  
  # Create vectors and matrices
  mu	<- matrix(NA,N,1) # Mean vector
  phi<- diag(N+1) # Vector Autoregressive coeffient matrix VAR(1)	
  H	<- diag(W) # Variance matrix of residuals
  Q	<- diag(N+1) # Transition covariance matrix of residuals
  
  # Loading matrix
  B	<- factor_loadings(l,m) 
  
  # Variance matrix of residuals
  start_ <- 11
  end_ <- 11 + W - 1
  H <- diag(para[start_:end_])^2
  
  # Vector autoregressive coeffient matrix: VAR(1)
  phi[1,1] <- para[2]
  phi[2,2] <- para[3]
  phi[3,3] <- para[4]
  phi[4,4] <- 0
  
  # Mean vector
  mu[1]<-para[5]
  mu[2]<-para[6]
  mu[3]<-para[7]
  
  # Transition covariance matrix of residuals
  Q[1,1] <- para[8]
  Q[2,2] <- para[9]
  Q[3,3] <- para[10]
  Q[4,4] <- 0
  
  # Take Square to ensure positivity
  Q <- Q %*% t(Q) 
  
  # Gamma vector
  start_ <- 11 + W
  end_ <- 11 + W -1
  Gamma <- abs(para[start_:end_])
  
  B <- cbind(B, Gamma)
  
  gamma0 <- 0.0001 #fixed for identifiability
  end_ <- length(para)
  
  gamma1 <- para[end_-1]
  gamma2 <- para[end_]
  
  # Set up arrays to store recursive variables
  a.tt <- matrix(NA, T, N+1)
  a.t  <- matrix(NA, (T+1), N+1)
  P.tt <- array(NA, c(T, N+1, N+1))
  P.t  <- array(NA, c((T+1), N+1, N+1))
  v1   <- matrix(NA,T,W)			  
  v2   <- matrix(NA,T,W)
  
  # Start state vector and variance matrix
  a.t[1, 1:3] <- mu
  a.t[1, 4] <- 0
  
  # INitialize estimation error covariance matrix
  lyapunov<-function(N,phi,Q){
    matrix(solve(diag(N^2) - kronecker(phi,phi)) %*% matrix(Q,(N^2),1),N,N)
  }  
  
  tryCatch(
    {
      cannot_be_done <- 0
      P.t[1,1:3,1:3] <-lyapunov(N=N,phi=phi[1:3,1:3],Q=Q[1:3,1:3])
    },
    warning = function(w){
      print(w)
      print("Difficult to evaluate lyapunov!")
      cannot_be_done <- 1
      logLik<- -1000000000000000
      return(1000000000000000)
      
    },
    error=function(e){
      print(e)
      print("Lyapunov could not be evaluated!")
      cannot_be_done <- 1
      logLik<- -1000000000000000
      return(1000000000000000)
    },
    finally={
      if(cannot_be_done==1){
        break
      }
    }
    
  )
  
  
  P.t[1, 4, 4] <- gamma0/(1-gamma1-gamma2)
  
  P.t[1,4, 0:3 ] <- rep(0, 3)
  P.t[1,0:3, 4 ] <- rep(0, 3)
  
  h_t <- rep(T, 0)
  
  # Initial log-likelihood	
  logLik <- 0
  
  for (t in 1:T) 
  {
    # Restrictions on Gamma
    if(gamma0 <= 0 || gamma1<=0 || gamma2<=0 || gamma1+gamma2>=1){
      print("Gamma error")
      logLik<- -1000000000000000;
      return(1000000000000000)
      break
    }
    v <- (as.numeric(Y[t, ])) - B %*% a.t[t, ] # prediction error vector
    
    F <- B %*% P.t[t, ,] %*% t(B) + H # prediction error variance matrix
    
    detF <- det(F)
    # Ensure invertibility of F
    if(is.na(detF) || is.nan(detF) || is.infinite(detF) || abs(detF)<1e-190 || detF<0){
      print("detF cannot be determined")
      logLik<- -1000000000000000
      return(1000000000000000)
      break
    }
    else{
      tryCatch(
        {
          cannot_be_done <<- 0
          F.inv  <- solve(F)
        },
        warning = function(w){
          print(w)
          print("Difficult to invert F!")
          cannot_be_done <<- 1
          logLik<- -1000000000000000
          return(1000000000000000)
          
        },
        error=function(e){
          print(e)
          print("F could not be inverted!")
          cannot_be_done <<- 1
          logLik<- -1000000000000000
          return(1000000000000000)
        },
        finally={
          if(cannot_be_done==1){
            return(1000000000000000)
            break
          }
        }
      )
      
      # Log-likelihood
      logLik <- logLik-0.5*(length(v))*log(2*pi)-0.5*log(detF)-0.5*t(v)%*%F.inv%*%v
      
      # Updating the state vector and its variance matrix
      a.tt[t, ]   <- a.t[t, ] +  P.t[t, , ] %*% t(B) %*% F.inv %*% v
      P.tt[t, , ] <- P.t[t, , ] - P.t[t, , ] %*% t(B) %*% F.inv %*% B %*% P.t[t, , ]
      
      # Filtered errors
      v1[t, ]	<- B %*% a.tt[t, ]
      v2[t, ] <- (as.numeric(Y[t, ])) - B %*% a.tt[t, ]
      
      # Predicting the state vector and its variance matrix
      a.t[t + 1, ]  <- phi %*% a.tt[t, ] + (diag(N+1) - phi) %*% c(mu, 0)
      
      h_t[t] <- Q[4,4]
      
      # GARCH process for common volatility term
      Q[4,4] <- gamma0 + gamma1*(a.tt[t, 4]^2+P.tt[t, 4, 4]) + gamma2 * h_t[t]
      
      P.t[t + 1, ,] <- phi %*% P.tt[t, ,] %*% t(phi) + diag(4) %*% Q %*% diag(4)
    }
  }
  # Forecasting
  y_t_h <- matrix(NA, h, W)
  if(t > T-1 & h>0){
    recurstion_term <- diag(N+1)
    for(i in 1:h){
      recurstion_term <- recurstion_term + ifelse(i>1, 1, 0)*phi**(i-1)
      a_t_pred <- recurstion_term%*%(diag(N+1)-phi)%*%c(mu, 0) + phi**(i)%*%a.tt[T, ]
      y_t_h[i, ] <- B%*%a_t_pred
    }
  }
  
  if(lik)
  {
    as.numeric(-logLik)
  }else{
    return(list(a.tt=a.tt,a.t=a.t,P.tt=P.tt,P.t=P.t,h_t=h_t,v2=v2,v1=v1, y_t_h=y_t_h, negLogLik=-logLik))
  }
}

##################################################################
##################################################################
##################################################################

para_init_TVV <- c( # derived from two-step method
  0.36,
  0.98, 0.97, 0.93,
  3.2, -1.6, -2.4,
  0.26, 0.31, 0.5,
  
  0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
  0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
  
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  
  0.4, 0.4)

# Automatically optimize model in recursive fashion (Iterative Nelder-Mead)
myDirectory = "/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/optimized_files/RollingWindowForcasts"
automated_rollingWindow(para_init=para_init_TVV, model=independent_DNS_TVV, maxiter=15, maxFuncEvals=2000, data=data, directoryPath=myDirectory, method_="BFGS")

para_optim <- as.matrix(read.csv("/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/optimized_files/RollingWindowForcastsTVV_1_final.csv"))
params_ <- para_optim[1,5:56]
TVV_model_object <- independent_DNS_TVV(para=params_, Y=data, lik=FALSE, h=12)

# Plot common volatility variance h_t
plot(TVV_model_object$h_t[10:191], type="l", main=expression(paste("Common Volatility ", h[t])))

# plot Gamma weights
plot(m, params_[31:50],
     type="b",
     main=expression(paste("6 Common Volatility Weights ", Gamma[epsilon])))

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
    write.csv(predictions, paste0(myDirectory,"predictions_window_", i, "TVV.csv"), row.names=TRUE)
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
    fileName <- paste0(myDirectory,"TVV_intermediate.csv")
    intermediate_result <- list("optim_statistics"=optim_statistics, "model_parameters"=model_parameters, "filtered_error_RMSE"=filtered_error_RMSE)
    write.csv(intermediate_result, file=fileName, row.names = TRUE)
  }
  fileName_finalOutput <- paste0(myDirectory,"TVV_", window_nr, "_final.csv")
  result <- list("optim_statistics"=optim_statistics, "model_parameters"=model_parameters, "filtered_error_RMSE"=filtered_error_RMSE)
  write.csv(result, fileName_finalOutput, row.names=TRUE)
  return(optim_values$par)
}
# Convert Log Likelihood of MLE estimates between models with differently scaled data
para_scaler <- function(parameters){
  parameters_Scaled <- parameters
  parameters_Scaled[5:50] <- parameters_Scaled[5:50]/100
  parameters_Scaled
}
