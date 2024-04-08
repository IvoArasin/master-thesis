library(readxl)
data <- read_excel("/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/zero_rates08til23.xlsx")
m <- c(1/365, 7/365, 14/365, 1/12, 2/12, 3/12, 6/12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15)
N <- 3
scaler <- 1
data <- as.matrix(subset(data, select=-c(1, 21, 22)))/scaler

factor_loadings <- function(l,m)
{
  column1 <- rep.int(1,length(m))
  column2 <- (1 - exp(-l * m))/(l * m)
  column3 <- (1 - exp(-l * m) - l*m*exp(-l*m))/(l*m)
  
  B <- cbind(column1,column2,column3)
  B
} 

print_param_counter <<- 0
total_iterations <<- 0
independent_DNS_TVV <- function(para,Y,lik, h=0) {
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
  H <- diag(para[11:30])^2
  
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
  
  Q <- Q %*% t(Q) 
  
  v1   <- matrix(NA,T,W)			  
  v2   <- matrix(NA,T,W) # Filtered errors: are defined as the difference between the observed yield curve and its filtered estimate from KF
  
  # Gamma vector
  Gamma <- abs(para[31:50])
  
  B <- cbind(B, Gamma)
  
  gamma0 <- 0.0001 #fixed for identifiability
  gamma1 <- para[51]
  gamma2 <- para[52]
  
  a.tt <- matrix(NA, T, N+1)
  a.t  <- matrix(NA, (T+1), N+1)
  P.tt <- array(NA, c(T, N+1, N+1))
  P.t  <- array(NA, c((T+1), N+1, N+1))
  
  # Start state vector and variance matrix
  a.t[1, 1:3] <- mu
  a.t[1, 4] <- 0
  
  # Start variance matrix
  lyapunov<-function(N,phi,Q){
    matrix(solve(diag(N^2) - kronecker(phi,phi)) %*% matrix(Q,(N^2),1),N,N)
  }  
  
  tryCatch(
    {
      cannot_be_done <- 0
      P.t[1,1:3,1:3] <-lyapunov(N=N,phi=phi[1:3,1:3],Q=Q[1:3,1:3]) # Start variance matrix. Pt0
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
  logLik <- 1
  
  for (t in 1:T) 
  {
    # Restrictions
    if(gamma0 <= 0 || gamma1<=0 || gamma2<=0 || gamma1+gamma2>=1){
      print("Gamma error")
      logLik<- -1000000000000000;
      return(1000000000000000)
      break
    }
    v <- (as.numeric(Y[t, ])) - B %*% a.t[t, ] # prediction error vector
    
    F <- B %*% P.t[t, ,] %*% t(B) + H # prediction error variance matrix
    
    detF <- det(F)
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


para_init <- c( # derived from two-step method
  0.36,
  0.98, 0.97, 0.93,
  3.2, -1.6, -2.4,
  0.26, 0.31, 0.5,
  
  0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
  0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
  
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  
  0.4, 0.4)

automated_rollingWindow(para_init=para_init, model=independent_DNS_TVV, start_=25, end_=c(106, 111, 116, 121, 126, 131), maxiter=15)

automated_rollingWindow <- function(para_init, model, start_, end_, maxiter=5){
  starting_values <- para_init
  start <- start_
  for(i in 1:length(end_)){
    end <- end_[i]
    window_data <- data[start:end,]
    starting_values <- automatic_optimization(starting_values, model, window_data, maxiter=maxiter, window_nr=i)
    preds <- model(para=starting_values, Y=window_data, lik=FALSE,h=12)
    preds <- preds$y_t_h
    preds_h12 <- preds[12,]
    preds_h6 <- preds[6,]
    preds_h1 <- preds[1,]
    predictions <- rbind(preds_h12, preds_h6, preds_h1)
    myDirectory = "/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/optimized_files/CalmTimeRollingWindow/"
    write.csv(predictions, paste0(myDirectory,"predictions_window_", i, "TVV.csv"), row.names=TRUE)
  }
}

automatic_optimization <- function(para_init, model, input_data, maxiter=5, window_nr=0){
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
    optim_values<-optim(parameters,model,control = list(trace=1, maxit=100000), Y=input_data, lik=TRUE, h=0)
    print("HERE")
    print(nrow(input_data))
    #convergence
    optim_statistics[i,1] <- optim_values$convergence
    
    #iterations
    optim_statistics[i,2] <- total_iterations
    
    # NegLogLike
    prevLogLike <- logLike
    logLike <- model(para=parameters, Y=input_data, lik=TRUE, h=0)
    optim_statistics[i,3] <- logLike
    
    # model parameters
    model_parameters[i, ] <- optim_values$par
    parameters <- optim_values$par
    
    # filtered RMSE
    model_output <- model(para=parameters, Y=input_data, lik=FALSE, h=0)
    
    filtered_error_RMSE[i, ] <- sqrt(colMeans((model_output$v2)^2))
    
    total_iterations <<- 0
    
    myDirectory = "/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/optimized_files/CalmTimeRollingWindow/"
    fileName <- paste0(myDirectory,"TVV_calmTime_window_intermediate.csv")
    intermediate_result <- list("optim_statistics"=optim_statistics, "model_parameters"=model_parameters, "filtered_error_RMSE"=filtered_error_RMSE)
    write.csv(intermediate_result, file=fileName, row.names = TRUE)
  }
  fileName_finalOutput <- paste0(myDirectory,"TVV_calmTime_window_", window_nr, "_final.csv")
  #result <- list("model_parameters"=model_parameters)
  result <- list("optim_statistics"=optim_statistics, "model_parameters"=model_parameters, "filtered_error_RMSE"=filtered_error_RMSE)
  write.csv(result, fileName_finalOutput, row.names=TRUE)
  return(optim_values$par)
}

parameters <- read.csv2("/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/optimized_files/RollingWindowForcasts/full_sample_parameters/TVV_fullSample.csv")
para_fullSamle_optim <- as.numeric(parameters[8,5:56])

# Convert Log Likelihood of MLE estimates between models with differently scaled data
para_scaler <- function(parameters){
  parameters_Scaled <- parameters
  parameters_Scaled[5:50] <- parameters_Scaled[5:50]/100
  parameters_Scaled
}

# apply model
time_ <- 131
k_obj_TVV <- independent_DNS_TVV(DNS_TVV_window_3_bfgs_normal$par, Y=train_data[25:time_,], lik=FALSE, h=12)
k_obj_TVV$negLogLik
plot(k_obj_TVV$h_t, type="l", main=expression(paste("Common Volatility ", h[t])))

plot(m, data[(time_+12),], type="l")
lines(m, data[(time_),], type="l", col="red", lty=2)
lines(m,k_obj_TVV$y_t_h[12,], type="l", lty=2, col="blue")

# Observe Predictive Tendency for 12mo ahead forecast
for(i in 1:12){
  if(i==1){plot(m, k_obj_TVV$y_t_h[i,], type="l",lwd=2)
    lines(m, data[(time_),], type="l", col="red", lty=2)
    lines(m, data[(time_+12),], type="l", col="blue", lty=2)
    }
  else{
    lines(m, k_obj_TVV$y_t_h[i,], type="l")
  }
}

# Reconstruct Yield Curve with Model
time_ <- 75
recon <- factor_loadings(para_init[1], m)%*%k_obj_TVV$a.tt[time_,1:3]+para_init[31:50]*k_obj_TVV$a.tt[time_,4]
plot(m, data[(time_),], type="l")
lines(m, recon, type="l", lty=2, col="blue")

# plot Gamma weights
plot(c(1/365, 7/365, 14/365, 1/12, 2/12, 3/12, 6/12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15), abs(DNS_TVV_window_1_bfgs_normal$par[31:50]), type="b", main=expression(paste("6 Common Volatility Weights ", Gamma[epsilon])))
# plot volatility variance
plot(k_obj_TVV$h_t, type="l", main=expression(paste("Common Volatility ", h[t])))
