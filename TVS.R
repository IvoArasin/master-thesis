library(readxl)
data <- read_excel("/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/zero_rates08til23.xlsx")
m <- c(1/365, 7/365, 14/365, 1/12, 2/12, 3/12, 6/12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15)
N <- 4
scaler <- 1
data <- as.matrix(subset(data, select=-c(1, 21, 22)))/scaler

print_param_counter <<- 0
total_iterations <<- 0

# Kalman Filter for DNS model with time-varying shape parameter
independent_DNS_TVS <- function(para,Y,lik, h=0) {
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
  phi<- diag(N) # Vector Autoregressive coeffient matrix VAR(1)	
  H	<- diag(W) # Variance matrix of residuals
  Q	<- diag(N) # Transition covariance matrix of residuals
  
  # Vector autoregressive coeffient matrix: VAR(1)
  phi[1,1] <- para[1]
  phi[2,2] <- para[2]
  phi[3,3] <- para[3]
  phi[4,4] <- para[4]
  
  diagA <- eigen(phi)
  diagA_eigenvalues <- diagA$values
  
  if(max(abs(Re(diagA_eigenvalues)))>=1){
    print("Phi has Eigenvalue >1")
    return(10000000000000)
    #break
  }
  
  # Mean vector
  mu[1:3]<-para[5:7]
  mu[4] <- para[8]
  
  # Transition covariance matrix of residuals
  Q[1,1] <- para[9]
  Q[2,2] <- para[10]
  Q[3,3] <- para[11]
  Q[4,4] <- para[12]
  
  Q <- Q %*% t(Q) 
  
  # Variance matrix of residuals
  H <- diag(para[13:32])^2
  
  v1   <- matrix(NA,T,W)			  
  v2   <- matrix(NA,T,W) # Filtered errors: are defined as the difference between the observed yield curve and its filtered estimate from KF
  
  a.tt <- matrix(NA, T, N)
  a.t  <- matrix(NA, (T+1), N)
  P.tt <- array(NA, c(T, N, N))
  P.t  <- array(NA, c((T+1), N, N))
  
  # Start state vector and variance matrix
  a.t[1, ]  <- mu # Start state vector: at0
  
  # Start variance matrix
  lyapunov<-function(N,phi,Q){ #here, Ivo changed "A" to the current "phi"
    matrix(solve(diag(N^2) - kronecker(phi,phi)) %*% matrix(Q,(N^2),1),N,N)
  }  
  P.t[1, ,] <-lyapunov(N=N,phi=phi,Q=Q) # Start variance matrix. Pt0
  # Initial log-likelihood	
  logLik <- 1
  
  # Kalman Filter and log-likelihood
  for (t in 1:T) 
  { 
    source("/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/factor_loadings.R")
    NS_factor_term <- factor_loadings(exp(a.t[t,4]), m) %*% a.t[t, 1:3]
    v <- (as.numeric(Y[t, ])) - NS_factor_term # prediction error vector
    
    source("~/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/jacobian.R")
    jac_w_params <- jacobian(a.t[t, 2], a.t[t, 3], exp(a.t[t, 4]), m)
    #print(jac_w_params)
    F <-  jac_w_params%*% P.t[t, ,] %*% t(jac_w_params) + H # prediction error variance matrix
    
    detF <- det(F)
    if(is.na(detF) || is.nan(detF) || is.infinite(detF) || abs(detF)<1e-160 || detF<0){
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
      
      logLik <- logLik-0.5*(length(v))*log(2*pi)-0.5*log(detF)-0.5*t(v)%*%F.inv%*%v # constructed via the prediction error decomposition
    
      # Kalman Gain
      Kalman_Gain <- P.t[t, , ] %*% t(jac_w_params) %*% F.inv
      
      # Updating the state vector and its variance matrix
      a.tt[t, ]   <- a.t[t, ] +  Kalman_Gain %*% v
      P.tt[t, , ] <- P.t[t, , ] - Kalman_Gain %*% jac_w_params %*% P.t[t, , ]
      
      source("/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/factor_loadings.R")
      v1[t, ]	<- factor_loadings(exp(a.tt[t,4]), m) %*% a.tt[t, 1:3] # Filtered prediction, i.e. estimate
      v2[t, ] <- (as.numeric(Y[t, ])) - v1[t, ] # Filtered errors
      
      # Predicting the state vector and its variance matrix
      a.t[t + 1, ]  <- phi %*% a.tt[t, ] + (diag(N) - phi) %*% mu  
      P.t[t + 1, ,] <- phi %*% P.tt[t, ,] %*% t(phi) + Q
    }
  }
  
  # Forecasting
  y_t_h <- matrix(NA, h, W)
  if(t > T-1 & h>0){
    recurstion_term <- diag(N)
    last_a.tt <- a.tt[T, ]
    for(i in 1:h){
      a_t_pred  <- phi %*% last_a.tt + (diag(N) - phi) %*% mu 
      #recurstion_term <- recurstion_term + ifelse(i>1, 1, 0)*phi**(i-1)
      #a_t_pred <- recurstion_term%*%(diag(N)-phi)%*%mu + phi**(i)%*%a.tt[T, ]
      
      source("/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/factor_loadings.R")
      y_t_h[i, ] <- factor_loadings(exp(a_t_pred[4]), m)%*%a_t_pred[1:3]
      last_a.tt <- a_t_pred
    }
    print(a_t_pred)
  }
  
  if(lik)
  {
    as.numeric(-logLik)
  }else{
    return(list(a.tt=a.tt,a.t=a.t,P.tt=P.tt,P.t=P.t,v2=v2,v1=v1, y_t_h=y_t_h, negLogLik=-logLik))
  }
}


# To convert parameter estimates between models with differently scaled data
para_scler <- function(parameters){
parameters_scaledDown <- parameters
parameters_scaledDown[5:7] <- parameters_scaledDown[5:7]/100
parameters_scaledDown[9:11] <- parameters_scaledDown[9:11]/100
parameters_scaledDown[13:32] <- parameters_scaledDown[13:32]/100
parameters_scaledDown
}

# Para optimized from init para from two-step method
para_twoStep_optim_DNS_TVS <- c(
  0.949920437604552,0.993759943955409,0.817536981785338,0.917059248925496,2.14302091461833,-1.06820668534872,-3.09498495150364,-0.34066824960347,0.199527379927803,0.254413090090596,0.548678786898038,0.192775614587237,0.160116703805882,0.0976644969006015,0.057284886089234,0.0223957982910576,-0.00203582774570009,0.0163050007603775,0.0385699155159099,0.0601552092331301,0.0465520008069969,0.0345834759533765,0.0217959728176574,0.0101969792545626,0.00435844250858639,0.00438558949513976,0.00496133082381787,0.0042722225255378,0.00794724309333226,0.0259698847891625,0.0348217551357714,0.0719224049208639
)


# Para from original optimization WITHOUT twoStep method init
para_noTwoStep_optim_DNS_TVS <- c(
  0.99807717,
  0.9716746,
  0.53807089,
  0.32368666,
  4.95666377,
  -1.133603,
  -3.0559929,
  -0.6719975,
  0.24183509,
  0.28121741,
  0.65147548,
  0.36694022,
  0.16938064,
  0.10236043,
  0.06129074,
  0.02564866,
  -2.00E-05,
  0.01856056,
  0.04688426,
  0.06774986,
  0.04693074,
  0.0226544,
  0.0105471,
  0.006731,
  0.00410906,
  0.00267378,
  0.00386889,
  -0.0052009,
  -0.0050708,
  0.01318372,
  0.01923475,
  0.04456685
)
  
# optimize model
#para_init <- c(
#  0.97, 0.97, 0.93, 0.94,
#  2.9, -0.9, -2.5, 0.42,
#  0.26, 0.31, 0.5, 0.03,
#  
#  0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
#  0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
#)
#para_init <- c( # This one works well
#  0.98, 0.97, 0.93, 0.6,
#  3.2, -1.6, -2.4, 0,
#  0.26, 0.31, 0.5, 0.2,
#  
#  0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
#  0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
#)

para_init <- c( # This one works very well along with the bounds
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

independent_DNS_TVS_fullSample <- optim(para_init,independent_DNS_TVS,Y=data,lik=TRUE,control = list(REPORT=2), method="L-BFGS-B", lower=lower_, upper=upper_)
independent_DNS_TVS_fullSample$par
independent_DNS_TVS_fullSample$convergence
independent_DNS_TVS_fullSample$message
independent_DNS_TVS_fullSample$par

# apply model
optimized_TVS_parameters_fullSample <- c(9.698777e-01,
                                         9.837768e-01,
                                         9.995413e-01,
                                         6.442378e-01,
                                         2.900151e+00,
                                         5.961248e-01,
                                         -2.019678e+00,
                                         -3.866426e-01,
                                         2.191424e-01,
                                         2.616273e-01,
                                         5.803580e-01,
                                         2.584501e-01,
                                         1.495545e-01,
                                         9.353353e-02,
                                         5.736415e-02,
                                         2.643208e-02,
                                         6.274542e-05,
                                         1.755474e-02,
                                         4.175561e-02,
                                         7.536156e-02,
                                         5.561795e-02,
                                         3.667176e-02,
                                         1.650445e-02,
                                         8.661341e-03,
                                         9.887177e-03,
                                         6.535739e-03,
                                         3.856203e-03,
                                         -1.176929e-02,
                                         2.257124e-02,
                                         3.709290e-02,
                                         5.180030e-02,
                                         8.466528e-02)

independent_TVS_model <- independent_DNS_TVS(para=optimized_TVS_parameters_fullSample, Y=data, lik=FALSE, h=12)
independent_TVS_model$negLogLik

# Reconstruct Yield Curve from Model
time_ <- 100
recon <- factor_loadings(exp(independent_TVS_model$a.tt[time_,4]), m)%*%independent_TVS_model$a.tt[time_,1:3]
plot(m, data[time_,], type="l") # actual
lines(m, recon, type="l", lty=2, col="blue") # model reconstruction

# Preds
plot(m, data[(time_+12),], type="l") # actuals
lines(m, independent_TVS_model$y_t_h[12,], type="l", lty=2, col="blue") # model predictions
lines(m, data[(time_),], type="l") # random walk

# Plot Level, Slope and Curvature and Shape
plot(exp(independent_TVS_model$a.tt[, 4]), type="l", col="blue")
plot(independent_TVS_model$a.tt[, 1], type="l", col="red")
plot(independent_TVS_model$a.tt[, 2], type="l", col="darkgreen")
plot(independent_TVS_model$a.tt[, 3], type="l", col="orange")

###############################################################
###############################################################
###############################################################
# Automatic Window Estimation
TVS_twoStepFullSample_inits <- c(
                             0.97, 0.98, 0.94, 0.95,
                             
                             1.9, -1.5, -2.3, log(0.35), # Shape parameter is log-bounded. So it's logarithm is estimated
                             
                             0.2, 0.2, 0.2, 0.4,
                             
                             0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                             0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
automated_rollingWindow(TVS_twoStepFullSample_inits, independent_DNS_TVS, 25, c(106, 111, 116, 121, 126, 131))

automated_rollingWindow <- function(para_init, model, start_, end_){
  # six expanding windows
  # 25:106
  # 25:111
  # 25:116
  # 25:121
  # 25:126
  # 25:131
  starting_values <- para_init
  start <- start_
  for(i in 1:length(end_)){
    end <- end_[i]
    window_data <- data[start:end,]
    starting_values <- automatic_optimization(starting_values, model, window_data, maxiter=15, window_nr=i)
    #function(para,Y=train_data,lik=TRUE, forecast_horizon=0) {
    preds <- model(para=starting_values, Y=window_data, lik=FALSE, h=12)
    preds <- preds$y_t_h
    preds_h12 <- preds[12,]
    preds_h6 <- preds[6,]
    preds_h1 <- preds[1,]
    predictions <- rbind(preds_h12, preds_h6, preds_h1)
    myDirectory = "/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/optimized_files/CalmTimeRollingWindow/"
    write.csv(predictions, paste0(myDirectory,"TVS_predictions_window_", i, ".csv"), row.names=TRUE)
  }
}

automatic_optimization <- function(para_init, model, input_data, maxiter=10, window_nr=0){
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
    plot(exp(model_output$a.tt[,4]), type="l", main=i)
    
    if(optim_values$convergence!=0){
      print("POLYTOPE ERROR")
      logLike_add <- 10
      prevLogLike <- 0
    }else{
      logLike_add <- 0
    }
    
    logLike <- logLike + logLike_add
    optim_statistics[i,3] <- logLike
    
    filtered_error_RMSE[i, ] <- sqrt(colMeans((model_output$v2)^2))
    
    total_iterations <<- 0
    
    myDirectory = "/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/optimized_files/CalmTimeRollingWindow/"
    fileName <- paste0(myDirectory,"TVS_window_intermediate.csv")
    intermediate_result <- list("optim_statistics"=optim_statistics, "model_parameters"=model_parameters, "filtered_error_RMSE"=filtered_error_RMSE)
    write.csv(intermediate_result, file=fileName, row.names = TRUE)
  }
  fileName_finalOutput <- paste0(myDirectory,"TVS_window_", window_nr, "_final.csv")
  #result <- list("model_parameters"=model_parameters)
  result <- list("optim_statistics"=optim_statistics, "model_parameters"=model_parameters, "filtered_error_RMSE"=filtered_error_RMSE)
  write.csv(result, fileName_finalOutput, row.names=TRUE)
  return(optim_values$par)
}
