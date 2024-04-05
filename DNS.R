library(readxl)
data <- read_excel("/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/zero_rates08til23.xlsx")
m <- c(1/365, 7/365, 14/365, 1/12, 2/12, 3/12, 6/12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15)
N <- 3
scaler <- 1
data <- as.matrix(subset(data, select=-c(1, 21, 22)))/scaler

train_data <- data#[1:152,]

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
independent_DNS <- function(para,Y=train_data,lik=TRUE, forecast_horizon=0) {
  if(print_param_counter >=100){
    print(para)
    print_param_counter <<- 0
  }
  print_param_counter <<- print_param_counter +1
  total_iterations <<- total_iterations <<- 0 +1
  T <- nrow(Y)
  W <- ncol(Y)
  l<- para[1]
  pars<-list()

  # Create vectors and matrices
  mu	<- matrix(NA,N,1) # Mean vector
  phi<- diag(N) # Vector Autoregressive coeffient matrix VAR(1)	
  H	<- diag(W) # Variance matrix of residuals
  Q	<- diag(N) # Transition covariance matrix of residuals
  
  # Loading matrix
  B	<- factor_loadings(l,m) 
  
  # Variance matrix of residuals
  H <- diag(para[11:30])^2
  
  # Vector autoregressive coeffient matrix: VAR(1)
  phi[1,1] <- para[2]
  phi[2,2] <- para[3]
  phi[3,3] <- para[4]
  
  # Mean vector
  mu[1]<-para[5]
  mu[2]<-para[6]
  mu[3]<-para[7]
  
  # Transition covariance matrix of residuals
  Q[1,1] <- para[8]
  Q[2,2] <- para[9]
  Q[3,3] <- para[10]
  
  Q <- Q %*% t(Q) 
  
  v1   <- matrix(NA,T,W)			  
  v2   <- matrix(NA,T,W) # Filtered errors: are defined as the difference between the observed yield curve and its filtered estimate from KF
  
  #Set up dataframes for parameter collection
  a.tt <- matrix(NA, T, N)
  a.t  <- matrix(NA, (T+1), N)
  P.tt <- array(NA, c(T, N, N))
  P.t  <- array(NA, c((T+1), N, N))
  
  # Start state vector and variance matrix
  a.t[1, ]  <- mu # Start state vector: at0
  
  # Start variance matrix
  lyapunov<-function(N,phi,Q){
    matrix(solve(diag(N^2) - kronecker(phi,phi)) %*% matrix(Q,(N^2),1),N,N)
  }  
  
  P.t[1, ,] <-lyapunov(N=N,phi=phi,Q=Q) # Start variance matrix Pt0
  
  # Initial log-likelihood	
  logLik <- 0
  
  for (t in 1:T) 
  {
    v <- (as.numeric(Y[t, ])) - B %*% a.t[t, ] # prediciton error vector
    
    F <- B %*% P.t[t, ,] %*% t(B) + H # prediciton error variance matrix
    
    detF <- det(F)
    if(detF<=1e-150 || is.na(detF) || is.nan(detF) || is.infinite(detF) || detF<0){
      print("break reached in KFilter")
      logLik<- -1000000000000000; break
    }
    else{
      tryCatch(
        {
          cannot_be_done <- 0
          F.inv  <- solve(F)
        },
        warning = function(w){
          print(w)
          print("Difficult to invert F!")
          cannot_be_done <- 1
          logLik<- -1000000000000000
          
        },
        error=function(e){
          print(e)
          print("F could not be inverted!")
          cannot_be_done <- 1
          logLik<- -1000000000000000
        },
        finally={
          if(cannot_be_done==1){
            break
          }
        }
        
      )
      
      logLik <- logLik-0.5*(W)*log(2*pi)-0.5*log(detF)-0.5*t(v)%*%F.inv%*%v
      }
    # Updating the state vector and its variance matrix
    a.tt[t, ]   <- a.t[t, ] +  P.t[t, , ] %*% t(B) %*% F.inv %*% v
    P.tt[t, , ] <- P.t[t, , ] - P.t[t, , ] %*% t(B) %*% F.inv %*% B %*% P.t[t, , ]
    v1[t, ]	<- B %*% a.tt[t, ]
    v2[t, ] <- (as.numeric(Y[t, ])) - B %*% a.tt[t, ] # Filtered errors
    
    # Predicting the state vector and its variance matrix
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
# Two Step Guess for full Sample
#Mu: 1.913442 -1.506884 -2.345552
#AR: 0.977471  0.983737  0.946256
twoStepFullSample_inits <- c(0.72,
                     
                     0.97, 0.98, 0.94,
                     
                     1.9, -1.5, -2.3,
                     
                     0.2, 0.2, 0.2,
                     
                     0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                     0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
)

automatic_optimization_test <- automatic_optimization(twoStepFullSample_inits, independent_DNS, train_data, maxiter=30)

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
    optim_values<-optim(parameters,model,control = list(trace=1, maxit=100000))
    
    #convergence
    optim_statistics[i,1] <- optim_values$convergence
    
    #iterations
    optim_statistics[i,2] <- total_iterations
    
    # NegLogLike
    prevLogLike <- logLike
    logLike <- model(para=parameters, Y=input_data, lik=TRUE)
    optim_statistics[i,3] <- logLike
    
    # model parameters
    model_parameters[i, ] <- optim_values$par
    parameters <- optim_values$par
    
    # filtered RMSE
    model_output <- model(para=parameters, Y=input_data, lik=FALSE)
    
    filtered_error_RMSE[i, ] <- sqrt(colMeans((model_output$v2)^2))
    
    total_iterations <<- 0
    intermediate_result <- list("optim_statistics"=optim_statistics, "model_parameters"=model_parameters, "filtered_error_RMSE"=filtered_error_RMSE)
    write.csv(intermediate_result, "intermediate_result_DNS_fullSample_init.csv", row.names = TRUE)
  }
  
  result <- list("model_parameters"=model_parameters)
  result <- write.csv(result, "optimization_file_DNS_fullSample_init.csv", row.names=TRUE)
}

# To convert likelihood between differently scaled data
para_scaler <- function(parameters){
  parameters_scaledDown <- parameters
  parameters_scaledDown[5:30] <- parameters_scaledDown[5:30]/100
  parameters_scaledDown
}

# Empirical yields
# (1/365, 7/365, 14/365, 1/12, 2/12, 3/12, 6/12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15)
empirical_beta1 <- array(unlist(data[,20]))
empirical_beta2 <- array(unlist(data[,17] - data[,6]))
empirical_beta3 <- array(unlist(2*data[,10]-data[,6]-data[,20]))
