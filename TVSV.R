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

factor_loadings <- function(l,m)
{
  column1 <- rep.int(1,length(m))
  column2 <- (1 - exp(-l * m))/(l * m)
  column3 <- (1 - exp(-l * m) - l*m*exp(-l*m))/(l*m)
  
  B <- cbind(column1,column2,column3)
  B
}

print_param_counter <<- 0
iteration_counter <<- 0
KalmanFilterTVSV <- function(para,Y=train_data,lik=TRUE, h=0) {
  # print parameters every 100th iteration
  if(print_param_counter >=100){
    print("Paras printed because of counter")
    print(para)
    print_param_counter <<- 0
  }
  print_param_counter <<- print_param_counter +1
  iteration_counter <<- iteration_counter +1
  
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
  
  # Variance matrix of residuals
  H <- diag(para[13:32])^2
  
  # Vector autoregressive coeffient matrix: VAR(1)
  phi[1,1] <- para[1] # level autoregressive factor
  phi[2,2] <- para[2] # slope autoregressive factor
  phi[3,3] <- para[3] # curvature autoregressive factor
  phi[4,4] <- para[4] # shape parameter autoregressive factor
  phi[5,5] <- 0
  
  
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
  
  # Transition covariance matrix of residuals
  Q[1,1] <- para[9]  # level
  Q[2,2] <- para[10]  # slope
  Q[3,3] <- para[11] #curvature
  Q[4,4] <- para[12] # shape parameter
  Q[5,5] <- 0
  
  Q <- Q %*% t(Q) 
  
  v1   <- matrix(NA,T,W)			  
  v2   <- matrix(NA,T,W) # Filtered errors: are defined as the difference between the observed yield curve and its filtered estimate from KF
  
  # Gamma vector
  Gamma <- abs(para[33:52])
  
  gamma0 <- 0.0001 #fixed for identifiability (original)
  gamma1 <- para[53]
  gamma2 <- para[54]
  
  
  a.tt <- matrix(NA, T, N+1)
  a.t  <- matrix(NA, (T+1), N+1)
  P.tt <- array(NA, c(T, N+1, N+1))
  P.t  <- array(NA, c((T+1), N+1, N+1))
  
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
  
  P.t[1, 5, 5] <- gamma0/(1-gamma1-gamma2)
  
  P.t[1,5, 1:4 ] <- rep(0, 4)
  P.t[1,1:4, 5 ] <- rep(0, 4)
  
  h_t <- rep(T, 0)
  
  # Initial log-likelihood	
  logLik <- 0
  
  for (t in 1:T) 
  {
    # Restrictions
    if(gamma0 <= 0 || gamma1<=0 || gamma2<=0 || gamma1+gamma2>=1){
      print(paste("Gamma error: ", (gamma1+gamma2)))
      logLik<- -1000000000000000
      break
      
    }else{
      
      NS_factor_term <- factor_loadings(a.t[t,4], m) %*% a.t[t, 1:3] # Gamma isn't considered here because the
      # predicted value of the common volatility term is in zero in expectation, which would result in
      # Gamma always being multiplied with zero, meaning its omission is just more efficient
      
      if(T==1){
        Y_data <- as.numeric(Y)
      }else{
        Y_data <- (as.numeric(Y[t, ]))
      }
      v <- Y_data - NS_factor_term # prediction error vector
      
      jacobian_EKF <- function(beta2, beta3, l,m){
        column1 <- rep.int(1,length(m))
        column2 <- (1 - exp(-l * m))/(l * m)
        column3 <- (1 - exp(-l * m) - l*m*exp(-l*m))/(l*m)
        column4 <- (( beta2*exp(-m*l)*(m*l-exp(m*l)+1) + beta3*exp(-m*l)*(m**2*l**2+m*l-exp(m*l)+1))/(m*l**2))#*l
        jac_dataframe <- cbind(column1,column2,column3, column4)
        jac_dataframe
      }
      
      jac <- jacobian_EKF(a.t[t, 2], a.t[t, 3], a.t[t, 4], m)
      
      jac_w_Gamma <- cbind(jac, Gamma)
      F <- jac_w_Gamma %*% P.t[t, ,] %*% t(jac_w_Gamma) + H # prediction error variance matrix
      
      detF <- det(F)
      if(is.na(detF) || is.nan(detF) || is.infinite(detF) || abs(detF)<1e-160 || detF<0){
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
            #return(1000000000000000)
            
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
        
        # Predicting the state vector and its variance matrix
        a.t[t + 1, ]  <- phi %*% a.tt[t, ] + (diag(N+1) - phi) %*% c(mu, 0)
        
        h_t[t] <- Q[5,5]
        
        # GARCH process for common volatility term h_t
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

# Initial parameters based on DNS, TVS and TVV
para_init <- c(
  0.98, 0.97, 0.93, 0.64,
  3.2, -1.6, -2.4, 0.35, # Mean of shape parameter initialized by MLE estimates of standard DNS
  0.26, 0.31, 0.5, 0.2, # Shape variance initialized at generic value 0.2
  
  0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2, # Measurement noise variances initialized at generic values
  0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,
  
  2.86823415, # Gamma vector initialized at values found by MLE estimates of TVV model
  2.66683849,
  2.59439467,
  2.40663238,
  2.1185691,
  1.83628053,
  1.14769614,
  0.3314479,
  0.02681608,
  0.2158445,
  0.50877774,
  0.7766088,
  0.99584189,
  1.15609028,
  1.25567294,
  1.29949906,
  1.32561092,
  1.32664272,
  1.34158344,
  1.38736245,
  
  # Gamma coefficients initialized by MLE estimates of TVV model
  0.46, 0.52)

# optimize model
TVSV_replication <- optim(par=para_init, fn=KalmanFilterTVSV, control=list(REPORT=5), method="BFGS")
TVSV_replication$value
TVSV_replication$counts
TVSV_replication$convergence
TVSV_replication$par

# function to convert MLE parameters between differently scaled data
para_scaler <- function(para){
  scaled_para <- para
  scaled_para[5:7] <- para[5:7]/100
  scaled_para[9:11] <- para[9:11]/100
  scaled_para[13:32] <- para[13:32]/100
  scaled_para[33:52] <- para[33:52]/100
  return(scaled_para)
}

# apply model
model <- KalmanFilterTVSV(TVSV_BFGS_attempt1$par, Y=train_data, lik=FALSE, h=12)
