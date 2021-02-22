# Gendata 2:
library(tidyverse)
library(mvtnorm)
library(corrplot)
library(glmnet)
# Conpound Symmetric

Gen_CS <- function(sigma, sigma_1, size){
  re.cs <- diag(sigma^2 + sigma_1^2,size,size)
  for (i in 1:(size-1)){
    for (j in (i+1):size){
      re.cs[i,j] = sigma_1^2
      re.cs[j,i] = sigma_1^2
    }
  }
  return(re.cs)
}


# Ar(1)

Gen_AR <- function(ro, size){
  re.ar <- diag(1, size, size)
  for (i in 1:(size-1)){
    for (j in (i+1):size){
      re.ar[i,j] = ro^(j-i)
      re.ar[j,i] = ro^(j-i)
    }
  }
  return(re.ar)
}

# Matrix operation

join_matrix <- function(matrix.1, matrix.2){
  dim1 <- dim(matrix.1)[1]
  dim2 <- dim(matrix.2)[1]
  re.matrix <- diag(0,dim1+dim2, dim1+dim2)
  re.matrix[1:dim1,1:dim1] = matrix.1
  re.matrix[(dim1+1):(dim1+dim2),(dim1+1):(dim1+dim2)] = matrix.2
  return(re.matrix)
}


## Generate dataests:

Gen_Method1 <- function(S.1, S.2, S.3, P, N, co){ 
  
  #S.2 is a vector with same length as S.1
  if(length(S.2) != S.1){
    return("Length of S.2 must be same as value of S.1")
  }
  if(S.1 + sum(S.2) + S.3 > P){
    return("invalid total length!")
  }
  
  # Region S.1 and S.2, with pivot specified:
  
  X <- Gen_CS(sqrt(1-co),sqrt(co),1+S.2[1])
  pivot <- S.2   # For the position of strong signal
  pivot[1] <- 1
  
  for(i in 2:S.1){
    X_PLUS <- Gen_CS(sqrt(1-co),sqrt(co),1+S.2[i])
    X <- join_matrix(X,X_PLUS)
    pivot[i] <- 1 + sum(S.2[1:(i-1)]) + (i-1) 

  }

  # Region S.3
  X_3 <- Gen_CS(sqrt(1-co),sqrt(co),S.3)
  X <- join_matrix(X, X_3)

  # Noise and error term (i.i.d)
  
  X_Ad <- diag(1, P-S.1-sum(S.2)- S.3 + 1, P-S.1-sum(S.2) - S.3 +1)
  X <- join_matrix(X, X_Ad)
  # Final data:
  
  Data <- rmvnorm(N, sigma = X) 
  
  # Put pivots in first
  
  Data <- cbind(Data[,pivot], Data[,-pivot])
  
  return(Data)
}

Gen_Method2 <- function(S.1, S.2, S.3, P, N, ro){ 
  #S.2 is a vector with same length as S.1
  if(length(S.2) != S.1){
    return("Length of S.2 must be same as value of S.1")
  }
  if(S.1 + sum(S.2) + S.3 > P){
    return("invalid total length!")
  }
  
  # Region S.1 and S.2, with pivot specified:
  
  X <- Gen_AR(ro, 1+S.2[1])
  pivot <- S.2   # For the position of strong signal
  pivot[1] <- 1
  
  for(i in 2:S.1){
    X_PLUS <- Gen_AR(ro, 1+S.2[i])
    X <- join_matrix(X,X_PLUS)
    pivot[i] <- 1 + sum(S.2[1:(i-1)]) + (i-1) 
  }
  
  # Region S.3
  X_3 <- Gen_AR(ro,S.3)
  X <- join_matrix(X, X_3)
  
  # Noise and error term (i.i.d)
  
  X_Ad <- diag(1, P-S.1-sum(S.2)- S.3 + 1, P-S.1-sum(S.2) - S.3 +1)
  X <- join_matrix(X, X_Ad)

  # Final data:
  
  Data <- rmvnorm(N, sigma = X) 
  
  # Put pivots in first
  
  Data <- cbind(Data[,pivot], Data[,-pivot])
  
  return(Data)
}

Gen_Method3 <- function(S.1, S.2, S.3, P, N, co){ 
  #S.2 is a vector with same length as S.1
  if(length(S.2) != S.1){
    return("Length of S.2 must be same as value of S.1")
  }
  if(S.1 + sum(S.2) + S.3 > P){
    return("invalid total length!")
  }
  WN = P - S.1 - sum(S.2) - S.3 
  if(WN %% S.1 !=0){
    return("Non - compatible length of white noise!")
  }
  # Region S.1 and S.2 + White noise, with pivot specified:
  # In this case, #White noise = c * #Strong predictor 
  
  X <- Gen_CS(sqrt(1-co),sqrt(co),1+S.2[1]+ WN/S.1)
  pivot <- S.2   # For the position of strong signal
  pivot[1] <- 1
  pivot_WN <- c((S.1 + S.2[1] + 1):(S.1 + S.2[1] + WN/S.1) ) # For the position of white noise
  
  for(i in 2:S.1){
    X_PLUS <- Gen_CS(sqrt(1-co),sqrt(co),1 + S.2[i] + WN/S.1)
    X <- join_matrix(X,X_PLUS)
    pivot[i] <- 1 + sum(S.2[1:(i-1)]) + (i-1) + (i-1) * WN/S.1 
    pivot_WN <- c(pivot_WN, 
                  (S.1 + sum(S.2[1:i]) + 1 + (i-1) * WN/S.1)
                  :(S.1 + sum(S.2[1:i]) + i * WN/S.1))
    
  }
  
  # Region S.3
  X_3 <- Gen_CS(sqrt(1-co),sqrt(co),S.3)
  X <- join_matrix(X, X_3)
  
  # Error term
  X_ER <- diag(1,1,1)
  X <- join_matrix(X, X_ER)
  
  # Final data:
  
  Data <- rmvnorm(N, sigma = X) 
  
  # Put pivots in first
  
  Data <- cbind(Data[,pivot], Data[,-pivot])
  
  # Put pivots of white noise.

  Data <- cbind(Data[,-c(pivot_WN,P+1)],Data[,pivot_WN], Data[,(P+1)])
                     
  return(Data)
}



