# Generator.R is used to generate the data
# Used for the simulation study

###### Using function in different R scripts ########

simple_fun <- function(input)
{
  print(input)
}

####################################################
library(tidyverse)
library(mvtnorm)
library(corrplot)


####################################################
#        Generate the simulation data              #
####################################################
# Each observation have a mean 0 and variance 1
#
#                   Inputs:
#
# Predictor parameters:
# P: number of predictors
# S.1: number of strong signal
# S.2: number of weak but dependent signal
# S.3: number of weak and independent signal
# Remain signals are consider to be white noise
# 
# 
# N: number of observations
# Limit: threshold of beta as strong signal
# 
#                   Outputs:
#
# dataframe for generated X
# beta: vector of Beta
# The vector of Y (with error term)
Gendata <- function(P, S.1, S.2, S.3, N = 100, SP=0.8){
  
  cor.matrix <- GenMatrix(P, S.1, S.2, S.3, sp = SP)
  
  Df <- rmvnorm(N, sigma = cor.matrix) %>% as.data.frame()
  
  for (x in 1:length(Df)){
    Df[,x] = Df[,x]/sqrt(var(Df[,x]))
  }
  return(Df)
}




####################################################
#        Generate the correlation matrix           #
####################################################
# Correlation matrix have an extra line and column 
# to describe the error term (all cor will be 0)
#
#                   Inputs:
# Size: total number of predictors
# S1 to S, same as previous function (S.1 to S.3)
#
#                   Outputs:
# Cor.matrix: correlation matrix with Size + 1 rows
GenMatrix <- function(Size, S1, S2, S3, sp = 0.6)
{
  Re.matrix = diag(1, Size + 1, Size + 1)
  for (i in 1 : (S1-1)){
    for (j in ((i+1) : (S1))){
      Re.matrix[i,j] = runif(1) * (runif(1) > (1 - sp))
    }
  }
  for (k in (S1+1) : (S1+S2-1)){
    for (l in ((k+1) : (S1+S2))){
      Re.matrix[k,l] = runif(1) * (runif(1) > (1 - sp))
    }
  }
  for (m in (S1+S2+1) : (S1+S2+S3-1)){
    for (n in ((m+1) : (S1+S2+S3))){
      Re.matrix[m,n] = runif(1) * (runif(1) > (1 - sp))
    }
  }
  return(t(Re.matrix) %*% Re.matrix)
}



### Validation of the data-generating method:

### Test distribution of mean/var/cor under different situation


### Setup: Gendata(30,4,12,12,N=20)

### Col of interest
### Mean/Var: select 2(S.1),14(S.2),26(S.3),31 supposed to be 0/1
### Cor:  Select pair 1,2 (S.1) <>, 13,14 (S.2) <>, 18,19 (S.3) <>
###                   1,14  (S.1/2) <0>
###                   1,19 (S.1/3) <0>
###                  14,19 (S.2/3) <0>
###                   1,31 (S.1/e) <0>
###                  14,31 (S.2/e) <0>
###                  19,31 (S.3/e) <0>

### Generating 500 times of data to see what will happened:


ME = tibble(Col.1 = c(),
            Col.14 =c(),
            Col.19 =c(),
            Col.31 =c())
VR = tibble(Col.1 = c(),
            Col.14 =c(),
            Col.19 =c(),
            Col.31 =c())

COV = tibble(C1.2 =c(),
             C13.14 = c(),
             C18.19 = c(),
             C1.14 = c(),
             C1.19 = c(),
             C14.19 = c(),
             C1.31 =c(),
             C14.31= c(),
             C19.31 = c())


for(total in 1:500){
  data <- Gendata(30,4,12,12,N=20)
  
  ME.S <- tibble(Col.1 = mean(data[,1]),
                 Col.14 = mean(data[,14]),
                 Col.19 = mean(data[,19]),
                 Col.31 = mean(data[,31]))
  ME <- bind_rows(ME,ME.S)
  
  VR.S <- tibble(Col.1 = mean(data[,1]),
                 Col.14 = mean(data[,14]),
                 Col.19 = mean(data[,19]),
                 Col.31 = mean(data[,31]))
  VR <- bind_rows(VR,VR.S)
  
  COV.S <- tibble(C1.2 = cor(data[,1],data[,2]),
                  C13.14 = cor(data[,13],data[,14]),
                  C18.19 = cor(data[,18],data[,19]),
                  C1.14 = cor(data[,1],data[,14]),
                  C1.19 = cor(data[,1],data[,19]),
                  C14.19 = cor(data[,14],data[,19]),
                  C1.31 = cor(data[,1],data[,31]),
                  C14.31= cor(data[,14],data[,31]),
                  C19.31 = cor(data[,19],data[,31]))
  COV <- bind_rows(COV,COV.S)
}









