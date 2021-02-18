#####################################
#        AIC Criteria Selection     #                     
#        Forward Selection          #
#####################################

library(tidyverse)

set.seed(1234)

source("Generator2.R")

# define number of diff signals
n_s.1 = 3
n = 20
n_exp = 500
p = 30

data = Gen_Method1(n_s.1,c(2,2,2),6,p,n)

# get a vector of strong signal
s.1 = c()
for (i in 1:n_s.1){
  s.1 = c(rbind(s.1,paste0("x",i)))
}

# get the threshold

thre = sqrt(log(p)/n)

# get y 

beta = rep(c(3,0.1,0,1),c(n_s.1,12,p-n_s.1-12,1))

y =  data %*% as.matrix(beta)

# generate dataset for regression 

col_name = c()
for (i in 1:p){
  col_name = rbind(col_name,paste("x",i,sep = ""))
}
col_name = rbind(col_name,"y")

df = 
  as_tibble(data[,-1]) %>% 
  mutate(y = y) %>% 
  setNames(col_name) 

df%>% View

# AIC selection process 

for (i in 1:n_exp)
{
  fit.forward = step(object = lm(y ~ 1, data = df),
                      scope = formula(lm(y ~ ., data = df)), 
                      direction = "forward", k = 2, trace = 0) 
  
 selected <- attr(x = fit.forward$terms, which = "variables")
 n_true = 0
 for (x in selected){
   ifelse(x %in% c("list of strong signal"),n_true = n_true+1,n_true = n_true)
 }
  prop_true = n_true/length(selected)
}


