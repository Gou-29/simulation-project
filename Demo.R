# About the simulation:


# Add a csv file in your local dirctory (repo)


beta <- ...
n <- ...

for (i = 1: number_of_simulation) {
  set.seed(i)
  Df<- Gen_Method1(...)
  beta <- 
  Y <- beta %*% Df 
  
  aic.result <- step(lm(....))
  
  aic$....
  
  result <- tibble( Iteration. no = i,
                    TPR = ...,
                   TNR = ...)
  
  result_pre <- read.csv("your file")
  
  result_pre <- rbind(result_pre, result)
  
  write.csv("your file", result_pre)
  
  
  print(str_c("finish iteration #", i.))
}