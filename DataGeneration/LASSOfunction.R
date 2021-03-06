source("Generator.R")

LASSOMETHOD <- function(TuneGrid, Iteration){
  # define number of diff signals
  n_s.1 = TuneGrid$S1
  n_s.2 = c(TuneGrid$S21,TuneGrid$S22,TuneGrid$S23,TuneGrid$S24,TuneGrid$S25)
  p = TuneGrid$P   # total parameter size
  n_s.3 = TuneGrid$S3
  n = TuneGrid$N   #sample size
  cor = TuneGrid$corr
  
  # get a vector of strong, wbc and weak signal
  vec_s.1 = str_c("V",c(1:n_s.1))
  vec_s.2 = str_c("V",c(n_s.1+1:sum(n_s.2)))
  vec_s.3 = str_c("V",c(n_s.1+sum(n_s.2)+1:n_s.3))
  
  
  # col names and white noise
  col_name = str_c("V", c(1:p))
  vec_w.n  = col_name[-which(col_name %in% c(vec_s.1,vec_s.2,vec_s.3))]
  # get the threshold (set c = 10)
  
  thre = sqrt(log(p)/n) * 10
  
  # Vector for results:
  S.1_Selected = c()
  S.2_Selected = c()
  S.3_Selected = c()
  Noi_Selected = c()
  Seed = c()
  
  
  TuneGrid_final = TuneGrid %>% slice(-1)
  
  
  cv.lambda <- c()
  
  # Final Lasso process:
  
  for (i in 1:Iteration)
  {
    set.seed(i)
    Seed <- c(Seed, i)
    # get beta 
    
    if(TuneGrid$beta == "close") {
      beta = c(runif(n_s.1,      min = thre * 1.001, max = thre * 1.25),  #S.1
               runif(sum(n_s.2), min = thre * 0.75, max = thre * 0.999),  #S.2
               runif(n_s.3,      min = thre * 0.75, max = thre * 0.999),  #S.3
               rep(0, p - n_s.1 - sum(n_s.2) - n_s.3),                   #White noise
               1)                                                        #Error term = 1
    }else{
      beta = c(runif(n_s.1,      min = thre * 1.75, max = thre * 2),  #S.1
               runif(sum(n_s.2), min = thre * 0.0, max = thre * 0.25),  #S.2
               runif(n_s.3,      min = thre * 0.0, max = thre * 0.25),  #S.3
               rep(0, p - n_s.1 - sum(n_s.2) - n_s.3),                   #White noise
               1)                                                        #Error term = 1
    }
    
    # Get Data
    
    if(TuneGrid$method ==1){
      data = Gen_Method1(n_s.1,n_s.2,n_s.3,P=p,N=n,co=cor)
    }else{
      if(TuneGrid$method == 2){
        data = Gen_Method2(n_s.1,n_s.2,n_s.3,P=p,N=n,ro=cor)
      }
      else{
        data = Gen_Method3(n_s.1,n_s.2,n_s.3,P=p,N=n,co=cor)
      }
    }
    
    # get y 
    y =  as.vector(data %*% beta)
    
    # generate dataset for regression 
    
    data = data %>% 
      as.tibble() %>% 
      mutate(y = y) %>% 
      select(c(y, col_name))
    
    
    
    fit.lasso <- cv.glmnet(y = y, x = as.matrix(data[,-1]), 
                           alpha = 1,    
                           lambda = exp(seq(0, -40, length = 500)),  
                           type.measure = "mse") 
    cv.lambda <- c(cv.lambda, fit.lasso$lambda.min)    #Validation check
    coeff = predict(fit.lasso, 
                    s = "lambda.min", 
                    type = "coefficients")[-1,] %>% 
      as.tibble() %>% 
      mutate(number = 1:p) %>% 
      filter(value != 0) 
    selected = paste(c('V'),coeff$number,  sep = '')
    
    
    
    S.1_Selected = c(S.1_Selected, sum(as.character(selected) %in% vec_s.1))
    S.2_Selected = c(S.2_Selected, sum(as.character(selected) %in% vec_s.2))
    S.3_Selected = c(S.3_Selected, sum(as.character(selected) %in% vec_s.3))
    Noi_Selected = c(Noi_Selected, sum(as.character(selected) %in% vec_w.n))
    
    TuneGrid_final = bind_rows(TuneGrid_final,TuneGrid)
    
    print(str_c("@@@@  Lasso, Iteration No.",i,"  @@@@"))
    
    
  }
  
  Filename = str_c("SimResult1/Lasso/Lasso - TuneGrid-", TuneGrid$Number ,".csv")
  TuneGrid_final %>% 
    mutate(S1_Selected = S.1_Selected,
           S2_Selected = S.2_Selected,
           S3_Selected = S.3_Selected,
           No_Selected = Noi_Selected,
           cv.lambda   = cv.lambda) %>% 
    mutate(VselecMethod = rep("Lasso",Iteration)) %>% 
    write_csv(Filename)
  
  print(str_c("Finish Lasso TuneGrid No. ", TuneGrid$Number))
  
}