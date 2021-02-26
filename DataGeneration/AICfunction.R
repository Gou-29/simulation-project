source("Generator.R")

# Final function for AIC method:

AICMETHOD <- function(TuneGrid, Iteration){
  
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
      as_tibble() %>% 
      mutate(y = y) %>% 
      select(c(y, col_name))
    
    fit.forward = step(object = lm(y ~ 1, data = data),
                       scope = formula(lm(y ~ ., data = data)), 
                       direction = "forward", k = 2, trace = 0) 
    
    selected <- attr(x = fit.forward$terms, which = "variables")
    
    S.1_Selected = c(S.1_Selected, sum(as.character(selected) %in% vec_s.1))
    S.2_Selected = c(S.2_Selected, sum(as.character(selected) %in% vec_s.2))
    S.3_Selected = c(S.3_Selected, sum(as.character(selected) %in% vec_s.3))
    Noi_Selected = c(Noi_Selected, sum(as.character(selected) %in% vec_w.n))
    
    TuneGrid_final = bind_rows(TuneGrid_final,TuneGrid)
    
    print(str_c("@@@@  AIC, Iteration No.",i,"  @@@@"))
    
    
  }
  
  
  Filename = str_c("SimResult1/AIC/AIC - TuneGrid-", TuneGrid$Number ,".csv")
  TuneGrid_final %>% 
    mutate(S1_Selected = S.1_Selected,
           S2_Selected = S.2_Selected,
           S3_Selected = S.3_Selected,
           No_Selected = Noi_Selected) %>% 
    mutate(VselecMethod = rep("AIC",Iteration)) %>% 
    write_csv(Filename)
  
  print(str_c("Finish AIC TuneGrid No. ", TuneGrid$Number))
  
}
