# File for question 2:

source("Generator.R")


# Add S.2 -> S.3 -> W.N

Est_Method1 <- function(TuneGrid, Iteration){
  # define number of diff signals
  n_s.1 = TuneGrid$S1
  n_s.2 = c(TuneGrid$S21,TuneGrid$S22,TuneGrid$S23,TuneGrid$S24,TuneGrid$S25)
  p = TuneGrid$P   # total parameter size
  n_s.3 = TuneGrid$S3
  n = TuneGrid$N   #sample size
  cor = TuneGrid$corr
  
  # get names
  vec_s.1 = str_c("V",c(1:n_s.1))
  col_name = str_c("V", c(1:p))
  
  # get the threshold (set c = 10)
  
  thre = sqrt(log(p)/n) * 10
  
  # Result table:
  DFF <-
    matrix(data = 1:(length(vec_s.1)+2),ncol = length(vec_s.1)+2) %>% 
    as.data.frame()
  names(DFF) <- c(vec_s.1, "No_Predictors","MSE")
  DFF <- DFF %>% 
    as.tibble() %>%
    slice(-1)
  
  TuneGrid_final = TuneGrid %>% slice(-1)
  
  TrueBeta <- matrix(data = 1:length(vec_s.1),ncol = length(vec_s.1)) %>% 
    as.data.frame()
  names(TrueBeta) <- c(str_c(vec_s.1,"_TURE"))
  TrueBeta <- TrueBeta %>% as.tibble() %>% slice(-1)
  
  # Method 1:
  
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
    
    
    # Refill S.2 and S.3 and white noise:
    
    vec_s.2 = str_c("V",c(n_s.1+1:sum(n_s.2)))
    vec_s.3 = str_c("V",c(n_s.1+sum(n_s.2)+1:n_s.3))
    vec_w.n  = col_name[-which(col_name %in% c(vec_s.1,vec_s.2,vec_s.3))]
    
    
    # Null model formula:
    f_fit <- "y~V1"
    
    for (len in 2:length(vec_s.1)){f_fit <- str_c(f_fit, "+", vec_s.1[len])}
    
    # Fit null model:
    
    model <- lm(formula = formula(f_fit), data = data)
    coeff_est <- model$coefficients[2:(length(vec_s.1)+1)]   # fit coeff.
    mse <- mean( (y - model$fitted.values)^2)                # fit mse
    
    # Update the result:
    dfresult <- 
      t(coeff_est %>% as.data.frame()) %>% 
      as.tibble() %>% 
      mutate(No_Predictors = length(vec_s.1))%>% 
      mutate(MSE = mse)
    
    TuneGrid_final <- bind_rows(TuneGrid_final, TuneGrid)
    
    betasim <- beta[1:length(vec_s.1)] %>% 
      matrix(data = ., ncol = length(vec_s.1)) %>% as.data.frame()
    names(betasim) <- names(TrueBeta) <- c(str_c(vec_s.1,"_TURE"))
    betasim <- betasim %>% as.tibble() 
    
    TrueBeta <- bind_rows(TrueBeta, betasim)
    # Do a whole variable adding process in this generated data:
    
    for (j in (length(vec_s.1)+1):min(p,n)){
      
      No_Predictors <- j
      # fit the model
      model <- lm(formula = formula(f_fit), data = data)
      coeff_est <- model$coefficients[2:(length(vec_s.1)+1)]   # fit coeff.
      mse <- mean( (y - model$fitted.values)^2)                # fit mse
      
      dfresult_update <- 
        t(coeff_est %>% as.data.frame()) %>% 
        as.tibble() %>% 
        mutate(No_Predictors = No_Predictors) %>% 
        mutate(MSE = mse)
      
      No_Predictors <- No_Predictors + 1
      
      # add variable:
      # Random select a row in S.2 -> S.3 -> White noise
      
        if(length(vec_s.2) != 0){
          add_col_number <- sample(1:length(vec_s.2),1)
          add_col_name   <- vec_s.2[ add_col_number]
          vec_s.2        <- vec_s.2[-add_col_number]
        }else{
          if (length(vec_s.3) != 0 ){
            add_col_number <- sample(1:length(vec_s.3),1)
            add_col_name   <- vec_s.3[ add_col_number]
            vec_s.3        <- vec_s.3[-add_col_number]
          }else{
            add_col_number <- sample(1:length(vec_w.n),1)
            add_col_name   <- vec_w.n[ add_col_number]
            vec_w.n        <- vec_w.n[-add_col_number]
          }
        }
      
      
      #Fit model update
      f_fit <- str_c(f_fit, "+",add_col_name)
      
      # Update result
      TuneGrid_final <- bind_rows(TuneGrid_final, TuneGrid)
      TrueBeta <- bind_rows(TrueBeta, betasim)
      dfresult <- bind_rows(dfresult, dfresult_update)
    }
    
    # Save the result:
    
    DFF <- bind_rows(DFF,dfresult)
    
    print(str_c("@@@@  Method 1, Iteration No.",i,"  @@@@"))
    
  }
  
  Filename = str_c("Simresult2/Method1/Method1 - TuneGrid-", TuneGrid$Number ,".csv")
  bind_cols(TuneGrid_final, DFF, TrueBeta) %>% 
    mutate(VaddMethod = rep("Method1",nrow(DFF))) %>% 
    write_csv(Filename)
  
  print(str_c("Finish Method1 TuneGrid No. ", TuneGrid$Number))
  
}

# Add S.3 -> S.2 -> W.N
Est_Method2 <- function(TuneGrid, Iteration){
  # define number of diff signals
  n_s.1 = TuneGrid$S1
  n_s.2 = c(TuneGrid$S21,TuneGrid$S22,TuneGrid$S23,TuneGrid$S24,TuneGrid$S25)
  p = TuneGrid$P   # total parameter size
  n_s.3 = TuneGrid$S3
  n = TuneGrid$N   #sample size
  cor = TuneGrid$corr
  
  # get names
  vec_s.1 = str_c("V",c(1:n_s.1))
  col_name = str_c("V", c(1:p))
  
  # get the threshold (set c = 10)
  
  thre = sqrt(log(p)/n) * 10
  
  # Result table:
  DFF <-
    matrix(data = 1:(length(vec_s.1)+2),ncol = length(vec_s.1)+2) %>% 
    as.data.frame()
  names(DFF) <- c(vec_s.1, "No_Predictors","MSE")
  DFF <- DFF %>% 
    as.tibble() %>%
    slice(-1)
  
  TuneGrid_final = TuneGrid %>% slice(-1)
  
  TrueBeta <- matrix(data = 1:length(vec_s.1),ncol = length(vec_s.1)) %>% 
    as.data.frame()
  names(TrueBeta) <- c(str_c(vec_s.1,"_TURE"))
  TrueBeta <- TrueBeta %>% as.tibble() %>% slice(-1)
  
  # Method 1:
  
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
    
    
    # Refill S.2 and S.3 and white noise:
    
    vec_s.2 = str_c("V",c(n_s.1+1:sum(n_s.2)))
    vec_s.3 = str_c("V",c(n_s.1+sum(n_s.2)+1:n_s.3))
    vec_w.n  = col_name[-which(col_name %in% c(vec_s.1,vec_s.2,vec_s.3))]
    
    
    # Null model formula:
    f_fit <- "y~V1"
    
    for (len in 2:length(vec_s.1)){f_fit <- str_c(f_fit, "+", vec_s.1[len])}
    
    # Fit null model:
    
    model <- lm(formula = formula(f_fit), data = data)
    coeff_est <- model$coefficients[2:(length(vec_s.1)+1)]   # fit coeff.
    mse <- mean( (y - model$fitted.values)^2)                # fit mse
    
    # Update the result:
    dfresult <- 
      t(coeff_est %>% as.data.frame()) %>% 
      as.tibble() %>% 
      mutate(No_Predictors = length(vec_s.1))%>% 
      mutate(MSE = mse)
    
    TuneGrid_final <- bind_rows(TuneGrid_final, TuneGrid)
    
    betasim <- beta[1:length(vec_s.1)] %>% 
      matrix(data = ., ncol = length(vec_s.1)) %>% as.data.frame()
    names(betasim) <- names(TrueBeta) <- c(str_c(vec_s.1,"_TURE"))
    betasim <- betasim %>% as.tibble() 
    
    TrueBeta <- bind_rows(TrueBeta, betasim)
    # Do a whole variable adding process in this generated data:
    
    for (j in (length(vec_s.1)+1):min(p,n)){
      
      No_Predictors <- j
      # fit the model
      model <- lm(formula = formula(f_fit), data = data)
      coeff_est <- model$coefficients[2:(length(vec_s.1)+1)]   # fit coeff.
      mse <- mean( (y - model$fitted.values)^2)                # fit mse
      
      dfresult_update <- 
        t(coeff_est %>% as.data.frame()) %>% 
        as.tibble() %>% 
        mutate(No_Predictors = No_Predictors) %>% 
        mutate(MSE = mse)
      
      No_Predictors <- No_Predictors + 1
      
      # add variable:
      # Random select a row in S.3 -> S.2 -> White noise
      
        if (length(vec_s.3) != 0 ){
            add_col_number <- sample(1:length(vec_s.3),1)
            add_col_name   <- vec_s.3[ add_col_number]
            vec_s.3        <- vec_s.3[-add_col_number]
        }else{
          if(length(vec_s.2) != 0){
            add_col_number <- sample(1:length(vec_s.2),1)
            add_col_name   <- vec_s.2[ add_col_number]
            vec_s.2        <- vec_s.2[-add_col_number]
          }else{
            add_col_number <- sample(1:length(vec_w.n),1)
            add_col_name   <- vec_w.n[ add_col_number]
            vec_w.n        <- vec_w.n[-add_col_number]
          }
        }
      
      #Fit model update
      f_fit <- str_c(f_fit, "+",add_col_name)
      
      # Update result
      TuneGrid_final <- bind_rows(TuneGrid_final, TuneGrid)
      TrueBeta <- bind_rows(TrueBeta, betasim)
      dfresult <- bind_rows(dfresult, dfresult_update)
    }
    
    # Save the result:
    
    DFF <- bind_rows(DFF,dfresult)
    
    print(str_c("@@@@  Method 2, Iteration No.",i,"  @@@@"))
    
  }
  
  Filename = str_c("Simresult2/Method2/Method2 - TuneGrid-", TuneGrid$Number ,".csv")
  bind_cols(TuneGrid_final, DFF, TrueBeta) %>% 
    mutate(VaddMethod = rep("Method1",nrow(DFF))) %>% 
    write_csv(Filename)
  
  print(str_c("Finish Method2 TuneGrid No. ", TuneGrid$Number))
  
}

# Add by chance:

Est_Method3 <- function(TuneGrid, Iteration){
  # define number of diff signals
  n_s.1 = TuneGrid$S1
  n_s.2 = c(TuneGrid$S21,TuneGrid$S22,TuneGrid$S23,TuneGrid$S24,TuneGrid$S25)
  p = TuneGrid$P   # total parameter size
  n_s.3 = TuneGrid$S3
  n = TuneGrid$N   #sample size
  cor = TuneGrid$corr
  
  # get names
  vec_s.1 = str_c("V",c(1:n_s.1))
  col_name = str_c("V", c(1:p))
  
  # get the threshold (set c = 10)
  
  thre = sqrt(log(p)/n) * 10
  
  # Result table:
  DFF <-
    matrix(data = 1:(length(vec_s.1)+2),ncol = length(vec_s.1)+2) %>% 
    as.data.frame()
  names(DFF) <- c(vec_s.1, "No_Predictors","MSE")
  DFF <- DFF %>% 
    as.tibble() %>%
    slice(-1)
  
  TuneGrid_final = TuneGrid %>% slice(-1)
  
  TrueBeta <- matrix(data = 1:length(vec_s.1),ncol = length(vec_s.1)) %>% 
    as.data.frame()
  names(TrueBeta) <- c(str_c(vec_s.1,"_TURE"))
  TrueBeta <- TrueBeta %>% as.tibble() %>% slice(-1)
  
  # Method 1:
  
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
    
    
    # Refill S.2 and S.3 and white noise:
    
    vec_s.2 = str_c("V",c(n_s.1+1:sum(n_s.2)))
    vec_s.3 = str_c("V",c(n_s.1+sum(n_s.2)+1:n_s.3))
    vec_w.n  = col_name[-which(col_name %in% c(vec_s.1,vec_s.2,vec_s.3))]
    
    
    # Null model formula:
    f_fit <- "y~V1"
    
    for (len in 2:length(vec_s.1)){f_fit <- str_c(f_fit, "+", vec_s.1[len])}
    
    # Fit null model:
    
    model <- lm(formula = formula(f_fit), data = data)
    coeff_est <- model$coefficients[2:(length(vec_s.1)+1)]   # fit coeff.
    mse <- mean( (y - model$fitted.values)^2)                # fit mse
    
    # Update the result:
    dfresult <- 
      t(coeff_est %>% as.data.frame()) %>% 
      as.tibble() %>% 
      mutate(No_Predictors = length(vec_s.1))%>% 
      mutate(MSE = mse)
    
    TuneGrid_final <- bind_rows(TuneGrid_final, TuneGrid)
    
    betasim <- beta[1:length(vec_s.1)] %>% 
      matrix(data = ., ncol = length(vec_s.1)) %>% as.data.frame()
    names(betasim) <- names(TrueBeta) <- c(str_c(vec_s.1,"_TURE"))
    betasim <- betasim %>% as.tibble() 
    
    TrueBeta <- bind_rows(TrueBeta, betasim)
    # Do a whole variable adding process in this generated data:
    
    for (j in (length(vec_s.1)+1):min(p,n)){
      
      No_Predictors <- j
      # fit the model
      model <- lm(formula = formula(f_fit), data = data)
      coeff_est <- model$coefficients[2:(length(vec_s.1)+1)]   # fit coeff.
      mse <- mean( (y - model$fitted.values)^2)                # fit mse
      
      dfresult_update <- 
        t(coeff_est %>% as.data.frame()) %>% 
        as.tibble() %>% 
        mutate(No_Predictors = No_Predictors) %>% 
        mutate(MSE = mse)
      
      No_Predictors <- No_Predictors + 1
      
      # add variable by chance (S.1 - S.2):
      
      vec23 <- c(vec_s.2, vec_s.3)
      
      if (length(vec23) != 0 ){
        add_col_number <- sample(1:length(vec_s.23),1)
        add_col_name   <- vec_s.23[ add_col_number]
        vec_s.23       <- vec_s.23[-add_col_number]
      }else{
          add_col_number <- sample(1:length(vec_w.n),1)
          add_col_name   <- vec_w.n[ add_col_number]
          vec_w.n        <- vec_w.n[-add_col_number]
        }
      
      
      #Fit model update
      f_fit <- str_c(f_fit, "+",add_col_name)
      
      # Update result
      TuneGrid_final <- bind_rows(TuneGrid_final, TuneGrid)
      TrueBeta <- bind_rows(TrueBeta, betasim)
      dfresult <- bind_rows(dfresult, dfresult_update)
    }
    
    # Save the result:
    
    DFF <- bind_rows(DFF,dfresult)
    
    print(str_c("@@@@  Method 3, Iteration No.",i,"  @@@@"))
    
  }
  
  Filename = str_c("Simresult2/Method3/Method3 - TuneGrid-", TuneGrid$Number ,".csv")
  bind_cols(TuneGrid_final, DFF, TrueBeta) %>% 
    mutate(VaddMethod = rep("Method3",nrow(DFF))) %>% 
    write_csv(Filename)
  
  print(str_c("Finish Method3 TuneGrid No. ", TuneGrid$Number))
  
}

# Add by proportion solved by AIC/Lasso model:

Est_Method4 <- function(TuneGrid, Iteration, PropGrid){
  # define number of diff signals
  n_s.1 = TuneGrid$S1
  n_s.2 = c(TuneGrid$S21,TuneGrid$S22,TuneGrid$S23,TuneGrid$S24,TuneGrid$S25)
  p = TuneGrid$P   # total parameter size
  n_s.3 = TuneGrid$S3
  n = TuneGrid$N   #sample size
  cor = TuneGrid$corr
  
  
  # Data from PropGrid
  
  p.s2 <- PropGrid$S2
  p.s3 <- PropGrid$S3
  
  # get names
  vec_s.1 = str_c("V",c(1:n_s.1))
  col_name = str_c("V", c(1:p))
  
  # get the threshold (set c = 10)
  
  thre = sqrt(log(p)/n) * 10
  
  # Result table:
  DFF <-
    matrix(data = 1:(length(vec_s.1)+2),ncol = length(vec_s.1)+2) %>% 
    as.data.frame()
  names(DFF) <- c(vec_s.1, "No_Predictors","MSE")
  DFF <- DFF %>% 
    as.tibble() %>%
    slice(-1)
  
  TuneGrid <- TuneGrid %>% mutate(prop.s2 = p.s2, prop.s3 = p.s3)
  
  TuneGrid_final = TuneGrid %>% slice(-1)
  
  TrueBeta <- matrix(data = 1:length(vec_s.1),ncol = length(vec_s.1)) %>% 
    as.data.frame()
  names(TrueBeta) <- c(str_c(vec_s.1,"_TURE"))
  TrueBeta <- TrueBeta %>% as.tibble() %>% slice(-1)
  
  # Method 1:
  
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
    
    
    # Refill S.2 and S.3 and white noise:
    
    vec_s.2 = str_c("V",c(n_s.1+1:sum(n_s.2)))
    vec_s.3 = str_c("V",c(n_s.1+sum(n_s.2)+1:n_s.3))
    vec_w.n  = col_name[-which(col_name %in% c(vec_s.1,vec_s.2,vec_s.3))]
    
    
    # Null model formula:
    f_fit <- "y~V1"
    
    for (len in 2:length(vec_s.1)){f_fit <- str_c(f_fit, "+", vec_s.1[len])}
    
    # Fit null model:
    
    model <- lm(formula = formula(f_fit), data = data)
    coeff_est <- model$coefficients[2:(length(vec_s.1)+1)]   # fit coeff.
    mse <- mean( (y - model$fitted.values)^2)                # fit mse
    
    # Update the result:
    dfresult <- 
      t(coeff_est %>% as.data.frame()) %>% 
      as.tibble() %>% 
      mutate(No_Predictors = length(vec_s.1))%>% 
      mutate(MSE = mse)
    
    TuneGrid_final <- bind_rows(TuneGrid_final, TuneGrid)
    
    betasim <- beta[1:length(vec_s.1)] %>% 
      matrix(data = ., ncol = length(vec_s.1)) %>% as.data.frame()
    names(betasim) <- names(TrueBeta) <- c(str_c(vec_s.1,"_TURE"))
    betasim <- betasim %>% as.tibble() 
    
    TrueBeta <- bind_rows(TrueBeta, betasim)
    # Do a whole variable adding process in this generated data:
    
    for (j in (length(vec_s.1)+1):min(p,n)){
      
      No_Predictors <- j
      # fit the model
      model <- lm(formula = formula(f_fit), data = data)
      coeff_est <- model$coefficients[2:(length(vec_s.1)+1)]   # fit coeff.
      mse <- mean( (y - model$fitted.values)^2)                # fit mse
      
      dfresult_update <- 
        t(coeff_est %>% as.data.frame()) %>% 
        as.tibble() %>% 
        mutate(No_Predictors = No_Predictors) %>% 
        mutate(MSE = mse)
      
      No_Predictors <- No_Predictors + 1
      
      # add variable by chance (S.1 - S.2):
      
      vec23 <- c(sample(vec_s.2, size = length(vec_s.2) * p.s2), 
                 sample(vec_s.3, size = length(vec_s.3) * p.s3))
      
      if (length(vec23) != 0 ){
        add_col_number <- sample(1:length(vec_s.23),1)
        add_col_name   <- vec_s.23[ add_col_number]
        vec_s.23       <- vec_s.23[-add_col_number]
      }else{
        add_col_number <- sample(1:length(vec_w.n),1)
        add_col_name   <- vec_w.n[ add_col_number]
        vec_w.n        <- vec_w.n[-add_col_number]
      }
      
      
      #Fit model update
      f_fit <- str_c(f_fit, "+",add_col_name)
      
      # Update result
      TuneGrid_final <- bind_rows(TuneGrid_final, TuneGrid)
      TrueBeta <- bind_rows(TrueBeta, betasim)
      dfresult <- bind_rows(dfresult, dfresult_update)
    }
    
    # Save the result:
    
    DFF <- bind_rows(DFF,dfresult)
    
    print(str_c("@@@@  Method 4, Iteration No.",i,"  @@@@"))
    
  }
  
  Filename = str_c("Simresult2/Method4/Method4 - TuneGrid-", TuneGrid$Number ,".csv")
  bind_cols(TuneGrid_final, DFF, TrueBeta) %>% 
    mutate(VaddMethod = rep("Method4",nrow(DFF))) %>% 
    write_csv(Filename)
  
  print(str_c("Finish Method3 TuneGrid No. ", TuneGrid$Number))
  
}