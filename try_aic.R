#####################################
#        AIC Criteria Selection     #                     
#        Forward Selection          #
#####################################

set.seed(1234)

# define number of diff signals
n_s.1 = 3
n = 200
n_exp = 500
p = 300

data = Gen_Method1(n_s.1,c(20,20,20),60,p,n)

# get a vector of strong signal
s.1 = c()
for (i in 1:n_s.1){
  s.1 = c(rbind(s.1,paste("x",i,sep = "")))
}

# get the threshold

thre = sqrt(log(p,10)/n)

# get y

y = 

# AIC selection process 

for (i in 1:n_exp)
{
  fit.forward = step(object = lm(y ~ 1, data = df),
                      scope = formula(lm(y ~ ., data = df)), 
                      direction = "forward", k = 2, trace = 0) 
  
 slected <- attr(x = fit.forward$terms, which = "variables")
 n_true = 0
 for (x in selected){
   ifelse(x %in% c("list of strong signal"),n_true = n_true+1,n_true = n_true)
 }
  prop_true = n_true/length(selected)
}


