#### AIC Criteria Selection ###

###### Forward Selection ######
n.strong = numeric(n) ## n i the number of iteration
n.weak = numeric(n)
for (i in 1:n)
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


# A sample codes for generating a data with a combinaton of true predictors and null predictors.
n <- 1000
p <- 50
X <- matrix(rnorm(n * p), n, p)
b.true <- rnorm(p) * (runif(p) < 0.2)
cat("True non-zero effects:", which(b.true != 0), "\n")
Y <- 1 + X %*% b.true + rnorm(n)
df <- data.frame(cbind(X, Y))
names(df)[p + 1] <- "y"
# Forward Selection
fit.forward <- step(object = lm(y ~ 1, data = df),
                    scope = formula(lm(y ~ ., data = df)), direction = "forward", k = 2, trace = 0) # A
summary(fit.forward)