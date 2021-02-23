# simulation .R
source("AICfunction.R")
source("LASSOfunction.R")
source("Estimation.R")

# Following codes are for simulations:

# Setup:
SimGrid <- read_csv("SimGrid.csv")

# Simulation for question 1:

#for(SIM in 1:dim(SimGrid)[1]){
#  Grid <- SimGrid[SIM,]
#  AICMETHOD(Grid, Iteration = 100)
#  LASSOMETHOD(Grid, Iteration = 100)
#}

# Simulation for question 2:

for(SIM in 1:144){
  Grid <- SimGrid[SIM,]
  Est_Method1(Grid, Iteration = 100)
  #Est_Method2(Grid, Iteration = 100)
  #Est_Method3(Grid, Iteration = 100)
}
