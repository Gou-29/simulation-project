# simulation .R
source("Generator2.R")
source("AICfunction.R")
source("LASSOfunction.R")

# Following codes are for simulations:

# Setup:
SimGrid <- read_csv("SimGrid.csv")


for(SIM in 1:dim(SimGrid)[1]){
  Grid <- SimGrid[SIM,]
  
  # Simulation of AIC
  #AICMETHOD(Grid, Iteration = 100)
  LASSOMETHOD(Grid, Iteration = 20)
}
