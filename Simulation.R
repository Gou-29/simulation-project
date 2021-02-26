# simulation .R
source("DataGeneration/AICfunction.R")
source("DataGeneration/LASSOfunction.R")
source("DataGeneration/Estimation.R")

# Following codes are for simulations:

# Setup:
SimGrid <- read_csv("SimGrid.csv")

# Simulation for question 1:

for(SIM in 1:dim(SimGrid)[1]){
  Grid <- SimGrid[SIM,]
  AICMETHOD(Grid, Iteration = 100)
  LASSOMETHOD(Grid, Iteration = 100)
}

# Simulation for question 2:

for(SIM in 1:2){
  Grid <- SimGrid[SIM,]
  Est_Method1(Grid, Iteration = 100)
  Est_Method2(Grid, Iteration = 100)
  Est_Method3(Grid, Iteration = 100)
  Est_Method4(Grid, Iteration = 100, p.s2 = 0.799, p.s3 = 0.806)
  Est_Method4(Grid, Iteration = 100, p.s2 = 0.867, p.s3 = 0.898)
}
