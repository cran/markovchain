## ----ctmcInit, echo = TRUE, message=FALSE, warning=FALSE-----------------
library(markovchain)
energyStates <- c("sigma", "sigma_star")
byRow <- TRUE
gen <- matrix(data = c(-3, 3,
                       1, -1), nrow = 2,
              byrow = byRow, dimnames = list(energyStates, energyStates))
molecularCTMC <- new("ctmc", states = energyStates, 
                 byrow = byRow, generator = gen, 
                 name = "Molecular Transition Model")      

## ----ctmcRandom0, echo = TRUE, message=FALSE, warning=FALSE--------------
statesDist <- c(0.8, 0.2)
rctmc(n = 3, ctmc = molecularCTMC, initDist = statesDist, out.type = "df", include.T0 = FALSE)

## ----ctmcRandom1, echo = TRUE, message=FALSE, warning=FALSE--------------
statesDist <- c(0.8, 0.2)
rctmc(n = Inf, ctmc = molecularCTMC, initDist = statesDist, T = 2)

## ----ctmcSteadyStates, echo = TRUE, message=FALSE, warning=FALSE---------
steadyStates(molecularCTMC)

## ----ctmcFitting, echo = TRUE, message=FALSE, warning=FALSE--------------
data <- list(c("a", "b", "c", "a", "b", "a", "c", "b", "c"), c(0, 0.8, 2.1, 2.4, 4, 5, 5.9, 8.2, 9))
ctmcFit(data)

