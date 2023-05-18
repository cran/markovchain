## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(fig.width=8.5, fig.height=6, out.width = "70%")
set.seed(123)
library(knitr)
hook_output <- knit_hooks$get("output")
knit_hooks$set(output = function(x, options) {
   lines <- options$output.lines
   if (is.null(lines)) {
     return(hook_output(x, options))  # pass to default hook
   }
   x <- unlist(strsplit(x, "\n"))
   more <- "..."
   if (length(lines)==1) {        # first n lines
     if (length(x) > lines) {
       # truncate the output, but add ....
       x <- c(head(x, lines), more)
     }
   } else {
     x <- c(more, x[lines], more)
   }
   # paste these lines together
   x <- paste(c(x, ""), collapse = "\n")
   hook_output(x, options)
 })

## ---- load, results='hide', message=FALSE-------------------------------------
library("markovchain")

## ---- showClass, echo=FALSE---------------------------------------------------
showClass("markovchain")
showClass("markovchainList")

## ----mcInitLong---------------------------------------------------------------
weatherStates <- c("sunny", "cloudy", "rain")
byRow <- TRUE
weatherMatrix <- matrix(data = c(0.70, 0.2, 0.1,
                       0.3, 0.4, 0.3,
                       0.2, 0.45, 0.35), byrow = byRow, nrow = 3,
                     dimnames = list(weatherStates, weatherStates))
mcWeather <- new("markovchain", states = weatherStates, byrow = byRow, 
               transitionMatrix = weatherMatrix, name = "Weather")

## ----mcInitShort--------------------------------------------------------------
mcWeather <- new("markovchain", states = c("sunny", "cloudy", "rain"),
                 transitionMatrix = matrix(data = c(0.70, 0.2, 0.1,
                       0.3, 0.4, 0.3,
                       0.2, 0.45, 0.35), byrow = byRow, nrow = 3), 
                 name = "Weather")

## ----defaultMc----------------------------------------------------------------
defaultMc <- new("markovchain")

## ----intromcList--------------------------------------------------------------
mcList <- new("markovchainList", markovchains = list(mcWeather, defaultMc), 
		          name = "A list of Markov chains")

## ----operations---------------------------------------------------------------
initialState <- c(0, 1, 0)
after2Days <- initialState * (mcWeather * mcWeather)
after7Days <- initialState * (mcWeather ^ 7)
after2Days
round(after7Days, 3)

## ----operations2--------------------------------------------------------------
initialState <- c(0, 1, 0)
after2Days <- (t(mcWeather) * t(mcWeather)) * initialState
after7Days <- (t(mcWeather) ^ 7) * initialState
after2Days
round(after7Days, 3)

## ----fval---------------------------------------------------------------------
fvals<-function(mchain,initialstate,n) {
  out<-data.frame()
  names(initialstate)<-names(mchain)
  for (i in 0:n)
  {
    iteration<-initialstate*mchain^(i)
    out<-rbind(out,iteration)
  }
  out<-cbind(out, i=seq(0,n))
  out<-out[,c(4,1:3)]
  return(out)
}
fvals(mchain=mcWeather,initialstate=c(90,5,5),n=4)

## ----otherMethods-------------------------------------------------------------
states(mcWeather)
names(mcWeather)
dim(mcWeather)

## ----otherMethods2------------------------------------------------------------
name(mcWeather)
name(mcWeather) <- "New Name"
name(mcWeather)

## ----sortMethod---------------------------------------------------------------
markovchain:::sort(mcWeather)

## ----transProb----------------------------------------------------------------
transitionProbability(mcWeather, "cloudy", "rain")
mcWeather[2,3]

## ----printAndShow-------------------------------------------------------------
print(mcWeather)
show(mcWeather)

## ----mcPlot, echo=FALSE, fig.cap="Weather example. Markov chain plot"---------
if (requireNamespace("igraph", quietly = TRUE)) {
  library(igraph)
  plot(mcWeather,layout = layout.fruchterman.reingold)
  } else {
  message("igraph unavailable")
  }

## ----mcPlotdiagram, echo=FALSE, fig.cap="Weather example. Markov chain plot with diagram"----
if (requireNamespace("diagram", quietly = TRUE)) {
  library(diagram)
  plot(mcWeather, package="diagram", box.size = 0.04)
  } else {
  message("diagram unavailable")
  }

## ----exportImport1------------------------------------------------------------
mcDf <- as(mcWeather, "data.frame")
mcNew <- as(mcDf, "markovchain")
mcDf
mcIgraph <- as(mcWeather, "igraph")

## ----exportImport2------------------------------------------------------------
if (requireNamespace("msm", quietly = TRUE)) {
require(msm)
Q <- rbind ( c(0, 0.25, 0, 0.25),
             c(0.166, 0, 0.166, 0.166),
             c(0, 0.25, 0, 0.25),
             c(0, 0, 0, 0) )
cavmsm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = Q, death = 4)
msmMc <- as(cavmsm, "markovchain")
msmMc
  } else {
  message("msm unavailable")
  }

## ----exporImport3-------------------------------------------------------------
if (requireNamespace("etm", quietly = TRUE)) {
library(etm)
data(sir.cont)
sir.cont <- sir.cont[order(sir.cont$id, sir.cont$time), ]
for (i in 2:nrow(sir.cont)) {
  if (sir.cont$id[i]==sir.cont$id[i-1]) {
    if (sir.cont$time[i]==sir.cont$time[i-1]) {
      sir.cont$time[i-1] <- sir.cont$time[i-1] - 0.5
    }
  }
}
tra <- matrix(ncol=3,nrow=3,FALSE)
tra[1, 2:3] <- TRUE
tra[2, c(1, 3)] <- TRUE
tr.prob <- etm::etm(sir.cont, c("0", "1", "2"), tra, "cens", 1)
tr.prob
etm2mc<-as(tr.prob, "markovchain")
etm2mc
  } else {
  message("etm unavailable")
}

## ----fromAndTo, echo=FALSE, fig.cap="The markovchain methods for import and export"----
library(igraph)
importExportGraph<-graph.formula(dataframe++markovchain,markovchain-+igraph,
                                 markovchain++matrix,table-+markovchain,msm-+markovchain,etm-+markovchain,
                                 markovchain++sparseMatrix)
plot(importExportGraph,main="Import - Export from and to markovchain objects")

## ----exportImport4------------------------------------------------------------
myMatr<-matrix(c(.1,.8,.1,.2,.6,.2,.3,.4,.3), byrow=TRUE, ncol=3)
myMc<-as(myMatr, "markovchain")
myMc

## ----cchcMcList---------------------------------------------------------------
stateNames = c("H", "I", "D")
Q0 <- new("markovchain", states = stateNames, 
        transitionMatrix =matrix(c(0.7, 0.2, 0.1,0.1, 0.6, 0.3,0, 0, 1), 
        byrow = TRUE, nrow = 3), name = "state t0")
Q1 <- new("markovchain", states = stateNames, 
        transitionMatrix = matrix(c(0.5, 0.3, 0.2,0, 0.4, 0.6,0, 0, 1), 
        byrow = TRUE, nrow = 3), name = "state t1")
Q2 <- new("markovchain", states = stateNames, 
        transitionMatrix = matrix(c(0.3, 0.2, 0.5,0, 0.2, 0.8,0, 0, 1), 
        byrow = TRUE,nrow = 3), name = "state t2")
Q3 <- new("markovchain", states = stateNames, 
          transitionMatrix = matrix(c(0, 0, 1, 0, 0, 1, 0, 0, 1), 
        byrow = TRUE, nrow = 3), name = "state t3")
mcCCRC <- new("markovchainList",markovchains = list(Q0,Q1,Q2,Q3), 
      name = "Continuous Care Health Community")
print(mcCCRC)

## ----cchcMcList2--------------------------------------------------------------
mcCCRC[[1]]
dim(mcCCRC)

## ----conditionalDistr---------------------------------------------------------
conditionalDistribution(mcWeather, "sunny")

## ----steadyStates-------------------------------------------------------------
steadyStates(mcWeather)

## ----gamblerRuin--------------------------------------------------------------
gamblerRuinMarkovChain <- function(moneyMax, prob = 0.5) {
  m <- markovchain:::zeros(moneyMax + 1)
  m[1,1] <- m[moneyMax + 1,moneyMax + 1] <- 1
  states <- as.character(0:moneyMax)
  rownames(m) <- colnames(m) <- states
  
  for(i in 2:moneyMax){ 
    m[i,i-1] <- 1 - prob
    m[i, i + 1] <- prob   
  }
  
  new("markovchain", transitionMatrix = m, 
      name = paste("Gambler ruin", moneyMax, "dim", sep = " "))
}

mcGR4 <- gamblerRuinMarkovChain(moneyMax = 4, prob = 0.5)
steadyStates(mcGR4)

## ----absorbingStates----------------------------------------------------------
absorbingStates(mcGR4)
absorbingStates(mcWeather)

## ----renaldoMatrix1-----------------------------------------------------------
P <- markovchain:::zeros(10)
P[1, c(1, 3)] <- 1/2;
P[2, 2] <- 1/3; P[2,7] <- 2/3;
P[3, 1] <- 1;
P[4, 5] <- 1;
P[5, c(4, 5, 9)] <- 1/3;
P[6, 6] <- 1;
P[7, 7] <- 1/4; P[7,9] <- 3/4;
P[8, c(3, 4, 8, 10)] <- 1/4;
P[9, 2] <- 1;
P[10, c(2, 5, 10)] <- 1/3;
rownames(P) <- letters[1:10] 
colnames(P) <- letters[1:10]
probMc <- new("markovchain", transitionMatrix = P, 
              name = "Probability MC")
summary(probMc)

## ----transientStates----------------------------------------------------------
transientStates(probMc)

## ----probMc2Canonic-----------------------------------------------------------
probMcCanonic <- canonicForm(probMc)
probMc
probMcCanonic

## ----isAccessible-------------------------------------------------------------
is.accessible(object = probMc, from = "a", to = "c")
is.accessible(object = probMc, from = "g", to = "c")

## ----periodicity--------------------------------------------------------------
E <- matrix(0, nrow = 4, ncol = 4)
E[1, 2] <- 1
E[2, 1] <- 1/3; E[2, 3] <- 2/3
E[3,2] <- 1/4; E[3, 4] <- 3/4
E[4, 3] <- 1

mcE <- new("markovchain", states = c("a", "b", "c", "d"), 
		transitionMatrix = E, 
		name = "E")
is.irreducible(mcE)
period(mcE)

## ----mathematica9Mc-----------------------------------------------------------
mathematicaMatr <- markovchain:::zeros(5)
mathematicaMatr[1,] <- c(0, 1/3, 0, 2/3, 0)
mathematicaMatr[2,] <- c(1/2, 0, 0, 0, 1/2)
mathematicaMatr[3,] <- c(0, 0, 1/2, 1/2, 0)
mathematicaMatr[4,] <- c(0, 0, 1/2, 1/2, 0)
mathematicaMatr[5,] <- c(0, 0, 0, 0, 1)
statesNames <- letters[1:5]
mathematicaMc <- new("markovchain", transitionMatrix = mathematicaMatr,
                   name = "Mathematica MC", states = statesNames)

## ----mcMathematics, fig=TRUE, echo=FALSE, fig.align='center', fig.cap="Mathematica 9 example. Markov chain plot."----
plot(mathematicaMc, layout = layout.fruchterman.reingold)

## ----mathematica9MC, echo=FALSE-----------------------------------------------
summary(mathematicaMc)

## ----fpTime1, eval=FALSE------------------------------------------------------
#  .firstpassageKernel <- function(P, i, n){
#    G <- P
#    H <- P[i,]
#    E <- 1 - diag(size(P)[2])
#    for (m in 2:n) {
#      G <- P %*% (G * E)
#      H <- rbind(H, G[i,])
#    }
#    return(H)
#  }

## ----fpTime2------------------------------------------------------------------
firstPassagePdF <- firstPassage(object = mcWeather, state = "sunny", 
                                n = 10)
firstPassagePdF[3, 3]

## ----mfpt1--------------------------------------------------------------------
meanFirstPassageTime(mcWeather)

## ----mfpt2--------------------------------------------------------------------
meanFirstPassageTime(mcWeather,"rain")

## ----mfpt3--------------------------------------------------------------------
firstPassagePdF.long <- firstPassage(object = mcWeather, state = "sunny",  n = 100)
sum(firstPassagePdF.long[,"rain"] * 1:100)

## ----mrt-weather--------------------------------------------------------------
meanRecurrenceTime(mcWeather)

## ----mrt-probMc---------------------------------------------------------------
recurrentStates(probMc)
meanRecurrenceTime(probMc)

## ----data-drunkard------------------------------------------------------------
drunkProbs <- markovchain:::zeros(5)
drunkProbs[1,1] <- drunkProbs[5,5] <- 1
drunkProbs[2,1] <- drunkProbs[2,3] <- 1/2
drunkProbs[3,2] <- drunkProbs[3,4] <- 1/2
drunkProbs[4,3] <- drunkProbs[4,5] <- 1/2

drunkMc <- new("markovchain", transitionMatrix = drunkProbs)
drunkMc

## ----rs-drunkard--------------------------------------------------------------
recurrentStates(drunkMc)

## ----ts-drunkard--------------------------------------------------------------
transientStates(drunkMc)

## ----ap-drunkard--------------------------------------------------------------
absorptionProbabilities(drunkMc)

## ----at-drunkard--------------------------------------------------------------
meanAbsorptionTime(drunkMc)

## -----------------------------------------------------------------------------
committorAB(mcWeather,3,1)

## ----hitting-data-------------------------------------------------------------
M <- markovchain:::zeros(5)
M[1,1] <- M[5,5] <- 1
M[2,1] <- M[2,3] <- 1/2
M[3,2] <- M[3,4] <- 1/2
M[4,2] <- M[4,5] <- 1/2

hittingTest <- new("markovchain", transitionMatrix = M)
hittingProbabilities(hittingTest)

## ----hitting-probabilities----------------------------------------------------
hittingProbabilities(hittingTest)

## ----hitting-weather----------------------------------------------------------
hittingProbabilities(mcWeather)

## ----simulatingAMarkovChain---------------------------------------------------
weathersOfDays <- rmarkovchain(n = 365, object = mcWeather, t0 = "sunny")
weathersOfDays[1:30]

## ----simulatingAListOfMarkovChain---------------------------------------------
patientStates <- rmarkovchain(n = 5, object = mcCCRC, t0 = "H", 
                              include.t0 = TRUE)
patientStates[1:10,]

## ----fitMcbyMLE2--------------------------------------------------------------
weatherFittedMLE <- markovchainFit(data = weathersOfDays, method = "mle",name = "Weather MLE")
weatherFittedMLE$estimate
weatherFittedMLE$standardError

## ----fitMcbyLAPLACE-----------------------------------------------------------
weatherFittedLAPLACE <- markovchainFit(data = weathersOfDays, 
                                    method = "laplace", laplacian = 0.01,
                                    name = "Weather LAPLACE")
weatherFittedLAPLACE$estimate

## ----fitSequenceMatrix--------------------------------------------------------
createSequenceMatrix(stringchar = weathersOfDays)

## ----fitSequenceMatrix2-------------------------------------------------------
myMatr<-matrix(c("a","b","b","a","a","b","b","b","b","a","a","a","b","a"),ncol=2)
createSequenceMatrix(stringchar = myMatr,toRowProbs = TRUE)

## ----fitMcbyBootStrap1--------------------------------------------------------
weatherFittedBOOT <- markovchainFit(data = weathersOfDays, 
                                    method = "bootstrap", nboot = 20)
weatherFittedBOOT$estimate
weatherFittedBOOT$standardError

## ----fitMcbyBootStrap2, eval=FALSE--------------------------------------------
#  weatherFittedBOOTParallel <- markovchainFit(data = weathersOfDays,
#                                      method = "bootstrap", nboot = 200,
#                                      parallel = TRUE)
#  weatherFittedBOOTParallel$estimate
#  weatherFittedBOOTParallel$standardError

## ----fitMcbyBootStrap3, eval=FALSE--------------------------------------------
#  RcppParallel::setNumThreads(2)

## ----fitMcbyMLE1--------------------------------------------------------------
weatherFittedMLE$logLikelihood
weatherFittedBOOT$logLikelihood

## ----confint------------------------------------------------------------------
weatherFittedMLE$confidenceInterval
weatherFittedBOOT$confidenceInterval

## ----multinomial--------------------------------------------------------------
multinomialConfidenceIntervals(transitionMatrix = 
        weatherFittedMLE$estimate@transitionMatrix, 
        countsTransitionMatrix = createSequenceMatrix(weathersOfDays))

## ----fitMclists---------------------------------------------------------------
data(holson)
singleMc<-markovchainFit(data=holson[,2:12],name="holson")

## ----fitMclistsFit1, output.lines=20------------------------------------------
mcListFit<-markovchainListFit(data=holson[,2:6],name="holson")
mcListFit$estimate

## ----fitMclistsFit2-----------------------------------------------------------
c1<-c("a","b","a","a","c","c","a")
c2<-c("b")
c3<-c("c","a","a","c")
c4<-c("b","a","b","a","a","c","b")
c5<-c("a","a","c",NA)
c6<-c("b","c","b","c","a")
mylist<-list(c1,c2,c3,c4,c5,c6)
mylistMc<-markovchainFit(data=mylist)
mylistMc

## ----fitAMarkovChainListfromAlist, output.lines=15----------------------------
markovchainListFit(data=mylist)

## ----markovchainPredict-------------------------------------------------------
predict(object = weatherFittedMLE$estimate, newdata = c("cloudy", "sunny"),
        n.ahead = 3)

## ----markovchainListPredict---------------------------------------------------
predict(mcCCRC, newdata = c("H", "H"), n.ahead = 5)

## ----markovchainListPredict2--------------------------------------------------
predict(mcCCRC, newdata = c("H", "H"), n.ahead = 5, continue = TRUE)

## ----test1--------------------------------------------------------------------
sample_sequence<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", 
                   "b", "a", "a", "b", "b", "b", "a")
verifyMarkovProperty(sample_sequence)

## ----test2--------------------------------------------------------------------
data(rain)
assessOrder(rain$rain)

## ----test3--------------------------------------------------------------------
assessStationarity(rain$rain, 10)

## ----divergence1--------------------------------------------------------------
sequence<-c(0,1,2,2,1,0,0,0,0,0,0,1,2,2,2,1,0,0,1,0,0,0,0,0,0,1,1,
2,0,0,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,2,1,0,
0,2,1,0,0,0,0,0,0,1,1,1,2,2,0,0,2,1,1,1,1,2,1,1,1,1,1,1,1,1,1,0,2,
0,1,1,0,0,0,1,2,2,0,0,0,0,0,0,2,2,2,1,1,1,1,0,1,1,1,1,0,0,2,1,1,
0,0,0,0,0,2,2,1,1,1,1,1,2,1,2,0,0,0,1,2,2,2,0,0,0,1,1)
mc=matrix(c(5/8,1/4,1/8,1/4,1/2,1/4,1/4,3/8,3/8),byrow=TRUE, nrow=3)
rownames(mc)<-colnames(mc)<-0:2; theoreticalMc<-as(mc, "markovchain")
verifyEmpiricalToTheoretical(data=sequence,object=theoreticalMc)

## ----divergence2--------------------------------------------------------------
data(kullback)
verifyHomogeneity(inputList=kullback,verbose=TRUE)

## ----rCtmcInit----------------------------------------------------------------
energyStates <- c("sigma", "sigma_star")
byRow <- TRUE
gen <- matrix(data = c(-3, 3,
                       1, -1), nrow = 2,
              byrow = byRow, dimnames = list(energyStates, energyStates))
molecularCTMC <- new("ctmc", states = energyStates, 
                 byrow = byRow, generator = gen, 
                 name = "Molecular Transition Model")      

## ----rctmcRandom0-------------------------------------------------------------
statesDist <- c(0.8, 0.2)
rctmc(n = 3, ctmc = molecularCTMC, initDist = statesDist, out.type = "df", include.T0 = FALSE)

## ----ctmcRandom1--------------------------------------------------------------
statesDist <- c(0.8, 0.2)
rctmc(n = Inf, ctmc = molecularCTMC, initDist = statesDist, T = 2)

## ----rctmcSteadyStates--------------------------------------------------------
steadyStates(molecularCTMC)

## ----rctmcFitting-------------------------------------------------------------
data <- list(c("a", "b", "c", "a", "b", "a", "c", "b", "c"), 
             c(0, 0.8, 2.1, 2.4, 4, 5, 5.9, 8.2, 9))
ctmcFit(data)

## ----mcWeatherQ---------------------------------------------------------------
mcWeatherQ <- expm::logm(mcWeather@transitionMatrix,method='Eigen')
mcWeatherQ

## ----mcWeatherHalfDay---------------------------------------------------------
mcWeatherHalfDayTM <- expm::expm(mcWeatherQ*.5)
mcWeatherHalfDay <- new("markovchain",transitionMatrix=mcWeatherHalfDayTM,name="Half Day Weather Transition Matrix")
mcWeatherHalfDay

## ----ctmcd1-------------------------------------------------------------------
if(requireNamespace(package='ctmcd', quietly = TRUE)) {
require(ctmcd)
require(expm)
#defines a function to transform a GM into a TM
gm_to_markovchain<-function(object, t=1) {
  if(!(class(object) %in% c("gm","matrix","Matrix")))
    stop("Error! Expecting either a matrix or a gm object")
  if ( class(object) %in% c("matrix","Matrix")) generator_matrix<-object else generator_matrix<-as.matrix(object[["par"]])
  #must add importClassesFrom("markovchain",markovchain) in the NAMESPACE
  #must add importFrom(expm, "expm")
  transitionMatrix<-expm(generator_matrix*t)
  out<-as(transitionMatrix,"markovchain")
  return(out)
}
#loading ctmcd dataset
data(tm_abs)
gm0=matrix(1,8,8) #initializing
diag(gm0)=0
diag(gm0)=-rowSums(gm0)
gm0[8,]=0
gmem=gm(tm_abs,te=1,method="EM",gmguess=gm0) #estimating GM
mc_at_2=gm_to_markovchain(object=gmem, t=2) #converting to TM at time 2
} else {
  warning('package ctmcd unavailable')
}

## ----pseudobayes--------------------------------------------------------------
pseudoBayesEstimator <- function(raw, apriori){
  v_i <- rowSums(raw) 
  K_i <- numeric(nrow(raw))
  sumSquaredY <- rowSums(raw^2)
  #get numerator
  K_i_num <- v_i^2-sumSquaredY
  #get denominator
  VQ <- matrix(0,nrow= nrow(apriori),ncol=ncol(apriori))
  for (i in 1:nrow(VQ)) {
    VQ[i,]<-v_i[i]*apriori[i,]
  }
  
  K_i_den<-rowSums((raw - VQ)^2)
  
  K_i <- K_i_num/K_i_den
  
  #get the alpha vector
  alpha <- K_i / (v_i+K_i)
  
  #empirical transition matrix
  Emp<-raw/rowSums(raw)
  
  #get the estimate
  out<-matrix(0, nrow= nrow(raw),ncol=ncol(raw))
  for (i in 1:nrow(out)) {
    out[i,]<-alpha[i]*apriori[i,]+(1-alpha[i])*Emp[i,]
  }
  return(out)
}

## ----pseudobayes2-------------------------------------------------------------
trueMc<-as(matrix(c(0.1, .9,.7,.3),nrow = 2, byrow = 2),"markovchain")
aprioriMc<-as(matrix(c(0.5, .5,.5,.5),nrow = 2, byrow = 2),"markovchain")

smallSample<-rmarkovchain(n=20,object = trueMc)
smallSampleRawTransitions<-createSequenceMatrix(stringchar = smallSample)
pseudoBayesEstimator(
  raw = smallSampleRawTransitions, 
  apriori = aprioriMc@transitionMatrix
) - trueMc@transitionMatrix

biggerSample<-rmarkovchain(n=100,object = trueMc)
biggerSampleRawTransitions<-createSequenceMatrix(stringchar = biggerSample)
pseudoBayesEstimator(
  raw = biggerSampleRawTransitions,
  apriori = aprioriMc@transitionMatrix
) - trueMc@transitionMatrix

bigSample<-rmarkovchain(n=1000,object = trueMc)
bigSampleRawTransitions<-createSequenceMatrix(stringchar = bigSample)
pseudoBayesEstimator(
  raw = bigSampleRawTransitions,
  apriori = aprioriMc@transitionMatrix
) - trueMc@transitionMatrix

## ----loadAndDoExample---------------------------------------------------------
weatherStates <- c("sunny", "cloudy", "rain")
byRow <- TRUE
weatherMatrix <- matrix(data = c(0.7, 0.2, 0.1, 
                                 0.3, 0.4, 0.3, 
                                 0.2, 0.4, 0.4), 
                        byrow = byRow, nrow = 3, 
                        dimnames = list(weatherStates, weatherStates))
mcWeather <- new("markovchain", states = weatherStates, 
                 byrow = byRow, transitionMatrix = weatherMatrix, 
                 name = "Weather")      
weathersOfDays <- rmarkovchain(n = 365, object = mcWeather, t0 = "sunny")

## ----MAPFit-------------------------------------------------------------------
hyperMatrix<-matrix(c(1, 1, 2, 
                      3, 2, 1,
                      2, 2, 3), 
                    nrow = 3, byrow = TRUE,
                    dimnames = list(weatherStates,weatherStates))
markovchainFit(weathersOfDays[1:200], method = "map", 
               confidencelevel = 0.92, hyperparam = hyperMatrix)
predictiveDistribution(weathersOfDays[1:200], 
                       weathersOfDays[201:365],hyperparam = hyperMatrix) 

## ----MAPFit2------------------------------------------------------------------
hyperMatrix2<- hyperMatrix[c(2,3,1), c(2,3,1)]
markovchainFit(weathersOfDays[1:200], method = "map", 
               confidencelevel = 0.92, hyperparam = hyperMatrix2)
predictiveDistribution(weathersOfDays[1:200], 
                       weathersOfDays[201:365],hyperparam = hyperMatrix2)

## ----inferHyperparam----------------------------------------------------------
inferHyperparam(transMatr = weatherMatrix, scale = c(10, 10, 10))

## ----inferHyperparam2---------------------------------------------------------
inferHyperparam(data = weathersOfDays[1:15])

## ----inferHyperparam3---------------------------------------------------------
hyperMatrix3 <- inferHyperparam(transMatr = weatherMatrix, 
                                scale = c(10, 10, 10))
hyperMatrix3 <- hyperMatrix3$scaledInference
hyperMatrix4 <- inferHyperparam(data = weathersOfDays[1:15])
hyperMatrix4 <- hyperMatrix4$dataInference

## ----MAPandMLE----------------------------------------------------------------
data(preproglucacon)
preproglucacon <- preproglucacon[[2]]
MLEest <- markovchainFit(preproglucacon, method = "mle")
MAPest <- markovchainFit(preproglucacon, method = "map")
MLEest$estimate
MAPest$estimate

## ----weatPred1----------------------------------------------------------------

mcWP <- new("markovchain", states = c("rainy", "nice", "snowy"),
         transitionMatrix = matrix(c(0.5, 0.25, 0.25,
                                   0.5, 0, 0.5,
                                   0.25,0.25,0.5), byrow = T, nrow = 3))

## ----weatPred2----------------------------------------------------------------
W0 <- t(as.matrix(c(0, 1, 0)))
W1 <- W0 * mcWP; W1

W2 <- W0 * (mcWP ^ 2); W2

W3 <- W0 * (mcWP ^ 3); W3

## ----weatPred3----------------------------------------------------------------
W7 <- W0 * (mcWP ^ 7)
W7

## ----weatPred4----------------------------------------------------------------
q <- steadyStates(mcWP)
q

## ----weatPred5----------------------------------------------------------------
R0 <- t(as.matrix(c(1, 0, 0)))
R7 <- R0 * (mcWP ^ 7); R7

S0 <- t(as.matrix(c(0, 0, 1)))
S7 <- S0 * (mcWP ^ 7); S7

## ----Alofi1-------------------------------------------------------------------
data("rain", package = "markovchain")
table(rain$rain)

## ----Alofi2-------------------------------------------------------------------
mcAlofi <- markovchainFit(data = rain$rain, name = "Alofi MC")$estimate
mcAlofi

## ----Alofi3-------------------------------------------------------------------
steadyStates(mcAlofi)

## ----ratings1-----------------------------------------------------------------
rc <- c("AAA", "AA", "A", "BBB", "BB", "B", "CCC", "D")
creditMatrix <- matrix(
  c(90.81, 8.33, 0.68, 0.06, 0.08, 0.02, 0.01, 0.01,
    0.70, 90.65, 7.79, 0.64, 0.06, 0.13, 0.02, 0.01,
    0.09, 2.27, 91.05, 5.52, 0.74, 0.26, 0.01, 0.06,
    0.02, 0.33, 5.95, 85.93, 5.30, 1.17, 1.12, 0.18,
    0.03, 0.14, 0.67, 7.73, 80.53, 8.84, 1.00, 1.06,
    0.01, 0.11, 0.24, 0.43, 6.48, 83.46, 4.07, 5.20,
    0.21, 0, 0.22, 1.30, 2.38, 11.24, 64.86, 19.79,
    0, 0, 0, 0, 0, 0, 0, 100
   )/100, 8, 8, dimnames = list(rc, rc), byrow = TRUE)

## ----ratings2-----------------------------------------------------------------
creditMc <- new("markovchain", transitionMatrix = creditMatrix, 
                name = "S&P Matrix")
absorbingStates(creditMc)

## ----economicAnalysis1--------------------------------------------------------
statesNames <- c("customer", "non customer")
P <- markovchain:::zeros(2); P[1, 1] <- .9; P[1, 2] <- .1; P[2, 2] <- .95; P[2, 1] <- .05;
rownames(P) <- statesNames; colnames(P) <- statesNames
mcP <- new("markovchain", transitionMatrix = P, name = "Telephone company")
M <- markovchain:::zeros(2); M[1, 1] <- -20; M[1, 2] <- -30; M[2, 1] <- -40; M[2, 2] <- 0

## ----economicAnalysis2--------------------------------------------------------
c1 <- 100 + conditionalDistribution(mcP, state = "customer") %*% M[1,]
c2 <- 0 + conditionalDistribution(mcP, state = "non customer") %*% M[2,]

## ----economicAnalysis3--------------------------------------------------------
as.numeric((c(1, 0)* mcP ^ 5) %*% (as.vector(c(c1, c2))))

## ----bonusMalus1--------------------------------------------------------------

getBonusMalusMarkovChain <- function(lambda) {
	bmMatr <- markovchain:::zeros(5)
	bmMatr[1, 1] <- dpois(x = 0, lambda)
	bmMatr[1, 3] <- dpois(x = 1, lambda)
	bmMatr[1, 5] <- 1 - ppois(q = 1, lambda)
	
	bmMatr[2, 1] <- dpois(x = 0, lambda)
	bmMatr[2, 4] <- dpois(x = 1, lambda)
	bmMatr[2, 5] <- 1 - ppois(q = 1, lambda)
	
	bmMatr[3, 2] <- dpois(x = 0, lambda)
	bmMatr[3, 5] <- 1 - dpois(x=0, lambda)
 
	bmMatr[4, 3] <- dpois(x = 0, lambda)
	bmMatr[4, 5] <- 1 - dpois(x = 0, lambda)
  
	bmMatr[5, 4] <- dpois(x = 0, lambda)
	bmMatr[5, 5] <- 1 - dpois(x = 0, lambda)
	stateNames <- as.character(1:5)
	out <- new("markovchain", transitionMatrix = bmMatr, 
             states = stateNames, name = "BM Matrix")
	return(out)
}

## ----bonusMalus2--------------------------------------------------------------
bmMc <- getBonusMalusMarkovChain(0.05)
as.numeric(steadyStates(bmMc))

## ----bonusMalus3--------------------------------------------------------------
sum(as.numeric(steadyStates(bmMc)) * c(0.5, 0.7, 0.9, 1, 1.25))

## ----healthIns6---------------------------------------------------------------
ltcDemoPath<-system.file("extdata", "ltdItaData.txt", 
                         package = "markovchain")
ltcDemo<-read.table(file = ltcDemoPath, header=TRUE, 
                    sep = ";", dec = ".")
head(ltcDemo)

## ----healthIns7---------------------------------------------------------------
ltcDemo<-transform(ltcDemo,
                   pIA=0,
                   pII=1-pID,
                   pDD=1,
                   pDA=0,
                   pDI=0)

## ----healthIns8---------------------------------------------------------------
possibleStates<-c("A","I","D")
getMc4Age<-function(age) {
  transitionsAtAge<-ltcDemo[ltcDemo$age==age,]
  
  myTransMatr<-matrix(0, nrow=3,ncol = 3,
                      dimnames = list(possibleStates, possibleStates))
  myTransMatr[1,1]<-transitionsAtAge$pAA[1]
  myTransMatr[1,2]<-transitionsAtAge$pAI[1]
  myTransMatr[1,3]<-transitionsAtAge$pAD[1]
  myTransMatr[2,2]<-transitionsAtAge$pII[1]
  myTransMatr[2,3]<-transitionsAtAge$pID[1]
  myTransMatr[3,3]<-1
  
  myMc<-new("markovchain", transitionMatrix = myTransMatr,
            states = possibleStates,
            name = paste("Age",age,"transition matrix"))
  
  return(myMc)

}

## ----healthIns8-prob----------------------------------------------------------
getFullTransitionTable<-function(age){
  ageSequence<-seq(from=age, to=120)
  k=1
  myList=list()
  for ( i in ageSequence) {
    mc_age_i<-getMc4Age(age = i)
    myList[[k]]<-mc_age_i
    k=k+1
  }
  myMarkovChainList<-new("markovchainList", markovchains = myList,
                         name = paste("TransitionsSinceAge", age, sep = ""))
  return(myMarkovChainList)
}
transitionsSince100<-getFullTransitionTable(age=100)

## ----healthIns9---------------------------------------------------------------
rmarkovchain(n = 10, object = transitionsSince100,
             what = "matrix", t0 = "A", include.t0 = TRUE)

## ----healthIns10--------------------------------------------------------------
transitionsSince80<-getFullTransitionTable(age=80)
lifeTrajectories<-rmarkovchain(n=1e3, object=transitionsSince80,
                               what="matrix",t0="A",include.t0=TRUE)
temp<-matrix(0,nrow=nrow(lifeTrajectories),ncol = ncol(lifeTrajectories))
temp[lifeTrajectories=="I"]<-1
expected_period_disabled<-mean(rowSums((temp)))
expected_period_disabled

## ----healthIns11--------------------------------------------------------------
mean(rowMeans(12000*temp%*%( matrix((1+0.02)^-seq(from=0, to=ncol(temp)-1)))))

## ----blandenEtAlii------------------------------------------------------------
data("blanden")
mobilityMc <- as(blanden, "markovchain")
mobilityMc

## ----mobility, fig=TRUE, echo=FALSE, fig.align='center', fig.cap="1970 UK cohort mobility data."----
plot(mobilityMc, main = '1970 mobility',vertex.label.cex = 2,
		layout = layout.fruchterman.reingold)

## ----blandenEtAlii3-----------------------------------------------------------
round(steadyStates(mobilityMc), 2)

## ----preproglucacon1----------------------------------------------------------
data("preproglucacon", package = "markovchain")

## ----preproglucacon2----------------------------------------------------------
mcProtein <- markovchainFit(preproglucacon$preproglucacon, 
                          name = "Preproglucacon MC")$estimate
mcProtein

## ----epid1--------------------------------------------------------------------
craigSendiMatr <- matrix(c(682, 33, 25,
              154, 64, 47,
              19, 19, 43), byrow = T, nrow = 3)
hivStates <- c("0-49", "50-74", "75-UP")
rownames(craigSendiMatr) <- hivStates
colnames(craigSendiMatr) <- hivStates
craigSendiTable <- as.table(craigSendiMatr)
mcM6 <- as(craigSendiTable, "markovchain")
mcM6@name <- "Zero-Six month CD4 cells transition"
mcM6

## ----epid2--------------------------------------------------------------------
eig <- eigen(mcM6@transitionMatrix)
D <- diag(eig$values)

## ----epid3--------------------------------------------------------------------
V <- eig$vectors 
V %*% D %*% solve(V)
d <- D ^ (1/6)
M <- V %*% d %*% solve(V)
mcM1 <- new("markovchain", transitionMatrix = M, states = hivStates)

