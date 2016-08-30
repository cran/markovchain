### R code from vignette source 'an_introduction_to_markovchain_package.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: setup
###################################################
	options(prompt = "R> ", continue = "+  ", 
			width = 70, useFancyQuotes = FALSE)
	set.seed(123)


###################################################
### code chunk number 2: load
###################################################
library("markovchain")


###################################################
### code chunk number 3: showClass
###################################################
showClass("markovchain")
showClass("markovchainList")


###################################################
### code chunk number 4: mcInitLong
###################################################
weatherStates <- c("sunny", "cloudy", "rain")
byRow <- TRUE
weatherMatrix <- matrix(data = c(0.70, 0.2, 0.1,
                       0.3, 0.4, 0.3,
                       0.2, 0.45, 0.35), byrow = byRow, nrow = 3,
                     dimnames = list(weatherStates, weatherStates))
mcWeather <- new("markovchain", states = weatherStates, byrow = byRow, 
               transitionMatrix = weatherMatrix, name = "Weather")


###################################################
### code chunk number 5: mcInitLong
###################################################
mcWeather <- new("markovchain", states = c("sunny", "cloudy", "rain"),
                 transitionMatrix = matrix(data = c(0.70, 0.2, 0.1,
                       0.3, 0.4, 0.3,
                       0.2, 0.45, 0.35), byrow = byRow, nrow = 3), 
                 name = "Weather")


###################################################
### code chunk number 6: defaultMc
###################################################
defaultMc <- new("markovchain")


###################################################
### code chunk number 7: intromcList
###################################################
mcList <- new("markovchainList", markovchains = list(mcWeather, defaultMc), 
		name = "A list of Markov chains")


###################################################
### code chunk number 8: operations
###################################################
initialState <- c(0, 1, 0)
after2Days <- initialState * (mcWeather * mcWeather)
after7Days <- initialState * (mcWeather ^ 7)
after2Days
round(after7Days, 3)


###################################################
### code chunk number 9: operations2
###################################################
initialState <- c(0, 1, 0)
after2Days <- (t(mcWeather) * t(mcWeather)) * initialState
after7Days <- (t(mcWeather) ^ 7) * initialState
after2Days
round(after7Days, 3)


###################################################
### code chunk number 10: fval
###################################################
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


###################################################
### code chunk number 11: otherMethods
###################################################
states(mcWeather)
names(mcWeather)
dim(mcWeather)


###################################################
### code chunk number 12: otherMethods
###################################################
name(mcWeather)
name(mcWeather) <- "New Name"
name(mcWeather)


###################################################
### code chunk number 13: transProb
###################################################
transitionProbability(mcWeather, "cloudy", "rain")
mcWeather[2,3]


###################################################
### code chunk number 14: printAndShow
###################################################
print(mcWeather)
show(mcWeather)


###################################################
### code chunk number 15: mcPlot
###################################################
library("igraph")
plot(mcWeather,layout = layout.fruchterman.reingold,main="Weather transition matrix")


###################################################
### code chunk number 16: mcPlotdiagram
###################################################
plot(mcWeather, package="diagram", box.size = 0.04)


###################################################
### code chunk number 17: exportImport1
###################################################
mcDf <- as(mcWeather, "data.frame")
mcNew <- as(mcDf, "markovchain")
mcDf
mcIgraph <- as(mcWeather, "igraph")


###################################################
### code chunk number 18: exporImport2
###################################################
library(msm)
data(cav)
Q <- rbind ( c(0, 0.25, 0, 0.25),
             c(0.166, 0, 0.166, 0.166),
             c(0, 0.25, 0, 0.25),
             c(0, 0, 0, 0) )
cavmsm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = Q, death = 4)
msmMc <- as(cavmsm, "markovchain")
msmMc


###################################################
### code chunk number 19: exporImport3
###################################################
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
tr.prob <- etm(sir.cont, c("0", "1", "2"), tra, "cens", 1)
tr.prob
etm2mc<-as(tr.prob, "markovchain")
etm2mc


###################################################
### code chunk number 20: importExportPlot
###################################################
library(igraph)
importExportGraph<-graph.formula(dataframe++markovchain,markovchain-+igraph,markovchain++matrix,table-+markovchain,msm-+markovchain,etm-+markovchain)
plot(importExportGraph,main="Import - Export from and to markovchain objects")


###################################################
### code chunk number 21: exportImport2
###################################################
myMatr<-matrix(c(.1,.8,.1,.2,.6,.2,.3,.4,.3), byrow=TRUE, ncol=3)
myMc<-as(myMatr, "markovchain")
myMc


###################################################
### code chunk number 22: cchcMcList
###################################################
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


###################################################
### code chunk number 23: cchcMcList2
###################################################
mcCCRC[[1]]
dim(mcCCRC)


###################################################
### code chunk number 24: conditionalDistr
###################################################
conditionalDistribution(mcWeather, "sunny")


###################################################
### code chunk number 25: steadyStates
###################################################
steadyStates(mcWeather)


###################################################
### code chunk number 26: gamblerRuin
###################################################
gamblerRuinMarkovChain <- function(moneyMax, prob = 0.5) {
  require(matlab)
  matr <- zeros(moneyMax + 1)
  states <- as.character(seq(from = 0, to = moneyMax, by = 1))
  rownames(matr) = states; colnames(matr) = states
  matr[1,1] = 1; matr[moneyMax + 1,moneyMax + 1] = 1
  for(i in 2:moneyMax)
  { matr[i,i-1] = 1 - prob; matr[i, i + 1] = prob   }
  out <- new("markovchain",  
           transitionMatrix = matr, 
           name = paste("Gambler ruin", moneyMax, "dim", sep = " ")
           )
  return(out)
}

mcGR4 <- gamblerRuinMarkovChain(moneyMax = 4, prob = 0.5)
steadyStates(mcGR4)


###################################################
### code chunk number 27: absorbingStates
###################################################
absorbingStates(mcGR4)
absorbingStates(mcWeather)


###################################################
### code chunk number 28: commclassKernel
###################################################
.commclassesKernel <- function(P){
  m <- ncol(P)
	stateNames <- rownames(P)
	T <- zeros(m) 
	i <- 1
	while (i <= m) { 
		a <- i 
		b <- zeros(1,m)
		b[1,i] <- 1
		old <- 1
		new <- 0
		while (old != new) {
			old <- sum(find(b > 0))
			n <- size(a)[2]
			matr <- matrix(as.numeric(P[a,]), ncol = m, 
                     nrow = n)
			c <- colSums(matr)
			d <- find(c)
			n <- size(d)[2]
			b[1,d] <- ones(1,n)
			new <- sum(find(b>0))
			a <- d
		}
		T[i,] <- b
		i <- i+1 }
	F <- t(T)  
	C <- (T > 0)&(F > 0)
	v <- (apply(t(C) == t(T), 2, sum) == m)
	colnames(C) <- stateNames
	rownames(C) <- stateNames
	names(v) <- stateNames
	out <- list(C = C, v = v)
	return(out)
}


###################################################
### code chunk number 29: renaldoMatrix1
###################################################
P <- matlab::zeros(10)
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
.commclassesKernel(P)
summary(probMc)


###################################################
### code chunk number 30: transientStates
###################################################
transientStates(probMc)


###################################################
### code chunk number 31: probMc2Canonic
###################################################
probMcCanonic <- canonicForm(probMc)
probMc
probMcCanonic


###################################################
### code chunk number 32: isAccessible
###################################################
is.accessible(object = probMc, from = "a", to = "c")
is.accessible(object = probMc, from = "g", to = "c")


###################################################
### code chunk number 33: periodicity
###################################################

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


###################################################
### code chunk number 34: mathematica9Mc
###################################################
require(matlab)
mathematicaMatr <- zeros(5)
mathematicaMatr[1,] <- c(0, 1/3, 0, 2/3, 0)
mathematicaMatr[2,] <- c(1/2, 0, 0, 0, 1/2)
mathematicaMatr[3,] <- c(0, 0, 1/2, 1/2, 0)
mathematicaMatr[4,] <- c(0, 0, 1/2, 1/2, 0)
mathematicaMatr[5,] <- c(0, 0, 0, 0, 1)
statesNames <- letters[1:5]
mathematicaMc <- new("markovchain", transitionMatrix = mathematicaMatr,
                   name = "Mathematica MC", states = statesNames)


###################################################
### code chunk number 35: mathematica9McFig
###################################################
plot(mathematicaMc, layout = layout.fruchterman.reingold)


###################################################
### code chunk number 36: mathematica9MC
###################################################
summary(mathematicaMc)


###################################################
### code chunk number 37: fpTime1 (eval = FALSE)
###################################################
## .firstpassageKernel <- function(P, i, n){
##   G <- P
##   H <- P[i,]
##   E <- 1 - diag(size(P)[2])
##   for (m in 2:n) {
##     G <- P %*% (G * E)
##     H <- rbind(H, G[i,])
##   }
##   return(H)
## }


###################################################
### code chunk number 38: fpTime2
###################################################
firstPassagePdF <- firstPassage(object = mcWeather, state = "sunny", 
                                n = 10)
firstPassagePdF[3, 3]


###################################################
### code chunk number 39: simulatingAMarkovChain
###################################################
weathersOfDays <- rmarkovchain(n = 365, object = mcWeather, t0 = "sunny")
weathersOfDays[1:30]


###################################################
### code chunk number 40: simulatingAListOfMarkovChain
###################################################
patientStates <- rmarkovchain(n = 5, object = mcCCRC, t0 = "H", 
                              include.t0 = TRUE)
patientStates[1:10,]


###################################################
### code chunk number 41: fitMcbyMLE
###################################################
weatherFittedMLE <- markovchainFit(data = weathersOfDays, method = "mle",
                                 name = "Weather MLE")
weatherFittedMLE$estimate
weatherFittedMLE$standardError


###################################################
### code chunk number 42: fitMcbyLAPLACE
###################################################
weatherFittedLAPLACE <- markovchainFit(data = weathersOfDays, 
                                    method = "laplace", laplacian = 0.01,
                                    name = "Weather LAPLACE")
weatherFittedLAPLACE$estimate


###################################################
### code chunk number 43: fitSequenceMatrix
###################################################
createSequenceMatrix(stringchar = weathersOfDays)


###################################################
### code chunk number 44: fitSequenceMatrix2
###################################################
myMatr<-matrix(c("a","b","b","a","a","b","b","b","b","a","a","a","b","a"),ncol=2)
createSequenceMatrix(stringchar = myMatr,toRowProbs = TRUE)


###################################################
### code chunk number 45: fitMcbyBootStrap
###################################################
weatherFittedBOOT <- markovchainFit(data = weathersOfDays, 
                                    method = "bootstrap", nboot = 20)
weatherFittedBOOT$estimate
weatherFittedBOOT$standardError


###################################################
### code chunk number 46: fitMcbyBootStrap (eval = FALSE)
###################################################
## weatherFittedBOOTParallel <- markovchainFit(data = weathersOfDays, 
##                                     method = "bootstrap", nboot = 20, 
##                                     parallel = TRUE)
## weatherFittedBOOTParallel$estimate
## weatherFittedBOOTParallel$standardError


###################################################
### code chunk number 47: fitMcbyBootStrap (eval = FALSE)
###################################################
## RcppParallel::setNumThreads(2)


###################################################
### code chunk number 48: fitMcbyMLE
###################################################
weatherFittedMLE$logLikelihood
weatherFittedBOOT$logLikelihood


###################################################
### code chunk number 49: confint
###################################################
weatherFittedMLE$confidenceInterval
weatherFittedBOOT$confidenceInterval


###################################################
### code chunk number 50: multinomial
###################################################
multinomialConfidenceIntervals(transitionMatrix = 
        weatherFittedMLE$estimate@transitionMatrix, 
        countsTransitionMatrix = createSequenceMatrix(weathersOfDays))


###################################################
### code chunk number 51: fitMclists
###################################################
data(holson)
singleMc<-markovchainFit(data=holson[,2:12],name="holson")
mcListFit<-markovchainListFit(data=holson[,2:6],name="holson")
mcListFit$estimate


###################################################
### code chunk number 52: markovchainPredict
###################################################
predict(object = weatherFittedMLE$estimate, newdata = c("cloudy", "sunny"),
        n.ahead = 3)


###################################################
### code chunk number 53: markovchainListPredict
###################################################
predict(mcCCRC, newdata = c("H", "H"), n.ahead = 5)


###################################################
### code chunk number 54: markovchainListPredict2
###################################################
predict(mcCCRC, newdata = c("H", "H"), n.ahead = 5, continue = TRUE)


###################################################
### code chunk number 55: tests
###################################################
sequence<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", 
            "b", "a", "a", "b", "b", "b", "a")


###################################################
### code chunk number 56: test1
###################################################
verifyMarkovProperty(sequence)


###################################################
### code chunk number 57: test2
###################################################
data(rain)
assessOrder(rain$rain)


###################################################
### code chunk number 58: test3
###################################################
assessStationarity(rain$rain, 10)


###################################################
### code chunk number 59: test4
###################################################
mcFit<-markovchainFit(data=sequence,byrow=FALSE)
divergenceTest(sequence, mcFit$estimate@transitionMatrix)


###################################################
### code chunk number 60: rCtmcInit
###################################################
energyStates <- c("sigma", "sigma_star")
byRow <- TRUE
gen <- matrix(data = c(-3, 3,
                       1, -1), nrow = 2,
              byrow = byRow, dimnames = list(energyStates, energyStates))
molecularCTMC <- new("ctmc", states = energyStates, 
                 byrow = byRow, generator = gen, 
                 name = "Molecular Transition Model")      


###################################################
### code chunk number 61: rctmcRandom0
###################################################
statesDist <- c(0.8, 0.2)
rctmc(n = 3, ctmc = molecularCTMC, initDist = statesDist, out.type = "df", include.T0 = FALSE)


###################################################
### code chunk number 62: ctmcRandom1
###################################################
statesDist <- c(0.8, 0.2)
rctmc(n = Inf, ctmc = molecularCTMC, initDist = statesDist, T = 2)


###################################################
### code chunk number 63: rctmcSteadyStates
###################################################
steadyStates(molecularCTMC)


###################################################
### code chunk number 64: rctmcFitting
###################################################
data <- list(c("a", "b", "c", "a", "b", "a", "c", "b", "c"), 
             c(0, 0.8, 2.1, 2.4, 4, 5, 5.9, 8.2, 9))
ctmcFit(data)


###################################################
### code chunk number 65: mcWeatherQ
###################################################
mcWeatherQ <- expm::logm(mcWeather@transitionMatrix,method='Eigen')
mcWeatherQ


###################################################
### code chunk number 66: mcWeatherHalfDay
###################################################
mcWeatherHalfDayTM <- expm::expm(mcWeatherQ*.5)
mcWeatherHalfDay <- new("markovchain",transitionMatrix=mcWeatherHalfDayTM,name="Half Day Weather Transition Matrix")
mcWeatherHalfDay


###################################################
### code chunk number 67: loadAndDoExample
###################################################

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


###################################################
### code chunk number 68: MAPFit
###################################################
hyperMatrix<-matrix(c(1, 1, 2, 
                      3, 2, 1,
                      2, 2, 3), 
                    nrow = 3, byrow = TRUE,
                    dimnames = list(weatherStates,weatherStates))
markovchainFit(weathersOfDays[1:200], method = "map", 
               confidencelevel = 0.92, hyperparam = hyperMatrix)
predictiveDistribution(weathersOfDays[1:200], 
                       weathersOfDays[201:365],hyperparam = hyperMatrix)


###################################################
### code chunk number 69: MAPFit2
###################################################
hyperMatrix2<- hyperMatrix[c(2,3,1), c(2,3,1)]
markovchainFit(weathersOfDays[1:200], method = "map", 
               confidencelevel = 0.92, hyperparam = hyperMatrix2)
predictiveDistribution(weathersOfDays[1:200], 
                       weathersOfDays[201:365],hyperparam = hyperMatrix2)



###################################################
### code chunk number 70: inferHyperparam
###################################################
inferHyperparam(transMatr = weatherMatrix, scale = c(10, 10, 10))


###################################################
### code chunk number 71: inferHyperparam2
###################################################
inferHyperparam(data = weathersOfDays[1:15])


###################################################
### code chunk number 72: inferHyperparam3
###################################################
hyperMatrix3 <- inferHyperparam(transMatr = weatherMatrix, scale = c(10, 10, 10))
hyperMatrix3 <- hyperMatrix3$scaledInference
hyperMatrix4 <- inferHyperparam(data = weathersOfDays[1:15])
hyperMatrix4 <- hyperMatrix4$dataInference


###################################################
### code chunk number 73: MAPandMLE
###################################################
data(preproglucacon)
preproglucacon <- preproglucacon[[2]]
MLEest <- markovchainFit(preproglucacon, method = "mle")
MAPest <- markovchainFit(preproglucacon, method = "map")
MLEest$estimate
MAPest$estimate


###################################################
### code chunk number 74: higherOrder
###################################################
library(Rsolnp)
data(rain)
fitHigherOrder(rain$rain, 2)
fitHigherOrder(rain$rain, 3)


###################################################
### code chunk number 75: hommcObject
###################################################
showClass("hommc")


###################################################
### code chunk number 76: hommcCreate
###################################################
states <- c('a', 'b')
P <- array(dim = c(2, 2, 4), dimnames = list(states, states))
P[ , , 1] <- matrix(c(1/3, 2/3, 1, 0), byrow = FALSE, nrow = 2, ncol = 2)

P[ , , 2] <- matrix(c(0, 1, 1, 0), byrow = FALSE, nrow = 2, ncol = 2)

P[ , , 3] <- matrix(c(2/3, 1/3, 0, 1), byrow = FALSE, nrow = 2, ncol = 2)

P[ , , 4] <- matrix(c(1/2, 1/2, 1/2, 1/2), byrow = FALSE, nrow = 2, ncol = 2)

Lambda <- c(.8, .2, .3, .7)

hob <- new("hommc", order = 1, Lambda = Lambda, P = P, states = states, 
           byrow = FALSE, name = "FOMMC")
hob


###################################################
### code chunk number 77: weatPred1
###################################################

mcWP <- new("markovchain", states = c("rainy", "nice", "snowy"),
         transitionMatrix = matrix(c(0.5, 0.25, 0.25,
                                   0.5, 0, 0.5,
                                   0.25,0.25,0.5), byrow = T, nrow = 3))


###################################################
### code chunk number 78: weatPred2
###################################################
W0 <- t(as.matrix(c(0, 1, 0)))
W1 <- W0 * mcWP; W1

W2 <- W0 * (mcWP ^ 2); W2

W3 <- W0 * (mcWP ^ 3); W3


###################################################
### code chunk number 79: weatPred3
###################################################
W7 <- W0 * (mcWP ^ 7)
W7


###################################################
### code chunk number 80: weatPred4
###################################################
q <- steadyStates(mcWP)
q


###################################################
### code chunk number 81: weatPred5
###################################################
R0 <- t(as.matrix(c(1, 0, 0)))
R7 <- R0 * (mcWP ^ 7); R7

S0 <- t(as.matrix(c(0, 0, 1)))
S7 <- S0 * (mcWP ^ 7); S7


###################################################
### code chunk number 82: Alofi1
###################################################
data("rain", package = "markovchain")
table(rain$rain)


###################################################
### code chunk number 83: Alofi2
###################################################
mcAlofi <- markovchainFit(data = rain$rain, name = "Alofi MC")$estimate
mcAlofi


###################################################
### code chunk number 84: Alofi3
###################################################
steadyStates(mcAlofi)


###################################################
### code chunk number 85: ratings1
###################################################

rc <- c("AAA", "AA", "A", "BBB", "BB", "B", "CCC", "D")
creditMatrix <- matrix(c(90.81, 8.33, 0.68, 0.06, 0.08, 0.02, 0.01, 0.01,
0.70, 90.65, 7.79, 0.64, 0.06, 0.13, 0.02, 0.01,
0.09, 2.27, 91.05, 5.52, 0.74, 0.26, 0.01, 0.06,
0.02, 0.33, 5.95, 85.93, 5.30, 1.17, 1.12, 0.18,
0.03, 0.14, 0.67, 7.73, 80.53, 8.84, 1.00, 1.06,
0.01, 0.11, 0.24, 0.43, 6.48, 83.46, 4.07, 5.20,
0.21, 0, 0.22, 1.30, 2.38, 11.24, 64.86, 19.79,
0, 0, 0, 0, 0, 0, 0, 100
)/100, 8, 8, dimnames = list(rc, rc), byrow = TRUE)


###################################################
### code chunk number 86: ratings2
###################################################
creditMc <- new("markovchain", transitionMatrix = creditMatrix, 
                name = "S&P Matrix")
absorbingStates(creditMc)


###################################################
### code chunk number 87: economicAnalysis1
###################################################
statesNames <- c("customer", "non customer")
P <- zeros(2); P[1, 1] <- .9; P[1, 2] <- .1; P[2, 2] <- .95; P[2, 1] <- .05;
rownames(P) <- statesNames; colnames(P) <- statesNames
mcP <- new("markovchain", transitionMatrix = P, name = "Telephone company")
M <- zeros(2); M[1, 1] <- -20; M[1, 2] <- -30; M[2, 1] <- -40; M[2, 2] <- 0


###################################################
### code chunk number 88: economicAnalysis2
###################################################
c1 <- 100 + conditionalDistribution(mcP, state = "customer") %*% M[1,]
c2 <- 0 + conditionalDistribution(mcP, state = "non customer") %*% M[2,]


###################################################
### code chunk number 89: economicAnalysis3
###################################################
as.numeric((c(1, 0)* mcP ^ 5) %*% (as.vector(c(c1, c2))))


###################################################
### code chunk number 90: hommsales
###################################################
data(sales)
head(sales)


###################################################
### code chunk number 91: hommcFit
###################################################
# fit 8th order multivariate markov chain
object <- fitHighOrderMultivarMC(sales, order = 8, Norm = 2)


###################################################
### code chunk number 92: result
###################################################
i <- c(1, 2, 2, 3, 4, 4, 4, 5, 5, 5)
j <- c(2, 2, 2, 5, 2, 5, 5, 2, 4, 5)
k <- c(1, 1, 3, 1, 8, 1, 2, 8, 1, 2)

if(object@byrow == TRUE) {
    direction <- "(by rows)" 
} else {
    direction <- "(by cols)" 
}

cat("Order of multivariate markov chain =", object@order, "\n")
cat("states =", object@states, "\n")

cat("\n")
cat("List of Lambda's and the corresponding transition matrix", direction,":\n")

for(p in 1:10) {
  t <- 8*5*(i[p]-1) + (j[p]-1)*8
  cat("Lambda", k[p], "(", i[p], ",", j[p], ") : ", object@Lambda[t+k[p]],"\n", sep = "")
  cat("P", k[p], "(", i[p], ",", j[p], ") : \n", sep = "")
  print(object@P[, , t+k[p]])
  cat("\n")
}  
  


###################################################
### code chunk number 93: bonusMalus1
###################################################

getBonusMalusMarkovChain <- function(lambda)
{
	bmMatr <- zeros(5)
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



###################################################
### code chunk number 94: bonusMalus2
###################################################
bmMc <- getBonusMalusMarkovChain(0.05)
as.numeric(steadyStates(bmMc))


###################################################
### code chunk number 95: bonusMalus3
###################################################
sum(as.numeric(steadyStates(bmMc)) * c(0.5, 0.7, 0.9, 1, 1.25))


###################################################
### code chunk number 96: healthIns1
###################################################

mcHI <- new("markovchain", states = c("active", "disable", "withdrawn", 
                                      "death"),
         transitionMatrix = matrix(c(0.5, .25, .15, .1,
                                   0.4, 0.4, 0.0, 0.2,
                                   0, 0, 1, 0,
                                   0, 0, 0, 1), byrow = TRUE, nrow = 4))
benefitVector <- as.matrix(c(0, 0, 500, 1000))


###################################################
### code chunk number 97: healthIns2
###################################################
T0 <- t(as.matrix(c(1, 0, 0, 0)))
T1 <- T0 * mcHI
T2 <- T1 * mcHI
T3 <- T2 * mcHI


###################################################
### code chunk number 98: healthIns3
###################################################
PVFB <- T0 %*% benefitVector * 1.05 ^ -0 + 
  T1 %*% benefitVector * 1.05 ^ -1+
  T2 %*% benefitVector * 1.05 ^ -2 + T3 %*% benefitVector * 1.05 ^ -3


###################################################
### code chunk number 99: healthIns4
###################################################
P <- PVFB / (T0[1] * 1.05 ^- 0 + T1[1] * 1.05 ^ -1 + T2[1] * 1.05 ^ -2)


###################################################
### code chunk number 100: healthIns5
###################################################
PVFB <- T2 %*% benefitVector * 1.05 ^ -1 + T3 %*% benefitVector * 1.05 ^ -2
PVFP <- P*(T1[1] * 1.05 ^ -0 + T2[1] * 1.05 ^ -1)
V <- PVFB - PVFP
V


###################################################
### code chunk number 101: blandenEtAlii
###################################################
data("blanden")
mobilityMc <- as(blanden, "markovchain")
mobilityMc


###################################################
### code chunk number 102: blandenEtAlii2
###################################################
plot(mobilityMc, main = '1970 mobility',vertex.label.cex = 2,
		layout = layout.fruchterman.reingold)


###################################################
### code chunk number 103: blandenEtAlii3
###################################################
round(steadyStates(mobilityMc), 2)


###################################################
### code chunk number 104: preproglucacon1
###################################################
data("preproglucacon", package = "markovchain")


###################################################
### code chunk number 105: preproglucacon2
###################################################
mcProtein <- markovchainFit(preproglucacon$preproglucacon, 
                          name = "Preproglucacon MC")$estimate
mcProtein


###################################################
### code chunk number 106: epid1
###################################################
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


###################################################
### code chunk number 107: epid2
###################################################
eig <- eigen(mcM6@transitionMatrix)
D <- diag(eig$values)


###################################################
### code chunk number 108: epid3
###################################################
V <- eig$vectors 
V %*% D %*% solve(V)
d <- D ^ (1/6)
M <- V %*% d %*% solve(V)
mcM1 <- new("markovchain", transitionMatrix = M, states = hivStates)


