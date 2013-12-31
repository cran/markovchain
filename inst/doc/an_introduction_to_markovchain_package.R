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
weatherStates<-c("sunny", "cloudy", "rain")
byRow<-TRUE
weatherMatrix<-matrix(data=c(0.70, 0.2,0.1,
                       0.3,0.4, 0.3,
                       0.2,0.45,0.35),byrow=byRow, nrow=3,
                     dimnames=list(weatherStates, weatherStates))
mcWeather<-new("markovchain",states=weatherStates, byrow=byRow, 
               transitionMatrix=weatherMatrix, name="Weather")


###################################################
### code chunk number 5: mcInitLong
###################################################
mcWeather<-new("markovchain", states=c("sunny", "cloudy", "rain"), transitionMatrix=matrix(data=c(0.70, 0.2,0.1,
                       0.3,0.4, 0.3,
                       0.2,0.45,0.35),byrow=byRow, nrow=3), name="Weather")



###################################################
### code chunk number 6: defaultMc
###################################################
defaultMc<-new("markovchain")


###################################################
### code chunk number 7: intromcList
###################################################
mcList<-new("markovchainList",markovchains=list(mcWeather, defaultMc), name="A list of Markov chains")


###################################################
### code chunk number 8: operations
###################################################
initialState<-c(0,1,0)
after2Days<-initialState*(mcWeather*mcWeather)
after7Days<-initialState*(mcWeather^7)
after2Days
after7Days


###################################################
### code chunk number 9: operations2
###################################################
initialState<-c(0,1,0)
mcWeatherTransposed<-t(mcWeather)
after2Days<-(mcWeatherTransposed*mcWeatherTransposed)*initialState
after7Days<-(mcWeather^7)*initialState
after2Days
after7Days


###################################################
### code chunk number 10: otherMethods
###################################################
states(mcWeather)
dim(mcWeather)


###################################################
### code chunk number 11: transProb
###################################################
transitionProbability(mcWeather, "cloudy","rain")
mcWeather[2,3]


###################################################
### code chunk number 12: printAndShow
###################################################
print(mcWeather)
show(mcWeather)


###################################################
### code chunk number 13: mcPlot
###################################################
plotMc(mcWeather)


###################################################
### code chunk number 14: exportImport
###################################################
mcDf<-as(mcWeather, "data.frame")
mcNew<-as(mcDf, "markovchain")


###################################################
### code chunk number 15: cchcMcList
###################################################
stateNames=c("H","I","D")
Q0<-new("markovchain", states=stateNames, 
        transitionMatrix=matrix(c(0.7, 0.2, 0.1,0.1, 0.6, 0.3,0, 0, 1),byrow=TRUE, nrow=3), name="state t0")
Q1<-new("markovchain", states=stateNames, 
        transitionMatrix=matrix(c(0.5, 0.3, 0.2,0, 0.4, 0.6,0, 0, 1),byrow=TRUE, nrow=3), name="state t1")
Q2<-new("markovchain", states=stateNames, 
        transitionMatrix=matrix(c(0.3, 0.2, 0.5,0, 0.2, 0.8,0, 0, 1),byrow=TRUE,nrow=3), name="state t2")
Q3<-new("markovchain", states=stateNames, transitionMatrix=matrix(c(0, 0, 1,0, 0, 1,0, 0, 1),byrow=TRUE, nrow=3), name="state t3")
mcCCRC<-new("markovchainList",markovchains=list(Q0,Q1,Q2,Q3), name="Continuous Care Health Community")
print(mcCCRC)


###################################################
### code chunk number 16: cchcMcList2
###################################################
mcCCRC[[1]]
dim(mcCCRC)


###################################################
### code chunk number 17: conditionalDistr
###################################################
conditionalDistribution(mcWeather, "sunny")


###################################################
### code chunk number 18: steadyStates
###################################################
steadyStates(mcWeather)


###################################################
### code chunk number 19: gamblerRuin
###################################################
gamblerRuinMarkovChain<-function(moneyMax, prob=0.5) {
  require(matlab)
  matr<-zeros(moneyMax+1)
  states<-as.character(seq(from=0, to=moneyMax, by=1))
  rownames(matr)=states; colnames(matr)=states
  matr[1,1]=1;matr[moneyMax+1,moneyMax+1]=1
  for(i in 2:moneyMax)
  {
    matr[i,i-1]=1-prob;matr[i,i+1]=prob
  }
  out<-new("markovchain",  
           transitionMatrix=matr, 
           name=paste("Gambler ruin",moneyMax,"dim",sep=" ")
           )
  return(out)
}

mcGR4<-gamblerRuinMarkovChain(moneyMax=4, prob=0.5)
steadyStates(mcGR4)


###################################################
### code chunk number 20: absorbingStates
###################################################
absorbingStates(mcGR4)
absorbingStates(mcWeather)


###################################################
### code chunk number 21: simulatingAMarkovChain
###################################################
weathersOfDays<-rmarkovchain(n=365,object=mcWeather,t0="sunny")
weathersOfDays[1:30]


###################################################
### code chunk number 22: simulatingAListOfMarkovChain
###################################################
patientStates<-rmarkovchain(n=5, object=mcCCRC,t0="H",include.t0=TRUE)
patientStates[1:10,]


###################################################
### code chunk number 23: fitMcbyMLE
###################################################
weatherFittedMLE<-markovchainFit(data=weathersOfDays, method="mle",name="Weather MLE")
weatherFittedMLE$estimate


###################################################
### code chunk number 24: fitMcbyLAPLACE
###################################################
weatherFittedLAPLACE<-markovchainFit(data=weathersOfDays, method="laplace",laplacian=0.01,name="Weather LAPLACE")
weatherFittedLAPLACE$estimate


###################################################
### code chunk number 25: fitSequenceMatrix
###################################################
createSequenceMatrix(stringchar = weathersOfDays)


###################################################
### code chunk number 26: fitMcbyBootStrap
###################################################
weatherFittedBOOT<-markovchainFit(data=weathersOfDays, method="bootstrap",nboot=100)
weatherFittedBOOT$estimate
weatherFittedBOOT$standardError


###################################################
### code chunk number 27: markovchainPredict
###################################################
predict(object=weatherFittedMLE$estimate,newdata=c("cloudy","sunny"),n.ahead=3)


###################################################
### code chunk number 28: markovchainListPredict
###################################################
predict(mcCCRC,newdata=c("H","H"),n.ahead=5)


###################################################
### code chunk number 29: markovchainListPredict2
###################################################
predict(mcCCRC,newdata=c("H","H"),n.ahead=5, continue=TRUE)


###################################################
### code chunk number 30: healthIns1
###################################################

mcHI=new("markovchain", states=c("active", "disable", "withdrawn", "death"),
         transitionMatrix=matrix(c(0.5,.25,.15,.1,
                                   0.4,0.4,0.0,.2,
                                   0,0,1,0,
                                   0,0,0,1), byrow=TRUE, nrow=4))
         

benefitVector=as.matrix(c(0,0,500,1000))



###################################################
### code chunk number 31: healthIns2
###################################################
T0=t(as.matrix(c(1,0,0,0)))
T1=T0*mcHI
T2=T1*mcHI
T3=T2*mcHI


###################################################
### code chunk number 32: healthIns3
###################################################
PVFB=T0%*%benefitVector*1.05^-0+T1%*%benefitVector*1.05^-1+
		T2%*%benefitVector*1.05^-2+T3%*%benefitVector*1.05^-3


###################################################
### code chunk number 33: healthIns4
###################################################
P=PVFB/(T0[1]*1.05^-0+T1[1]*1.05^-1+T2[1]*1.05^-2)


###################################################
### code chunk number 34: healthIns5
###################################################
PVFB=(T2%*%benefitVector*1.05^-1+T3%*%benefitVector*1.05^-2)
PVFP=P*(T1[1]*1.05^-0+T2[1]*1.05^-1)
V=PVFB-PVFP
V


###################################################
### code chunk number 35: weatPred1
###################################################

mcWP=new("markovchain", states=c("rainy", "nice", "snowy"),
         transitionMatrix=matrix(c(0.5, 0.25, 0.25,
                                   0.5, 0, 0.5,
                                   0.25,0.25,0.5), byrow=TRUE, nrow=3))


###################################################
### code chunk number 36: weatPred2
###################################################
W0=t(as.matrix(c(0,1,0)))
W1=W0*mcWP
W1

W2=W0*(mcWP^2)
W2

W3=W0*(mcWP^3)
W3


###################################################
### code chunk number 37: weatPred3
###################################################
W7=W0*(mcWP^7)
W7


###################################################
### code chunk number 38: weatPred4
###################################################
q=steadyStates(mcWP)
q


###################################################
### code chunk number 39: weatPred5
###################################################
R0=t(as.matrix(c(1,0,0)))
R7=W0*(mcWP^7)
R7

S0=t(as.matrix(c(0,0,1)))
R7=W0*(mcWP^7)
R7


###################################################
### code chunk number 40: Alofi1
###################################################
data(rain, package="markovchain")
table(rain$rain)


###################################################
### code chunk number 41: Alofi2
###################################################
mcAlofi<-markovchainFit(data=rain$rain, name="Alofi MC")$estimate
mcAlofi


###################################################
### code chunk number 42: Alofi3
###################################################
steadyStates(mcAlofi)


###################################################
### code chunk number 43: preproglucacon1
###################################################
data(preproglucacon, package="markovchain")


###################################################
### code chunk number 44: preproglucacon2
###################################################
mcProtein<-markovchainFit(preproglucacon$preproglucacon, name="Preproglucacon MC")$estimate


###################################################
### code chunk number 45: epid1
###################################################
craigSendiMatr<-matrix(c(682,33,25,
              154,64,47,
              19,19,43), byrow=T,nrow=3)
hivStates<-c("0-49", "50-74", "75-UP")
rownames(craigSendiMatr)<-hivStates
colnames(craigSendiMatr)<-hivStates
craigSendiTable<-as.table(craigSendiMatr)
mcM6<-as(craigSendiTable,"markovchain")
mcM6@name="Zero-Six month CD4 cells transition"
mcM6


###################################################
### code chunk number 46: epid2
###################################################
autov=eigen(mcM6@transitionMatrix)
D=diag(autov$values)


###################################################
### code chunk number 47: epid3
###################################################
P=autov$vectors 
P%*%D%*%solve(P)
d=D^(1/6)
M=P%*%d%*%solve(P)
mcM1<-new("markovchain",transitionMatrix=M,states=hivStates)


