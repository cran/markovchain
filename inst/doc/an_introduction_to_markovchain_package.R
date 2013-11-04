### R code from vignette source 'an_introduction_to_markovchain_package.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: setup
###################################################
	options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
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
### code chunk number 8: showClassesAndMethods
###################################################
showMethods(class="markovchain")
showMethods(class="markovchainList")


###################################################
### code chunk number 9: operations
###################################################
initialState<-c(0,1,0)
after2Days<-initialState*(mcWeather*mcWeather)
after7Days<-initialState*(mcWeather^7)
after2Days
after7Days


###################################################
### code chunk number 10: operations2
###################################################
initialState<-c(0,1,0)
mcWeatherTransposed<-t(mcWeather)
after2Days<-(mcWeatherTransposed*mcWeatherTransposed)*initialState
after7Days<-(mcWeather^7)*initialState
after2Days
after7Days


###################################################
### code chunk number 11: otherMethods
###################################################
states(mcWeather)
dim(mcWeather)


###################################################
### code chunk number 12: transProb
###################################################
transitionProbability(mcWeather, "cloudy","rain")
mcWeather[2,3]


###################################################
### code chunk number 13: printAndShow
###################################################
print(mcWeather)
show(mcWeather)


###################################################
### code chunk number 14: mcPlot
###################################################
plotMc(mcWeather)


###################################################
### code chunk number 15: exportImport
###################################################
mcDf<-as(mcWeather, "data.frame")
mcNew<-as(mcDf, "markovchain")


###################################################
### code chunk number 16: cchcMcList
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


###################################################
### code chunk number 17: cchcMcList2
###################################################
mcCCRC[[1]]
dim(mcCCRC)


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
weatherFittedMLE<-markovchainFit(data=weathersOfDays, method="mle")
weatherFittedMLE$estimate


###################################################
### code chunk number 24: fitMcbyBootStrap
###################################################
weatherFittedBOOT<-markovchainFit(data=weathersOfDays, method="bootstrap",nboot=50)
weatherFittedBOOT$estimate
weatherFittedBOOT$standardError


###################################################
### code chunk number 25: healthIns1
###################################################

mcHI=new("markovchain", states=c("active", "disable", "withdrawn", "death"),
         transitionMatrix=matrix(c(0.5,.25,.15,.1,
                                   0.4,0.4,0.0,.2,
                                   0,0,1,0,
                                   0,0,0,1), byrow=TRUE, nrow=4))
         

benefitVector=as.matrix(c(0,0,500,1000))



###################################################
### code chunk number 26: healthIns2
###################################################
T0=t(as.matrix(c(1,0,0,0)))
T1=T0*mcHI
T2=T1*mcHI
T3=T2*mcHI


###################################################
### code chunk number 27: healthIns3
###################################################
PVFB=T0%*%benefitVector*1.05^-0+T1%*%benefitVector*1.05^-1+T2%*%benefitVector*1.05^-2+T3%*%benefitVector*1.05^-3


###################################################
### code chunk number 28: healthIns4
###################################################
P=PVFB/(T0[1]*1.05^-0+T1[1]*1.05^-1+T2[1]*1.05^-2)


###################################################
### code chunk number 29: healthIns5
###################################################
PVFB=(T2%*%benefitVector*1.05^-1+T3%*%benefitVector*1.05^-2)
PVFP=P*(T1[1]*1.05^-0+T2[1]*1.05^-1)
V=PVFB-PVFP
V


