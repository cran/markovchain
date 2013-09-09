\name{markovchain-package}
\alias{markovchain-package}
\alias{markovchain}
\docType{package}
\title{
A package for easily handling discrete Markov chains
}
\description{
It contains S4 classes and methods to create and operates with Markov chains}
\details{
\tabular{ll}{
Package: \tab markovchain\cr
Type: \tab Package\cr
Version: \tab 0.0.2\cr
Date: \tab 2013-08-31\cr
License: \tab GPL-2\cr
Depends: \tab R (>= 2.14), methods, expm, matlab, igraph, Matrix\cr
}
}
\author{
Giorgio Alfredo Spedicato
Maintainer: Giorgio Alfredo Spedicato <spedicato_giorgio@yahoo.it>

}
\references{
Introduction to Probability Models, Sheldon Ross
}
\keyword{ package }

\examples{
#create some markov chains
statesNames=c("a","b")
 mcA<-new("markovchain", transitionMatrix=matrix(c(0.7,0.3,0.1,0.9),byrow=TRUE,
 nrow=2, dimnames=list(statesNames,statesNames)))

statesNames=c("a","b","c")
mcB<-new("markovchain", states=statesNames, transitionMatrix=
          matrix(c(0.2,0.5,0.3,
                   0,1,0,
                   0.1,0.8,0.1),nrow=3, byrow=TRUE, dimnames=list(statesNames,
				   statesNames)
                 ))

statesNames=c("a","b","c","d")
matrice<-matrix(c(0.25,0.75,0,0,0.4,0.6,0,0,0,0,0.1,0.9,0,0,0.7,0.3), 
nrow=4, byrow=TRUE)
mcC<-new("markovchain", states=statesNames, transitionMatrix=matrice)
mcD<-new("markovchain", transitionMatrix=matrix(c(0,1,0,1), nrow=2,byrow=TRUE))



#operations with S4 methods

mcA^2
steadyStates(mcB)
absorbingStates(mcB)
markovchainSequence(n=20, markovchain=mcC, include=TRUE)

}
