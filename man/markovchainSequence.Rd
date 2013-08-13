\name{markovchainSequence}
\alias{markovchainSequence}
\title{
Function to generate a sequence of states from a given discrete markov chain.
}
\description{
Provided a \code{markovchain} object, it returns a sequence of states coming from the underlying 
stationary distribution. 
}
\usage{
markovchainSequence(n, markovchain, t0 = sample(markovchain@states, 1), 
include.t0 = FALSE)
}
\arguments{
  \item{n}{
length of the sequence
}
  \item{markovchain}{
markov chain
}
  \item{t0}{
initial state
}
  \item{include.t0}{
logical, tells if t0 shall be added as initial value.
}
}
\details{
It is assumed the process to arise from a single steady state markov chain.
}
\value{a vector of characters}
\author{
Giorgio Spedicato
}
\note{
It is assumed the process to arise from a single steady state markov chain.
}


\seealso{
\code{\link{markovchainFit}}
}
\examples{
#define the markov chain
statesNames=c("a","b","c")
mcB<-new("markovchain", states=statesNames, transitionMatrix=
          matrix(c(0.2,0.5,0.3,
                   0,0.2,0.8,
                   0.1,0.8,0.1),nrow=3, byrow=TRUE, dimnames=list(statesNames,
				   statesNames)
                 ))

#show the sequence
outs<-markovchainSequence(n=100,markovchain=mcB, t0="a")
}

