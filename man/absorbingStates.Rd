\name{absorbingStates}
\alias{absorbingStates}
\alias{transientStates}
\alias{canonicForm}
\title{
Function to get the absorbing, transient states of a Markov chain and to put the matrix in canonic form
}
\description{
These functions returns absorbing and transient states of the \code{markovchain} object.
}
\usage{
absorbingStates(object)
transientStates(object)
canonicForm(object)
}
\arguments{
  \item{object}{A \code{markovchain} object.}
}
\value{A matrix}
\references{
Feres, Matlab listing for markov chain.
}
\author{
Giorgio Spedicato
}
\seealso{
\code{\linkS4class{markovchain}}
}
\examples{
statesNames=c("a","b","c")
markovB<-new("markovchain", states=statesNames, transitionMatrix=
          matrix(c(0.2,0.5,0.3,
                   0,1,0,
                   0.1,0.8,0.1),nrow=3, byrow=TRUE, dimnames=list(statesNames,statesNames)
                 ))
absorbingStates(markovB)
transientStates(markovB)
canonicForm(markovB)

}

