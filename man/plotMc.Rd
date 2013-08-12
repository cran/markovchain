\name{plotMc}
\alias{plotMc}
\title{
Function to plot a discrete Markov chain
}
\description{
This function returns absorbing states in matricial form.
}
\usage{
plotMc(object,...)
}
\arguments{
  \item{object}{A discrete Markov chain.}
  \item{...}{Additional parameters for plot method.}
}
\value{A matrix}
\references{
	Sheldon Ross, Probability Models
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
plotMc(markovB)
}

