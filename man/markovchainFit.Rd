\name{markovchainFit}
\alias{markovchainFit}
\title{
Function to fit a discrete Markov chain
}
\description{
Given a sequence of states arising from a stationary state, it fits the underlying markov chain distribution using a frequentist approach.
}
\usage{
markovchainFit(data, method = "mle", byrow = TRUE, nboot = 10)
}

\arguments{
  \item{data}{
A character list.
}
\item{method}{
Method used to estimate the markov chain. Either "mle" or "bootstrap"
}
  \item{byrow}{
it tells whether the output markov chain should show the transition probabilities by row.
}
  \item{nboot}{
Number of bootstrap replicates in case "bootstrap" is used.
}
}
\value{
A \code{markovchain} object}
\references{
Sheldon Ross, Probability Models.
}
\author{
Giorgio Spedicato
}
\note{
The lists contains one entry: estimate. It is expected standard errors to 
be added in the future.
}

\section{Warning }{This function will probably be modified in the future}

\seealso{
\code{\link{markovchainSequence}}
}
\examples{
sequence<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", 
"b", "b", "b", "a")
mcFitMLE<-markovchainFit(data=sequence)
mcFitBSP<-markovchainFit(data=sequence,method="bootstrap",nboot=5)
}

