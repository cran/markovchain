\name{markovchainFit}
\alias{markovchainFit}
\title{
Function to fit a discrete Markov chain
}
\description{
Given a sequence of states arising from a stationary state, it fits the underlying markov chain distribution using a frequentist approach.
}
\usage{
markovchainFit(data, byrow = TRUE)
}

\arguments{
  \item{data}{
A character list.
}
  \item{byrow}{
it tells whether the output markov chain should show the transition probabilities by row.
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
mcFit<-markovchainFit(data=sequence)
}

