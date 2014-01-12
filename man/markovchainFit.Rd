\name{markovchainFit}
\alias{markovchainFit}
\alias{createSequenceMatrix}
\title{
Function to fit a discrete Markov chain
}
\description{
Given a sequence of states arising from a stationary state, 
it fits the underlying Markov chain distribution using either MLE (also using a Laplacian smoother) or bootstrap.
}
\usage{
markovchainFit(data, method = "mle", byrow = TRUE, nboot = 10,laplacian=0,name)
createSequenceMatrix(stringchar, toRowProbs = FALSE, sanitize = TRUE)
}

\arguments{
  \item{data}{
A character list.
}
\item{method}{
Method used to estimate the markov chain. Either "mle" or "bootstrap" or "laplace"
}
  \item{byrow}{
it tells whether the output markov chain should show the transition probabilities by row.
}
  \item{nboot}{
Number of bootstrap replicates in case "bootstrap" is used.
}
  \item{laplacian}{
Laplacian smoothing parameter, default zero. It is only used when "laplace" method is chosen.
}
  \item{name}{
Optional character for name slot.
}
  \item{stringchar}{
Equivalent to data
}
  \item{toRowProbs}{
converts a sequence matrix into a probability matrix
}
  \item{sanitize}{
put 1 in all rows having rowSum equal to zero
}
}
\value{
A list containing an estimate and, when "bootstrap" method is used, a matrix 
of standards deviations and the bootstrap samples.
}
\references{
Sheldon Ross, Probability Models.
}
\author{
Giorgio Spedicato
}
\note{
When MLE method is called, the lists contains one entry: estimate. 
Bootstrap method has been defined "euristically".
}

\section{Warning }{"mle" method calls \code{createSequenceMatrix} function using \code{sanitize} parameter set to TRUE.}

\seealso{
\code{\link{markovchainSequence}}
}
\examples{
sequence<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", 
"b", "b", "b", "a")
sequenceMatr<-createSequenceMatrix(sequence,sanitize=FALSE)
mcFitMLE<-markovchainFit(data=sequence)
mcFitBSP<-markovchainFit(data=sequence,method="bootstrap",nboot=5, name="Bootstrap Mc")
}

