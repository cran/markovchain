% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/probabilistic.R
\name{firstPassageMultiple}
\alias{firstPassageMultiple}
\title{function to calculate first passage probabilities}
\usage{
firstPassageMultiple(object, state, set, n)
}
\arguments{
\item{object}{a markovchain-class object}

\item{state}{intital state of the process (charactervector)}

\item{set}{set of states A, first passage of which is to be calculated}

\item{n}{Number of rows on which compute the distribution}
}
\value{
A vector of size n showing the first time proabilities
}
\description{
The function calculates first passage probability for a subset of
states given an initial state.
}
\examples{
statesNames <- c("a", "b", "c")
markovB <- new("markovchain", states = statesNames, transitionMatrix =
matrix(c(0.2, 0.5, 0.3,
         0, 1, 0,
         0.1, 0.8, 0.1), nrow = 3, byrow = TRUE,
       dimnames = list(statesNames, statesNames)
))
firstPassageMultiple(markovB,"a",c("b","c"),4)  

}
\references{
Renaldo Feres, Notes for Math 450 Matlab listings for Markov chains;
MIT OCW, course - 6.262, Discrete Stochastic Processes, course-notes, chap -05
}
\seealso{
\code{\link{firstPassage}}
}
\author{
Vandit Jain
}
