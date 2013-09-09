\name{markovchain-class}
\Rdversion{1.1}
\docType{class}
\alias{markovchain-class}
\alias{*,markovchain,markovchain-method}
\alias{*,markovchain,matrix-method}
\alias{*,markovchain,numeric-method}
\alias{*,matrix,markovchain-method}
\alias{*,numeric,markovchain-method}
\alias{^,markovchain,numeric-method}
\alias{==,markovchain,markovchain-method}
\alias{absorbingStates,markovchain-method}
\alias{coerce,data.frame,markovchain-method}
\alias{coerce,markovchain,data.frame-method}
\alias{dim,markovchain-method}
\alias{initialize,markovchain-method}
\alias{plotMc,markovchain-method}
\alias{print,markovchain-method}
\alias{show,markovchain-method}
\alias{states,markovchain-method}
\alias{steadyStates,markovchain-method}
\alias{t,markovchain-method}
\alias{transitionProbability,markovchain-method}
\title{Class \code{"markovchain"}}
\description{
The S4 class that describes \code{markovchain} objects.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("markovchain", states, byrow, transitionMatrix, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{states}:}{Name of the states. Must be the same of \code{colnames} and \code{rownames} of the transition matrix}
    \item{\code{byrow}:}{Binary flag.}
    \item{\code{transitionMatrix}:}{Square transition matrix}
    \item{\code{name}:}{Optional character name of the Markov chain}
  }
}
\section{Methods}{
  \describe{
    \item{*}{\code{signature(e1 = "markovchain", e2 = "markovchain")}: multiply two \code{markovchain} objects}
    \item{*}{\code{signature(e1 = "markovchain", e2 = "matrix")}: markovchain by matrix multiplication}
    \item{*}{\code{signature(e1 = "markovchain", e2 = "numeric")}: markovchain by numeric vector multiplication }
    \item{*}{\code{signature(e1 = "matrix", e2 = "markovchain")}: matrix by markov chain}
    \item{*}{\code{signature(e1 = "numeric", e2 = "markovchain")}: numeric vector by \code{markovchain} multiplication   }
    \item{^}{\code{signature(e1 = "markovchain", e2 = "numeric")}: power of a \code{markovchain} object}
    \item{==}{\code{signature(e1 = "markovchain", e2 = "markovchain")}: equality of two \code{markovchain} object}
    \item{absorbingStates}{\code{signature(object = "markovchain")}: method to get absorbing states }
    \item{coerce}{\code{signature(from = "markovchain", to = "data.frame")}: coerce method from markovchain to \code{data.frame}}
    \item{coerce}{\code{signature(from = "data.frame", to = "markovchain")}: coerce method from \code{data.frame} to \code{markovchain}}
    \item{dim}{\code{signature(x = "markovchain")}: method to get the size}
    \item{initialize}{\code{signature(.Object = "markovchain")}: initialize method }
    \item{plotMc}{\code{signature(object = "markovchain")}: plot method }
    \item{print}{\code{signature(x = "markovchain")}: print method. }
    \item{show}{\code{signature(object = "markovchain")}: show method. }
    \item{states}{\code{signature(object = "markovchain")}: states method. }
    \item{steadyStates}{\code{signature(object = "markovchain")}: method to get the steady vector. }
    \item{t}{\code{signature(x = "markovchain")}: transpose matrix }
    \item{transitionProbability}{\code{signature(object = "markovchain")}: transition probability }
	 }
}
\references{
Introduction to Probability Models, Sheldon Ross
}
\author{
Giorgio Spedicato
}
\note{
It is an S4 class.
}
\section{Warning }{Validation method is used to assess whether either columns or rows totals to one. Rounding is used up to 5th decimal}

\seealso{
\code{\link{markovchainSequence}},\code{\link{markovchainFit}}
}
\examples{
#show markov chain definition
showClass("markovchain")
#create a simple markov chains
simpleMc<-new("markovchain", states=c("a","b"),
transitionMatrix=matrix(c(0.4,0.6,.3,.7),nrow=2,byrow=TRUE))
#power
simpleMc^4
#some methods
steadyStates(simpleMc)
}
\keyword{classes}
