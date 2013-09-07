\name{markovchainList-class}
\Rdversion{1.1}
\docType{class}
\alias{markovchainList-class}

\title{Class \code{"markovchainList"}}
\description{
A class to handle non - homogeneous markov chains
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("markovchainList", ...)}.
Each item in the list is a \code{markovchain} object.
}
\section{Slots}{
  \describe{
    \item{\code{markovchains}:}{Object of class \code{"list"}: a list of markovchains}
    \item{\code{name}:}{Object of class \code{"character"}: optional name of the class }
  }
}
\section{Methods}{
No methods defined with class "markovchainList" in the signature.
}
\references{
Sheldon Ross
}
\author{
Giorgio Spedicato
}
\note{
It can help to deal with non - homogeneous markov chains.
}
\seealso{
\code{\linkS4class{markovchain}}
}
\examples{
showClass("markovchainList")
}
\keyword{classes}
