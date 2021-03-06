% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{preproglucacon}
\alias{preproglucacon}
\title{Preprogluccacon DNA protein bases sequences}
\format{
A data frame with 1572 observations on the following 2 variables.
\describe{
  \item{\code{V1}}{a numeric vector, showing original coding}
  \item{\code{preproglucacon}}{a character vector, showing initial of DNA bases (Adenine, Cytosine, Guanine, Thymine)}
}
}
\source{
Avery Henderson
}
\usage{
data(preproglucacon)
}
\description{
Sequence of bases for preproglucacon DNA protein
}
\examples{
data(preproglucacon)
preproglucaconMc<-markovchainFit(data=preproglucacon$preproglucacon)
}
\references{
Averuy Henderson, Fitting markov chain models on discrete time series such as DNA sequences
}
\keyword{datasets}
