% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{blanden}
\alias{blanden}
\title{Mobility between income quartiles}
\format{
An object of class \code{table} with 4 rows and 4 columns.
}
\source{
Personal reworking
}
\usage{
data(blanden)
}
\description{
This table show mobility between income quartiles for father and sons for the 1970 cohort born
}
\details{
The rows represent fathers' income quartile when the son is aged 16, whilst the columns represent sons' income quartiles when he is aged 30 (in 2000).
}
\examples{
data(blanden)
mobilityMc<-as(blanden, "markovchain")
}
\references{
Jo Blanden, Paul Gregg and Stephen Machin, Intergenerational Mobility in Europe and North America, Center for Economic Performances (2005)
}
\keyword{datasets}
