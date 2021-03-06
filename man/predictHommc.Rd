% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hommc.R
\name{predictHommc}
\alias{predictHommc}
\title{Simulate a higher order multivariate markovchain}
\usage{
predictHommc(hommc,t,init)
}
\arguments{
\item{hommc}{a hommc-class object}

\item{t}{no of iterations to predict}

\item{init}{matrix of previous states size of which depends on hommc}
}
\value{
The function returns a matrix of size s X t displaying t predicted states 
in each row coressponding to every categorical sequence.
}
\description{
This function provides a prediction of states for a higher order 
multivariate markovchain object
}
\details{
The user is required to provide a matrix of giving n previous coressponding
every categorical sequence. Dimensions of the init are s X n, where s is 
number of categorical sequences and n is order of the homc.
}
\author{
Vandit Jain
}
