## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(fig.width=8.5, fig.height=6, out.width = "70%")
set.seed(123)

## ----load, results='hide', warning=FALSE, message=FALSE-----------------------
require(markovchain)


## ----higherOrder--------------------------------------------------------------
if (requireNamespace("Rsolnp", quietly = TRUE)) {
library(Rsolnp)
data(rain)
fitHigherOrder(rain$rain, 2)
fitHigherOrder(rain$rain, 3)
}

## ----hommcObject--------------------------------------------------------------
showClass("hommc")

## ----hommcCreate--------------------------------------------------------------
states <- c('a', 'b')
P <- array(dim = c(2, 2, 4), dimnames = list(states, states))
P[ , , 1] <- matrix(c(1/3, 2/3, 1, 0), byrow = FALSE, nrow = 2, ncol = 2)

P[ , , 2] <- matrix(c(0, 1, 1, 0), byrow = FALSE, nrow = 2, ncol = 2)

P[ , , 3] <- matrix(c(2/3, 1/3, 0, 1), byrow = FALSE, nrow = 2, ncol = 2)

P[ , , 4] <- matrix(c(1/2, 1/2, 1/2, 1/2), byrow = FALSE, nrow = 2, ncol = 2)

Lambda <- c(.8, .2, .3, .7)

hob <- new("hommc", order = 1, Lambda = Lambda, P = P, states = states, 
           byrow = FALSE, name = "FOMMC")
hob

## ----hommsales----------------------------------------------------------------
data(sales)
head(sales)

## ----hommcFit, warning = FALSE, message = FALSE-------------------------------
# fit 8th order multivariate markov chain
if (requireNamespace("Rsolnp", quietly = TRUE)) {
object <- fitHighOrderMultivarMC(sales, order = 8, Norm = 2)
}

## ----result, echo = FALSE-----------------------------------------------------
if (requireNamespace("Rsolnp", quietly = TRUE)) {
i <- c(1, 2, 2, 3, 4, 4, 4, 5, 5, 5)
j <- c(2, 2, 2, 5, 2, 5, 5, 2, 4, 5)
k <- c(1, 1, 3, 1, 8, 1, 2, 8, 1, 2)

if(object@byrow == TRUE) {
    direction <- "(by rows)" 
} else {
    direction <- "(by cols)" 
}

cat("Order of multivariate markov chain =", object@order, "\n")
cat("states =", object@states, "\n")

cat("\n")
cat("List of Lambda's and the corresponding transition matrix", direction,":\n")

for(p in 1:10) {
  t <- 8*5*(i[p]-1) + (j[p]-1)*8
  cat("Lambda", k[p], "(", i[p], ",", j[p], ") : ", object@Lambda[t+k[p]],"\n", sep = "")
  cat("P", k[p], "(", i[p], ",", j[p], ") : \n", sep = "")
  print(object@P[, , t+k[p]])
  cat("\n")
}
} else {
  print("package Rsolnp unavailable")
}

