## -----------------------------------------------------------------------------
#| label: global_options
#| include: false
knitr::opts_chunk$set(fig.width=8.5, fig.height=6, out.width = "70%")
set.seed(123)
# ANSI colour codes (U+001B) break the LaTeX build, so keep the output plain
options(crayon.enabled = FALSE, cli.num_colors = 1, cli.hyperlink = FALSE)
Sys.setenv(NO_COLOR = "1")


## -----------------------------------------------------------------------------
#| label: setup
#| include: false
knitr::opts_chunk$set(echo = TRUE,
                      collapse = TRUE,
                      comment = "#>")


## -----------------------------------------------------------------------------
#| label: setup_2
#| include: false
#| message: false
#| echo: false
require(markovchain)


## -----------------------------------------------------------------------------
#| label: higherOrder
#| message: false
#| warning: false
if (requireNamespace("Rsolnp", quietly = TRUE)) {
  data(rain)
  rain_small <- rain$rain[1:150]
  fitHigherOrder(rain_small, 2)
}


## -----------------------------------------------------------------------------
#| label: compareOrders
#| message: false
#| warning: false
if (requireNamespace("Rsolnp", quietly = TRUE)) {
  compareOrders <- function(sequence, orders = 1:3) {
    byMethod <- lapply(c("lsq", "mle"), function(method) {
      fits <- lapply(orders, function(k) fitHigherOrder(sequence, k, method = method))
      sapply(fits, function(f)
        unlist(higherOrderLogLik(sequence, f, start = max(orders) + 1)[c("logLik", "BIC")]))
    })
    out <- rbind(byMethod[[1]], byMethod[[2]])
    rownames(out) <- c("logLik (lsq)", "BIC (lsq)", "logLik (mle)", "BIC (mle)")
    colnames(out) <- paste("order", orders)
    round(out, 1)
  }
  data(rain)
  print(compareOrders(rain$rain))
  data(preproglucacon)
  print(compareOrders(preproglucacon$preproglucacon))
}


## -----------------------------------------------------------------------------
#| label: selectOrder
data(rain)
rainOrder <- selectOrder(rain$rain, maxOrder = 3)
rainOrder$order
print(rainOrder$table, digits = 4, row.names = FALSE)
selectOrder(rain$rain, maxOrder = 3, criterion = "AIC")$order


## -----------------------------------------------------------------------------
#| label: berchtoldRaftery
koeberg <- as.character(read.csv(system.file("extdata", "koeberg_wind.csv",
                                             package = "markovchain"))$state)
start <- 15
n_eff <- length(koeberg) - (start - 1)

# first-order Markov chain, estimated on the transitions that enter the likelihood
Q1 <- seq2matHigh(koeberg[(start - 1):length(koeberg)], 1)
mc1 <- higherOrderLogLik(koeberg, list(lambda = 1, Q = list(Q1)), start = start)

# MTD(2) with lambda and Q as printed in Section 1.3 of the paper
Qpaper <- matrix(c(0.8301, 0.0689, 0.0077, 0.0933,
                   0.0369, 0.9012, 0.0619, 0.0000,
                   0.0155, 0.1553, 0.8070, 0.0222,
                   0.0779, 0.0000, 0.0528, 0.8693), 4, 4, byrow = TRUE)
Q <- t(Qpaper)
dimnames(Q) <- list(as.character(1:4), as.character(1:4))
mtd2 <- higherOrderLogLik(koeberg, list(lambda = c(0.7569, 0.2431), Q = list(Q, Q)),
                          start = start)

# parameters not forced to zero: 11 for the chain (one empty transition),
# 4 * 3 - 2 + (2 - 1) = 11 for the MTD(2) (two structural zeros in Q)
comparison <- data.frame(
  model = c("Markov chain, order 1", "MTD, order 2"),
  logLik = c(mc1$logLik, mtd2$logLik),
  BIC = c(-2 * mc1$logLik + 11 * log(n_eff), -2 * mtd2$logLik + 11 * log(n_eff)),
  logLik_published = c(-413.3, -393.4),
  BIC_published = c(899.1, 859.3))
print(comparison, digits = 4, row.names = FALSE)


## -----------------------------------------------------------------------------
#| label: selectOrderKoeberg
print(selectOrder(koeberg, maxOrder = 3, start = 15, parameters = "observed")$table[, 1:5],
      digits = 5, row.names = FALSE)


## -----------------------------------------------------------------------------
#| label: fitMTD
readSeries <- function(file, column)
  read.csv(system.file("extdata", file, package = "markovchain"))[[column]]
series <- list(Koeberg = readSeries("koeberg_wind.csv", "state"),
               seizures = readSeries("epileptic_seizures.csv", "seizure"))
published <- data.frame(series = rep(c("Koeberg", "seizures"), each = 2),
                        order = c(2, 3, 2, 3),
                        logLik_published = c(-393.4, -393.2, -119.5, -117.7),
                        BIC_published = c(859.3, 865.6, 254.7, 256.4))
fits <- Map(function(s, k) fitMTD(series[[s]], order = k, start = 15),
            published$series, published$order)
names(fits) <- paste(published$series, published$order)
estimated <- t(sapply(fits, function(fit) {
  zeros <- sum(fit$estimate@transitionMatrix < 1e-8)
  c(logLik = fit$logLikelihood,
    BIC = -2 * fit$logLikelihood + (fit$npar - zeros) * log(fit$nobs))
}))
print(cbind(published, round(estimated, 1)), row.names = FALSE)


## -----------------------------------------------------------------------------
#| label: fitMTDestimates
fits[["Koeberg 2"]]$lambda
round(fits[["Koeberg 2"]]$estimate@transitionMatrix, 4)


## -----------------------------------------------------------------------------
#| label: higherOrderPredict
fit <- fits[["Koeberg 2"]]
# next direction after two hours from directions 1 and then 2, and after 2 and then 2
round(higherOrderPredict(fit, rbind("1 then 2" = c(1, 2), "2 then 2" = c(2, 2))), 3)
set.seed(123)
higherOrderSimulate(24, fit, t0 = c(2, 2))


## -----------------------------------------------------------------------------
#| label: hommcObject
showClass("hommc")


## -----------------------------------------------------------------------------
#| label: hommcCreate
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


## -----------------------------------------------------------------------------
#| label: hommsales
data(sales)
head(sales)


## -----------------------------------------------------------------------------
#| label: hommcFit
#| warning: false
#| message: false
#| eval: false

# 
# # fit 8th order multivariate markov chain
# if (requireNamespace("Rsolnp", quietly = TRUE)) {
# object <- fitHighOrderMultivarMC(sales, order = 8, Norm = 2)
# }


## -----------------------------------------------------------------------------
#| label: result
#| echo: false
#| eval: false

# 
# if (requireNamespace("Rsolnp", quietly = TRUE)) {
# i <- c(1, 2, 2, 3, 4, 4, 4, 5, 5, 5)
# j <- c(2, 2, 2, 5, 2, 5, 5, 2, 4, 5)
# k <- c(1, 1, 3, 1, 8, 1, 2, 8, 1, 2)
# 
# if(object@byrow == TRUE) {
#     direction <- "(by rows)"
# } else {
#     direction <- "(by cols)"
# }
# 
# cat("Order of multivariate markov chain =", object@order, "\n")
# cat("states =", object@states, "\n")
# 
# cat("\n")
# cat("List of Lambda's and the corresponding transition matrix", direction,":\n")
# 
# for(p in 1:10) {
#     t <- 8*5*(i[p]-1) + (j[p]-1)*8
#     cat("Lambda", k[p], "(", i[p], ",", j[p], ") : ", object@Lambda[t+k[p]],"\n", sep = "")
#     cat("P", k[p], "(", i[p], ",", j[p], ") : \n", sep = "")
#     print(object@P[, , t+k[p]])
#     cat("\n")
# }
# } else {
#   print("package Rsolnp unavailable")
# }

