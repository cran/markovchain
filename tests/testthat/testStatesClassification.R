context("Checking classification of states: recurrentStates, transientStates, absorbingStates")

A <- structure(c(0, 0, 0, 0, 0, 0, 0, 0, 
                 1, 0, 0, 0, 0, 0, 0, 0, 
                 0, 0.5, 0.3, 0.3, 0, 0, 0, 0, 
                 0, 0.5, 0.7, 0.7, 0, 0, 0, 0, 
                 0, 0, 0, 0, 0, 0.3, 0.4, 0, 
                 0, 0, 0, 0, 0.4, 0, 0.5, 0, 
                 0, 0, 0, 0, 0.6, 0.7, 0, 0, 
                 0, 0, 0, 0, 0, 0, 0.1, 1), .Dim = c(8L, 8L), 
                .Dimnames = list(c("1", "2", "3", "4", "5", "6", "7", "8"), 
                                 c("1", "2", "3", "4", "5", "6", "7", "8")))
mchain <- new("markovchain", transitionMatrix = A)


mc1Matrix <- matrix(c(0, 0, 1/2, 1/2,
                      1, 0, 0, 0,
                      0, 1, 0, 0,
                      0, 1, 0, 0), 
                    ncol = 4, byrow=TRUE)
mc1 <- as(mc1Matrix, "markovchain")


mc2Matrix <- matrix(c(0, 1, 0, 0, 0, 0,
                      0.4, 0.6, 0, 0, 0, 0,
                      0.3, 0, 0.4, 0.2, 0.1, 0,
                      0, 0, 0, 0.3, 0.7, 0,
                      0, 0, 0, 0.5, 0, 0.5,
                      0, 0, 0, 0.3, 0, 0.7), 
                    nrow = 6, byrow=TRUE)
mc2 <- as(mc2Matrix,"markovchain")


mc3Matrix <- matlab::zeros(5)
mc3Matrix[1:2,1:2] <- 0.5*matlab::ones(2)
mc3Matrix[5,1] <- 1
mc3Matrix[3,3] <- 1
mc3Matrix[4,3:4] <- 0.5
mc3 <- as(mc3Matrix,"markovchain")


test_that("Test recurrent / transient / absorbing states for known Markov chains", {
  expect_equal(recurrentClasses(mc1), list(c("s1","s2","s3","s4")))
  expect_equal(recurrentClasses(mc2), list(c("s1","s2"),c("s4","s5","s6") ))
  expect_equal(transientStates(mc2), "s3")
  expect_equal(recurrentClasses(mc3), list(c("s1","s2"),c("s3")))
  expect_equal(absorbingStates(mc3), "s3")
  expect_equal(transientStates(mc3), c("s4","s5"))
  expect_equal(recurrentClasses(mchain), list(c("3", "4"), c("8")))
  expect_equal(transientStates(mchain), c("1", "2", "5", "6", "7"))
  expect_equal(absorbingStates(mchain), "8")
})


test_that("A state is absorbing iff it is singleton recurrent class", {
  
  for (mc in allMCs) {
    states <- mc$states
    classes <- mc$recurrentClasses
    absorbing <- mc$absorbingStates
    
    expect_true(.testthatAbsorbingAreRecurrentClassRcpp(absorbing, classes))
  }
})


test_that("Recurrent states and transient states are a partition of states", {
  
  for (mc in allMCs) {
    states <- mc$states
    recurrentStates <- mc$recurrentStates
    transientStates <- mc$transientStates
    states <- mc$states
    statesUnion <- sort(unique(append(recurrentStates, transientStates)))
    
    expect_equal(statesUnion, sort(states))
  }
})


test_that("hittingProb(i,i) < 1 for i a transient state", {
  
  for (mc in allMCs) {
    transStates <- mc$transientStates
    hitting <- mc$hittingProbabilities
    transientHittingLessOne <- all(sapply(transStates, function(s){ hitting[s, s] < 1}))
    expect_true(transientHittingLessOne)
  }
})


test_that("All states are recurrent in a identity Markov chain", {
  
  for (mc in allDiagonalMCs) {
    states <- mc$states
    recurrentStates <- mc$recurrentStates
    expect_true(setequal(recurrentStates, states))
  }
})

test_that("If Markov chain is irreducible then all states are recurrent", {
    
  for (mc in allMCs) {
    states <- mc$states
    recurrent <- mc$recurrentStates
    irreducible <- mc$irreducible
    allRecurrent <- setequal(states, recurrent)
    
    if (irreducible)
      expect_true(allRecurrent)
  }
})


test_that("If there are transient states then Markov chain is not irreducible", {
  
  for (mc in allMCs) {
    states <- mc$states
    transient  <- mc$transientStates
    irreducible <- mc$irreducible
    
    if (length(transient) > 0)
      expect_false(irreducible)
  }
})


context("Checking recurrentClasses method")

P <- matlab::zeros(10)
P[1, c(1, 3)] <- 1/2;
P[2, 2] <- 1/3; P[2,7] <- 2/3;
P[3, 1] <- 1;
P[4, 5] <- 1;
P[5, c(4, 5, 9)] <- 1/3;
P[6, 6] <- 1;
P[7, 7] <- 1/4; P[7,9] <- 3/4;
P[8, c(3, 4, 8, 10)] <- 1/4;
P[9, 2] <- 1;
P[10, c(2, 5, 10)] <- 1/3;
rownames(P) <- letters[1:10]
colnames(P) <- letters[1:10]
probMc <- new("markovchain", transitionMatrix = P,
              name = "Probability MC")

test_that("Checking recurrent classes for known Markov chains", {
  expect_equal(recurrentClasses(probMc), list(c("a", "c")
                                              , c("b", "g", "i")
                                              , c("f")))
})


test_that("hittingProb(i,j) = 1 for i, j in same recurrent class, hittingProb(i, k) = 0 for k otherwise", {
  
  for (mc in allMCs) {
    byrow <- mc$byrow
    states <- mc$states
    hitting <- mc$hittingProbabilities
    recurrentClasses <- mc$recurrentClasses
    expect_true(.testthatRecurrentHittingRcpp(recurrentClasses, hitting, states, byrow))
  }
})


test_that("Union of recurrentClasses is recurrentStates", {
  
  for (mc in allMCs) {
    recClasses <- mc$recurrentClasses
    recStates <- as.character(sort(unlist(recClasses)))
    target <- sort(mc$recurrentStates)
    
    expect_equal(recStates, target)
  }
})

test_that("Recurrent classes are disjoint", {
  
  for (mc in allMCs) {
    recClasses <- mc$recurrentClasses
    lengthRecClasses <- sapply(recClasses, function(c){ length(c) })
    numRecurrentStates <- ifelse(length(recClasses) > 0, sum(lengthRecClasses), 0)
    numUnion <- length(unique(unlist(recClasses)))
    
    expect_equal(numRecurrentStates, numUnion)
  }
})

context("Checking transientClasses method")

test_that("Checking recurrent classes for known Markov chains", {
  expect_equal(transientClasses(probMc), list(c("d", "e"), c("h"), c("j")))
})


test_that("Union of transientClases is transientStates", {
  
  for (mc in allMCs) {
    transClasses <- mc$transientClasses
    # as.character forces the result to be a char vector when it is empty
    transStates <- as.character(sort(unlist(transClasses)))
    target <- sort(mc$transientStates)
    
    expect_equal(transStates, target)
  }
})

test_that("Transient classes are disjoint", {
  
  for (mc in allMCs) {
    transClasses <- mc$transientClasses
    lengthTransClasses <- sapply(transClasses, function(c){ length(c) })
    numTransientStates <- ifelse(length(transClasses) > 0, sum(lengthTransClasses), 0)
    numUnion <- length(unique(unlist(transClasses)))
    
    expect_equal(numTransientStates, numUnion)
  }
})