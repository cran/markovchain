# given a markovchain object is it possible to reach goal state from 
# a given state

#' @name is.accessible
#' @title Verify if a state j is reachable from state i.
#' @description This function verifies if a state is reachable from another, i.e., 
#'              if there exists a path that leads to state j leaving from state i with 
#'              positive probability
#'              
#' @param object A \code{markovchain} object.
#' @param from The name of state "i" (beginning state).
#' @param to The name of state "j" (ending state).
#' 
#' @details It wraps an internal function named \code{reachabilityMatrix}.
#' @return A boolean value.
#' 
#' @references James Montgomery, University of Madison
#' 
#' @author Giorgio Spedicato, Ignacio Cordón
#' @seealso \code{is.irreducible}
#' 
#' @examples 
#' statesNames <- c("a", "b", "c")
#' markovB <- new("markovchain", states = statesNames, 
#'                transitionMatrix = matrix(c(0.2, 0.5, 0.3,
#'                                              0,   1,   0,
#'                                            0.1, 0.8, 0.1), nrow = 3, byrow = TRUE, 
#'                                          dimnames = list(statesNames, statesNames)
#'                                         )
#'                )
#' is.accessible(markovB, "a", "c")
#' 
#' @exportMethod is.accessible
setGeneric("is.accessible", function(object, from, to) standardGeneric("is.accessible"))

setMethod("is.accessible", c("markovchain", "character", "character"), 
  function(object, from, to) {
    # O(n²) procedure to see if to state is reachable starting at from state
    return(.isAccessibleRcpp(object, from, to))
  }
)

setMethod("is.accessible", c("markovchain", "missing", "missing"), 
  function(object, from, to) {
    .reachabilityMatrixRcpp(object)
  }
)

# a markov chain is irreducible if it is composed of only one communicating class

#' @name is.irreducible
#' @title Function to check if a Markov chain is irreducible (i.e. ergodic)
#' @description This function verifies whether a \code{markovchain} object transition matrix 
#'              is composed by only one communicating class.
#' @param object A \code{markovchain} object
#' 
#' @details It is based on \code{.communicatingClasses} internal function.
#' @return A boolean values.
#' 
#' @references Feres, Matlab listings for Markov Chains.
#' @author Giorgio Spedicato
#' 
#' @seealso \code{\link{summary}}
#' 
#' @examples 
#' statesNames <- c("a", "b")
#' mcA <- new("markovchain", transitionMatrix = matrix(c(0.7,0.3,0.1,0.9),
#'                                              byrow = TRUE, nrow = 2, 
#'                                              dimnames = list(statesNames, statesNames)
#'            ))
#' is.irreducible(mcA)
#' 
#' @exportMethod is.irreducible
setGeneric("is.irreducible", function(object) standardGeneric("is.irreducible"))

setMethod("is.irreducible", "markovchain", function(object) {
  .isIrreducibleRcpp(object)
})


# what this function will do?
# It calculates the probability to go from given state
# to all other states in k steps
# k varies from 1 to n

#' @name firstPassage
#' @title First passage across states
#' @description This function compute the first passage probability in states
#' 
#' @param object A \code{markovchain} object
#' @param state Initial state
#' @param n Number of rows on which compute the distribution
#' 
#' @details Based on Feres' Matlab listings
#' @return A matrix of size 1:n x number of states showing the probability of the 
#'         first time of passage in states to be exactly the number in the row.
#'
#' @references Renaldo Feres, Notes for Math 450 Matlab listings for Markov chains
#' 
#' @author Giorgio Spedicato
#' @seealso \code{\link{conditionalDistribution}}
#' 
#' @examples 
#' simpleMc <- new("markovchain", states = c("a", "b"),
#'                  transitionMatrix = matrix(c(0.4, 0.6, .3, .7), 
#'                                     nrow = 2, byrow = TRUE))
#' firstPassage(simpleMc, "b", 20)
#'
#' @export
firstPassage <- function(object, state, n) {
  P <- object@transitionMatrix
  stateNames <- states(object)
  
  # row number
  i <- which(stateNames == state)

  
  outMatr <- .firstpassageKernelRcpp(P = P, i = i, n = n)
  colnames(outMatr) <- stateNames
  rownames(outMatr) <- 1:n
  return(outMatr)
}




#' function to calculate first passage probabilities
#' 
#' @description The function calculates first passage probability for a subset of
#' states given an initial state.
#' 
#' @param object a markovchain-class object
#' @param state intital state of the process (charactervector)
#' @param set set of states A, first passage of which is to be calculated
#' @param n Number of rows on which compute the distribution
#' 
#' @return A vector of size n showing the first time proabilities
#' @references
#' Renaldo Feres, Notes for Math 450 Matlab listings for Markov chains;
#' MIT OCW, course - 6.262, Discrete Stochastic Processes, course-notes, chap -05
#' 
#' @author Vandit Jain
#' 
#' @seealso \code{\link{firstPassage}}
#' @examples 
#' statesNames <- c("a", "b", "c")
#' markovB <- new("markovchain", states = statesNames, transitionMatrix =
#' matrix(c(0.2, 0.5, 0.3,
#'          0, 1, 0,
#'          0.1, 0.8, 0.1), nrow = 3, byrow = TRUE,
#'        dimnames = list(statesNames, statesNames)
#' ))
#' firstPassageMultiple(markovB,"a",c("b","c"),4)  
#' 
#' @export 
firstPassageMultiple <- function(object,state,set, n){
  
  # gets the transition matrix
  P <- object@transitionMatrix
  
  # character vector of states of the markovchain
  stateNames <- states(object)
  
  k <- -1
  k <- which(stateNames == state)
  if(k==-1)
    stop("please provide a valid initial state")
  
  # gets the set in numeric vector
  setno <- rep(0,length(set))
  for(i in 1:length(set))
  {
    setno[i] = which(set[i] == stateNames)
    if(setno[i] == 0)
      stop("please provide proper set of states")
  }
  
  # calls Rcpp implementation
  outMatr <- .firstPassageMultipleRCpp(P,k,setno,n)
  
  #sets column and row names of output
  colnames(outMatr) <- "set"
  rownames(outMatr) <- 1:n
  return(outMatr)
}


#' @name communicatingClasses
#' @rdname structuralAnalysis
#' @aliases transientStates recurrentStates absorbingStates communicatingClasses
#'   transientClasses recurrentClasses
#' @title Various function to perform structural analysis of DTMC
#' @description These functions return absorbing and transient states of the \code{markovchain} objects.
#' 
#' @param object A \code{markovchain} object.
#' 
#' @return
#' \describe{
#'   \item{\code{period}}{returns a integer number corresponding to the periodicity of the Markov 
#'     chain (if it is irreducible)}
#'   \item{\code{absorbingStates}}{returns a character vector with the names of the absorbing 
#'     states in the Markov chain}
#'   \item{\code{communicatingClasses}}{returns a list in which each slot contains the names of
#'     the states that are in that communicating class}
#'   \item{\code{recurrentClasses}}{analogously to \code{communicatingClasses}, but with 
#'     recurrent classes}
#'   \item{\code{transientClasses}}{analogously to \code{communicatingClasses}, but with 
#'     transient classes}
#'   \item{\code{transientStates}}{returns a character vector with all the transient states
#'     for the Markov chain}
#'   \item{\code{recurrentStates}}{returns a character vector with all the recurrent states 
#'     for the Markov chain}
#'   \item{\code{canonicForm}}{returns the Markov chain reordered by a permutation of states 
#'     so that we have blocks submatrices for each of the recurrent classes and a collection 
#'     of rows in the end for the transient states}
#' }
#' 
#' @references Feres, Matlab listing for markov chain.
#' 
#' @author Giorgio Alfredo Spedicato, Ignacio Cordón
#' 
#' @seealso \code{\linkS4class{markovchain}}
#' 
#' @examples 
#' statesNames <- c("a", "b", "c")
#' mc <- new("markovchain", states = statesNames, transitionMatrix =
#'           matrix(c(0.2, 0.5, 0.3,
#'                    0,   1,   0,
#'                    0.1, 0.8, 0.1), nrow = 3, byrow = TRUE,
#'                  dimnames = list(statesNames, statesNames))
#'          )
#' 
#' communicatingClasses(mc)
#' recurrentClasses(mc)
#' recurrentClasses(mc)
#' absorbingStates(mc)
#' transientStates(mc)
#' recurrentStates(mc)
#' canonicForm(mc)
#' 
#' # periodicity analysis
#' A <- matrix(c(0, 1, 0, 0, 0.5, 0, 0.5, 0, 0, 0.5, 0, 0.5, 0, 0, 1, 0), 
#'             nrow = 4, ncol = 4, byrow = TRUE)
#' mcA <- new("markovchain", states = c("a", "b", "c", "d"), 
#'           transitionMatrix = A,
#'           name = "A")
#'
#' is.irreducible(mcA) #true
#' period(mcA) #2
#'
#' # periodicity analysis
#' B <- matrix(c(0, 0, 1/2, 1/4, 1/4, 0, 0,
#'                    0, 0, 1/3, 0, 2/3, 0, 0,
#'                    0, 0, 0, 0, 0, 1/3, 2/3,
#'                    0, 0, 0, 0, 0, 1/2, 1/2,
#'                    0, 0, 0, 0, 0, 3/4, 1/4,
#'                    1/2, 1/2, 0, 0, 0, 0, 0,
#'                    1/4, 3/4, 0, 0, 0, 0, 0), byrow = TRUE, ncol = 7)
#' mcB <- new("markovchain", transitionMatrix = B)
#' period(mcB)
#' 
#' @exportMethod communicatingClasses
setGeneric("communicatingClasses", function(object) standardGeneric("communicatingClasses"))

setMethod("communicatingClasses", "markovchain", function(object) {
  return(.communicatingClassesRcpp(object))
})


# A communicating class will be a recurrent class if 
# there is no outgoing edge from this class
# Recurrent classes are subset of communicating classes
#' @rdname structuralAnalysis
#' 
#' @exportMethod recurrentClasses
setGeneric("recurrentClasses", function(object) standardGeneric("recurrentClasses"))

setMethod("recurrentClasses", "markovchain", function(object) {
  return(.recurrentClassesRcpp(object))
})


# A communicating class will be a transient class iff
# there is an outgoing edge from this class to an state
# outside of the class
# Transient classes are subset of communicating classes
#' @rdname structuralAnalysis
#' 
#' @exportMethod transientClasses
setGeneric("transientClasses", function(object) standardGeneric("transientClasses"))

setMethod("transientClasses", "markovchain", function(object) {
  return(.transientClassesRcpp(object))
})


#' @rdname structuralAnalysis
#' 
#' @exportMethod transientStates
setGeneric("transientStates", function(object) standardGeneric("transientStates"))


setMethod("transientStates", "markovchain", function(object) {
    .transientStatesRcpp(object)
  }
)


#' @rdname structuralAnalysis
#' 
#' @exportMethod recurrentStates
setGeneric("recurrentStates", function(object) standardGeneric("recurrentStates"))


setMethod("recurrentStates", "markovchain", function(object) {
    .recurrentStatesRcpp(object)
  }
)

# generic function to extract absorbing states

#' @rdname structuralAnalysis
#' 
#' @exportMethod absorbingStates
setGeneric("absorbingStates", function(object) standardGeneric("absorbingStates"))

setMethod("absorbingStates", "markovchain", function(object) {
    .absorbingStatesRcpp(object)
  }
)


#' @rdname structuralAnalysis
#' 
#' @exportMethod canonicForm
setGeneric("canonicForm", function(object) standardGeneric("canonicForm"))

setMethod("canonicForm", "markovchain", function(object) {
    .canonicFormRcpp(object)
  }
)


#' @title Calculates committor of a markovchain object with respect to set A, B
#' 
#' @description Returns the probability of hitting states rom set A before set B 
#' with different initial states
#' 
#' @usage committorAB(object,A,B,p)
#' 
#' @param object a markovchain class object
#' @param A a set of states
#' @param B a set of states
#' @param p initial state (default value : 1)
#' 
#' @details The function solves a system of linear equations to calculate probaility that the process hits
#' a state from set A before any state from set B
#' 
#' @return Return a vector of probabilities in case initial state is not provided else returns a number
#' 
#' @examples 
#' transMatr <- matrix(c(0,0,0,1,0.5,
#'                       0.5,0,0,0,0,
#'                       0.5,0,0,0,0,
#'                       0,0.2,0.4,0,0,
#'                       0,0.8,0.6,0,0.5),
#'                       nrow = 5)
#' object <- new("markovchain", states=c("a","b","c","d","e"),transitionMatrix=transMatr)
#' committorAB(object,c(5),c(3))
#' 
#' @export
committorAB <- function(object,A,B,p=1) {
  
  if(!is(object,"markovchain"))
    stop("please provide a valid markovchain object")
  
  matrix <- object@transitionMatrix
  
  noofstates <- length(object@states)
  
  for(i in length(A))
  {
    if(A[i] <= 0 || A[i] > noofstates)
      stop("please provide a valid set A")
  }
  
  for(i in length(B))
  {
    if(B[i] <= 0 || B[i] > noofstates)
      stop("please provide a valid set B")
  }
  
  for(i in 1:noofstates)
  {
    if(i %in% A && i %in% B)
      stop("intersection of set A and B in not null")
  }
  
  if(p <=0 || p > noofstates)
    stop("please provide a valid initial state")
  
  I <- diag(noofstates)
  
  matrix <- matrix - I
  
  A_size = length(A)
  B_size = length(B)
  
  # sets the matrix according to the provided states
  for(i in 1:A_size)
  {
    for(j in 1:noofstates)
    {
      if(A[i]==j)
        matrix[A[i],j] = 1
      else
        matrix[A[i],j] = 0
    }
  }
  
  # sets the matrix according to the provided states
  for(i in 1:B_size)
  {
    for(j in 1:noofstates)
    {
      if(B[i]==j)
        matrix[B[i],j] = 1
      else
        matrix[B[i],j] = 0
    }
  }
  
  # initialises b in the equation the system of equation AX =b
  b <- rep(0,noofstates)
  
  
  for(i in 1:A_size)
  {
    b[A[i]] = 1
  }
  
  # solve AX = b according using solve function from base package
  out <- solve(matrix,b)
  
  
  if(missing(p))
    return(out)
  else
    return(out[p])
}


#' Expected Rewards for a markovchain
#' 
#' @description Given a markovchain object and reward values for every state,
#' function calculates expected reward value after n steps.
#' 
#' @usage expectedRewards(markovchain,n,rewards)
#' 
#' @param markovchain the markovchain-class object
#' @param n no of steps of the process
#' @param rewards vector depicting rewards coressponding to states
#' 
#' @details the function uses a dynamic programming approach to solve a 
#' recursive equation described in reference.
#' 
#' @return
#' returns a vector of expected rewards for different initial states
#' 
#' @author Vandit Jain
#' 
#' @references Stochastic Processes: Theory for Applications, Robert G. Gallager,
#' Cambridge University Press
#' 
#' @examples 
#' transMatr<-matrix(c(0.99,0.01,0.01,0.99),nrow=2,byrow=TRUE)
#' simpleMc<-new("markovchain", states=c("a","b"),
#'              transitionMatrix=transMatr)
#' expectedRewards(simpleMc,1,c(0,1))
#' @export
expectedRewards <- function(markovchain, n, rewards) {
  
  # gets the transition matrix
  matrix <- markovchain@transitionMatrix
  
  # Rcpp implementation of the function
  out <- .expectedRewardsRCpp(matrix,n, rewards)
  
  noofStates <- length(states(markovchain))
  
  result <- rep(0,noofStates)
  
  for(i in 1:noofStates)
    result[i] = out[i]
  
  #names(result) <- states(markovchain)
  return(result)
}

#' Expected first passage Rewards for a set of states in a markovchain
#' 
#' @description Given a markovchain object and reward values for every state,
#' function calculates expected reward value for a set A of states after n 
#' steps. 
#'  
#' @usage expectedRewardsBeforeHittingA(markovchain, A, state, rewards, n)
#'  
#' @param markovchain the markovchain-class object
#' @param A set of states for first passage expected reward
#' @param state initial state
#' @param rewards vector depicting rewards coressponding to states
#' @param n no of steps of the process
#'  
#' @details The function returns the value of expected first passage 
#' rewards given rewards coressponding to every state, an initial state
#' and number of steps.
#'  
#' @return returns a expected reward (numerical value) as described above
#'  
#' @author Sai Bhargav Yalamanchi, Vandit Jain
#'  
#' @export
expectedRewardsBeforeHittingA <- function(markovchain, A, state, rewards, n) {
  
  ## gets the markovchain matrix
  matrix <- markovchain@transitionMatrix
  
  # gets the names of states
  stateNames <- states(markovchain)
  
  # no of states
  S <- length(stateNames)
  
  # vectors for states in S-A
  SAno <- rep(0,S-length(A))
  rewardsSA <- rep(0,S-length(A))
  
  # for initialisation for set S-A 
  i=1
  ini = -1
  for(j in 1:length(stateNames))
  {
    if(!(stateNames[j] %in% A)){
      SAno[i] = j
      rewardsSA[i] = rewards[j]
      if(stateNames[j] == state)
        ini = i
      i = i+1
    }
  }
  
  ## get the matrix coressponding to S-A
  matrix <- matrix[SAno,SAno]
  
  ## cals the cpp implementation
  out <- .expectedRewardsBeforeHittingARCpp(matrix, ini, rewardsSA, n)
  
  return(out)
  
}




#' Mean First Passage Time for irreducible Markov chains
#'
#' @description Given an irreducible (ergodic) markovchain object, this function
#'   calculates the expected number of steps to reach other states
#'
#' @param object the markovchain object
#' @param destination a character vector representing the states respect to
#'   which we want to compute the mean first passage time. Empty by default
#'
#' @details For an ergodic Markov chain it computes: 
#' \itemize{ 
#'   \item If destination is empty, the average first time (in steps) that takes
#'   the Markov chain to go from initial state i to j. (i, j) represents that 
#'   value in case the Markov chain is given row-wise, (j, i) in case it is given
#'   col-wise. 
#'   \item If destination is not empty, the average time it takes us from the 
#'   remaining states to reach the states in \code{destination} 
#' }
#'
#' @return a Matrix of the same size with the average first passage times if
#'   destination is empty, a vector if destination is not
#'
#' @author Toni Giorgino, Ignacio Cordón
#'
#' @references C. M. Grinstead and J. L. Snell. Introduction to Probability.
#' American Mathematical Soc., 2012.
#'
#' @examples
#' m <- matrix(1 / 10 * c(6,3,1,
#'                        2,3,5,
#'                        4,1,5), ncol = 3, byrow = TRUE)
#' mc <- new("markovchain", states = c("s","c","r"), transitionMatrix = m)
#' meanFirstPassageTime(mc, "r")
#'
#'
#' # Grinstead and Snell's "Oz weather" worked out example
#' mOz <- matrix(c(2,1,1,
#'                 2,0,2,
#'                 1,1,2)/4, ncol = 3, byrow = TRUE)
#'
#' mcOz <- new("markovchain", states = c("s", "c", "r"), transitionMatrix = mOz)
#' meanFirstPassageTime(mcOz)
#'
#' @export meanFirstPassageTime
setGeneric("meanFirstPassageTime", function(object, destination) {
  standardGeneric("meanFirstPassageTime")
})


setMethod("meanFirstPassageTime",  signature("markovchain", "missing"),
  function(object, destination) {
    destination = character()
    .meanFirstPassageTimeRcpp(object, destination)
  }
)

setMethod("meanFirstPassageTime",  signature("markovchain", "character"),
  function(object, destination) {
    states <- object@states
    incorrectStates <- setdiff(destination, states)
    
    if (length(incorrectStates) > 0)
      stop("Some of the states you provided in destination do not match states from the markovchain")

    result <- .meanFirstPassageTimeRcpp(object, destination)
    asVector <- as.vector(result)
    names(asVector) <- colnames(result)
    
    asVector
  }
)

#' Mean recurrence time
#'
#' @description Computes the expected time to return to a recurrent state
#'   in case the Markov chain starts there
#'
#' @usage meanRecurrenceTime(object)
#'
#' @param object the markovchain object
#'
#' @return For a Markov chain it outputs is a named vector with the expected 
#'   time to first return to a state when the chain starts there.
#'   States present in the vector are only the recurrent ones. If the matrix
#'   is ergodic (i.e. irreducible), then all states are present in the output
#'   and order is the same as states order for the Markov chain
#'
#' @author Ignacio Cordón
#'
#' @references C. M. Grinstead and J. L. Snell. Introduction to Probability.
#' American Mathematical Soc., 2012.
#'
#' @examples
#' m <- matrix(1 / 10 * c(6,3,1,
#'                        2,3,5,
#'                        4,1,5), ncol = 3, byrow = TRUE)
#' mc <- new("markovchain", states = c("s","c","r"), transitionMatrix = m)
#' meanRecurrenceTime(mc)
#'
#' @export meanRecurrenceTime
setGeneric("meanRecurrenceTime", function(object) {
  standardGeneric("meanRecurrenceTime")
})

setMethod("meanRecurrenceTime", "markovchain", function(object) {
  .meanRecurrenceTimeRcpp(object)
})


#' Mean absorption time
#'
#' @description Computes the expected number of steps to go from any of the
#'   transient states to any of the recurrent states. The Markov chain should
#'   have at least one transient state for this method to work
#'
#' @usage meanAbsorptionTime(object)
#'
#' @param object the markovchain object
#'
#' @return A named vector with the expected number of steps to go from a
#'   transient state to any of the recurrent ones
#'
#' @author Ignacio Cordón
#'
#' @references C. M. Grinstead and J. L. Snell. Introduction to Probability.
#' American Mathematical Soc., 2012.
#'
#' @examples
#' m <- matrix(c(1/2, 1/2, 0,
#'               1/2, 1/2, 0,
#'                 0, 1/2, 1/2), ncol = 3, byrow = TRUE)
#' mc <- new("markovchain", states = letters[1:3], transitionMatrix = m)
#' times <- meanAbsorptionTime(mc)
#'
#' @export meanAbsorptionTime
setGeneric("meanAbsorptionTime", function(object) {
  standardGeneric("meanAbsorptionTime")
})

setMethod("meanAbsorptionTime",  "markovchain", function(object) {
  .meanAbsorptionTimeRcpp(object)
})

#' Absorption probabilities
#'
#' @description Computes the absorption probability from each transient
#'   state to each recurrent one (i.e. the (i, j) entry or (j, i), in a 
#'   stochastic matrix by columns, represents the probability that the
#'   first not transient state we can go from the transient state i is j
#'   (and therefore we are going to be absorbed in the communicating
#'   recurrent class of j)
#'
#' @usage absorptionProbabilities(object)
#'
#' @param object the markovchain object
#'
#' @return A named vector with the expected number of steps to go from a
#'   transient state to any of the recurrent ones
#'
#' @author Ignacio Cordón
#'
#' @references C. M. Grinstead and J. L. Snell. Introduction to Probability.
#' American Mathematical Soc., 2012.
#'
#' @examples
#' m <- matrix(c(1/2, 1/2, 0,
#'               1/2, 1/2, 0,
#'                 0, 1/2, 1/2), ncol = 3, byrow = TRUE)
#' mc <- new("markovchain", states = letters[1:3], transitionMatrix = m)
#' absorptionProbabilities(mc)
#'
#' @export absorptionProbabilities
setGeneric("absorptionProbabilities", function(object) {
  standardGeneric("absorptionProbabilities")
})

setMethod("absorptionProbabilities",  "markovchain", function(object) {
  .absorptionProbabilitiesRcpp(object)
})


#' @title Check if a DTMC is regular
#' 
#' @description Function to check wether a DTCM is regular
# 
#' @details A Markov chain is regular if some of the powers of its matrix has all elements 
#'   strictly positive
#' 
#' @param object a markovchain object
#'
#' @return A boolean value
#'
#' @author Ignacio Cordón
#' @references Matrix Analysis. Roger A.Horn, Charles R.Johnson. 2nd edition. 
#'   Corollary 8.5.8, Theorem 8.5.9
#'
#' 
#' @examples 
#' P <- matrix(c(0.5,  0.25, 0.25,
#'               0.5,     0, 0.5,
#'               0.25, 0.25, 0.5), nrow = 3)
#' colnames(P) <- rownames(P) <- c("R","N","S")
#' ciao <- as(P, "markovchain")
#' is.regular(ciao)
#' 
#' @seealso \code{\link{is.irreducible}}
#' 
#' @exportMethod is.regular
setGeneric("is.regular", function(object) standardGeneric("is.regular"))

setMethod("is.regular", "markovchain", function(object) {
  .isRegularRcpp(object)
})


#' Hitting probabilities for markovchain
#' 
#' @description Given a markovchain object,
#' this function calculates the probability of ever arriving from state i to j
#' 
#' @usage hittingProbabilities(object)
#' 
#' @param object the markovchain-class object
#' 
#' @return a matrix of hitting probabilities
#' 
#' @author Ignacio Cordón
#' 
#' @references R. Vélez, T. Prieto, Procesos Estocásticos, Librería UNED, 2013
#' 
#' @examples
#' M <- markovchain:::zeros(5)
#' M[1,1] <- M[5,5] <- 1
#' M[2,1] <- M[2,3] <- 1/2
#' M[3,2] <- M[3,4] <- 1/2
#' M[4,2] <- M[4,5] <- 1/2
#' 
#' mc <- new("markovchain", transitionMatrix = M)
#' hittingProbabilities(mc)
#' 
#' @exportMethod hittingProbabilities
setGeneric("hittingProbabilities", function(object) standardGeneric("hittingProbabilities"))

setMethod("hittingProbabilities", "markovchain", function(object) {
  .hittingProbabilitiesRcpp(object)
})



#' Mean num of visits for markovchain, starting at each state
#' 
#' @description Given a markovchain object, this function calculates 
#' a matrix where the element (i, j) represents the expect number of visits
#' to the state j if the chain starts at i (in a Markov chain by columns it
#' would be the element (j, i) instead)
#' 
#' @usage meanNumVisits(object)
#' 
#' @param object the markovchain-class object
#' 
#' @return a matrix with the expect number of visits to each state
#' 
#' @author Ignacio Cordón
#' 
#' @references R. Vélez, T. Prieto, Procesos Estocásticos, Librería UNED, 2013
#' 
#' @examples
#' M <- markovchain:::zeros(5)
#' M[1,1] <- M[5,5] <- 1
#' M[2,1] <- M[2,3] <- 1/2
#' M[3,2] <- M[3,4] <- 1/2
#' M[4,2] <- M[4,5] <- 1/2
#' 
#' mc <- new("markovchain", transitionMatrix = M)
#' meanNumVisits(mc)
#' 
#' @exportMethod meanNumVisits
setGeneric("meanNumVisits", function(object) standardGeneric("meanNumVisits"))

setMethod("meanNumVisits", "markovchain", function(object) {
  .minNumVisitsRcpp(object)
})


setMethod(
  "steadyStates",
  "markovchain", 
  function(object) {
    .steadyStatesRcpp(object)
  }
)


#' @exportMethod summary
setGeneric("summary")

# summary method for markovchain class
# lists: closed, transient classes, irreducibility, absorbint, transient states
setMethod("summary", signature(object = "markovchain"),
  function(object){
    
    # list of closed, recurrent and transient classes
    outs <- .summaryKernelRcpp(object)
    
    # display name of the markovchain object
    cat(object@name," Markov chain that is composed by:", "\n")
    
    # number of closed classes
    check <- length(outs$closedClasses)
    
    cat("Closed classes:","\n")
    
    # display closed classes
    if(check == 0) cat("NONE", "\n") else {
      for(i in 1:check) cat(outs$closedClasses[[i]], "\n")
    }
    
    # number of recurrent classes
    check <- length(outs$recurrentClasses)
    
    cat("Recurrent classes:", "\n")
    
    # display recurrent classes
    if(check == 0) cat("NONE", "\n") else {
      cat("{")
      cat(outs$recurrentClasses[[1]], sep = ",")
      cat("}")
      if(check > 1) {
        for(i in 2:check) {
          cat(",{")
          cat(outs$recurrentClasses[[i]], sep = ",")
          cat("}")
        }
      }
      cat("\n")
    }
    
    # number of transient classes
    check <- length(outs$transientClasses)
    
    cat("Transient classes:","\n")
    
    # display transient classes
    if(check == 0) cat("NONE", "\n") else {
      cat("{")
      cat(outs$transientClasses[[1]], sep = ",")
      cat("}")
      if(check > 1) { 
        for(i in 2:check) {
          cat(",{")
          cat(outs$transientClasses[[i]], sep = ",")
          cat("}")
        }
      }
      cat("\n")
    }
    
    # bool to say about irreducibility of markovchain
    irreducibility <- is.irreducible(object)
    
    if(irreducibility) 
      cat("The Markov chain is irreducible", "\n") 
    else cat("The Markov chain is not irreducible", "\n")
    
    # display absorbing states
    check <- absorbingStates(object)
    if(length(check) == 0) check <- "NONE"
    cat("The absorbing states are:", check )
    cat("\n")
    
    # return outs
    # useful when user will assign the value returned
    invisible(outs) 
  }
)
