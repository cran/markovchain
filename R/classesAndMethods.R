#define Markov Chain class

setClass("markovchain", #classe lifetable
         representation(states="character",byrow="logical",transitionMatrix="matrix"),
         prototype(states=c("a","b"), byrow=TRUE, transitionMatrix=matrix(data=c(0,1,1,0),
                                           nrow=2,
                                           byrow=TRUE, 
                                           dimnames=list(c("a","b"), c("a","b"))
                                                              )
         ))

		 
setMethod("initialize",
    signature(.Object = "markovchain"),
    function (.Object, states, byrow, transitionMatrix,...) 
    {
		if(missing(transitionMatrix)) transitionMatrix=matrix(data=c(0,1,1,0), #create a naive matrix
                                           nrow=2,
                                           byrow=TRUE, 
                                           dimnames=list(c("a","b"), c("a","b")))
	
		#check names of transition matrix
		if(all(is.null(rownames(transitionMatrix)), is.null(colnames(transitionMatrix)))==TRUE) { #if all names are missing it initializes them to "1", "2",...
			if(missing(states)) {
			nr=nrow(transitionMatrix)
			stateNames<-as.character(seq(1:nr))
			} else stateNames=states
			
			rownames(transitionMatrix)=stateNames
			colnames(transitionMatrix)=stateNames
		} else if(is.null(rownames(transitionMatrix))) { #fix when rownames null
		  rownames(transitionMatrix)=colnames(transitionMatrix)
		} else if(is.null(colnames(transitionMatrix))) { #fix when colnames null
			colnames(transitionMatrix)=rownames(transitionMatrix)
		} else if(!setequal(rownames(transitionMatrix),colnames(transitionMatrix)))  colnames(transitionMatrix)=rownames(transitionMatrix) #fix when different
		if(missing(states)) states=rownames(transitionMatrix) #assign
		if(missing(byrow)) byrow=TRUE #set byrow as true by default
		callNextMethod(.Object, states = states, byrow = byrow, transitionMatrix=transitionMatrix,...)
    }
)
		 
# #test
# stateNames=c("a","b")
# ciao<-new("markovchain", states=stateNames, transitionMatrix=matrix(c(1,0,0,1),byrow=TRUE, nrow=2, 
#                                                  dimnames=list(stateNames,stateNames)
#                                                  ))

#function to check if a value can be deemed a probability

.isProb<-function(prob)
{
	if(class(prob)!="numeric") return(FALSE)
	if(prob<0 | prob >1) return(FALSE)
	return(TRUE)
}


#generic function to print out states

setGeneric("states", function(object) standardGeneric("states"))
setMethod("states","markovchain", 
          function(object) {
            out<-object@states
            return(out)
          }
)

#generic function to get the dim of a markovchain

setMethod("dim","markovchain", 
		function(x) {
			out<-nrow(x@transitionMatrix)
			return(out)
		}
)

setValidity("markovchain",
            function(object) {
              check<-NULL
            if(any(sapply(as.numeric(object@transitionMatrix),.isProb))==FALSE) check<-"Error! Some elements are not probabilities" #checks if probability
            if(object@byrow==TRUE) {if(any(round(rowSums(object@transitionMatrix),5)!=1)) check<-"Error! Row sums not equal to one"} else {if(any(round(colSums(object@transitionMatrix),5)!=1)) check<-"Error! Col sums not equal to one"} #checks if col sums not equal to one
            if(nrow(object@transitionMatrix)!=ncol(object@transitionMatrix)) check<-"Error! Not squared matrix" #check if squalre matrix
            if(!setequal(colnames(object@transitionMatrix),object@states)) check<-"Error! Colnames <> states" #checks if 
            if(!setequal(rownames(object@transitionMatrix),object@states)) check<-"Error! Rownames <> states"
              if(is.null(check)) return(TRUE) else 
                return(check)
            }
)


#internal function to extract steadyStates 

.mcEigen<-function(matr, transpose=TRUE)
{
  if(transpose) tMatr<-t(matr) else tMatr<-matr #need to transpose
  eigenResults<-eigen(x=tMatr,symmetric=FALSE) #perform the eigenvalue extraction
  onesIndex<-which(round(eigenResults$values,3)==1) #takes the one eigenvalue
  #do the following: 1:get eigenvectors whose eigenvalues==1
	#2: normalize
  if(transpose==TRUE)
  {
	  eigenTake<-as.matrix(t(eigenResults$vectors[,onesIndex])) 
	  out<-eigenTake/rowSums(eigenTake)
	  
  } else {
	  eigenTake<-as.matrix(eigenResults$vectors[,onesIndex]) 
	  out<-eigenTake/colSums(eigenTake)
	  
  }
   #subset the eigenvectors
  #normalize
 
  return(out)
}

#
setGeneric("steadyStates", function(object) standardGeneric("steadyStates"))
setMethod("steadyStates","markovchain", 
          function(object) {
			transposeYN=FALSE
			if(object@byrow==TRUE) transposeYN=TRUE
            out<-.mcEigen(matr=object@transitionMatrix, transpose=transposeYN) #wrapper for .mcEigen
			if(transposeYN==TRUE) colnames(out)<-object@states else rownames(out)<-object@states
            #if(nrow(out)==1) out<-as.numeric(out)
            return(out)
          }
)


#generic function to extract absorbing states

setGeneric("absorbingStates", function(object) standardGeneric("absorbingStates"))
setMethod("absorbingStates","markovchain", 
          function(object) {
          out<-character()
          matr<-object@transitionMatrix #extract the byrow transition matrix
		  transposeYN=FALSE
		  if(object@byrow==TRUE) transposeYN=TRUE
          steadyStates<-.mcEigen(matr, transpose=transposeYN) #checkk
          if(transposeYN==TRUE) maxCols<-apply(steadyStates, 2, "max") else maxCols<-apply(steadyStates, 1, "max")  
          index<-which(maxCols==1)
          if(length(index)>0) out<-object@states[index]
            return(out)
          }
)


#generic function to extract transition probability
setGeneric("transitionProbability", function(object, t0, t1) standardGeneric("transitionProbability"))
setMethod("transitionProbability","markovchain", 
	function(object, t0, t1) {
		out<-numeric(1)
		fromState=which(object@states==t0)
		toState=which(object@states==t1)
		if(object@byrow==TRUE) out<-object@transitionMatrix[fromState, toState] else out<-object@transitionMatrix[toState, fromState] 
		return(out)
	}
)




#print plot show methods

.showInt<-function(object, verbose=TRUE)
{
	if(object@byrow==TRUE) direction="(by rows)" else direction="(by cols)"
	if(verbose==TRUE) cat("A",dim(object),"- dimensional discrete Markov Chain with following states \n",states(object), "\n and following ", direction," transition probabilities \n")
	print(object@transitionMatrix)
	cat("\n")
}


#show method 
setMethod("show","markovchain", #metodo show
          function(object){
          .showInt(object)
          }
)

setMethod("print","markovchain", #metodo print
          function(x){
           object<-x
		   .showInt(object,verbose=FALSE)
          }
)

.getNet<-function(object)
{
  require(igraph)
  net<-graph.adjacency(adjmatrix=object@transitionMatrix*10,
                       weighted=TRUE,mode="directed")
  return(net)
}


setGeneric("plotMc", function(object,...) standardGeneric("plotMc"))
setMethod("plotMc","markovchain", #metodo print
          function(object,...){
            netMc<-.getNet(object)
            plot.igraph(x=netMc, ...)
          }
)

###converting from and to df.

.mc2Df<-function(from)
{
	nr<-nrow(from@transitionMatrix)
	for(i in 1:nr){
		for(j in 1:nr){
				t0<-from@states[i]
				t1<-from@states[j]
				prob<-transitionProbability(object=from, t0=t0, t1=t1)
				rowDf<-data.frame(t0=t0, t1=t1, prob=prob)
							if(exists("outDf")) outDf=rbind(outDf, rowDf) else outDf=rowDf
			}
	}
	return(outDf)
}

setAs(from="markovchain", to="data.frame", def=.mc2Df)

.whichColProb<-function(df)
{
	out=0
	if(ncol(df)>3) warning("Warning! More than three column. Only the first three will be used")
	if(ncol(df)<3) stop("Error! Three column needed")
	
	for(i in 1:ncol(df))
		{
			if((class(df[,i])=="numeric")&(all(sapply(df[,i], .isProb)==TRUE))) #when found the first numeric and probability col
				{
					out=i
					break
				}
		}
	return(out)
}

# .df2Mc<-function(df)
# {
	# statesNames=unique(df[,1])
	# colProb<-.whichColProb(df)
	

# }

#aritmethics


#transposing method for markov chain

setMethod("t", "markovchain", 
	function(x) {
		out=new("markovchain", byrow=!x@byrow, transitionMatrix=t(x@transitionMatrix))
		return(out)
	} 
)





#multiplicationMethod
#by a markov chain (like a power by 2)
setMethod("*", c("markovchain", "markovchain"),
function(e1, e2) {
	 if(!setequal(e1@states,e1@states)) warning("Warning! Different states")
	 if(!setequal(dim(e1@transitionMatrix),dim(e2@transitionMatrix))) stop("Error! Different size")
	 newStates<-e1@states
	 newTransMatr<-e1@transitionMatrix%*%e2@transitionMatrix
	 out<-new("markovchain", states=newStates, transitionMatrix=newTransMatr)
	 return(out)
}
)

setMethod("*", c("matrix", "markovchain"),
		function(e1, e2) {
			out<-e1%*%e2@transitionMatrix
			return(out)
		}
)

setMethod("*", c("markovchain","matrix"),
		function(e1, e2) {
			out<-e1@transitionMatrix%*%e2
			return(out)
		}
)

setMethod("*", c("numeric","markovchain"),
		function(e1, e2) {
			if(length(e1)!=dim(e2)) stop("Error! Uncompatible dimensions") else out<-e1%*%e2@transitionMatrix
			return(out)
		}
)

setMethod("*", c("markovchain","numeric"),
		function(e1, e2) {
			if(length(e2)!=dim(e1)) stop("Error! Uncompatible dimensions") else out<-e1@transitionMatrix%*%e2
			return(out)
		}
)

#power method

# setMethod("^", c("markovchain", "numeric"),
# function(e1, e2) {
	 # if(e2==1) return(e1)
	 # newStates<-e1@states
	 # if(e2<1) stop("Error. Power must be greater or equal than one")
	 # newTransMatr<-e1@transitionMatrix
	 # for(i in 1:(e2-1))
		# {
			# newTransMatr<-newTransMatr%*%e1@transitionMatrix
		# }
	 # out<-new("markovchain", states=newStates, transitionMatrix=newTransMatr)
	 # return(out)
# }
# )


# library(microbenchmark)


setMethod("^", c("markovchain", "numeric"),
function(e1, e2) {
	require(expm)
	 out<-new("markovchain", states=e1@states, transitionMatrix=e1@transitionMatrix%^%e2)
	 return(out)
}
)
