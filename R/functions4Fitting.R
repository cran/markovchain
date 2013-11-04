#sampler for univariate markov chains
markovchainSequence<-function(n,markovchain, t0=sample(markovchain@states,1),include.t0=FALSE)
{
  if(!(t0%in%markovchain@states)) stop("Error! Initial state not defined")
  chain=character()
  state=t0
  for(i in 1:n) {
    rowProbs<-markovchain@transitionMatrix[which(markovchain@states==state),]
    outstate<-sample(size=1, x=markovchain@states, prob=rowProbs)
    chain=c(chain, outstate)
    state=outstate
  }
  if(include.t0) out<-c(t0, chain) else out<-chain
  return(out)
}
################
#random sampler#
################

#check if the subsequent states are included in the previous ones

.checkSequence<-function(object)
{
  out<-TRUE
  if(dim(object)==1) return(out) #if the size of the list is one do
  for(i in 2:dim(object))
  {
    statesNm1<-states(object[[i-1]]) #evalutate mc n.1
    statesN<-states(object[[i]]) #evaluate mc n
    intersection<-intersect(statesNm1,statesN) #check the ibntersection
    if(setequal(intersection, statesNm1)==FALSE) { #the states at n-1 
      out<-FALSE
      break
    }
  }
  return(out)
}


#function to perform random sampling
rmarkovchain<-function(n,object,...)
{
  if(class(object)=="markovchain") out<-markovchainSequence(n=n, markovchain=object,...)
  if(class(object)=="markovchainList")
  {
    verify<-.checkSequence(object=object)
    if(!verify) warning("Warning: some states in the markovchain sequences are not contained in the following states!")
    iteration<-numeric()
    values<-character()
    for(i in 1:n) #number of replicates
    {
      #the first iteration may include initial state
      sampledValues<-markovchainSequence(n=1,markovchain=object[[1]],...)
      outIter<-rep(i,length(sampledValues))
      if(dim(object)>1)
        {for(j in 2:dim(object))
        {
          pos2take<-length(sampledValues)
          newVals<-markovchainSequence(n=1,markovchain=object[[j]],t0=sampledValues[pos2take]) #the initial state is in the ending position of the mclist
          outIter<-c(outIter,i)
          sampledValues<-c(sampledValues,newVals)
        }
      }
      iteration<-c(iteration, outIter)
      values<-c(values, sampledValues)
    }
    out<-data.frame(iteration=iteration, values=values)
  }
  return(out)
}

#core function to get sequence matrix

.createSequenceMatrix<-function(stringchar, toRowProbs=FALSE)
{
  elements<-sort(unique(stringchar))
  freqMatrix<-zeros(length(elements))
  rownames(freqMatrix)<-elements
  colnames(freqMatrix)<-elements
  for(i in 1:(length(stringchar)-1))
  {
    posFrom<-which(rownames(freqMatrix)==stringchar[i])
    posTo<-which(rownames(freqMatrix)==stringchar[i+1])
    freqMatrix[posFrom,posTo]=freqMatrix[posFrom,posTo]+1
  }
  if(toRowProbs==TRUE)
  {
    freqMatrix<-freqMatrix/rowSums(freqMatrix)
  }
  return(freqMatrix)
}

#functon to fit a Markov chain by MLE

.mcFitMle<-function(stringchar,byrow)
{
  initialMatr<-.createSequenceMatrix(stringchar=stringchar,toRowProbs=TRUE)
  outMc<-new("markovchain", transitionMatrix=initialMatr,name="MLE Fit")
  if(byrow==FALSE) outMc<-t(outMc)
  out<-list(estimate=outMc)
  return(out)
}


#example

#data=markovchainSequence(10000,markovA,t0="a")
#ciao<-markovchainFit(data=data)


#given a sting of characters, returns the associate one step transition matrix
.bootstrapCharacterSequences<-function(stringchar, n, size=length(stringchar))
{
  contingencyMatrix<-.createSequenceMatrix(stringchar=stringchar)
  samples<-list()
  itemset<-rownames(contingencyMatrix)
  for(i in 1:n) #cicle to fill the samples
  {
    charseq<-character()
    char<-sample(x=itemset,size=1)
    charseq<-c(charseq,char)
    for(j in 2:size) #cicle to define the element in a list
    {
      probsVector<-contingencyMatrix[which(rownames(contingencyMatrix)==char),]
      char<-sample(x=itemset,size=1, replace=TRUE,prob=probsVector)
      charseq<-c(charseq,char)
    }
    samples[[length(samples)+1]]<-charseq #increase the list
  }
  return(samples)
}


.fromBoot2Estimate<-function(listMatr)
{
  sampleSize<-length(listMatr)
  matrDim<-nrow(listMatr[[1]])
  #create the estimate output
  matrMean<-zeros(matrDim)
  matrSd<-zeros(matrDim)
  #create the sample output
  for(i in 1:matrDim) #move by row
  {
    for(j in 1:matrDim) #move by cols
    {
      probsEstimated<-numeric()
      #fill the probs
      for(k in 1:sampleSize) probsEstimated<-c(probsEstimated,listMatr[[k]][i,j])
      muCell<-mean(probsEstimated)
      sdCell<-sd(probsEstimated)
      matrMean[i,j]<-muCell
      matrSd[i,j]<-sdCell
    }
  }
out<-list(estMu=matrMean, estSigma=matrSd)
  return(out)
}


.mcFitBootStrap<-function(data, nboot=10,byrow=TRUE)
{
  #create the list of bootstrap sequence sample
  theList<-.bootstrapCharacterSequences(stringchar=data, n=nboot)
  #convert the list in a probability matrix
  pmsBootStrapped<-lapply(X=theList, FUN=.createSequenceMatrix, toRowProbs=TRUE)
  estimateList<-.fromBoot2Estimate(listMatr=pmsBootStrapped)
  #from raw to estimate
  temp<-estimateList$estMu
  transMatr<-sweep(temp, 1, rowSums(temp), FUN="/")
  estimate<-new("markovchain",transitionMatrix=transMatr, byrow=byrow, name="BootStrap Estimate")
  out<-list(estimate=estimate, standardError=estimateList$estSigma,bootStrapSamples=pmsBootStrapped)
  return(out)
}

#fit

markovchainFit<-function(data,method="mle", byrow=TRUE,nboot=10)
{
  #MLE FIT
  if(method=="mle") out<-.mcFitMle(stringchar=data,byrow=byrow)
  if(method=="bootstrap") out<-.mcFitBootStrap(data=data,byrow=byrow,nboot=nboot)
  return(out)
}




