
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

#functon to fit a Markov chain

.mcFitMle<-function(stringchar,byrow)
{
  states<-unique(stringchar)
  initialMatr=zeros(length(states))
  rownames(initialMatr)<-states
  colnames(initialMatr)<-states
  for(i in 1:(length(stringchar)-1))
  {
    state1<-stringchar[i];rowIndex<-which(states==state1)
    state2<-stringchar[i+1];colIndex<-which(states==state2)
    initialMatr[rowIndex, colIndex]<-initialMatr[rowIndex, colIndex]+1
  }
  initialMatr<-initialMatr/rowSums(initialMatr)
  outMc<-new("markovchain", transitionMatrix=initialMatr)
  if(byrow==FALSE) outMc<-t(outMc)
  out<-list(estimate=outMc)
  return(out)
}

markovchainFit<-function(data,method="mle", byrow=TRUE)
{
	#frequency fit
	if(method=="mle") out<-.mcFitMle(stringchar=data,byrow=byrow)
  return(out)
}

#example

#data=markovchainSequence(10000,markovA,t0="a")
#ciao<-markovchainFit(data=data)