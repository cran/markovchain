

.commclassesKernel <- function(P){
	m=ncol(P)
	stateNames<-rownames(P)
	T=zeros(m) 
	i=1
	while (i<=m) { 
		a=i 
		b<-zeros(1,m)
		b[1,i]<-1
		old<-1
		new<-0
		while (old != new) {
			old=sum(find(b>0))
			n=size(a)[2]
			matr<-matrix(as.numeric(P[a,]),ncol=m,nrow=n) #fix
			c=colSums(matr)
			d=find(c)
			n=size(d)[2]
			b[1,d]<-ones(1,n)
			new<-sum(find(b>0))
			a<-d
		}
		T[i,]=b
		i=i+1 }
	F=t(T)  
	C=(T>0)&(F>0)
	v=(apply(t(C)==t(T),2,sum)==m)
	colnames(C)=stateNames
	rownames(C)=stateNames
	names(v)=stateNames
	out<-list(C=C,v=v)
	return(out)
}

#returns the underlying communicating classes
.communicatingClasses<-function(adjMatr)
{
  len<-dim(adjMatr)[1]
  classesList<-list()
  for(i in 1:len)
  {
    row2Check<-adjMatr[i,]
    proposedCommClass<-names(which(row2Check==TRUE))
    if(i>1) 
    {
      for(j in 1:(length(classesList)))
      {
        check<-FALSE
        check<-setequal(classesList[[j]],proposedCommClass)
        if(check==TRUE) {proposedCommClass<-NULL; break}
      }
    }
    if(!is.null(proposedCommClass)) classesList[[length(classesList)+1]]<-proposedCommClass
  }
  return(classesList)
}

#communicating states
.commStatesFinder<-function(matr)
{
  #Reachability matrix
  dimMatr<-dim(matr)[1]
  Z<-sign(matr)
  temp<-(eye(dimMatr)+Z)%^%(dimMatr-1)
  R<-sign(temp)
  return(R)
}

is.accessible<-function(object, from, to)
{
  out<-FALSE
  statesNames<-states(object)
  fromPos<-which(statesNames==from)
  toPos<-which(statesNames==to)
  R<-.commStatesFinder(object@transitionMatrix)
  if(R[fromPos,toPos]==TRUE) out<-TRUE
  return(out)
}

#a markov chain is irreducible if is composed by only one communicating class
is.irreducible<-function(object)
{
  out<-FALSE
  tocheck<-.communicatingClasses(.commclassesKernel(object@transitionMatrix)$C)
  if(length(tocheck)==1) out<-TRUE
  return(out)
}

.summaryKernel<-function(object)
{
  matr<-object@transitionMatrix
  temp<-.commclassesKernel(matr)
  communicatingClassList<-.communicatingClasses(temp$C)
  transientStates<-names(which(temp$v==FALSE))
  closedClasses<-list()
  transientClasses<-list()
  for(i in 1:length(communicatingClassList))
  {
    class2Test<-communicatingClassList[[i]]
    if(length(intersect(class2Test,transientStates))>0) transientClasses[[length(transientClasses)+1]]<-class2Test else closedClasses[[length(closedClasses)+1]]<-class2Test
  }
  summaryMc<-list(closedClasses=closedClasses, 
                  transientClasses=transientClasses)
  return(summaryMc)
}


#qui la funzione firstpassage
.firstpassageKernel<-function(P,i,n){
  G<-P
  H<-P[i,]
  E<-1-diag(size(P)[2])
  for (m in 2:n) {
    G<-P%*%(G*E)
    H<-rbind(H,G[i,])
  }
  return(H)
}

firstPassage<-function(object,state,n)
{
  P<-object@transitionMatrix
  stateNames<-states(object)
  i<-which(stateNames==state)
  outMatr<-.firstpassageKernel(P=P,i=i,n=n)
  colnames(outMatr)<-stateNames
  rownames(outMatr)<-1:n
  return(outMatr)
}

