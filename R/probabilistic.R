

.commclassesKernel <- function(P){
	m=ncol(P)
	stateNames<-rownames(P)

	T=zeros(m) #crea una matrice di zeri
	i=1
	while (i<=m) { 
		a=i # differenza con MatLab: a e b qui sono VETTORI
		b<-zeros(1,m)
		b[1,i]<-1
		old<-1
		new<-0
		while (old != new) {
			old=sum(find(b>0))
			#[ignore,n]=size(a);
			#ignore=1
			n=size(a)[2]
			#c=sum(P(a,:),1);
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

