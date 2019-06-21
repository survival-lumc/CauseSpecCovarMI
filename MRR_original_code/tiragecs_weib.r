tiragecs_we<-function(n=100,h10=1,h20=1,a1=1,a2=1,b1X=1,b2X=1,b1Z=1,b2Z=1,
fXZ=function(nn){data.frame(X=rbinom(nn,1,0.5),Z=rnorm(nn,mean=0,sd=1))},
censure=list(tirc=function(nn,hc){rexp(nn,hc)},varc=1),
miss=function(nn,zz){rbinom(nn,1,0.5)})
{
### Tirage covariable
baz<-fXZ(n)

### Tirage temps
baz$t1<-rweibull(n,shape=a1,scale=h10*exp(-(b1X*baz$X+b1Z*baz$Z)/a1))
baz$t2<-rweibull(n,shape=a2,scale=h20*exp(-(b2X*baz$X+b2Z*baz$Z)/a2))

### Tirage Censure
baz$tc<-censure$tirc(n,censure$varc)

### Formalisation des donn?es de survie
baz$tt<-pmin(baz$t1,baz$t2,baz$tc)
baz$ee<-rep(1,n)
baz$ee[baz$tt==baz$t2]<-2
baz$ee[baz$tt==baz$tc]<-0

### M?canime des donn?es manquantes
baz$missind<-miss(n,zz=baz$Z)
baz$Xm<-baz$X
baz$Xm[baz$missind==1]<-NA

return(baz)
}


