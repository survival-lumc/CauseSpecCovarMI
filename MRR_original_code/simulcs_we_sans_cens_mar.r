rm(list=ls())

require(survival)
require(cmprsk) # not available for R 4.1
require(MASS)
require(mice)

source("tiragecs_weib.r")

repet<-100
set.seed(17042000)

ref.1<-ref.2<-NULL
cc.1<-cc.2<-NULL
impc.1<-impc.2<-NULL
imp1.1<-imp1.2<-NULL
imp2.1<-imp2.2<-NULL

beta1X<- 1
beta2X<- -0.5
beta1Z<- 0.5 # -0.5 change
beta2Z<- 0.5
reference<-c(beta1X,beta2X,beta1Z,beta2Z)
names(reference)<-c("beta1X","beta2X","beta1Z","beta2Z")

maxx<-1

tmp_temps<-Sys.time()

for (i in 1:repet) {
  print(paste(i,round(difftime(Sys.time(),tmp_temps,units="min")*(repet-i++1)/i,2)))
  # Tirage
  baz<-tiragecs_we(n=200,h10=1,h20=0.5,a1=0.3,a2=1.7,b1X=beta1X,b2X=beta2X,
                   b1Z=beta1Z,b2Z=beta2Z,
                fXZ=function(nn){
                  tmp<-mvrnorm(n=nn,mu=c(1,1),Sigma=matrix(c(1,0,0,1),nrow=2))
                  data.frame(X=tmp[,1],Z=tmp[,2])
                #},censure=list(tirc=function(nn,bb){runif(nn,0,bb)},varc=1))
                },
                miss=function(nn,zz){
                  tmp.nn<-table(zz>0)
                  ifelse(zz>0,rbinom(tmp.nn[2],1,0.1),rbinom(tmp.nn[1],1,0.9))},
                censure=list(tirc=function(nn,hc){rexp(nn,hc)},varc=0.0000001))
  
  baz$ee<-as.factor(baz$ee)
  baz$ee1<-(baz$ee=="1")*1
  baz$ee2<-(baz$ee=="2")*1
  baz$cens<-(baz$ee1+baz$ee2==0)*1

  #baz$XX1<-baz$XX*as.numeric(as.character(baz$ee1))
  #baz$XX2<-baz$XX*as.numeric(as.character(baz$ee2))
  
  # Calcul des cumulative hazard
  #tmp.surv<-summary(survfit(Surv(tt,I(ee=="1"))~1,type="fleming-harrington",data=baz),
  #                  times=baz$tt,censored=TRUE,extend=T)
  #tmp.ind<-match(baz$tt,tmp.surv$time)
  baz$H01<-nelsonaalen(baz,tt,ee1)
  
  ## v?rification
  #verif_tmp<-nelsonaalen(baz,tt,ee1)
  #cbind(verif_tmp,baz$H01)
  
#   tmp.surv<-summary(survfit(Surv(tt,I(ee=="2"))~1,type="fleming-harrington",data=baz),
#                     times=baz$tt,censored=TRUE,extend=T)
#   tmp.ind<-match(baz$tt,tmp.surv$time)
  baz$H02<-nelsonaalen(baz,tt,ee2)

  baz$H0cens<-nelsonaalen(baz,tt,cens)

  baz$H01Z<-baz$H01*baz$Z
  baz$H02Z<-baz$H02*baz$Z
  baz$H0censZ<-baz$H0cens*baz$Z

  # Analyse reference
  res.ref.1<-summary(coxph(Surv(tt,ee1)~X+Z,data=baz))
  #res.ref.2<-summary(coxph(Surv(tt,I(ee=="2"))~X+Z,data=baz))
  ref.1<-rbind(ref.1,c(res.ref.1$coefficients["X",c("coef","se(coef)","Pr(>|z|)")],
                       log(res.ref.1$conf.int["X",c("lower .95","upper .95")]),
                       res.ref.1$coefficients["Z",c("coef","se(coef)","Pr(>|z|)")],
                       log(res.ref.1$conf.int["Z",c("lower .95","upper .95")])))
#   ref.2<-rbind(ref.2,c(res.ref.2$coefficients["X",c("coef","se(coef)","Pr(>|z|)")],
#                        log(res.ref.2$conf.int["X",c("lower .95","upper .95")]),
#                        res.ref.2$coefficients["Z",c("coef","se(coef)","Pr(>|z|)")],
#                       log(res.ref.2$conf.int["Z",c("lower .95","upper .95")])))
  
  # Analyse complete case
  res.cc.1<-summary(coxph(Surv(tt,ee1)~Xm+Z,data=baz))
  #res.cc.2<-summary(coxph(Surv(tt,I(ee=="2"))~Xm+Z,data=baz))
  cc.1<-rbind(cc.1,c(res.cc.1$coefficients["Xm",c("coef","se(coef)","Pr(>|z|)")],
                     log(res.cc.1$conf.int["Xm",c("lower .95","upper .95")]),
                       res.cc.1$coefficients["Z",c("coef","se(coef)","Pr(>|z|)")],
                         log(res.cc.1$conf.int["Z",c("lower .95","upper .95")])))
#   cc.2<-rbind(cc.2,c(res.cc.2$coefficients["Xm",c("coef","se(coef)","Pr(>|z|)")],
#                      log(res.cc.2$conf.int["Xm",c("lower .95","upper .95")]),
#                        res.cc.2$coefficients["Z",c("coef","se(coef)","Pr(>|z|)")],
#                          log(res.cc.2$conf.int["Z",c("lower .95","upper .95")])))
#   
  # Imputation avec H01 et H02 new method CH12
  # basic 1
  # avec interaction 2
  tmp.imp<-mice(baz,m=5,maxit=0,print=FALSE)
  tmp.imp$pred[,c("X","t1","t2","tc","tt","missind","ee1","ee2","H02Z","H01Z")]<-0
  tmp.imp<-mice(baz,m=10,maxit=maxx,predictorMatrix=tmp.imp$pred,print=FALSE)
  tmp.ana.imp<-with(tmp.imp,expr=coxph(Surv(tt,ee1)~Xm+Z))
  res.impc.1<-summary(pool(tmp.ana.imp), conf.int = T)

  tmp.imp<-mice(baz,m=5,maxit=0,print=FALSE)
  tmp.imp$pred[,c("X","t1","t2","tc","tt","missind","ee1","ee2")]<-0
  tmp.imp<-mice(baz,m=10,maxit=maxx,predictorMatrix=tmp.imp$pred,print=FALSE)
  tmp.ana.imp<-with(tmp.imp,expr=coxph(Surv(tt,ee1)~Xm+Z))
  res.impc.2<-summary(pool(tmp.ana.imp), conf.int = T)

  impc.1<-rbind(impc.1,c(res.impc.1["Xm",c("estimate","std.error","p.value","2.5 %","97.5 %")],
                         res.impc.1["Z",c("estimate","std.error","p.value","2.5 %","97.5 %")]))
  impc.2<-rbind(impc.2,c(res.impc.2["Xm",c("estimate","std.error","p.value","2.5 %","97.5 %")],
                         res.impc.2["Z",c("estimate","std.error","p.value","2.5 %","97.5 %")]))
  
  #Imputation avec que H01
  tmp.imp<-mice(baz,m=5,maxit=0,print=FALSE)
  tmp.imp$pred[,c("X","t1","t2","tc","tt","missind","H02","H02Z","ee","ee2","H01Z")]<-0
  tmp.imp$pred["Xm","ee1"]<-1
  tmp.imp<-mice(baz,m=10,maxit=1,predictorMatrix=tmp.imp$pred,print=FALSE)
  tmp.ana.imp<-with(tmp.imp,expr=coxph(Surv(tt,ee1)~Xm+Z))
  res.imp1.1<-summary(pool(tmp.ana.imp),  conf.int = T)
#   tmp.ana.imp<-with(tmp.imp,expr=coxph(Surv(tt,I(ee=="2"))~Xm+Z))
#   res.imp1.2<-summary(pool(tmp.ana.imp))
#   tmp.imp<-mice(baz,m=5,maxit=0,print=FALSE)
#   tmp.imp$pred[,c("X","t1","t2","tc","tt","missind","H01","H01Z","ee","ee1","H02Z","H01Z")]<-0
#   tmp.imp<-mice(baz,m=10,maxit=1,predictorMatrix=tmp.imp$pred,print=FALSE)
#   tmp.ana.imp<-with(tmp.imp,expr=coxph(Surv(tt,I(ee=="1"))~Xm+Z))
#   res.imp2.1<-summary(pool(tmp.ana.imp))
#   tmp.ana.imp<-with(tmp.imp,expr=coxph(Surv(tt,I(ee=="2"))~Xm+Z))
#   res.imp2.2<-summary(pool(tmp.ana.imp))
  imp1.1<-rbind(imp1.1,c(res.imp1.1["Xm",c("estimate","std.error","p.value","2.5 %","97.5 %")],
                         res.imp1.1["Z",c("estimate","std.error","p.value","2.5 %","97.5 %")]))
#   imp1.2<-rbind(imp1.2,c(res.imp1.2["Xm",c("est","se","Pr(>|t|)","lo 95","hi 95")],
#                          res.imp1.2["Z",c("est","se","Pr(>|t|)","lo 95","hi 95")]))
#   imp2.1<-rbind(imp2.1,c(res.imp2.1["Xm",c("est","se","Pr(>|t|)","lo 95","hi 95")],
#                          res.imp2.1["Z",c("est","se","Pr(>|t|)","lo 95","hi 95")]))
#   imp2.2<-rbind(imp2.2,c(res.imp2.2["Xm",c("est","se","Pr(>|t|)","lo 95","hi 95")],
#                          res.imp2.2["Z",c("est","se","Pr(>|t|)","lo 95","hi 95")]))
  
}

ref.1<-data.frame(ref.1)
#ref.2<-data.frame(ref.2)
cc.1<-data.frame(cc.1)
#cc.2<-data.frame(cc.2)
impc.1<-data.frame(impc.1)
impc.2<-data.frame(impc.2)
imp1.1<-data.frame(imp1.1)
#imp1.2<-data.frame(imp1.2)
#imp2.1<-data.frame(imp2.1)
#imp2.2<-data.frame(imp2.2)

names(ref.1)<-c("X.est","X.se","X.p","X.lo","X.hi","Z.est","Z.se","Z.p","Z.lo","Z.hi")
names(cc.1)<-names(impc.1)<-names(impc.2)<-names(imp1.1)<-names(ref.1)

write.table(ref.1,file=paste(file_result,"/ref.1.txt",sep=""),sep="\t",row.names=FALSE)
#write.table(ref.2,file=paste(file_result,"/ref.2.txt",sep=""),sep="\t",row.names=FALSE)
write.table(cc.1,file=paste(file_result,"/cc.1.txt",sep=""),sep="\t",row.names=FALSE)
#write.table(cc.2,file=paste(file_result,"/cc.2.txt",sep=""),sep="\t",row.names=FALSE)
write.table(impc.1,file=paste(file_result,"/impc.1.txt",sep=""),sep="\t",row.names=FALSE)
write.table(impc.2,file=paste(file_result,"/impc.2.txt",sep=""),sep="\t",row.names=FALSE)
write.table(imp1.1,file=paste(file_result,"/imp1.1.txt",sep=""),sep="\t",row.names=FALSE)
#write.table(imp1.2,file=paste(file_result,"/imp1.2.txt",sep=""),sep="\t",row.names=FALSE)
#write.table(imp2.1,file=paste(file_result,"/imp2.1.txt",sep=""),sep="\t",row.names=FALSE)
#write.table(imp2.2,file=paste(file_result,"/imp2.2.txt",sep=""),sep="\t",row.names=FALSE)

write.table(reference,file=paste(file_result,"/reference.txt",sep=""),sep="\t",row.names=TRUE)

















