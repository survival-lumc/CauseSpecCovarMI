##******************************##
## smcfcs with time-fixed coxph ##
##******************************##


# Cite Bartlett here

# Simply a version from https://github.com/jwb133/smcfcs , with all
# cox models adding the control = coxph.control(timefix = FALSE)

smcfcs_timefix <- function(originaldata,
                           smtype,
                           smformula,
                           method,
                           predictorMatrix = NULL,
                           m=5,
                           numit=10,
                           rjlimit=1000,
                           noisy=FALSE,
                           errorProneMatrix=NULL) {
  
  #call core  smcfcs function, passing through arguments
  smcfcs.core_timefix(originaldata,smtype,smformula,method,predictorMatrix,
              m,numit,rjlimit,noisy,errorProneMatrix=errorProneMatrix)
}



#this is the core of the smcfcs function, called by wrapper functions for certain different substantive models
smcfcs.core_timefix <- function(originaldata,
                                smtype,
                                smformula,
                                method,
                                predictorMatrix=NULL,
                                m=5,
                                numit=10,
                                rjlimit=1000,
                                noisy=FALSE,
                                errorProneMatrix=NULL, ...) {
  
  #get extra arguments passed in ...
  extraArgs <- list(...)
  
  stopifnot(is.data.frame(originaldata))
  if (ncol(originaldata)!=length(method)) stop("Method argument must have the same length as the number of columns in the data frame.")
  
  n <- dim(originaldata)[1]
  
  #create matrix of response indicators
  r <- 1*(is.na(originaldata)==0)
  
  if ((smtype %in% c("lm", "logistic", "poisson", "coxph", "compet", "casecohort","nestedcc", "weibull"))==FALSE)
    stop(paste("Substantive model type ",smtype," not recognised.",sep=""))
  
  #find column numbers of partially observed, fully observed variables, and outcome
  if (smtype=="coxph") {
    
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    d <- originaldata[,dCol]
    
    nullMod <- survival::coxph(Surv(originaldata[,timeCol],originaldata[,dCol])~1,
                               control = survival::coxph.control(timefix = FALSE))
    basehaz <- survival::basehaz(nullMod)
    H0indices <- match(originaldata[,timeCol], basehaz[,2])
    rm(nullMod)
  } else if (smtype=="weibull") {
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    d <- originaldata[,dCol]
  } else if (smtype=="compet") {
    
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula[[1]])[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula[[1]])[[2]][[3]][[2]])]
    outcomeCol <- c(timeCol, dCol)
    d <- originaldata[,dCol]
    numCauses <- length(smformula)
    H0 <- vector("list", numCauses)
    H0indices <- vector("list", numCauses)
    outcomeModBeta <- vector("list", numCauses)
    linpred <- vector("list", numCauses)
    for (cause in 1:numCauses) {
      nullMod <- survival::coxph(as.formula(paste(strsplit(smformula[[cause]],"~")[[1]][1],"~1")), 
                                 originaldata,
                                 control = survival::coxph.control(timefix = FALSE))
      basehaz <- survival::basehaz(nullMod)
      H0[[cause]] <- basehaz[,1]
      H0indices[[cause]] <- match(originaldata[,timeCol], basehaz[,2])
      linpred[[cause]] <- as.formula(smformula[[cause]])
    }
    rm(nullMod)
  } else if (smtype=="casecohort") {
    
    subcoCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% extraArgs$in.subco]
    #subcoMembers is a vector of row numbers of those in the subcohort
    subcoMembers <- which(originaldata[,subcoCol]==1)
    
    #generate weights for use in later analysis which we use to obtain baseline cumulative hazard
    #assign a weight of /samp.frac to individuals in the subcohort and 0 to those outside the subcohort
    subco.weight<-ifelse(originaldata[,subcoCol]==1,1/extraArgs$sampfrac,0)
    
    
    entertimeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[4]])]
    outcomeCol <- c(entertimeCol,timeCol, dCol)
    d <- originaldata[,dCol]
    
    #list of unique event times - used in calculation of baseline cumulative hazard
    list.times=sort(unique(originaldata[,timeCol][originaldata[,dCol]==1]))  #RUTH 21/03/17: ADDED THIS LINE
    
    smcfcsid <- 1:n
    smformula2<-paste(smformula,"+cluster(smcfcsid)",sep="")
  }
  else if (smtype=="nestedcc") {
    
    timeCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[2]])]
    dCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(as.formula(smformula)[[2]][[3]])]
    outcomeCol <- c(timeCol, dCol)
    setCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(extraArgs$set)]
    nriskCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(extraArgs$nrisk)]
    eventCol <- (1:dim(originaldata)[2])[colnames(originaldata) %in% toString(extraArgs$event)] #this is distinct from the dCol
    
    #the below command creates "Surv(t,case)~x" from "Surv(t,case)~x+strata(setno)" (for example) (i.e. it removes the strata part of the formula)
    #this is used when obtaining outmodxb
    #note this is done slightly oddly, but this is because as.formula does not work well for long formulas as it splits across lines
    exp1=as.formula(paste(smformula))[[2]]
    exp2=as.formula(smformula)[[3]][[2]]
    smformula2<-paste(deparse(exp1),"~",deparse(exp2,width.cutoff = 500L))
    
    #This is the indicator of whether an individual ever has the event (regardless of whether they are sometimes used as a control and sometimes (one) as a case)
    d <- originaldata[,eventCol]
    
    #noncases is a vector of row numbers of those who never have the event (which is a subset of the controls)
    noncases <- which(originaldata[,eventCol]==0)
    
    #number of individuals in each sampled risk set (matched set)
    num.sampriskset<-ave(rep(1,dim(originaldata)[1]), originaldata[,setCol], FUN = function(x) sum(x))
    
  }
  else {
    outcomeCol <- which(colnames(originaldata)==as.formula(smformula)[[2]])
  }
  
  if (smtype=="logistic") {
    if (is.numeric(originaldata[,outcomeCol])==FALSE) {
      stop("For logistic substantive models the outcome variable must be numeric 0/1.")
    } else {
      if (all.equal(unique(originaldata[,outcomeCol]),c(0,1))==FALSE) {
        stop("For logistic substantive models the outcome variable must be coded 0/1.")
      }
    }
  }
  
  if (smtype=="compet") {
    smcovnames <- attr(terms(as.formula(smformula[[1]])), "term.labels")
    for (cause in 2:numCauses) {
      smcovnames <- c(smcovnames, attr(terms(as.formula(smformula[[cause]])), "term.labels"))
    }
    smcovnames <- unique(smcovnames)
  } else  if (smtype=="nestedcc") {
    smcovnames <- attr(terms(as.formula(smformula)), "term.labels")[-length(attr(terms(as.formula(smformula)), "term.labels"))]
  }
  else {
    smcovnames <- attr(terms(as.formula(smformula)), "term.labels")
  }
  smcovcols <- (1:ncol(originaldata))[colnames(originaldata) %in% smcovnames]
  
  #partial vars are those variables for which an imputation method has been specified among the available regression types
  partialVars <- which((method=="norm") | (method=="latnorm") | (method=="logreg") | (method=="poisson") | (method=="podds") | (method=="mlogit"))
  
  if (length(partialVars)==0) stop("You have not specified any valid imputation methods in the method argument.")
  
  #check that methods are given for each partially observed column, and not given for fully observed columns
  for (colnum in 1:ncol(originaldata)) {
    if (method[colnum]!="") {
      #an imputation method has been specified
      if (colnum %in% outcomeCol) {
        stop(paste("An imputation method has been specified for ",colnames(originaldata)[colnum],
                   ". Elements of the method argument corresponding to the outcome variable(s) should be empty.",sep=""))
      }
      else {
        if (sum(r[,colnum])==n) {
          stop(paste("An imputation method has been specified for ",colnames(originaldata)[colnum],
                     ", but it appears to be fully observed.",sep=""))
        }
      }
    }
    else {
      #no imputation method has been specified
      if (sum(r[,colnum])<n) {
        #some values are missing
        if (((colnum %in% outcomeCol)==FALSE) & (sum(errorProneMatrix[,colnum])==0)) {
          stop(paste("Variable ",colnames(originaldata)[colnum], " does not have an imputation method specified, yet appears to have missing values.",sep=""))
        }
      }
    }
    
  }
  
  #check that if any variables have latnorm specified as method, that at least two error-prone measures
  #are specified for it in errorProneMatrix
  if ("latnorm" %in% method) {
    if (is.null(errorProneMatrix)==TRUE) {
      stop("If you specify method latnorm you must specify the errorProneMatrix argument.")
    } else {
      #check errorProneMatrix is of correct dimensional
      if (identical(dim(errorProneMatrix), c(length(originaldata),length(originaldata)))==FALSE) {
        stop("The errorProneMatrix should be a square matrix with number of rows equal to the number of variables in the dataset.")
      }
      #check entries of errorProneMatrix only consists of 0s and 1s
      if (identical(sort(unique(as.vector(errorProneMatrix))),c(0,1))==FALSE) {
        stop("The errorProneMatrix should only consist of 0s and 1s.")
      }
      #check each latnorm variable has at least 2 error-prones
      for (varNum in 1:length(method)) {
        if (method[varNum]=="latnorm") {
          if (sum(errorProneMatrix[varNum,])<2) {
            stop("Each latnorm variable must have two or more error prone measurements specified in the errorProneMatrix argument.")
          }
        }
      }
      #check no error-prone measurement is allocated to more than one latnorm
      if (sum(colSums(errorProneMatrix)>1)>0) {
        stop("Each error-prone measurement should be allocated to exactly one latnorm variable.")
      }
    }
  }
  if (is.null(errorProneMatrix)==FALSE) {
    if (("latnorm" %in% method)==FALSE) {
      stop("If you specify errorProneMatrix then at least one variable must be imputed using latnorm.")
    }
  }
  
  #fully observed vars are those that are fully observed and are covariates in the substantive model
  fullObsVars <- which((colSums(r)==n) & (colnames(originaldata) %in% smcovnames))
  
  #passive variables
  passiveVars <- which((method!="") & (method!="norm") & (method!="logreg") & (method!="poisson") & (method!="podds") & (method!="mlogit") & (method!="latnorm"))
  
  print(paste("Outcome variable(s):", paste(colnames(originaldata)[outcomeCol],collapse=',')))
  print(paste("Passive variables:", paste(colnames(originaldata)[passiveVars],collapse=',')))
  print(paste("Partially obs. variables:", paste(colnames(originaldata)[partialVars],collapse=',')))
  print(paste("Fully obs. substantive model variables:", paste(colnames(originaldata)[fullObsVars],collapse=',')))
  
  imputations <- list()
  for (imp in 1:m) {
    imputations[[imp]] <- originaldata
  }
  
  rjFailCount <- 0
  
  for (imp in 1:m) {
    
    print(paste("Imputation ",imp))
    
    #initial imputation of each partially observed variable based on observed values
    for (var in 1:length(partialVars)) {
      targetCol <- partialVars[var]
      if (method[targetCol]=="latnorm") {
        #first impute any missing replicate error-prone measurements of this variable by a randomly chosen observed value
        errorProneCols <- which(errorProneMatrix[targetCol,]==1)
        for (measure in 1:length(errorProneCols)) {
          if (sum(r[,errorProneCols[measure]])<n) {
            imputations[[imp]][r[,errorProneCols[measure]]==0,errorProneCols[measure]] <- sample(imputations[[imp]][r[,errorProneCols[measure]]==1,
                                                                                                                    errorProneCols[measure]], size=sum(r[,errorProneCols[measure]]==0), replace=TRUE)
          }
        }
        
        #initialize latent predictors with mean of their error-prone measurements
        imputations[[imp]][,targetCol] <- apply(imputations[[imp]][,errorProneCols], 1, mean)
      }
      else {
        imputations[[imp]][r[,targetCol]==0,targetCol] <- sample(imputations[[imp]][r[,targetCol]==1,targetCol], size=sum(r[,targetCol]==0), replace=TRUE)
      }
    }
    
    #initial imputations of missing outcomes, if present (using improper imputation)
    if ((smtype=="lm") | (smtype=="logistic") | (smtype=="poisson")) {
      if (sum(r[,outcomeCol])<n) {
        if (imp==1) {
          print("Imputing missing outcomes using specified substantive model.")
        }
        #update passive variable(s)
        imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
        
        imputationNeeded <- (1:n)[r[,outcomeCol]==0]
        #estimate parameters of substantive model
        if (smtype=="lm") {
          ymod <- stats::lm(as.formula(smformula),imputations[[imp]])
          beta <- ymod$coef
          sigmasq <- summary(ymod)$sigma^2
          #fill out missing values so that model.matrix works for all rows
          imputations[[imp]][imputationNeeded,outcomeCol] <- 0
          outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% beta
          imputations[[imp]][imputationNeeded,outcomeCol] <- rnorm(length(imputationNeeded),outmodxb[imputationNeeded], sigmasq^0.5)
        }
        else if (smtype=="logistic") {
          ymod <- glm(as.formula(smformula),family="binomial",imputations[[imp]])
          beta <- ymod$coef
          imputations[[imp]][imputationNeeded,outcomeCol] <- 0
          outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% beta
          prob <- expit(outmodxb[imputationNeeded])
          imputations[[imp]][imputationNeeded,outcomeCol] <- rbinom(length(imputationNeeded),1,prob)
        }
        else if (smtype=="poisson") {
          ymod <- glm(as.formula(smformula),family="poisson",imputations[[imp]])
          beta <- ymod$coef
          imputations[[imp]][imputationNeeded,outcomeCol] <- 0
          outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% beta
          imputations[[imp]][imputationNeeded,outcomeCol] <- rpois(length(imputationNeeded),exp(outmodxb[imputationNeeded]))
        }
      }
    }
    
    for (cyclenum in 1:numit) {
      
      if (noisy==TRUE) {
        print(paste("Iteration ",cyclenum))
      }
      #update passive variable(s)
      imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
      
      for (var in 1:length(partialVars)) {
        targetCol <- partialVars[var]
        if (is.null(predictorMatrix)) {
          predictorCols <- c(partialVars[! partialVars %in% targetCol], fullObsVars)
        }
        else {
          predictorCols <- which(predictorMatrix[targetCol,]==1)
          #ensure that user has not included outcome variable(s) here
          predictorCols <- predictorCols[! predictorCols %in% outcomeCol]
        }
        if ((imp==1) & (cyclenum==1)) {
          if (method[targetCol]=="latnorm") {
            print(paste("Imputing: ",colnames(imputations[[imp]])[targetCol]," using ",paste(colnames(imputations[[imp]])[c(predictorCols,which(errorProneMatrix[targetCol,]==1))],collapse=',')," plus outcome",collapse=','))
          }
          else {
            print(paste("Imputing: ",colnames(imputations[[imp]])[targetCol]," using ",paste(colnames(imputations[[imp]])[predictorCols],collapse=',')," plus outcome",collapse=','))
          }
        }
        if (length(predictorCols)>0) {
          xmodformula <- as.formula(paste(colnames(imputations[[imp]])[targetCol], "~", paste(colnames(imputations[[imp]])[predictorCols], collapse="+"),sep=""))
        }
        else {
          xmodformula <- as.formula(paste(colnames(imputations[[imp]])[targetCol], "~1",sep=""))
        }
        if (smtype=="casecohort") {
          xmoddata <- imputations[[imp]][subcoMembers,]
        } else if (smtype=="nestedcc"){
          xmoddata <- imputations[[imp]][noncases,]
        } else {
          xmoddata <- imputations[[imp]]
        }
        if (method[targetCol]=="norm") {
          #estimate parameters of covariate model
          xmod <- lm(xmodformula, data=xmoddata)
          #take draw from posterior of covariate model parameters
          beta <- xmod$coef
          sigmasq <- summary(xmod)$sigma^2
          newsigmasq <- (sigmasq*xmod$df) / rchisq(1,xmod$df)
          covariance <- (newsigmasq/sigmasq)*vcov(xmod)
          newbeta = beta + MASS::mvrnorm(1, mu=rep(0,ncol(covariance)), Sigma=covariance)
          #calculate fitted values
          if ((smtype=="casecohort")|(smtype=="nestedcc")) {
            xfitted <- model.matrix(xmodformula, data=imputations[[imp]]) %*% newbeta
          } else {
            xfitted <- model.matrix(xmod) %*% newbeta
          }
        } else if  (method[targetCol]=="latnorm") {
          #estimate parameters of covariate model
          xmod <- lm(xmodformula, data=xmoddata)
          #take draw from posterior of covariate model parameters
          beta <- xmod$coef
          sigmasq <- summary(xmod)$sigma^2
          #draw from sigmasq posterior based on proper inverse gamma prior for sigmasq
          #prior equivalent to 1 observation and guess of sigmasq=1
          newsigmasq <- 1/rgamma(1,shape=((n+1)/2), rate=((n*sigmasq+1)/2))
          covariance <- (newsigmasq/sigmasq)*vcov(xmod)
          newbeta = beta + MASS::mvrnorm(1, mu=rep(0,ncol(covariance)), Sigma=covariance)
          #calculate fitted values
          if ((smtype=="casecohort")|(smtype=="nestedcc")) {
            xfitted <- model.matrix(xmodformula, data=imputations[[imp]]) %*% newbeta
          } else {
            xfitted <- model.matrix(xmod) %*% newbeta
          }
          
          #estimate error variance and draw new value of error variance
          errorProneCols <- which(errorProneMatrix[targetCol,]==1)
          xmat <- matrix(imputations[[imp]][,targetCol], nrow=nrow(imputations[[imp]]), ncol=length(errorProneCols))
          uVec <- c(as.matrix(imputations[[imp]][,errorProneCols] - xmat))
          sigmausq <- mean(uVec^2)
          #take draw from posterior of error variance, using proper inverse gamma prior
          sum_ni <- n*length(errorProneCols)
          sigmausq <- 1/rgamma(1,shape=((sum_ni+1)/2), rate=((sum_ni*sigmausq+1)/2))
          
          #re-impute any originally missing error-prone measurements, based on classical error model assumption
          for (measure in 1:length(errorProneCols)) {
            nToImpute <- n-sum(r[,errorProneCols[measure]])
            if (nToImpute>0) {
              #then some values need imputing
              imputations[[imp]][r[,errorProneCols[measure]]==0,errorProneCols[measure]] <- imputations[[imp]][r[,errorProneCols[measure]]==0,targetCol] +
                rnorm(nToImpute, 0, sd=sqrt(sigmausq))
            }
          }
          
          #calculate conditional mean and variance of X|everything else except outcome
          wmean <- rowMeans(imputations[[imp]][,errorProneCols])
          lambda <- newsigmasq/(newsigmasq+sigmausq/length(errorProneCols))
          xfitted <- xfitted + lambda * (wmean - xfitted)
          newsigmasq <- rep(newsigmasq*(1-lambda), n)
          
        } else if (method[targetCol]=="logreg") {
          xmod <- glm(xmodformula, family="binomial",data=xmoddata)
          newbeta = modPostDraw(xmod)
          if ((smtype=="casecohort")|(smtype=="nestedcc")) {
            xfitted <- expit(model.matrix(xmodformula, data=imputations[[imp]]) %*% newbeta)
          } else {
            xfitted <- expit(model.matrix(xmod) %*% newbeta)
          }
        } else if (method[targetCol]=="poisson") {
          xmod <- glm(xmodformula, family="poisson", data=xmoddata)
          newbeta = modPostDraw(xmod)
          if ((smtype=="casecohort")|(smtype=="nestedcc")) {
            xfitted <- exp(model.matrix(xmodformula, data=imputations[[imp]]) %*% newbeta)
          } else {
            xfitted <- exp(model.matrix(xmod) %*% newbeta)
          }
        } else if (method[targetCol]=="podds") {
          if (is.ordered(imputations[[imp]][,targetCol])==FALSE) stop("Variables to be imputed using method podds must be stored as ordered factors.")
          xmod <- VGAM::vglm(xmodformula, VGAM::propodds, data=xmoddata)
          xmod.dummy <- VGAM::vglm(xmodformula, VGAM::propodds, data=imputations[[imp]])
          newbeta <- VGAM::coef(xmod) + MASS::mvrnorm(1, mu=rep(0,ncol(VGAM::vcov(xmod))), Sigma=VGAM::vcov(xmod))
          linpreds <- matrix((VGAM::model.matrix(xmod.dummy)) %*% newbeta, byrow=TRUE, ncol=(nlevels(imputations[[imp]][,targetCol])-1))
          cumprobs <- cbind(1/(1+exp(linpreds)), rep(1,nrow(linpreds)))
          xfitted <- cbind(cumprobs[,1] ,cumprobs[,2:ncol(cumprobs)] - cumprobs[,1:(ncol(cumprobs)-1)])
        } else if (method[targetCol]=="mlogit") {
          if (is.factor(imputations[[imp]][,targetCol])==FALSE) stop("Variables to be imputed using method modds must be stored as factors.")
          xmod <- VGAM::vglm(xmodformula, VGAM::multinomial(refLevel=1), data=xmoddata)
          xmod.dummy <- VGAM::vglm(xmodformula, VGAM::multinomial(refLevel=1), data=imputations[[imp]])
          newbeta <- VGAM::coef(xmod) + MASS::mvrnorm(1, mu=rep(0,ncol(VGAM::vcov(xmod))), Sigma=VGAM::vcov(xmod))
          linpreds <- matrix((VGAM::model.matrix(xmod.dummy)) %*% newbeta, byrow=TRUE, ncol=(nlevels(imputations[[imp]][,targetCol])-1))
          denom <- 1+rowSums(exp(linpreds))
          xfitted <-cbind(1/denom, exp(linpreds) / denom)
        }
        if (noisy==TRUE) {
          print(summary(xmod))
        }
        
        #estimate parameters of substantive model
        if (smtype=="lm") {
          ymod <- lm(as.formula(smformula),imputations[[imp]])
          beta <- ymod$coef
          sigmasq <- summary(ymod)$sigma^2
          varcov <- vcov(ymod)
          outcomeModResVar <- (sigmasq*ymod$df) / rchisq(1,ymod$df)
          covariance <- (outcomeModResVar/sigmasq)*vcov(ymod)
          outcomeModBeta = beta + MASS::mvrnorm(1, mu=rep(0,ncol(covariance)), Sigma=covariance)
          if (noisy==TRUE) {
            print(summary(ymod))
          }
        }
        else if (smtype=="logistic") {
          ymod <- glm(as.formula(smformula),family="binomial",imputations[[imp]])
          outcomeModBeta = modPostDraw(ymod)
          if (noisy==TRUE) {
            print(summary(ymod))
          }
        }
        else if (smtype=="poisson") {
          ymod <- glm(as.formula(smformula),family="poisson",imputations[[imp]])
          outcomeModBeta = modPostDraw(ymod)
          if (noisy==TRUE) {
            print(summary(ymod))
          }
        }
        else if (smtype=="coxph") {
          ymod <- survival::coxph(as.formula(smformula), imputations[[imp]],
                                  control = survival::coxph.control(timefix = FALSE))
          outcomeModBeta <- modPostDraw(ymod)
          ymod$coefficients <- outcomeModBeta
          basehaz <- survival::basehaz(ymod, centered=FALSE)[,1]
          H0 <- basehaz[H0indices]
          if (noisy==TRUE) {
            print(summary(ymod))
          }
        }
        else if (smtype=="weibull") {
          ymod <- survival::survreg(as.formula(smformula), data=imputations[[imp]], dist="weibull")
          outcomeModBeta <-  c(coef(ymod), log(ymod$scale)) +
            MASS::mvrnorm(1, mu=rep(0,ncol(vcov(ymod))), Sigma=vcov(ymod))
          weibullScale <- exp(utils::tail(outcomeModBeta,1))
          outcomeModBeta <- utils::head(outcomeModBeta, -1)
          if (noisy==TRUE) {
            print(summary(ymod))
          }
        }
        else if (smtype=="compet") {
          for (cause in 1:numCauses) {
            ymod <- survival::coxph(as.formula(smformula[[cause]]), imputations[[imp]],
                                    control = survival::coxph.control(timefix = FALSE))
            outcomeModBeta[[cause]] <- modPostDraw(ymod)
            ymod$coefficients <- outcomeModBeta[[cause]]
            basehaz <- survival::basehaz(ymod, centered=FALSE)[,1]
            H0[[cause]] <- basehaz[H0indices[[cause]]]
            if (noisy==TRUE) {
              print(summary(ymod))
            }
          }
        }
        else if (smtype=="casecohort") {
          ymod <- survival::coxph(as.formula(smformula2), imputations[[imp]],
                                  control = survival::coxph.control(timefix = FALSE))
          outcomeModBeta <- modPostDraw(ymod)
          
          cumhaz.denom.elements=exp(model.matrix(as.formula(smformula),imputations[[imp]])[,-1] %*% outcomeModBeta)
          cumhaz.denom=sapply(list.times,function(x){sum(cumhaz.denom.elements[which(originaldata[,timeCol]>=x)]*subco.weight[which(originaldata[,timeCol]>=x)])})
          exp.func.denom=cumsum(1/cumhaz.denom)
          H0.fun=stepfun(list.times,c(0,exp.func.denom))
          H0=H0.fun(originaldata[,timeCol])
          
          if (noisy==TRUE) {
            print(summary(ymod))
          }
        }
        else if (smtype=="nestedcc") {
          ymod <- survival::coxph(as.formula(smformula), imputations[[imp]],
                                  control = survival::coxph.control(timefix = FALSE))
          outcomeModBeta <- modPostDraw(ymod)
          
          explan.matrix<-model.matrix(ymod)
          cumbasehaz.denom<-exp(matrix(outcomeModBeta,nrow=1)%*%t(explan.matrix))*originaldata[,nriskCol]/num.sampriskset
          cumbasehaz.denom<-ave(cumbasehaz.denom, originaldata[,setCol], FUN = sum)[originaldata[,dCol]==1] #this is the denominator of the contribution to the cumulative baseline hazard at each event time
          cumbasehaz.t<-originaldata[,timeCol][originaldata[,dCol]==1] #times to which the baseline cumulative hazards refer
          H0<-unlist(lapply(originaldata[,timeCol],function(x) {sum((1/cumbasehaz.denom)[cumbasehaz.t<=x])}))
          
          if (noisy==TRUE) {
            print(summary(ymod))
          }
        }
        
        if ((imp==1) & (cyclenum==1) & (var==1)) {
          if (smtype=="compet") {
            totalCoefVec <- outcomeModBeta[[1]]
            for (cause in 2:numCauses) {
              totalCoefVec <- c(totalCoefVec, outcomeModBeta[[cause]])
            }
            smCoefIter <- array(0, dim=c(m, length(totalCoefVec), numit))
          }
          else {
            smCoefIter <- array(0, dim=c(m, length(outcomeModBeta), numit))
          }
        }
        
        if (var==length(partialVars)) {
          #then we have reached end of a cycle
          if (smtype=="compet") {
            totalCoefVec <- outcomeModBeta[[1]]
            for (cause in 2:numCauses) {
              totalCoefVec <- c(totalCoefVec, outcomeModBeta[[cause]])
            }
            smCoefIter[imp,,cyclenum] <- totalCoefVec
          }
          else {
            smCoefIter[imp,,cyclenum] <- outcomeModBeta
          }
        }
        
        #impute x, either directly where possibly, or using rejection sampling otherwise
        imputationNeeded <- (1:n)[r[,targetCol]==0]
        
        if ((method[targetCol]=="logreg") | (method[targetCol]=="podds") | (method[targetCol]=="mlogit")) {
          #directly sample
          if (method[targetCol]=="logreg") {
            numberOutcomes <- 2
            fittedMean <- cbind(1-xfitted, xfitted)
          }
          else {
            numberOutcomes <- nlevels(imputations[[imp]][,targetCol])
            fittedMean <- xfitted
          }
          
          outcomeDensCovDens = array(dim=c(length(imputationNeeded),numberOutcomes),0)
          
          for (xMisVal in 1:numberOutcomes) {
            if (method[targetCol]=="logreg") {
              if (is.factor(imputations[[imp]][,targetCol])==TRUE) {
                valToImpute <- levels(imputations[[imp]][,targetCol])[xMisVal]
              }
              else {
                valToImpute <- xMisVal-1
              }
            }
            else {
              valToImpute <- levels(imputations[[imp]][,targetCol])[xMisVal]
            }
            imputations[[imp]][imputationNeeded,targetCol] <- valToImpute
            
            #update passive variables
            imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
            
            if (smtype=="lm") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
              deviation <- imputations[[imp]][imputationNeeded,outcomeCol] - outmodxb[imputationNeeded]
              outcomeDens <- dnorm(deviation, mean=0, sd=outcomeModResVar^0.5)
            }
            else if (smtype=="logistic") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
              prob <- expit(outmodxb[imputationNeeded])
              outcomeDens <- prob*imputations[[imp]][imputationNeeded,outcomeCol] + (1-prob)*(1-imputations[[imp]][imputationNeeded,outcomeCol])
            }
            else if (smtype=="poisson") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
              outcomeDens <- dpois(imputations[[imp]][imputationNeeded,outcomeCol], exp(outmodxb[imputationNeeded]))
            }
            else if (smtype=="weibull") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
              #weibull survival function
              outcomeDens <- (1-d[imputationNeeded])*(1-survival::psurvreg(imputations[[imp]][imputationNeeded,timeCol], mean=outmodxb[imputationNeeded], scale=weibullScale))+
                d[imputationNeeded]*(survival::dsurvreg(imputations[[imp]][imputationNeeded,timeCol], mean=outmodxb[imputationNeeded], scale=weibullScale))
            }
            else if ((smtype=="coxph") | (smtype=="casecohort")) {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]])
              outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
              outcomeDens <- exp(-H0[imputationNeeded] * exp(outmodxb[imputationNeeded]))* (exp(outmodxb[imputationNeeded])^d[imputationNeeded])
            }
            else if (smtype=="nestedcc") {
              outmodxb <-  model.matrix(as.formula(smformula2),imputations[[imp]])
              outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
              outcomeDens <- exp(-H0[imputationNeeded] * exp(outmodxb[imputationNeeded]))* (exp(outmodxb[imputationNeeded])^d[imputationNeeded])
            }
            else if (smtype=="compet") {
              outcomeDens <- rep(1,length(imputationNeeded))
              for (cause in 1:numCauses) {
                outmodxb <-  model.matrix(linpred[[cause]],imputations[[imp]])
                outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta[[cause]])
                outcomeDens <- outcomeDens * exp(-H0[[cause]][imputationNeeded] * exp(outmodxb[imputationNeeded]))* (exp(outmodxb[imputationNeeded])^(d[imputationNeeded]==cause))
              }
            }
            outcomeDensCovDens[,xMisVal] <- outcomeDens * fittedMean[imputationNeeded,xMisVal]
          }
          directImpProbs = outcomeDensCovDens / rowSums(outcomeDensCovDens)
          
          if (method[targetCol]=="logreg") {
            directImpProbs = directImpProbs[,2]
            if (is.factor(imputations[[imp]][,targetCol])==TRUE) {
              imputations[[imp]][imputationNeeded,targetCol] <- levels(imputations[[imp]][,targetCol])[1]
              imputations[[imp]][imputationNeeded,targetCol][rbinom(length(imputationNeeded),1,directImpProbs)==1] <- levels(imputations[[imp]][,targetCol])[2]
            }
            else {
              imputations[[imp]][imputationNeeded,targetCol] <- rbinom(length(imputationNeeded),1,directImpProbs)
            }
          }
          else {
            imputations[[imp]][imputationNeeded,targetCol] <- levels(imputations[[imp]][,targetCol])[apply(directImpProbs, 1, catdraw)]
          }
          
          #update passive variables
          imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
        }
        else {
          #use rejection sampling
          #first draw for all subjects who need imputing, using a small number of attempts
          firstTryLimit <- 25
          j <- 1
          
          while ((length(imputationNeeded)>0) & (j<firstTryLimit)) {
            #sample from covariate model
            if ((method[targetCol]=="norm") | (method[targetCol]=="latnorm")) {
              imputations[[imp]][imputationNeeded,targetCol] <- rnorm(length(imputationNeeded),xfitted[imputationNeeded],newsigmasq^0.5)
            }
            else if (method[targetCol]=="poisson") {
              imputations[[imp]][imputationNeeded,targetCol] <- rpois(length(imputationNeeded),xfitted[imputationNeeded])
            }
            
            #update passive variables
            imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
            
            #accept/reject
            uDraw <- runif(length(imputationNeeded))
            if (smtype=="lm") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
              deviation <- imputations[[imp]][imputationNeeded,outcomeCol] - outmodxb[imputationNeeded]
              reject = 1*(log(uDraw) > -(deviation^2) / (2*array(outcomeModResVar,dim=c(length(imputationNeeded),1))))
            }
            else if (smtype=="logistic") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
              prob = expit(outmodxb[imputationNeeded])
              prob = prob*imputations[[imp]][imputationNeeded,outcomeCol] + (1-prob)*(1-imputations[[imp]][imputationNeeded,outcomeCol])
              reject = 1*(uDraw>prob)
            }
            else if (smtype=="poisson") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
              prob = dpois(imputations[[imp]][imputationNeeded,outcomeCol], exp(outmodxb[imputationNeeded]))
              reject = 1*(uDraw>prob)
            } else if (smtype=="weibull") {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
              s_t <- 1-survival::psurvreg(imputations[[imp]][imputationNeeded,timeCol], mean=outmodxb[imputationNeeded], scale=weibullScale)
              prob <- -exp(1)*log(s_t)*s_t
              prob <- d[imputationNeeded]*prob + (1-d[imputationNeeded])*s_t
              reject <- 1*(uDraw>prob)
            } else if ((smtype=="coxph") | (smtype=="casecohort")) {
              outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]])
              outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
              s_t = exp(-H0[imputationNeeded]* exp(outmodxb[imputationNeeded]))
              prob = exp(1 + outmodxb[imputationNeeded] - (H0[imputationNeeded]* exp(outmodxb[imputationNeeded])) ) * H0[imputationNeeded]
              prob = d[imputationNeeded]*prob + (1-d[imputationNeeded])*s_t
              reject = 1*(uDraw > prob )
            }
            else if (smtype=="nestedcc") {
              outmodxb <-  model.matrix(as.formula(smformula2),imputations[[imp]])
              outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
              s_t = exp(-H0[imputationNeeded]* exp(outmodxb[imputationNeeded]))
              prob = exp(1 + outmodxb[imputationNeeded] - (H0[imputationNeeded]* exp(outmodxb[imputationNeeded])) ) * H0[imputationNeeded]
              prob = d[imputationNeeded]*prob + (1-d[imputationNeeded])*s_t
              reject = 1*(uDraw > prob )
            }
            else if (smtype=="compet") {
              prob <- rep(1,length(imputationNeeded))
              for (cause in 1:numCauses) {
                outmodxb <-  model.matrix(linpred[[cause]],imputations[[imp]])
                outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta[[cause]])
                prob = prob * exp(-H0[[cause]][imputationNeeded] * exp(outmodxb[imputationNeeded]))* (H0[[cause]][imputationNeeded]*exp(1+outmodxb[imputationNeeded]))^(d[imputationNeeded]==cause)
              }
              reject = 1*(uDraw > prob )
            }
            imputationNeeded <- imputationNeeded[reject==1]
            
            j <- j+1
          }
          
          #now, for those remaining, who must have low acceptance probabilities, sample by subject
          for (i in imputationNeeded) {
            
            tempData <- imputations[[imp]][i,]
            tempData <- tempData[rep(1,rjlimit),]
            if (method[targetCol]=="norm") {
              tempData[,targetCol] <- rnorm(rjlimit,xfitted[i],newsigmasq^0.5)
            }
            else if (method[targetCol]=="logreg") {
              tempData[,targetCol] <- rbinom(rjlimit,size=1,xfitted[i])
            }
            else if (method[targetCol]=="poisson") {
              tempData[,targetCol] <- rpois(rjlimit,xfitted[i])
            }
            else if (method[targetCol]=="latnorm") {
              tempData[,targetCol] <- rnorm(rjlimit,xfitted[i],newsigmasq[i]^0.5)
            }
            
            #passively impute
            tempData <- updatePassiveVars(tempData, method, passiveVars)
            
            #accept reject
            uDraw <- runif(rjlimit)
            if (smtype=="lm") {
              outmodxb <-  model.matrix(as.formula(smformula),tempData) %*% outcomeModBeta
              deviation <- tempData[,outcomeCol] - outmodxb
              reject = 1*(log(uDraw) > -(deviation^2) / (2*array(outcomeModResVar,dim=c(rjlimit,1))))
            }
            else if (smtype=="logistic") {
              outmodxb <-  model.matrix(as.formula(smformula),tempData) %*% outcomeModBeta
              prob = expit(outmodxb)
              prob = prob*tempData[,outcomeCol] + (1-prob)*(1-tempData[,outcomeCol])
              reject = 1*(uDraw>prob)
            }
            else if (smtype=="poisson") {
              outmodxb <-  model.matrix(as.formula(smformula),tempData) %*% outcomeModBeta
              prob = dpois(tempData[,outcomeCol], exp(outmodxb))
              reject = 1*(uDraw>prob)
            }
            else if (smtype=="weibull") {
              outmodxb <-  model.matrix(as.formula(smformula),tempData) %*% outcomeModBeta
              s_t <- 1-survival::psurvreg(tempData[,timeCol], mean=outmodxb, scale=weibullScale)
              if (d[i]==1) {
                prob <- -exp(1)*log(s_t)*s_t
                #the following line fixes a numerical error which occurs if s_t=0 for any draws
                prob[is.na(prob)] <- 0
              } else {
                prob <- s_t
              }
              reject = 1*(uDraw>prob)
            }
            else if ((smtype=="coxph") | (smtype=="casecohort")) {
              outmodxb <-  model.matrix(as.formula(smformula),tempData)
              outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
              s_t = exp(-H0[i]* exp(outmodxb))
              prob = exp(1 + outmodxb - (H0[i]* exp(outmodxb)) ) * H0[i]
              prob = d[i]*prob + (1-d[i])*s_t
              reject = 1*(uDraw > prob )
            }
            else if (smtype=="nestedcc") {
              outmodxb <-  model.matrix(as.formula(smformula2),tempData)
              outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta)
              s_t = exp(-H0[i]* exp(outmodxb))
              prob = exp(1 + outmodxb - (H0[i]* exp(outmodxb)) ) * H0[i]
              prob = d[i]*prob + (1-d[i])*s_t
              reject = 1*(uDraw > prob )
            }
            else if (smtype=="compet") {
              prob <- rep(1,rjlimit)
              for (cause in 1:numCauses) {
                outmodxb <-  model.matrix(linpred[[cause]],tempData)
                outmodxb <- as.matrix(outmodxb[,2:dim(outmodxb)[2]]) %*% as.matrix(outcomeModBeta[[cause]])
                prob = prob * exp(-H0[[cause]][i] * exp(outmodxb))* (H0[[cause]][i]*exp(1+outmodxb))^(d[i]==cause)
              }
              reject = 1*(uDraw > prob )
            }
            if (sum(reject)<rjlimit) {
              imputations[[imp]][i,targetCol] <- tempData[reject==0,targetCol][1]
            } else {
              rjFailCount <- rjFailCount + 1
            }
          }
          #update passive variables
          imputations[[imp]] <- updatePassiveVars(imputations[[imp]], method, passiveVars)
        }
      }
      
      #imputations of missing outcomes, if present (using proper imputation), for regression and logistic
      #substantive models
      if ((smtype=="lm") | (smtype=="logistic")) {
        if (sum(r[,outcomeCol])<n) {
          imputationNeeded <- (1:n)[r[,outcomeCol]==0]
          #estimate parameters of substantive model using those with outcomes observed
          if (smtype=="lm") {
            ymod <- lm(as.formula(smformula),imputations[[imp]][r[,outcomeCol]==1,])
            beta <- ymod$coef
            sigmasq <- summary(ymod)$sigma^2
            varcov <- vcov(ymod)
            outcomeModResVar <- (sigmasq*ymod$df) / rchisq(1,ymod$df)
            covariance <- (outcomeModResVar/sigmasq)*vcov(ymod)
            outcomeModBeta = beta + MASS::mvrnorm(1, mu=rep(0,ncol(covariance)), Sigma=covariance)
            outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
            imputations[[imp]][imputationNeeded,outcomeCol] <- rnorm(length(imputationNeeded),outmodxb[imputationNeeded], sigmasq^0.5)
          }
          else if (smtype=="logistic") {
            ymod <- glm(as.formula(smformula),family="binomial",imputations[[imp]][r[,outcomeCol]==1,])
            outcomeModBeta = modPostDraw(ymod)
            outmodxb <-  model.matrix(as.formula(smformula),imputations[[imp]]) %*% outcomeModBeta
            prob <- expit(outmodxb[imputationNeeded])
            imputations[[imp]][imputationNeeded,outcomeCol] <- rbinom(length(imputationNeeded),1,prob)
          }
        }
      }
    }
    
  }
  
  if (rjFailCount>0) {
    warning(paste("Rejection sampling failed ",rjFailCount," times (across all variables, iterations, and imputations). You may want to increase the rejection sampling limit.",sep=""))
  }
  
  list(impDatasets=imputations, smCoefIter=smCoefIter)
  
}

updatePassiveVars <- function(data, method, passivecols) {
  for (i in passivecols) {
    data[,i] <- with(data, eval(parse(text=method[i])))
  }
  data
}

expit <- function(x) {
  exp(x)/(1+exp(x))
}

sumna <- function(x) {
  sum(is.na(x)==FALSE)
}

#returns first non missing entry of x
firstnonna <- function(x) {
  x[is.na(x)==FALSE][1]
}

catdraw <- function(prob) {
  (1:length(prob))[rmultinom(1,size=1,prob=prob)==1]
}

modPostDraw <- function(modobj) {
  beta <- modobj$coef
  varcov <- vcov(modobj)
  beta + MASS::mvrnorm(1, mu=rep(0,ncol(varcov)), Sigma=varcov)
}

