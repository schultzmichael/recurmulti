############################
# Model Term Functions
############################
presence.matrix <- function(sequences){
  speakers <- sort(unique(unlist(sequences)))

  x <- matrix(NA,length(sequences),length(speakers))
  for(i in 1:nrow(x)){
    for(j in 1:length(speakers))
      if(speakers[j]%in%sequences[[i]])
        x[i,j] <- 1
      x[i,which(x[i,]==1)[1]] <- 0
  }

  x
}
sequence.matrix <- function(sequences,seq.padding){
  x <- matrix(0,ncol=length(sequences),nrow=sum(sapply(sequences,length))+(length(sequences)-1)*seq.padding)
  seq.start <- c(0,c(cumsum(sapply(sequences,length)) + c(1:length(sequences))*seq.padding))
  seq.start <- seq.start[1:length(sequences)]
  seq.end <- cumsum(sapply(sequences,length))+ c(0:(length(sequences)-1))*seq.padding
  for(i in 1:nrow(x)){
    col <- which(i>seq.start & i <= seq.end)
    x[i,col] <- 1
  }

  x
}
seqnum.term <- function(sequences,n,seq.padding){
  speakers <- sort(unique(unlist(sequences)))
  S <- length(speakers)

  x <- rep(0,length(sequences[[1]]))
  for(i in 2:length(sequences))
    x <- c(x,rep(NA,seq.padding),rep(1*(i==n),length(sequences[[i]])))

  x0 <- cbind(x)
  x <- matrix(0,nrow=nrow(x0),ncol=S)
  for(i in 1:length(S))
    x[,i] <- x0

  #colnames(x) <- paste0('seq',n)
  x
}
seqfrac.term <- function(sequences,seq.padding){
  speakers <- sort(unique(unlist(sequences)))
  S <- length(speakers)

  x <- seq(0,1,length.out=length(sequences[[1]]))
  if(length(sequences)>1)
    for(i in 2:length(sequences))
      x <- c(x, rep(NA,seq.padding),
             seq(0,1,length.out=length(sequences[[i]])))
  x0 <- cbind(x)
  x <- matrix(0,nrow=nrow(x0),ncol=S)
  for(i in 1:length(S))
    x[,i] <- x0

  #colnames(x) <- paste0('seq',n)
  x
}
seqfracsq.term <- function(sequences,seq.padding){
  speakers <- sort(unique(unlist(sequences)))
  S <- length(speakers)

  x <- (seq(0,1,length.out=length(sequences[[1]])))^2
  if(length(sequences)>1)
    for(i in 2:length(sequences))
      x <- c(x, rep(NA,seq.padding),
             (seq(0,1,length.out=length(sequences[[i]])))^2)
  x0 <- cbind(x)
  x <- matrix(0,nrow=nrow(x0),ncol=S)
  for(i in 1:length(S))
    x[,i] <- x0

  #colnames(x) <- paste0('seq',n)
  x
}
presence.term <- function(sequences,seq.padding){
  speakers <- sort(unique(unlist(sequences)))
  S <- length(speakers)
  x <- sapply(speakers,function(sp) rep(1*sp%in%sequences[[1]],length(sequences[[1]])))
  if(length(sequences)>1)
    for(i in 2:length(sequences))
      x <- rbind(x, matrix(NA,nrow=seq.padding,ncol=ncol(x)),
                 sapply(speakers,function(sp) rep(1*sp%in%sequences[[i]],length(sequences[[i]]))))
  x
}
absence.term <- function(sequences,seq.padding){
  speakers <- sort(unique(unlist(sequences)))
  S <- length(speakers)
  x <- sapply(speakers,function(sp) rep(1*!(sp%in%sequences[[1]]),length(sequences[[1]])))
  if(length(sequences)>1)
    for(i in 2:length(sequences))
      x <- rbind(x, matrix(NA,nrow=seq.padding,ncol=ncol(x)),
                 sapply(speakers,function(sp) rep(1*!(sp%in%sequences[[i]]),length(sequences[[i]]))))
  x
}

corr.term <- function(sequences, seq.padding,covar,name){
  s <- flatten.sequences(sequences,seq.padding)
  speakers <- sort(unique(s))
  #browser()
  if(is.function(covar)){
    x <- Reduce(cbind,lapply(speakers,function(a)
      flatten.sequences(
        lapply(1:length(sequences),function(i) sapply(1:length(sequences[[i]]),function(l)
          covar(sequences[[i]],i,l,a))),seq.padding)))
  }else{
    stop('Nonfunction correlation covariates not supported yet')
  }
  x
}
flatten.sequences <- function(seq,padding){
  N <- length(seq)
  s <- seq[[1]]
  if(N>1)
    for(i in 2:N)
      s <- c(s,rep(NA,padding),seq[[i]])
  s
}

#' @importFrom abind abind
prep.terms <- function(sequences, seq.padding,
                       prob.terms=c(),
                       corr.terms=c(),
                       corr.seq.covar=list(), ...){
  #Add dim names, gradient
  #Word length
  #holding the floor
  #Terms - Constant, Sequence , Sequence fraction, sequence covars, person covars, recency, recent frac

  cat(paste0("Flattening ",length(sequences)," sequences...\n"))

  if (!is.list(sequences))
    sequences <- list(sequences)
  N <- length(sequences)
  s <- flatten.sequences(sequences,seq.padding)

  L <- length(s)
  speakers <- sort(unique(s))
  S <- length(speakers)

  P <- presence.matrix(sequences)
  Q <- sequence.matrix(sequences, seq.padding)

  X <- array(1,dim=c(L,0,S))
  Z <- array(1,dim=c(L,0,S))


  if('constant' %in% prob.terms){
    X <- array(1,dim=c(L,1,S))
    dimnames(X) <- list(c(),'(Intercept)',speakers)
  }
  if('sequence' %in% prob.terms & N>1){
    dnx <- dimnames(X)[[2]]
    for(i in 2:length(sequences))
      X <- abind(X,matrix(seqnum.term(sequences,i,seq.padding),nrow=L,ncol=S),along=2)
    dimnames(X)[[2]] <- c(dnx,paste0('seq',2:length(sequences)))
  }
  if('sequence fraction' %in% prob.terms)
    X <- abind(X,frac=matrix(seqfrac.term(sequences,seq.padding),nrow=L,ncol=S),along=2)
  if('presence' %in% prob.terms)
    X <- abind(X,presence=presence.term(sequences,seq.padding),along=2)
  if('absence' %in% prob.terms)
    X <- abind(X,absence=absence.term(sequences,seq.padding),along=2)

  if('constant' %in% corr.terms){
    Z <- array(1,dim=c(L,1,S))
    dimnames(Z) <- list(c(),'(Intercept)',speakers)
  }
  if('sequence' %in% corr.terms & N>1){
    dnz <- dimnames(Z)[[2]]
    #browser()
    for(i in 2:length(sequences))
      Z <- abind(Z,seqnum.term(sequences,i,seq.padding),along=2)
    dimnames(Z)[[2]] <- c(dnz,paste0('seq',2:length(sequences)))
  }
  if('sequence fraction' %in% corr.terms)
    Z <- abind(Z,'frac'=seqfrac.term(sequences,seq.padding),along=2)
  if('sequence fraction squared' %in% corr.terms)
    Z <- abind(Z,'fracsq'=seqfracsq.term(sequences,seq.padding),along=2)

  if(length(corr.seq.covar)>0){
    dnz <- dimnames(Z)[[2]]
    for (i in 1:length(corr.seq.covar))
      Z<- abind(Z,
                corr.term(sequences,seq.padding,corr.seq.covar[[i]],names(corr.seq.covar)[i]),
                along=2)
    dimnames(Z)[[2]] <- c(dnz,names(corr.seq.covar))
  }

  list(s=s,X=X,Z=Z,Q=Q,P=P)
}

#############
# Top Level Functions
############
#' Print an RMM
#'
#' Summarizes a fitted RMM
#'
#' @param r A fitted rmm object estimated using \code{rmm}
#'
#' @export
print.rmm <- function(model,...){
  cat(paste0('Model of ',model$N,' sequences (L=',model$L,') with ',model$D,' lags\n'))
  cat(paste0('Elements: ',model$S,', ','Lambda: ',model$lambda,'\n'))
  cat('\n')
  # if('beta.se' %in% names(r)){
  #   coef <- cbind(beta=round(r$beta,4),
  #                 t=round(r$beta/r$beta.se,4))
  #   colnames(coef) <- c(colnames(r$beta),paste0(colnames(r$beta),'.t'))
  # }else{
  #   coef <- cbind(beta=round(r$beta,4))
  #   colnames(coef) <- colnames(r$beta)
  # }
  if(prod(dim(model$beta))>0){
    cat('Beta Coefficients:\n')
    if('beta.se' %in% names(model)){
      print.coef.matrix(model$beta,model$beta.se,...)
    }else{
      print.matrix(model$beta)
    }
  }
  cat('\n')
  cat(paste0('Log Likelihood: ',round(model$loglik,2),', '))
  cat(paste0('AIC: ',round(model$AIC,2),', '))
  cat(paste0('BIC: ',round(model$BIC,2),'\n'))
  cat(paste0('Convergence Code: ',model$opt$convergence,' - ',model$opt$message,'\n'))
}

#' Estimate a Recurrent Multinomial Model
#'
#' \code{RMM} estimates a recurrent mutlinomial model of specified length with
#'   covariate terms.
#'
#' This is the primary estimation function for fitting RMMs
#'
#' @param sequences A list of vectors with sequence data
#' @param prob.terms A list of predefined covariate terms to include in the X
#'   matrix predicting baseline element-specific probability. See details.
#' @param recur.terms A list of predefined covariate terms to include in the Z
#'   matrix predicting recurrence coefficients. See details.
#' @param recurr.seq.covar A list of functions for defining custom recurrence
#'   covariates
#' @param hessian Logical value indicating whether the hessian should be
#'   computed
#'
#' @param etastart,betastart,zetastart Starting values for optimization.
#'    Defaults to 0.
#' @param parallel Logical values indicating whether to use parallel computation
#' @param cores Number of cores to use in parallel computation
#' @param num.grad Logical value indicating whether to use a numerical gradient
#' @param lambda,lambda2 Regularization parameter for recurrence and baseline
#'   probability coefficients respectively
#'
#' @section Baseline Covariates:
#' These covariates provide an element specific estimate for the probability
#' of an element appearing, as with a normal multinomial model coefficient.
#' Prespecified baseline covariates include:
#' * \code{'constant'} - Intercept
#' * \code{'sequence'} - Sequence specific probability
#' * \code{'sequence fraction'} - Fraction of sequence complete
#' * \code{'presence'} - Whether the element occurs at all in the observed
#'   sequence
#'
#' @section Recurrence Covariates:
#' These covariates provide an lag specific estimate for recurrence of a
#' previously observed element reappearing.
#' Prespecified recurrence covariates include:
#' * \code{'constant'} - Intercept
#' * \code{'sequence'} - Sequence specific probability
#' * \code{'sequence fraction'} - Fraction of sequence complete
#' * \code{'sequence fraction squared'} - (Fraction of sequence complete)^2
#'
#' @section Custom Covariates:
#' Custom covariates can be defined by adding functions to the
#' \code{recur.seq.covar} list. These functions should take four parameters:
#' \code{s0}, \code{i}, \code{l}, and \code{a} where \code{s0} is a observed
#' element in at position \code{l} in the \code{i}th sequence and \code{a} is a
#' possible element.
#'
#' Possible uses for these functions include controlling for elements
#' belonging to a certain category, sequences of a certain type, timing
#' specific recurrences, or combinations of those effects.
#'
#' @return An fitted recurrent multinomial model of class 'rmm'
#' @export
rmm <- function(sequences,D,
                prob.terms=c(),
                recur.terms=c(),
                recur.seq.covar=list(),
                hessian=TRUE, ...){

  ps <- prep.terms(sequences,seq.padding=D,
                   prob.terms=prob.terms,
                   corr.terms=recur.terms,
                   corr.seq.covar=recur.seq.covar,...)

  out <- fit.seq(ps$s,D,ps$X,ps$Z,ps$Q,ps$P, ...)
  if(hessian)
    out <- infer.seq(out)


  class(out) <- 'rmm'
  out
}

#' @describeIn rmm Evaluate the probability of observed sequences on a prefit
#' rmm model
#' @export
predict.rmm <- function(model, newdata,
                        prob.terms=c(),
                        corr.terms=c(),
                        corr.seq.covar=list(), ...){
  ps <- prep.terms(newdata,seq.padding=model$D,
                   prob.terms=prob.terms,
                   corr.terms=corr.terms,
                   corr.seq.covar=corr.seq.covar,...)

  out <- eval.seq(s=ps$s,X=ps$X,Z=ps$Z,Q=ps$Q,P=ps$P,beta=model$beta, ...)
  class(out) <- 'rmm'
  out
}
