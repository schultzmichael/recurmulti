#####################
## General Functions
#####################
aselect <- function(x,i){
  matrix(t(sapply(1:dim(x)[1],function(j) x[j,,i[j]])),ncol=dim(x)[2])
}
row.select <- function(x,i)
  x[cbind(i,1:ncol(x))]

decode.params <- function(params,S,J,D,K,P){#holdouts here
  eta <- c(rep(0,J),params[1:((S-1)*J)])
  beta <- params[(S-1)*J + (1:(D*K))]
  #eta <- c(0,params[1:(S*J-1)])
  #beta <- params[(S*J-1) + (1:(D*K))]
  eta <- matrix(eta,nrow=S,ncol=J,byrow=TRUE)
  if(D>0){
    beta <- matrix(beta,nrow=D,ncol=K)
  }else{
    beta <- matrix(NA,nrow=0,ncol=K)
  }
  zeta <- matrix(NA,nrow=nrow(P),ncol=ncol(P))
  zeta[!is.na(P)] <- 0
  zeta[which(P==1)] <- params[((S-1)*J + D*K) + (1:sum(P,na.rm=TRUE))]
  list(eta=eta,beta=beta,zeta=zeta)
}
encode.params <- function(eta,beta,zeta,P){#holdouts here
  params <- c()
  # Add etas
  if(length(eta)>0){
    J <- ncol(eta)
    params <- as.vector(matrix(eta,nrow=1))
    params <- params[(J+1):length(params)]
    #params <- params[2:length(params)]
  }
  # Add betas
  if(length(beta)>0)
    params <- c(params,as.vector(matrix(beta,nrow=1)))

  #Add zetas
  params <- c(params, as.vector(zeta[which(P==1)]))
  params
}
decode.params.pred <- function(params,S,J,P){#holdouts here
  eta <- c(rep(0,J),params[1:((S-1)*J)])
  eta <- matrix(eta,nrow=S,ncol=J,byrow=TRUE)
  # zeta <- matrix(NA,nrow=nrow(P),ncol=ncol(P))
  # zeta[!is.na(P)] <- 0
  # zeta[which(P==1)] <- params[((S-1)*J) + (1:sum(P,na.rm=TRUE))]
  zeta <- matrix(0,nrow=nrow(P),ncol=ncol(P))
  #zeta[!is.na(P)] <- 0
  zeta[which(P==1)] <- params[((S-1)*J) + (1:sum(P,na.rm=TRUE))]
  list(eta=eta,zeta=zeta)
}
encode.params.pred <- function(eta,zeta,P){#holdouts here
  params <- c()

  # Add etas
  if(length(eta)>0){
    J <- ncol(eta)
    params <- as.vector(matrix(eta,nrow=1))
    params <- params[(J+1):length(params)]
    #params <- params[2:length(params)]
  }

  #Add zetas
  params <- c(params, as.vector(zeta[which(P==1)]))
  params
}

calc.prob <- function(s,X,Z,Q,eta,beta,zeta){

  D <- nrow(beta)
  S <- ncol(zeta)
  L <- length(s[!is.na(s)])

  pi <- matrix(0,nrow=dim(Q)[1],ncol=S)

  # eta/X part
  if(ncol(eta)>0){
    for(a in 1:S)
      pi[,a] <- X[,,a,drop=FALSE]%*%eta[a,,drop=FALSE]
  }

  # beta/Z part
  if(D>0){
    l <- dim(Z)[1]
    G <- matrix(0,nrow=l,ncol=nrow(beta))
    for(a in 1:S){
      for(d in 1:D){
        G[(d+1):l,d] <- Z[1:(l-d),,a]%*%t(beta[d,,drop=FALSE])
      }
      delta <- sapply(1:D,function(d){ # should be LxD
        s1 <- c(rep(NA,d),s[1:(length(s)-d)])
        eq<-(s1==a)*1
        eq[is.na(eq)] <- 0
        eq
      })
      pi[,a] <- pi[,a]+rowSums(G*delta,na.rm=TRUE)
    }
  }

  # zeta/Q part
  z <- zeta
  z[is.na(z)] <- 0
  Qz <- Q%*%z
  Qzna <- Q%*%(1*is.na(zeta))
  Qz[which(Qzna>0)] <- NA
  pi <- pi + Qz #r$Q%*%t(zeta) #new!

  epi <- exp(pi)
  prob <- t(epi/rowSums(epi,na.rm=TRUE))

  list(prob=prob,epi=epi)
}
#######################
## Evaluation
#######################

#' @importFrom optimParallel optimParallel
#' @importFrom parallel makeCluster
#' @importFrom parallel setDefaultCluster
eval.seq <- function(s,
                     X,Z,
                     beta=beta,
                     Q=NA,P=NA,
                     etastart=0,zetastart=0,
                     parallel=FALSE,cores=5, num.grad=FALSE,
                     lambda=0, lambda2=0, ...){
  eps <- 1e-6

  D <- nrow(beta)

  ns <- sort(unique(s))
  S <- length(unique(s[!is.na(s)]))
  L <- length(s)

  if(is.na(Q[1])) # Probability non-time varying covariates
    Q <- matrix(1,nrow=L,ncol=1)
  if(all(is.na(P))) # Structural zeros for probability
    P <- matrix(1,nrow=1,ncol=S)

  K <- dim(Z)[2]
  J <- dim(X)[2]
  N <- dim(Q)[2]

  if(!parallel)
    cores <- 1

  cat(paste0("Evaluating sequence of ",S," speakers and length ",L," with degree ",D," (cores = ",cores,")\n"))
  cat(paste0("Regularization constant: ",lambda," / ", lambda2, "\n"))
  cat(paste0("Correlation/Speaker covariates: ",K," / ",J,"\n"))

  s.transform <- 1:S
  names(s.transform) <- ns
  s <- s.transform[s]

  eta0 <- matrix(etastart,nrow=S,ncol=J)
  zeta0 <- matrix(zetastart,nrow=N,ncol=S)

  f <- function(param){

    row.select <- function(x,i)
      x[cbind(i,1:ncol(x))]

    decode.params.pred <- function(params,S,J,P){#holdouts here
      eta <- c(rep(0,J),params[1:((S-1)*J)])
      eta <- matrix(eta,nrow=S,ncol=J,byrow=TRUE)
      # zeta <- matrix(NA,nrow=nrow(P),ncol=ncol(P))
      # zeta[!is.na(P)] <- 0
      # zeta[which(P==1)] <- params[((S-1)*J) + (1:sum(P,na.rm=TRUE))]
      zeta <- matrix(0,nrow=nrow(P),ncol=ncol(P))
      #zeta[!is.na(P)] <- 0
      zeta[which(P==1)] <- params[((S-1)*J) + (1:sum(P,na.rm=TRUE))]
      list(eta=eta,zeta=zeta)
    }
    encode.params.pred <- function(eta,zeta,P){#holdouts here
      params <- c()
      # Add etas
      if(length(eta)>0){
        J <- ncol(eta)
        params <- as.vector(matrix(eta,nrow=1))
        params <- params[(J+1):length(params)]
        #params <- params[2:length(params)]
      }

      #Add zetas
      params <- c(params, as.vector(zeta[which(P==1)]))
      params
    }
    calc.prob <- function(s,X,Z,Q,eta,beta,zeta){

      D <- nrow(beta)
      S <- ncol(zeta)
      L <- length(s[!is.na(s)])

      pi <- matrix(0,nrow=dim(Q)[1],ncol=S)

      # eta/X part
      if(ncol(eta)>0){
        for(a in 1:S)
          pi[,a] <- X[,,a,drop=FALSE]%*%eta[a,,drop=FALSE]
      }

      # beta/Z part
      if(D>0){
        l <- dim(Z)[1]
        G <- matrix(0,nrow=l,ncol=nrow(beta))
        for(a in 1:S){
          for(d in 1:D){
            G[(d+1):l,d] <- Z[1:(l-d),,a]%*%t(beta[d,,drop=FALSE])
          }
          delta <- sapply(1:D,function(d){ # should be LxD
            s1 <- c(rep(NA,d),s[1:(length(s)-d)])
            eq<-(s1==a)*1
            eq[is.na(eq)] <- 0
            eq
          })
          pi[,a] <- pi[,a]+rowSums(G*delta,na.rm=TRUE)
        }
      }

      # zeta/Q part
      z <- zeta
      z[is.na(z)] <- 0
      Qz <- Q%*%z
      Qzna <- Q%*%(1*is.na(zeta))
      Qz[which(Qzna>0)] <- NA
      pi <- pi + Qz #r$Q%*%t(zeta) #new!

      epi <- exp(pi)
      prob <- t(epi/rowSums(epi,na.rm=TRUE))

      list(prob=prob,epi=epi)
    }
    log.likelihood <- function(samples, X,Z,Q,eta,beta,zeta,...){

      cp <- calc.prob(samples,X,Z,Q,eta,beta,zeta)
      cp$prob[cp$prob<1e-8] <- 1e-8

      ll <- log(row.select(cp$prob,samples))
      ll <- sum(ll,na.rm=TRUE)
      if(!is.finite(ll))
        browser()
      if(D>0)
        ll <- ll - lambda*L/D*sqrt(sum(beta^2))
      ll
    }
    LL <- function(params,s,X,Z,Q,P, ...){
      S <- dim(X)[3]
      K <- ncol(Z)
      J <- dim(X)[2]
      theta <- decode.params.pred(params=params,S,J,P)
      log.likelihood(s,X,Z,Q,eta=theta$eta,beta=beta,zeta=theta$zeta,...)
    }
    -LL(param,s,X,Z,Q,P)
  }

  # Gradient function for optimization
  g <- function(param){
    row.select <- function(x,i)
      x[cbind(i,1:ncol(x))]
    col.select <- function(x,i)
      x[cbind(1:nrow(x),i)]

    aselect <- function(x,i){
      matrix(t(sapply(1:dim(x)[1],function(j) x[j,,i[j]])),ncol=dim(x)[2])
    }

    decode.params.pred <- function(params,S,J,P){#holdouts here
      eta <- c(rep(0,J),params[1:((S-1)*J)])
      eta <- matrix(eta,nrow=S,ncol=J,byrow=TRUE)
      zeta <- matrix(NA,nrow=nrow(P),ncol=ncol(P))
      zeta[!is.na(P)] <- 0
      zeta[which(P==1)] <- params[((S-1)*J) + (1:sum(P,na.rm=TRUE))]
      list(eta=eta,zeta=zeta)
    }
    encode.params.pred <- function(eta,zeta,P){#holdouts here
      params <- c()
      # Add etas
      if(length(eta)>0){
        J <- ncol(eta)
        params <- as.vector(matrix(eta,nrow=1))
        params <- params[(J+1):length(params)]
        #params <- params[2:length(params)]
      }

      #Add zetas
      params <- c(params, as.vector(zeta[which(P==1)]))
      params
    }
    calc.prob <- function(s,X,Z,Q,eta,beta,zeta){

      D <- nrow(beta)
      S <- ncol(zeta)
      L <- length(s[!is.na(s)])

      pi <- matrix(0,nrow=dim(Q)[1],ncol=S)

      # eta/X part
      if(ncol(eta)>0){
        for(a in 1:S)
          pi[,a] <- X[,,a,drop=FALSE]%*%eta[a,,drop=FALSE]
      }

      # beta/Z part
      if(D>0){
        l <- dim(Z)[1]
        G <- matrix(0,nrow=l,ncol=nrow(beta))
        for(a in 1:S){
          for(d in 1:D){
            G[(d+1):l,d] <- Z[1:(l-d),,a]%*%t(beta[d,,drop=FALSE])
          }
          delta <- sapply(1:D,function(d){ # should be LxD
            s1 <- c(rep(NA,d),s[1:(length(s)-d)])
            eq<-(s1==a)*1
            eq[is.na(eq)] <- 0
            eq
          })
          pi[,a] <- pi[,a]+rowSums(G*delta,na.rm=TRUE)
        }
      }

      # zeta/Q part
      z <- zeta
      z[is.na(z)] <- 0
      Qz <- Q%*%z
      Qzna <- Q%*%(1*is.na(zeta))
      Qz[which(Qzna>0)] <- NA
      pi <- pi + Qz #r$Q%*%t(zeta) #new!

      epi <- exp(pi)
      prob <- t(epi/rowSums(epi,na.rm=TRUE))

      list(prob=prob,epi=epi)
    }
    log.lik.grad <- function(samples, X,Z,Q,eta,beta,zeta,...){

      s <- samples
      D <- nrow(beta)
      S <- ncol(zeta)
      L <- length(s[!is.na(s)])

      cp <- calc.prob(s,X,Z,Q,eta,beta,zeta)
      p <- row.select(cp$prob,s)
      epi <- cp$epi

      eg <- eta*0
      zg <- zeta*0
      bg <- beta*0

      if(dim(eta)[2]>0){
        for(i in 1:dim(eg)[1])
          for(j in 1:dim(eg)[2])
            eg[i,j] <- sum( X[,j,i]*(1*(s==i)) - X[,j,i]*epi[,i]/rowSums(epi,na.rm=TRUE),na.rm=TRUE)
      }

      if(dim(zeta)[2]>0){
        for(i in 1:dim(zg)[2])
          for(j in 1:dim(zg)[1])
            zg[j,i] <- sum( Q[,j]*(1*(s==i)) - Q[,j]*epi[,i]/rowSums(epi,na.rm=TRUE), na.rm=TRUE)
      }

      if(D>0){
        ZYi <- aselect(Z,s)
        for(d in 1:dim(bg)[1]){
          delta <- sapply(1:S,function(a){
            s1 <- c(rep(NA,d),s[1:(length(s)-d)])
            eq<-(s1==a)*1
            eq[is.na(eq)] <- 0
            eq
          })
          dYi <- col.select(delta,s)
          for(j in 1:dim(bg)[2]){
            bg[d,j] <- sum( ZYi[,j,drop=TRUE]*dYi - rowSums(Z[,j,]*delta *epi) / rowSums(epi,na.rm=TRUE), na.rm=TRUE)
          }
        }
      }

      list(beta=bg,eta=eg,zeta=zg)
    }
    LLgrad <- function(params,s,D,X,Z,Q,P, ...){
      S <- dim(X)[3]
      K <- ncol(Z)
      J <- dim(X)[2]
      theta <- decode.params.pred(params=params,S,J,P)
      grads <- log.lik.grad(s,X,Z,Q,eta=theta$eta,beta=beta,zeta=theta$zeta,...)
      encode.params.pred(grads$eta,grads$zeta,P)
    }
    -LLgrad(param,s,D,X,Z,Q,P)
  }

  params <- encode.params.pred(eta0,zeta0,P)

  gr <- g
  if(num.grad)
    gr <- NULL
  if(parallel & cores>1){
    cl <- makeCluster(cores)     # set the number of processor cores
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    opt <- optimParallel(params,f,gr=gr,hessian=FALSE,method='L-BFGS-B',
                         parallel=list(forward=TRUE))
    stopCluster(cl)
  }else{
    opt <- optim(params,f,gr=gr,hessian=FALSE,method='L-BFGS-B')
  }
  theta <- decode.params.pred(opt$par,S,J,P)

  eta<- theta$eta
  rownames(eta) <- ns
  colnames(eta) <- dimnames(X)[[2]]

  if(nrow(beta)>0)
    rownames(beta) <- paste0('Lag',1:nrow(beta))
  colnames(beta) <- colnames(Z)

  zeta <- t(theta$zeta)
  rownames(zeta) <- ns
  colnames(zeta) <- paste0('Seq',1:ncol(zeta))

  cp <- calc.prob(s,X,Z,Q,eta,beta,t(zeta))
  prob <- cp$prob

  fit <- list(s=s,
              X=X,Z=Z,Q=Q,P=P,
              S=S,D=D,J=J,K=K,N=N,
              L = length(s[!is.na(s)]),
              k = length(opt$par),

              deviance = 2*opt$value,
              AIC = 2*opt$value +2*length(opt$par),
              BIC = 2*opt$value +log(length(s[!is.na(s)]))*length(opt$par),

              loglik=-opt$value,

              eta=eta, beta=beta, zeta=zeta,
              prob = prob,
              params = opt$par,
              opt=opt,
              lambda=lambda,
              lambda2=lambda2
  )

  fit
}

###########################
## Fit
###########################

#' @importFrom optimParallel optimParallel
#' @importFrom parallel makeCluster
#' @importFrom parallel setDefaultCluster
fit.seq <- function(s, D,
                    X,Z,
                    Q=NA,P=NA,
                    etastart=0,betastart=0,zetastart=0,
                    parallel=FALSE,cores=5, num.grad=FALSE,
                    lambda=0, lambda2=0, ...){
  eps <- 1e-6

  ns <- sort(unique(s))
  S <- length(unique(s[!is.na(s)]))
  L <- length(s)

  if(is.na(Q[1])) # Probability non-time varying covariates
    Q <- matrix(1,nrow=L,ncol=1)
  if(all(is.na(P))) # Structural zeros for probability
    P <- matrix(1,nrow=1,ncol=S)

  K <- dim(Z)[2]
  J <- dim(X)[2]
  N <- dim(Q)[2]

  if(!parallel)
    cores <- 1

  cat(paste0("Modeling sequence of ",S," speakers and length ",L," with degree ",D," (cores = ",cores,")\n"))
  cat(paste0("Regularization constant: ",lambda," / ", lambda2, "\n"))
  cat(paste0("Correlation/Speaker covariate: ",K," / ",J,"\n"))

  s.transform <- 1:S
  names(s.transform) <- ns
  s <- s.transform[s]

  eta0 <- matrix(etastart,nrow=S,ncol=J)
  beta0 <- matrix(betastart,nrow=D,ncol=K)
  zeta0 <- matrix(zetastart,nrow=N,ncol=S)

  beta.bounds <- c(-Inf,Inf)
  if('beta.bounds' %in% names(list(...)))
    beta.bounds <- list(...)$beta.bounds


  # Helper functions are included within f so they are accessible in parallel computation
  # Likelihood function for optimization
  f <- function(param){
    row.select <- function(x,i)
      x[cbind(i,1:ncol(x))]

    decode.params <- function(params,S,J,D,K,P){#holdouts here
      eta <- c(rep(0,J),params[1:((S-1)*J)])
      beta <- params[(S-1)*J + (1:(D*K))]
      eta <- matrix(eta,nrow=S,ncol=J,byrow=TRUE)
      if(D>0){
        beta <- matrix(beta,nrow=D,ncol=K)
      }else{
        beta <- matrix(NA,nrow=0,ncol=K)
      }
      zeta <- matrix(0,nrow=nrow(P),ncol=ncol(P))
      #zeta[!is.na(P)] <- 0
      zeta[which(P==1)] <- params[((S-1)*J + D*K) + (1:sum(P,na.rm=TRUE))]
      list(eta=eta,beta=beta,zeta=zeta)
    }
    encode.params <- function(eta,beta,zeta,P){#holdouts here
      params <- c()
      # Add etas
      if(length(eta)>0){
        J <- ncol(eta)
        params <- as.vector(matrix(eta,nrow=1))
        params <- params[(J+1):length(params)]
      }
      # Add betas
      if(length(beta)>0)
        params <- c(params,as.vector(matrix(beta,nrow=1)))

      #Add zetas
      params <- c(params, as.vector(zeta[which(P==1)]))
      params
    }
    calc.prob <- function(s,X,Z,Q,eta,beta,zeta){

      D <- nrow(beta)
      S <- ncol(zeta)
      L <- length(s[!is.na(s)])

      pi <- matrix(0,nrow=dim(Q)[1],ncol=S)

      # eta/X part
      if(ncol(eta)>0){
        for(a in 1:S)
          pi[,a] <- X[,,a,drop=FALSE]%*%eta[a,,drop=FALSE]
      }

      # beta/Z part
      if(D>0){
        l <- dim(Z)[1]
        G <- matrix(0,nrow=l,ncol=nrow(beta))
        for(a in 1:S){
          for(d in 1:D){
            G[(d+1):l,d] <- Z[1:(l-d),,a]%*%t(beta[d,,drop=FALSE])
          }
          delta <- sapply(1:D,function(d){ # should be LxD
            s1 <- c(rep(NA,d),s[1:(length(s)-d)])
            eq<-(s1==a)*1
            eq[is.na(eq)] <- 0
            eq
          })
          pi[,a] <- pi[,a]+rowSums(G*delta,na.rm=TRUE)
        }
      }

      # zeta/Q part
      z <- zeta
      z[is.na(z)] <- 0
      Qz <- Q%*%z
      Qzna <- Q%*%(1*is.na(zeta))
      Qz[which(Qzna>0)] <- NA
      pi <- pi + Qz #r$Q%*%t(zeta) #new!

      epi <- exp(pi)
      prob <- t(epi/rowSums(epi,na.rm=TRUE))

      list(prob=prob,epi=epi)
    }


    log.likelihood <- function(samples, X,Z,Q,eta,beta,zeta, lambda...){
      #pprint(t(zeta)[1:min(ncol(zeta),18),])
      #Sys.sleep(.5)

      cp <- calc.prob(samples,X,Z,Q,eta,beta,zeta)
      cp$prob[cp$prob<1e-8] <- 1e-8

      ll <- log(row.select(cp$prob,samples))
      ll <- sum(ll,na.rm=TRUE)
      if(D>0)
        ll <- ll - lambda*L/D*sqrt(sum(beta^2))
      if(!is.finite(ll))
        browser()

      ll
    }
    LL <- function(params,s,D,X,Z,Q,P, ...){
      S <- dim(X)[3]
      #L <- length(unique(s[!is.na(s)]))
      K <- ncol(Z)
      J <- dim(X)[2]
      theta <- decode.params(params=params,S,J,D,K,P)
      log.likelihood(s,X,Z,Q,eta=theta$eta,beta=theta$beta,zeta=theta$zeta,...)
    }
    -LL(param,s,D,X,Z,Q,P)
  }

  # Gradient function for optimization
  g <- function(param){
    row.select <- function(x,i)
      x[cbind(i,1:ncol(x))]
    col.select <- function(x,i)
      x[cbind(1:nrow(x),i)]

    aselect <- function(x,i){
      matrix(t(sapply(1:dim(x)[1],function(j) x[j,,i[j]])),ncol=dim(x)[2])
    }

    decode.params <- function(params,S,J,D,K,P){#holdouts here
      eta <- c(rep(0,J),params[1:((S-1)*J)])
      beta <- params[(S-1)*J + (1:(D*K))]
      eta <- matrix(eta,nrow=S,ncol=J,byrow=TRUE)
      if(D>0){
        beta <- matrix(beta,nrow=D,ncol=K)
      }else{
        beta <- matrix(NA,nrow=0,ncol=K)
      }
      zeta <- matrix(0,nrow=nrow(P),ncol=ncol(P))
      #zeta[!is.na(P)] <- 0
      zeta[which(P==1)] <- params[((S-1)*J + D*K) + (1:sum(P,na.rm=TRUE))]
      list(eta=eta,beta=beta,zeta=zeta)
    }
    encode.params <- function(eta,beta,zeta,P){#holdouts here
      params <- c()
      # Add etas
      if(length(eta)>0){
        J <- ncol(eta)
        params <- as.vector(matrix(eta,nrow=1))
        params <- params[(J+1):length(params)]
      }
      # Add betas
      if(length(beta)>0)
        params <- c(params,as.vector(matrix(beta,nrow=1)))

      #Add zetas
      params <- c(params, as.vector(zeta[which(P==1)]))
      params
    }
    calc.prob <- function(s,X,Z,Q,eta,beta,zeta){

      D <- nrow(beta)
      S <- ncol(zeta)
      L <- length(s[!is.na(s)])

      pi <- matrix(0,nrow=dim(Q)[1],ncol=S)

      # eta/X part
      if(ncol(eta)>0){
        for(a in 1:S)
          pi[,a] <- X[,,a,drop=FALSE]%*%eta[a,,drop=FALSE]
      }

      # beta/Z part
      if(D>0){
        l <- dim(Z)[1]
        G <- matrix(0,nrow=l,ncol=nrow(beta))
        for(a in 1:S){
          for(d in 1:D){
            G[(d+1):l,d] <- Z[1:(l-d),,a]%*%t(beta[d,,drop=FALSE])
          }
          delta <- sapply(1:D,function(d){ # should be LxD
            s1 <- c(rep(NA,d),s[1:(length(s)-d)])
            eq<-(s1==a)*1
            eq[is.na(eq)] <- 0
            eq
          })
          pi[,a] <- pi[,a]+rowSums(G*delta,na.rm=TRUE)
        }
      }

      # zeta/Q part
      z <- zeta
      z[is.na(z)] <- 0
      Qz <- Q%*%z
      Qzna <- Q%*%(1*is.na(zeta))
      Qz[which(Qzna>0)] <- NA
      pi <- pi + Qz #r$Q%*%t(zeta) #new!

      epi <- exp(pi)
      prob <- t(epi/rowSums(epi,na.rm=TRUE))

      list(prob=prob,epi=epi)
    }
    log.lik.grad <- function(samples, X,Z,Q,eta,beta,zeta,
                             beta.bounds=c(-Inf,Inf), ...){

      s <- samples
      D <- nrow(beta)
      S <- ncol(zeta)
      L <- length(s[!is.na(s)])

      cp <- calc.prob(s,X,Z,Q,eta,beta,zeta)
      p <- row.select(cp$prob,s)
      epi <- cp$epi

      eg <- eta*0
      zg <- zeta*0
      bg <- beta*0

      if(dim(eta)[2]>0){
        for(i in 1:dim(eg)[1])
          for(j in 1:dim(eg)[2])
            eg[i,j] <- sum( X[,j,i]*(1*(s==i)) - X[,j,i]*epi[,i]/rowSums(epi,na.rm=TRUE),na.rm=TRUE)
      }

      if(dim(zeta)[2]>0){
        for(i in 1:dim(zg)[2])
          for(j in 1:dim(zg)[1])
            zg[j,i] <- sum( Q[,j]*(1*(s==i)) - Q[,j]*epi[,i]/rowSums(epi,na.rm=TRUE), na.rm=TRUE)
      }

      if(D>0){
        ZYi <- aselect(Z,s)
        for(d in 1:dim(bg)[1]){
          delta <- sapply(1:S,function(a){
            s1 <- c(rep(NA,d),s[1:(length(s)-d)])
            eq<-(s1==a)*1
            eq[is.na(eq)] <- 0
            eq
          })
          dYi <- col.select(delta,s)
          # browser()
          for(j in 1:dim(bg)[2]){
            bg[d,j] <- sum( ZYi[,j,drop=TRUE]*dYi - rowSums(Z[,j,]*delta *epi) / rowSums(epi,na.rm=TRUE), na.rm=TRUE)
            #Regularization
            if(sum(beta)^2>0){
              bg[d,j] <- bg[d,j] - 2*lambda*L/D*beta[d,j]/sqrt(sum(beta^2))
            }
            #Bounds
            #bg[beta<beta.bounds[1] & bg<0] <- 0
            #bg[beta>beta.bounds[2] & bg>0] <- 0
          }
        }

      }

      list(beta=bg,eta=eg,zeta=zg)
    }
    LLgrad <- function(params,s,D,X,Z,Q,P, ...){
      S <- dim(X)[3]
      K <- ncol(Z)
      J <- dim(X)[2]
      theta <- decode.params(params=params,S,J,D,K,P)
      grads <- log.lik.grad(s,X,Z,Q,eta=theta$eta,beta=theta$beta,zeta=theta$zeta, beta.bounds=beta.bounds, ...)
      encode.params(grads$eta,grads$beta,grads$zeta,P)
    }
    -LLgrad(param,s,D,X,Z,Q,P)
  }

  cat(paste0("Fitting..."))

  params <- encode.params(eta0,beta0,zeta0,P)
  gr <- g
  if(num.grad)
    gr <- NULL
  if(parallel & cores>1){
    cl <- makeCluster(cores)     # set the number of processor cores
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    opt <- optimParallel(params,f,gr=gr,hessian=FALSE,method='L-BFGS-B',
                         parallel=list(forward=TRUE))
    stopCluster(cl)
  }else{
    opt <- optim(params,f,gr=gr,hessian=FALSE,method='L-BFGS-B')
  }
  theta <- decode.params(opt$par,S,J,D,K,P)

  cat(paste0("Done\n"))

  eta<- theta$eta
  rownames(eta) <- ns
  colnames(eta) <- dimnames(X)[[2]]

  beta <- theta$beta
  if(nrow(beta)>0)
    rownames(beta) <- paste0('Lag',1:nrow(beta))
  colnames(beta) <- colnames(Z)

  zeta <- t(theta$zeta)
  rownames(zeta) <- ns
  colnames(zeta) <- paste0('Seq',1:ncol(zeta))

  if('boundtype' %in% names(list(...)))
    if(list(...)$boundtype=='hard'){
      beta[beta<beta.bounds[1]] <- beta.bounds[1]
      beta[beta>beta.bounds[2]] <- beta.bounds[2]
    }

  cp <- calc.prob(s,X,Z,Q,eta,beta,t(zeta))
  prob <- cp$prob

  fit <- list(s=s,
              X=X,Z=Z,Q=Q,P=P,
              S=S,D=D,J=J,K=K,N=N,
              L = length(s[!is.na(s)]),
              k = length(opt$par),

              deviance = 2*opt$value,
              AIC = 2*opt$value +2*length(opt$par),
              BIC = 2*opt$value +log(length(s[!is.na(s)]))*length(opt$par),

              loglik=-opt$value,

              eta=eta, beta=beta, zeta=zeta,
              prob = prob,
              lambda=lambda, lambda2=lambda2,
              params = opt$par, opt = opt
  )
  #sp2 <- colSums(prob^2)

  #fit$hessian <- optimHess(encode.params(eta,beta,zeta,P),f=f)#opt$hessian

  # try({
  #   theta.se <- decode.params(sqrt(diag(solve(opt$hessian))),S,J,D,K,P)
  #   fit$eta.se <- theta.se$eta
  #   fit$beta.se <- theta.se$beta
  #   fit$zeta.se <- t(theta.se$zeta)
  #   dimnames(fit$eta.se) <- dimnames(fit$eta)
  #   dimnames(fit$beta.se) <- dimnames(fit$beta)
  #   dimnames(fit$zeta.se) <- dimnames(fit$zeta)
  # })
  fit
}

seq.hessian <- function(r){

  p <- r$params
  k <- 1:length(p)
  hk <- decode.params(k,r$S,r$J,r$D,r$K,r$P)
  hk$zeta[hk$zeta==0] <- NA

  cp <- calc.prob(r$s,r$X,r$Z,r$Q,r$eta,r$beta,t(r$zeta))

  hessian <- matrix(0,nrow=length(r$params),ncol=length(r$params))
  #beta hessian
  if(r$D>0){
    for(d1 in 1:dim(r$beta)[1]){
      for(d2 in 1:dim(r$beta)[1]){
        s1 <- c(rep(NA,d1),r$s[1:(length(r$s)-d1)])
        s2 <- c(rep(NA,d2),r$s[1:(length(r$s)-d2)])
        delta<-(s1==s2)*1
        delta[is.na(delta)] <- 0
        #if(dim(r$Z)[2]==1){
        ZYd1 <- aselect(r$Z,s1)
        ZYd2 <- aselect(r$Z,s2)
        #}else{
        #  ZYd1 <- t(aselect(r$Z,s1))
        #  ZYd2 <- t(aselect(r$Z,s2))
        #}
        ps1 <- row.select(cp$prob,s1)
        ps2 <- row.select(cp$prob,s2)
        for(j in 1:dim(r$beta)[2]){
          for(k in 1:dim(r$beta)[2]){
            hessian[hk$beta[d1,j],hk$beta[d2,k]] <- sum( ZYd1[,j,drop=TRUE]*ZYd2[,k,drop=TRUE]*ps1*(delta - ps2), na.rm=TRUE)
          }
        }
      }
    }
  }
  #hessian[hk$beta,hk$beta]

  #eta hessian
  if(dim(r$eta)[2]>0){
    for(i1 in 1:dim(r$eta)[1])
      for(i2 in 1:dim(r$eta)[1])
        for(j1 in 1:dim(r$eta)[2])
          for(j2 in 1:dim(r$eta)[2])
            hessian[hk$eta[i1,j1],hk$eta[i2,j2]] <- sum( r$X[,j1,i1]*r$X[,j2,i2]*cp$prob[i1,]*(1*(i1==i2) - cp$prob[i2,]),na.rm=TRUE)
  }
  #image(hessian[hk$eta,hk$eta])

  #zeta hessian
  if(dim(r$zeta)[2]>0){
    for(i1 in 1:dim(r$zeta)[1])
      for(i2 in 1:dim(r$zeta)[1])
        for(j1 in 1:dim(r$zeta)[2])
          for(j2 in 1:dim(r$zeta)[2])
            hessian[hk$zeta[j1,i1],hk$zeta[j2,i2]] <- sum( r$Q[,j1]*r$Q[,j2]*cp$prob[i1,]*(1*(i1==i2) - cp$prob[i2,]),na.rm=TRUE)
  }
  #zk <- as.vector(hk$zeta)
  #zk <- zk[!is.na(zk)]
  #hessian[zk,zk]

  #beta-eta hessian
  if(r$D>0){
    if(dim(r$eta)[2]>0){
      for(d1 in 1:dim(r$beta)[1]){
        s1 <- c(rep(NA,d1),r$s[1:(length(r$s)-d1)])
        ZYd1 <- aselect(r$Z,s1)#t(aselect(r$Z,s1))
        ps1 <- row.select(cp$prob,s1)
        for(i1 in 1:dim(r$eta)[1]){
          for(j1 in 1:dim(r$eta)[2]){
            for(j in 1:dim(r$beta)[2]){
              hessian[hk$beta[d1,j],hk$eta[i1,j1]] <- sum( ZYd1[,j,drop=TRUE]*r$X[,j1,i1]*ps1*(1*(s1==i1) - cp$prob[i1,]), na.rm=TRUE)
              hessian[hk$eta[i1,j1],hk$beta[d1,j]] <- hessian[hk$beta[d1,j],hk$eta[i1,j1]]
            }
          }
        }
      }
    }
  }

  #beta-zeta hessian
  if(r$D>0){
    if(dim(r$zeta)[2]>0){
      for(d1 in 1:dim(r$beta)[1]){
        s1 <- c(rep(NA,d1),r$s[1:(length(r$s)-d1)])
        ZYd1 <- aselect(r$Z,s1)
        ps1 <- row.select(cp$prob,s1)
        for(i1 in 1:dim(r$zeta)[1]){
          for(j1 in 1:dim(r$zeta)[2]){
            for(j in 1:dim(r$beta)[2]){
              hessian[hk$beta[d1,j],hk$zeta[j1,i1]] <- sum( ZYd1[,j,drop=TRUE]*r$Q[,j1]*ps1*(1*(s1==i1) - cp$prob[i1,]), na.rm=TRUE)
              hessian[hk$zeta[j1,i1],hk$beta[d1,j]] <- hessian[hk$beta[d1,j],hk$zeta[j1,i1]]
            }
          }
        }
      }
    }
  }
  #eta-zeta hessian
  if(dim(r$eta)[2]>0){
    if(dim(r$zeta)[2]>0){
      for(i1 in 1:dim(r$eta)[1])
        for(i2 in 1:dim(r$zeta)[1])
          for(j1 in 1:dim(r$eta)[2])
            for(j2 in 1:dim(r$zeta)[2]){
              hessian[hk$eta[i1,j1],hk$zeta[j2,i2]] <- sum( r$X[,j1,i1]*r$Q[,j2]*cp$prob[i1,]*(1*(i1==i2) - cp$prob[i2,]),na.rm=TRUE)
              hessian[hk$zeta[j2,i2],hk$eta[i1,j1]] <- hessian[hk$eta[i1,j1],hk$zeta[j2,i2]]
            }
    }
  }

  zk <- as.vector(hk$zeta)
  zk <- zk[!is.na(zk)]
  hessian[is.na(hessian[,zk])] <- 0
  hessian[is.na(hessian[zk,])] <- 0

  hessian
}


infer.seq <- function(r){
  hessian <- seq.hessian(r)
  theta.se <- decode.params(sqrt(diag(solve(hessian))),r$S,r$J,r$D,r$K,r$P)
  dimnames(theta.se$eta) <- dimnames(r$eta)
  dimnames(theta.se$beta) <- dimnames(r$beta)
  dimnames(theta.se$zeta) <- dimnames(t(r$zeta))

  r$eta.se = theta.se$eta
  r$beta.se = theta.se$beta
  r$zeta.se = t(theta.se$zeta)
  r
}

##############
decode.hessian <- function(r,partial=FALSE){
  p <- r$params
  k <- 1:length(p)
  hk <- decode.params(k,r$S,r$J,r$D,r$K,r$P)
  if(partial){
    list(etabeta = r$hessian[c(hk$eta,hk$beta),c(hk$eta,hk$beta)])
  }else{
    list(eta = r$hessian[hk$eta,hk$eta],
         beta = r$hessian[hk$beta,hk$beta],
         zeta = r$hessian[hk$zeta,hk$zeta])
  }
}
diagnose.sequence.model <- function(r){
  pp <- pointwise.probability(r)
  list(zeros=sum(pp==0,na.rm=TRUE),
       actives=sum(pp>0 & pp<1,na.rm=TRUE),
       ones=sum(pp==1,na.rm=TRUE),
       nas= sum(is.na(pp)),
       exp.nas =(r$N-1)*r$D
  )
}

# Censor probability
pointwise.probability <- function(r,censor=FALSE){

  cp <- calc.prob(r$s,r$X,r$Z,r$Q,r$eta,r$beta,t(r$zeta))
  prob <- cp$prob
  p <- row.select(prob,r$s)
  # if(D>0){
  #   G <- sapply(1:D,function(delta)
  #     c(rep(1,delta),gamma.vec(r$s,delta,gamma[delta,],p0)))
  #   l <- l*apply(G,1,prod)
  #   if(censor)
  #     l[l>1] <- 1
  # }

  p
}
#################
