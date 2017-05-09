logpost <- function(sigma, p, r, Z, N, S, nonzero){
  lp <- Z*log(1 + (1/sigma-1)*(1-p)^r)
  lp <- lp + Z*log(sigma)
  lp <- lp + (N-Z)*log(1-sigma)
  lp <- lp + (N-Z)*r*log(1-p)
  lp <- lp + S*log(p)
  lp <- lp - log(r)/2
  lp <- lp - (N-Z)*lgamma(r)
  lp <- lp + sum(lgamma(r+nonzero))
  return(lp)
}

sample <- function(r.mean, p.mean, sigma.mean, r.var, p.var, sigma.var){
  sigma.star <- rtruncnorm(1, 0, 1, sigma.mean, sqrt(sigma.var))
  #rbeta(1, sigma.mean*sigma.N, (1-sigma.mean)*sigma.N)
  p.star <- rtruncnorm(1, 0, 1, p.mean, sqrt(p.var))
  #rbeta(1, p.mean*p.N, (1-p.mean)*p.N)
  r.star <- rtruncnorm(1, 0, Inf, r.mean, sqrt(r.var))
  return(c(sigma.star, p.star, r.star))
}

g <-function(sigma, p, r, 
             sigma.mean, p.mean, r.mean, 
             sigma.var, p.var, r.var){
    return(dtruncnorm(sigma, 0, 1, sigma.mean, sqrt(sigma.var))*
           dtruncnorm(p, 0, 1, p.mean, sqrt(p.var))*
           dtruncnorm(r, 0, Inf, r.mean, sqrt(r.var)))
}


zinb.sampler <- function(df, year, event_type, chainlen, 
                         r.var=1, p.var=1, sigma.var=1, 
                         burnin=0, thinning=1, nchains=1, random=TRUE){
  X <- df[df$YEAR>year&df$EVENT_TYPE==event_type,"DEATHS_DIRECT"]
  N <- length(X)
  Z <- sum(X==0)
  S <- sum(X)
  nonzero <- X[X!=0]
  B <- chainlen
  b <- burnin
  
  r.chains <- c()
  sigma.chains <- c()
  p.chains <- c()
  ar <- c()
  
  for(j in 1:nchains){
    r.array <- rep(0, B)
    sigma.array <- rep(0, B)
    p.array <- rep(0, B)
    ar.array <- rep(0, B)
    posts <- rep(0, B)
    
    if (random==FALSE){
      r.array[1] <- seq(0,3,length.out=(nchains+2))[j+1]
      p.array[1] <- seq(0,1,length.out=(nchains+2))[j+1]
      sigma.array[1] <- seq(0,1,length.out=(nchains+2))[j+1]
    }
    else {
      r.array[1] <- runif(1, 0, 1)
      p.array[1] <- runif(1, 0, .5)
      sigma.array[1] <- runif(1, .5, 1)
    }
    print(c(r.array[1], p.array[1], sigma.array[1]))
    posts[1] <- logpost(sigma.array[1], p.array[1], r.array[1], Z, N, S, nonzero)
  
    for(i in 2:chainlen){
      sampled <- sample(r.array[i-1], p.array[i-1], sigma.array[i-1], 
                        r.var, p.var, sigma.var)
      r.star <- sampled[3]
      p.star <- sampled[2]
      sigma.star <- sampled[1]
      ## update r
      lpost <-logpost(sigma.star, p.star, r.star, Z, N, S, nonzero)
      reject.prob <- exp(lpost-posts[i-1])
      mh <- g(sigma.array[i-1], p.array[i-1], r.array[i-1], 
              sigma.star, p.star, r.star, 
              sigma.var, p.var, r.var)/
            g(sigma.star, p.star, r.star,
              sigma.array[i-1], p.array[i-1], r.array[i-1], 
              sigma.var, p.var, r.var)
      #print(c(reject.prob, mh))
      u <- runif(1,0,1)
      if(u < min(reject.prob*mh, 1)){
        r.array[i] <- r.star
        p.array[i] <- p.star
        sigma.array[i] <- sigma.star
        ar.array[i] <- 1
        posts[i] <- lpost
      }
      else{
        r.array[i] <- r.array[i-1]
        p.array[i] <- p.array[i-1]
        sigma.array[i] <- sigma.array[i-1]  
        posts[i] <- posts[i-1]
      }
    }
    r.chains <- c(r.chains, r.array[seq(b+1, B, by=thinning)])
    sigma.chains <- c(sigma.chains, sigma.array[seq(b+1, B, by=thinning)])
    p.chains <- c(p.chains, p.array[seq(b+1, B, by=thinning)])
    ar <- c(ar, ar.array[seq(1, B)])
  }
  return(list(sigma=sigma.chains, 
              p=p.chains, 
              r=r.chains, 
              ar=ar))
}

flood.zinb <- zinb.sampler(storm_details, 1995, "Flash Flood", 200000, 
                           r.var=.00002, p.var=.00001, sigma.var=.001, 
                           burnin=0, thinning=1, nchains=1, random=TRUE)
thinning <- 100
B <- length(flood.zinb$p)
b <- 0
print(mean(flood.zinb$sigma[seq(b+1, B, by=thinning)]))
print(mean(flood.zinb$p[seq(b+1, B, by=thinning)]))
print(mean(flood.zinb$r[seq(b+1, B, by=thinning)]))
print(mean(flood.zinb$ar))


par(mfrow=c(1,3))
plot(density(flood.zinb$p[seq(b+1, B, by=thinning)]))
plot(density(flood.zinb$sigma[seq(b+1, B, by=thinning)]))
plot(density(flood.zinb$r[seq(b+1, B, by=thinning)]))
par(mfrow=c(3,1))
plot(cumsum(flood.zinb$p)/1:B, type='l')
plot(cumsum(flood.zinb$sigma)/1:B, type='l')
plot(cumsum(flood.zinb$r)/1:B, type='l')
plot(1:B, flood.zinb$p)
plot(1:B, flood.zinb$sigma)
plot(1:B, flood.zinb$r)



