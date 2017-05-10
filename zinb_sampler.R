# ---- zinbsampler ----
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

loglik <- function(sigma, p, r, Z, N, S, nonzero){
  lp <- Z*log(1 + (1/sigma-1)*(1-p)^r)
  lp <- lp + Z*log(sigma)
  lp <- lp + (N-Z)*log(1-sigma)
  lp <- lp + (N-Z)*r*log(1-p)
  lp <- lp + S*log(p)
  lp <- lp - sum(lgamma(nonzero+1))
  lp <- lp - (N-Z)*lgamma(r)
  lp <- lp + sum(lgamma(r+nonzero))
  return(lp)
}

g <-function(sigma, p, r, 
             sigma.mean, p.mean, r.mean, 
             sigma.var, p.var, r.var){
    return(dtruncnorm(sigma, 0, 1, sigma.mean, sqrt(sigma.var))*
           dtruncnorm(p, 0, 1, p.mean, sqrt(p.var))*
           dtruncnorm(r, 0, Inf, r.mean, sqrt(r.var)))
}


zinb.sampler <- function(df, event_type, chainlen, 
                         r.var=1, p.var=1, sigma.var=1, 
                         burnin=0, thinning=1){
  X <- df[df$EVENT_TYPE==event_type,"DEATHS_DIRECT"]
  N <- length(X)
  Z <- sum(X==0)
  S <- sum(X)
  nonzero <- X[X!=0]
  B <- chainlen
  b <- burnin
  r.sd <- sqrt(r.var)
  p.sd <- sqrt(p.var)
  sigma.sd <- sqrt(sigma.var)
  
  r.array <- rep(0, B)
  sigma.array <- rep(0, B)
  p.array <- rep(0, B)
  ar.array <- rep(0, B)
  posts <- rep(0, B)
    
  r.array[1] <- .5
  p.array[1] <- min(mean(X), .5)
  sigma.array[1] <- mean(X==0)/2
  #print(c(r.array[1], p.array[1], sigma.array[1]))
  posts[1] <- logpost(sigma.array[1], p.array[1], r.array[1], Z, N, S, nonzero)
  for(i in 2:chainlen){
    sigma.star <- rtruncnorm(1, 0, 1, sigma.array[i-1], sigma.sd)
    p.star <- rtruncnorm(1, 0, 1, p.array[i-1], p.sd)
    r.star <- rtruncnorm(1, 0, Inf, r.array[i-1], r.sd)
    lpost <-logpost(sigma.star, p.star, r.star, Z, N, S, nonzero)
    reject.prob <- exp(lpost-posts[i-1])*
                   dtruncnorm(sigma.array[i-1], 0, 1, 
                              sigma.star, sigma.sd)*
                   dtruncnorm(p.array[i-1], 0, 1, p.star, p.sd)*
                   dtruncnorm(r.array[i-1], 0, Inf, r.star, r.sd)/
                  (dtruncnorm(sigma.star, 0, 1, 
                              sigma.array[i-1], sigma.sd)*
                     dtruncnorm(p.star, 0, 1, p.array[i-1], p.sd)*
                     dtruncnorm(r.star, 0, Inf, r.array[i-1], r.sd))
    u <- runif(1,0,1)
    if(u < min(reject.prob, 1)){
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
  inds <- seq(b+1, B, by=thinning)
  return(list(sigma=sigma.array[inds], 
              p=p.array[inds], 
              r=r.array[inds], 
              ar=ar.array))
}

# ---- densityplot ----

density.plot <- function(theta, theta.density, col='black') {
  theta.mean <- mean(theta)
  theta.quantiles <- quantile(theta, probs=c(.025, .975))
  theta.lower <- min(which(theta.density$x >= theta.quantiles[1]))
  theta.upper <- max(which(theta.density$x < theta.quantiles[2]))
  theta.med <- max(which(theta.density$x<theta.mean))
  lines(c(theta.mean,theta.mean), 
        c(theta.density$y[theta.med], 0), lwd=.5, col=col)
  lines(c(theta.quantiles[1],theta.quantiles[1]), 
        c(theta.density$y[theta.lower], 0), lty=2, col=col)
  lines(c(theta.quantiles[2],theta.quantiles[2]), 
        c(theta.density$y[theta.upper], 0), lty=2, col=col)
}

double.density.plot<- function(theta1, theta2,
                               title, xlabel, ylabel, 
                               col=c('black', 'blue')) {
  theta1.density <- density(theta1)
  theta2.density <- density(theta2)
  xmax <- max(max(theta1), max(theta2))
  xmin <- min(min(theta1), min(theta2))
  ymax <- max(max(theta1.density$y), max(theta2.density$y))
  plot(theta1.density,
       col=col[1],
       main=title,
       xlab=xlabel,
       ylab=ylabel,
       xlim=c(xmin,xmax),
       ylim=c(0,ymax))
  lines(theta2.density, col=col[2])
  density.plot(theta1, theta1.density, col=col[1])
  density.plot(theta2, theta2.density, col=col[2])
}

# ---- evaluation ----
dic <- function(df, column, sigma, r, p){
  X <- df[df$EVENT_TYPE==column,"DEATHS_DIRECT"]
  N <- length(X)
  Z <- sum(X==0)
  S <- sum(X)
  B <- length(sigma)
  nonzero <- X[X!=0]
  l <- loglik(mean(sigma), mean(p), mean(r), Z, N, S, nonzero)
  d <- -2*l
  pdic2 <- 0 
  for(i in 1:B){
    pdic2 <- pdic2 + loglik(sigma[i], p[i], r[i], Z, N, S, nonzero)
  }
  return(-2*l + 2*(l - pdic2/B))
}

wiac <- function(df, column, sigma, r, p) {
  X <- as.matrix(df[df$EVENT_TYPE==column,"DEATHS_DIRECT"])
  l_mat <- rep(0, length(X)*length(sigma))
  dim(l_mat) <- c(length(X), length(sigma))
  #print(dim(l_mat))
  l_mat[X==0,] <- sigma + (1-sigma)*(1-p)^r
  l_mat[X!=0,] <- (1-sigma)*mapply(dnbinom, prob=p, size=r, 
                                   MoreArgs=list(x=X[X!=0]))
  lppd <- sum(log(rowMeans(l_mat)))
  pwaic <- sum(apply(log(l_mat),1,var))
  return(2*pwaic - 2*lppd)
}





