coupled_euler <- function(x1,x2,a,b,n,delta,N){
  
  for(k in 1:n){
    
    w <- rnorm(N)
    x1 <- x1 + a(x1)*delta + b(x1)*sqrt(delta)*w
    
    if(k %% 2 == 1){
      
      w_old <- w
      
    } else {
      
      x2 <- x2 + a(x2)*(2*delta) + b(x2)*sqrt(delta)*(w + w_old)
    }
  }
  
  return(list(x1=x1,x2=x2))
  
}

coupled_resampling3 <- function(w1, w2, N){
  
  wmin <- pmin(w1,w2)
  a <- sum(wmin)
  
  U <- runif(N)
  no_paired <- sum(U < a)
  no_independent <- N - no_paired
  
  I1 <- I2 <- integer(N)

  q <- wmin
  
  if(no_paired > 0){
    I1[1:no_paired] <- I2[1:no_paired] <- sample(N, no_paired, replace=TRUE, prob = q)
  }
  
  q1 <- w1 - wmin; q2 <- w2 - wmin
  if(no_independent > 0){
    I1[(no_paired+1):N] <- sample(N, no_independent, replace=TRUE, prob = q1)
    I2[(no_paired+1):N] <- sample(N, no_independent, replace=TRUE, prob = q2)
  }
  
  return(list(I1=I1,I2=I2))
}

test <- function(ell,id){
  
  set.seed(id)
  
  theta <- 1
  mu <- 0
  sigma <- 3
  
  a <- function(x){
    out <- theta*(mu-x)
    return(out)
  }
  
  b <- function(x){
    return(sigma)
  }
  
  N <- 50000
  
  x1 <- x2 <- rep(0,N)
  x_0 <- rep(0,N)
  w1 <- w2 <- x_0
  
  alpha <- 1
  
  V <- function(x){
    return(x)
  }
  
  V_0 <- V(x_0)
  
  G <- function(x,w){
    out <- exp(alpha*(V(x) - V(w)))
  }
  
  n <- 10
  h_ell <- 2^(-ell)
  
  r <- 10
  
  nc1 <- nc2 <- 1
  
  for(i in 1:n){
    G_1s <- G(x1,w1); G_2s <- G(x2,w2)
    Gbar_1 <- mean(G_1s); Gbar_2 <- mean(G_2s)  
    weights1 <- G_1s/(N*Gbar_1); weights2 <- G_2s/(N*Gbar_2)
    nc1 <- nc1*Gbar_1; nc2 <- nc2*Gbar_2
  
    new_indices <- coupled_resampling3(weights1, weights2, N)
    I1 <- new_indices$I1; I2 <- new_indices$I2
    
    x1 <- x1[I1]; x2 <- x2[I2]
    w1 <- w1[I1]; w2 <- w2[I2]
    
    xs <- coupled_euler(x1,x2,a,b,2^ell,h_ell,N)
    w1 <- x1; w2 <- x2
    x1 <- xs$x1; x2 <- xs$x2
  }
  
  reached1 <- (V(x1) > r); reached2 <- (V(x2) > r) 
  p1 <- nc1*mean(reached1*G(-w1,-x_0)); p2 <- nc2*mean(reached2*G(-w2,-x_0))
  
  out <- p1 - p2
  
  name <- paste(paste("beta-50-results/L",ell,"task",id,sep="_"),".RData",sep="")
  save(out, file=name)
  
  return(out)
}

