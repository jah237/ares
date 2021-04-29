ips <- function(N, lambda, init_weights, V, G, K, update_weights, 
                                                    A, n, estimator, alpha, V_0=0){
  X <- lambda(N); weights <- init_weights; norm_constants <- numeric(n)
  
  for(i in 1:n){
    
    G_N <- G(alpha, V, X, weights)
    eta <- mean(G_N); norm_constants[i] <- eta
    p <- G_N/(N*eta)
    new_indices <- sample(1:N, replace = TRUE, prob = p)

    X <- X[new_indices]; weights <- weights[new_indices]
    
    weights <- update_weights(alpha, X, weights, G, V)
    
    X <- K(X)
  }
  
  V_fin <- V(X)
  verdict <- (V_fin > A)
  nc <- prod(norm_constants)
  V_0s <- rep(V_0, N)

  out <- estimator(nc, verdict, weights, alpha, V, V_0s)
  return(out)
}

ips_with_histories <- function(N, lambda, init_weights, V, G, K, update_weights, 
                                                      A, n, estimator, alpha, V_0=0){
  
  X <- lambda(N); weights <- init_weights; norm_constants <- numeric(n)
  X_hist <- matrix(0,n,N); X_hist[1, ] <- X
  
  for(i in 1:n){
    
    G_N <- G(alpha, V, X, weights)
    eta <- mean(G_N); norm_constants[i] <- eta
    p <- G_N/(N*eta)
    new_indices <- sample(1:N, replace = TRUE, prob = p)
    
    X <- X[new_indices]; weights <- weights[new_indices]
    
    weights <- update_weights(alpha, X, weights, G, V)
    X <- K(X)
  }
  
  V_fin <- V(X)
  verdict <- (V_fin > A)
  nc <- prod(norm_constants)
  V_0s <- rep(V_0, N)
  
  out <- estimator(nc, verdict, weights, alpha, V, V_0s)
  return(out)
}

#Toy model from del Moral & Garnier, w/ Alg 1

ips_ex1 <- function(){
  
  N <- 20000
  A = 20
  n = 15
  beta = 0.15
  
  lambda <- function(N){
    out <- rep(0,N)
    return(out)
  }
  
  init_weights <- rep(1,N)
  
  V <- function(x){return(x)}
  
  G <- function(beta, V, x, weights){
    out <- exp(beta*V(x))
    return(out)
  }
  
  K <- function(x){
    out <- x + rnorm(N)
    return(out)
  }
  
  update_weights <- function(beta, X, weights, G, V){
    out <- weights/G(beta, V, X, weights)
    return(out)
  }
  
  estimator <- function(nc, verdict, weights, alpha, V, V_0s){
    out <- nc * mean(verdict*weights)
    return(out)
  }
  
  out <- ips(N, lambda, init_weights, V, G, K, update_weights, 
                                                      A, n, estimator, beta)
  return(out)
}

#w/ alg 2

ips_ex2 <- function(){
  
  N <- 20000
  A = 20
  n = 15
  alpha = 1
  
  x_0 <- 0
  
  lambda <- function(N){
    out <- rep(0,N)
    return(out)
  }
  
  V <- function(x){return(x)}
  init_weights <- rep(x_0, N)
  V_0 <- V(x_0)
  
  G <- function(alpha, V, x, weights){
    out <- exp(alpha*(V(x)-V(weights)))
    return(out)
  }
  
  K <- function(x){
    out <- x + rnorm(N)
    return(out)
  }
  
  update_weights <- function(alpha, X, weights, G, V){
    out <- X
    return(out)
  }
  
  estimator <- function(nc, verdict, weights, alpha, V, V_0s){
    out <- nc * mean(verdict*exp(-alpha*(V(weights) - V_0s)))
  }
  
  out <- ips(N, lambda, init_weights, V, G, K, update_weights, 
                                                        A, n, estimator, alpha)
  return(out)
}

ips_ex3 <- function(){
  
  N <- 20000
  A = 5
  n = 8
  alpha = 1
  x_0 <- 0
  
  lambda <- function(N){
    out <- rep(0,N)
    return(out)
  }
  
  V <- function(x){return(x)}
  init_weights <- rep(x_0, N)
  V_0 <- V(x_0)
  
  G <- function(alpha, V, x, weights){
    out <- exp(alpha*(V(x)-V(weights)))
    return(out)
  }
  
  K <- function(x){
    
    n_steps <- 2^3
    delta <- 1/n_steps
    sigma <- 3
    theta <- 1
    
    for(i in 1:n_steps){
      x <- x - theta*delta*x + rnorm(N,0,sd=sqrt(delta)*sigma)
    }
    return(x)
  }
  
  update_weights <- function(alpha, X, weights, G, V){
    out <- X
    return(out)
  }
  
  estimator <- function(nc, verdict, weights, alpha, V, V_0s){
    out <- nc * mean(verdict*exp(-alpha*(V(weights) - V_0s)))
  }
  
  out <- ips(N, lambda, init_weights, V, G, K, update_weights, 
                                                    A, n, estimator, alpha, V_0)
  return(out)
}

coupled_resampling <- function(Gx1, Gx2){
  
  N <- length(Gx1)
  
  w1 <- Gx1/sum(Gx1); w2 <- Gx2/sum(Gx2)
  w_min <- pmin(w1,w2)
  
  a <- sum(w_min)
  
  U <- runif(1)
  
  if(U < a){
    
    prob=w_min/a
    
    I1 <- sample(1:N, size=1, prob=prob)
    I2 <- I1
    
  } else {
    
    prob1 <- (w1-w_min)/sum(w1-w_min)
    prob2 <- (w2-w_min)/sum(w2-w_min)
    
    I1 <- sample(1:N, size=1, prob=prob1)
    I2 <- sample(1:N, size=1, prob=prob2)
  }
  
  out <- list(I1=I1, I2=I2)
  return(out)
}

ips_coupled <- function(N, lambda, init_weights, V, G, K_coupled, update_weights, 
                                                        A, n, estimator, alpha, V_0=0){
  
  
  X1 <- X2 <- lambda(N); 
  weights1 <- weights2 <- init_weights; 
  
  norm_constants1 <- norm_constants2 <- numeric(n)
  
  for(i in 1:n){
    
    GX1 <- G(alpha, V, X1, weights1); GX2 <- G(alpha, V, X2, weights2)
    norm_constants1[i] <- mean(GX1); norm_constants2[i] <- mean(GX2) 
    
    
    new_indices <- coupled_resampling(GX1, GX2)
    I1 <- new_indices$I1; I2 <- new_indices$I2
    
    X1 <- X1[I1]; weights1 <- weights1[I1]
    X2 <- X2[I2]; weights2 <- weights2[I2]
    
    weights1 <- update_weights(alpha, X1, weights1, G, V)
    weights2 <- update_weights(alpha, X2, weights2, G, V)
    
    Xs <- K_coupled(X1,X2)
    X1 <- Xs$X1; X2 <- Xs$X2
  }
  
  V_fin1 <- V(X1); V_fin2 <- V(X2)
  verdict1 <- (V_fin1 > A); verdict2 <- (V_fin2 > A)
  nc1 <- prod(norm_constants1); nc2 <- prod(norm_constants2)
  V_0s <- rep(V_0, N)
  
  p1 <- estimator(nc1, verdict1, weights1, alpha, V, V_0s)
  p2 <- estimator(nc2, verdict2, weights2, alpha, V, V_0s)
  out <- (p1 - p2)
  return(out)
}


ips_exact <- function(L, L_mass, N, lambda, init_weights, V, G, 
                                K, K_coupled, update_weights, A, n, estimator, alpha, V_0=0){
  
  if(L == 0){
    
    p_est <- ips(N, lambda, init_weights, V, G, K, update_weights, 
                                                    A, n, estimator, alpha, V_0=0)
    out <- p_est/L_mass(L)
    
  } else {
    
    #estimate using coupled scheme
    p_est <- ips_coupled(N, lambda, init_weights, V, G, K_coupled, update_weights, 
                                                              A, n, estimator, alpha, V_0=0)
    out <- p_est/L_mass(L)
  }
  
  return(list(p=out, L=L))
}

ips_ex4 <- function(){

  N <- 20000
  A <- 5
  n <- 8
  alpha <- 1
  x_0 <- 0

  L_gen <- function(){
    rgeom(1,0.5)
  }

  L_mass <- function(x){
    dgeom(x,0.5)
  }
  L <- L_gen()

  n_steps <- 2^(L+3)
  delta <- 1/n_steps
  
  lambda <- function(N){
    out <- rep(0,N)
    return(out)
  }

  V <- function(x){return(x)}
  init_weights <- rep(x_0, N)
  V_0 <- V(x_0)

  G <- function(alpha, V, x, weights){
    out <- exp(alpha*(V(x)-V(weights)))
    return(out)
  }

  sigma <- 3
  theta <- 1

  K <- function(x){
    
    for(i in 1:n_steps){
      x <- x - theta*delta*x + rnorm(N,0,sd=sqrt(delta)*sigma)
    }
    return(x)
  }

  K_coupled <- function(x1,x2){

    sgn <- 0

    for(i in 1:n_steps){
      sgn <- (sgn + 1)%%2
      w <- rnorm(N, mean = 0, sd = sqrt(delta)*sigma)
      x1 <- x1 - theta*delta*x1 + w
      
      if(sgn == 1){
        w_old <- w
      } else {
        x2 <- x2 - theta*delta*x1 + (w_old + w)
      }
    }

    return(list(X1=x1,X2=x2))
  }
  
  update_weights <- function(alpha, X, weights, G, V){
    out <- X
    return(out)
  }

  estimator <- function(nc, verdict, weights, alpha, V, V_0s){
    out <- nc * mean(verdict*exp(-alpha*(V(weights) - V_0s)))
  }

  out <- ips_exact(L, L_mass, N, lambda, init_weights, V, G, K, K_coupled,
                          update_weights, A, n, estimator, alpha, V_0)
  return(out)
}