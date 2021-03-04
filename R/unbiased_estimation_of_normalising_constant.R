#normalised importance sampling
importance_sample <- function(N, phi, m_gen, q_eval, m_eval){

  sum <- 0

  for(i in 1:N){
    X <- m_gen()
    importance_weight <- q_eval(X)/m_eval(X)
    sum <- sum + phi(X)*importance_weight
  }

  return(sum)
}

#auto-normalised importance sampling
AN_importance_sample <- function(N, phi, m_gen, q_eval, m_eval){

  num <- 0
  denom <- 0

  for(i in 1:N){
    X <- m_gen
    unnormalised_importance_weight <- q_eval(X)/m_eval(X)
    num <- num + phi(X)*unnormalised_importance_weight
    denom <- denom + unnormalised_importance_weight
  }

  out <- num/denom
  return(out)
}

#multinomial respampling

ordered_uniform_sample <- function(N){

  exps <- rexp(N+1)

  cumsum_exps <- cumsum(exps)
  sum_exps <- sum(exps)

  out <- cumsum_exps[1:N]/sum_exps
  return(out)
}

multinomial_sample <- function(weights){

  N <- length(weights)
  U <- ordered_uniform_sample(N)

  A <- numeric(N)

  s <- weights[1]
  m <- 1

  for(n in 1:N){
    while(s < U[n]){
      m <- m + 1
      s <- s + weights[m]
    }
    A[n] <- m
  }

  return(A)
}

#particle_filter
generic_pf <- function(T_steps, N, M_list, G_list){

  M_0 <- M_list[[1]]
  G_0 <- G_list[[1]]

  X <- numeric(N)
  w <- numeric(N)

  for(n in 1:N){
    X[n] <- M_0()
    w[n] <- G_0(X[n])
  }

  W <- w/sum(w)

  for(t in 1:T_steps){

    M <- M_list[[t+1]]
    G <- G_list[[t+1]]

    A <- multinomial_sample(W)
    X_r <- X[A]

    for(n in 1:N){

      X[n] <- M(X_r[n])
      w[n] <- G(X_r[n],X[n])
    }

    W <- w/sum(w)
  }

  return(list(X=X,W=W))
}

norm_const <- function(T_steps, N, M_list, G_list){

  M_0 <- M_list[[1]]
  G_0 <- G_list[[1]]

  X <- numeric(N)
  w <- numeric(N)

  L <- 1

  for(n in 1:N){
    X[n] <- M_0()
    w[n] <- G_0(X[n])
  }

  L <- L*mean(w)

  W <- w/sum(w)

  for(t in 1:T_steps){

    M <- M_list[[t+1]]
    G <- G_list[[t+1]]

    A <- multinomial_sample(W)
    X_r <- X[A]

    for(n in 1:N){

      X[n] <- M(X_r[n])
      w[n] <- G(X_r[n],X[n])
    }

    L <- L*mean(w)

    W <- w/sum(w)
  }

  return(L)
}


#Algorithm 1. Discretised MLS
cdm_euler2 <- function(x, delta, a, b, reac_coord, plot){

  d <- length(x)

  x_rec <- x

  while((a <= reac_coord(x)) & (reac_coord(x) <= b)){
    old_x <- x
    x <- x + rnorm(d, mean = 0, sd = sqrt(delta))

    #if(plot==TRUE){lines(rbind(old_x,x))}
  }

  if(reac_coord(x) < a){
    return(list(survived=0))
  } else {
    return(list(survived=1,x=x))
  }
}

msplitting_euler3 <- function(N, d, lambda, z_A, levels, reac_coord, delta, delta_scale, L, plot, save_seed, id){

  m <- length(levels)
  survivors <- matrix(lambda(N), ncol = d)
  n_surv <- integer(m)

  if(plot == TRUE){
    points(survivors, pch=4)
  }

  for (i in 1:m){
    new_survivors <- matrix(0,nrow=0,ncol=d)

    for (j in 1:N){
      
      if(!exists(".Random.seed")){
        set.seed(NULL)
      }

      if(save_seed==TRUE){
        saved_seed <- .Random.seed
        if((j==1) & (i>1)){
          unlink(paste(paste("data/task",id,sep="_"),".RData",sep=""))
        } else if(j>1){
          unlink(paste(paste("data/task",id,sep="_"),".RData",sep=""))
        }

        name <- paste(paste("data/task",id,sep="_"),".RData",sep="")
        save(list=ls(), file=name)
      }

      trial <- cdm_euler2(survivors[j,], delta, z_A, levels[i], reac_coord, plot)

      if(trial$survived==1){
        new_survivors <- rbind(new_survivors,trial$x)

        if(plot==TRUE){
          lines(rbind(survivors[j,],trial$x))
        }
      }
    }

    if(plot==TRUE){
      points(new_survivors,pch=4)
    }

    n_surv[i] <- nrow(new_survivors)

    if(n_surv[i] == 0){
      return(0)
    }

    if(i < m){
      indices <- sample(1:n_surv[i], N, replace=TRUE)
      survivors <- new_survivors[indices,]
      delta <- delta*delta_scale[i]
    }
  }
  out <- prod(n_surv) / N^m
  return(out)
}

#Algorithm 3
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

    type <- 0

  } else {

    prob1 <- (w1-w_min)/sum(w1-w_min)
    prob2 <- (w2-w_min)/sum(w2-w_min)

    I1 <- sample(1:N, size=1, prob=prob1)
    I2 <- sample(1:N, size=1, prob=prob2)

    type <- 1
  }

  out <- list(I1=I1, I2=I2, type = type)
  return(out)
}

#Algorithm 2
cd_euler_coupled <- function(x1, x2, delta, z_a, z_b, xi, plot){

  #add plotting option

  sd <- sqrt(delta)

  sgn <- 0

  res1 <- 0
  res2 <- 0
  xi_tot <- 0

  while(res1==0){
    sgn <- (sgn + 1)%%2
    w <- rnorm(2, mean = 0, sd = sd)
    x1 <- x1 + w

    if(sgn == 1){
      w_old <- w
    } else {
      x2 <- x2 + w_old + w
    }

    if(xi(x1) < z_a){
      res1 <- -1
    } else if (xi(x1) > z_b){
      res1 <- 1
    }
  }

  if((z_a < xi(x2)) & (xi(x2) < z_b)){

    w <- rnorm(2, mean = 0, sd = sd)
    x2 <- x2 + w_old + w
    delta <- 2*delta
    sd <- sqrt(delta)

    while((z_a <= xi(x2)) & (xi(x2) <= z_b)){
      x2 <- x2 + rnorm(2, mean = 0, sd = sd)
    }
  }

  if(xi(x2) < z_a){
    res2 <- -1
  } else if (xi(x2) > z_b){
    res2 <- 1
  }

  return(list(x1=x1,x2=x2,results=(c(res1,res2)+1)/2))
}

#replace multinomial with
#stratified sampling

#Algorithm 4
coupled_splitting <- function(N, d, lambda, z_a, levels, delta, xi, L,  plot, save_seed, id){

  #add plotting option
  if(plot==TRUE){
    cols <- rep(c("green", "purple", "orange", "yellow"),20)
  }

  m <- length(levels)
  init_vals <<- matrix(lambda(N), ncol = d)
  x1_vals <- x2_vals <- asplit(init_vals, 1)
  sd <- sqrt(delta)

  prop_ind <- numeric(m)

  if(plot==TRUE){
    points(init_vals, pch=4)
  }

  x1_sample <- vector("list",N); x2_sample <- vector("list", N)
  Gx1 <- numeric(N); Gx2 <- numeric(N)
  survivors1 <- numeric(m); survivors2 <- numeric(m)

  for (i in 1:m){

    print(c("i",i))
    z_b <- levels[i]

    for (j in 1:N){

      if(!exists(".Random.seed")){
        set.seed(NULL)
      }

      if(save_seed==TRUE){
        saved_seed <- .Random.seed
        if((j==1) & (i>1)){
          unlink(paste(paste("data/task",id,sep="_"),".RData",sep=""))
        } else if(j>1){
          unlink(paste(paste("data/task",id,sep="_"),".RData",sep=""))
        }

        name <- paste(paste("data/task",id,sep="_"),".RData",sep="")
        save(list=ls(), file=name)
      }

      x1 <- x1_vals[[j]]; x2 <- x2_vals[[j]]
      M_sample <- cd_euler_coupled(x1, x2, delta, z_a, z_b, xi, plot)
      x1_sample[[j]] <- M_sample$x1; x2_sample[[j]] <- M_sample$x2
      Gx1[j] <- M_sample$results[1]; Gx2[j] <- M_sample$results[2]

      if(plot == TRUE){
        lines(rbind(x1_vals[[j]],x1_sample[[j]]), col=cols[j])
        lines(rbind(x2_vals[[j]],x2_sample[[j]]), col=cols[j])
        if(identical(M_sample$results,c(1,1))){
          points(rbind(M_sample$x1, M_sample$x2), pch=4)
        } else if(identical(M_sample$results,c(0,1))){
          points(rbind(M_sample$x1, M_sample$x2), pch=4, col=c("red","black"))
        } else if(identical(M_sample$results,c(1,0))){
          points(rbind(M_sample$x1, M_sample$x2), pch=4, col=c("black","red"))
        } else {
          points(rbind(M_sample$x1, M_sample$x2), pch=4, col="red")
        }
      }
    }

    survivors1[i] <- sum(Gx1); survivors2[i] <- sum(Gx2)

    if((survivors1[i] == 0) | (survivors2[i] == 0)){

      return(0)
    } else if(i < m) {

      indices1 <- numeric(N); indices2 <- numeric(N)
      no_ind <- 0
      for(k in 1:N){
        new_indices <- coupled_resampling(Gx1, Gx2)
        indices1[k] <- new_indices$I1; indices2[k] <- new_indices$I2
        no_ind <- no_ind + new_indices$type
      }
      #print("indices1"); print(indices1); print("indices2"); print(indices2)
      x1_vals <- x1_sample[indices1]; x2_vals <- x2_sample[indices2]
    }
    prop_ind[i] <- no_ind/N
  }
  proportion_independent <<- prop_ind
  out <- (prod(survivors1) - prod(survivors2))/N^m
  return(out)
}

#Algorithm 5
mlpf <- function(L_gen, L_mass, N, d, lambda, xi, z_A, levels, plot=FALSE, save_seed, id){

  if(plot == TRUE){

    #rectangular plot
    #m <- length(levels)-1
    #plot(1, type="n", xlab="", ylab="", xlim=20*c(-1,1), ylim=c(0,levels[m+1]), asp=1)
    #abline(h = levels, col="blue")
    #abline(h=z_A, col = "red")

    #circular plot
    #m <- length(levels)-1
    #plot(1, type="n", xlab="", ylab="", xlim=(2^m+1)*c(-1,1), ylim=(2^m+1)*c(-1,1), asp=1)
    #draw.circle(0,0,(2^m+1)-2^(m:1), border="blue")
    #draw.circle(0,0,2^m+1, border="red")

    #L_shaped plot
    m <- length(levels)
    plot(1, type="n", xlab="", ylab="", xlim=1.3*2^((m+1)/2)*c(0,1), ylim=2^((m+1)/2)*c(0,1), asp=1)
    segments(levels,levels,levels,2^6, col="blue")
    segments(levels,levels,2^6,levels, col="blue")
    segments(0,0,0,2^12,col="red"); segments(0,0,2^12,0,col="red")

  }

  L <- L_gen()

  delta <- 2^(-(L+4))

  if(L == 0){

    #estimate using an ordinary particle filter
    out <- msplitting_euler3(N, d, lambda, z_A, levels, xi, delta,
                             delta_scale=rep(sqrt(2),length(levels)-1),L,plot,
                             save_seed, id)/L_mass(L)
  } else {

    #estimate using coupled scheme
    out <- coupled_splitting(N, d, lambda, z_A, levels, delta, xi, L, plot, save_seed, id)/L_mass(L)
  }

  return(list(p=out, L=L))
}
