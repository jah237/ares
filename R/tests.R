ex1_compare <- function( no_trials,no_particles,no_levels,save_seed,random_seed=1){

  n <- no_trials

 L_gen <- function(){
   rgeom(1,0.5)
 }

 L_dens <- function(x){
   dgeom(x,0.5)
 }

 #L_gen <- function(){return(L)}
 #L_dens <- function(L){return(1)}

  #L_gen <- function(){
  #  rgeom(1,0.7)
  #}
  #L_dens <- function(){
  #  rgeom(1, 0.7)
  #}
 

  lambda <- function(n){
    return(matrix(1.5,nrow=n,ncol=2))
  }

  xi <- function(x){return(min(x))}

  levels <- 2^(seq(1,length=no_levels,by=0.5))

  set.seed(random_seed)
  v <- matrix(0,nrow=no_trials,ncol=2)

  for(i in 1:no_trials){
    trial <-  mlpf(L_gen, L_dens, N=no_particles, d=2, lambda, xi=xi, z_A=0, levels, plot=FALSE,save_seed,random_seed)
    v[i,1] <- trial$p
    v[i,2] <- trial$L
  }

  name <- paste(paste("mls-results/task",random_seed,sep="_"),".RData",sep="")
  save(v, file=name)

  return(v)
}

gather_results <- function(folder_name){

  filenames <- list.files(folder_name, full.names=TRUE)
  L <- length(filenames)
  p_estimates <- numeric(L)
  L_samples <- numeric(L)

  for(i in 1:L){

    load(filenames[[i]])
    p_estimates[i] <- v[1]
    if(v[1]==0){print(filenames[[i]])}
    L_samples[i] <- v[2]
  }

  return(list(p=p_estimates, L=L_samples))
}

complete_filter <- function(id){

  task <- id

  data_file <- paste("data/task_", task, ".RData", sep="")
  load(data_file)

  print(survivors1); print(survivors2)

  new_filename <- paste("continued-data/task_", task, ".RData", sep="")

  .Random.seed <<- saved_seed

  jx <- j
  ix <- i

  print(c("L",L))
  print(c("ix",ix,"jx",jx))

  if(L == 0){

    print("L=0")

    for (j in jx:N){

      print(c("i",i,"j",j))

      if(save_seed==TRUE){

        saved_seed <- .Random.seed
        save(list=ls(), file=new_filename)
      }

      trial <- cdm_euler2(survivors[j,], delta, z_A, levels[i], reac_coord, plot)

      if(trial$survived==1){
        new_survivors <- rbind(new_survivors,trial$x)

        if(plot==TRUE){
          lines(rbind(survivors[j,],trial$x))
        }
      }
    }

    if(ix < m){
      for (i in ix:m){
        new_survivors <- matrix(0,nrow=0,ncol=d)

        for (j in 1:N){

          print(c("i",i,"j",j))

          if(save_seed==TRUE){

            saved_seed <- .Random.seed
            save(list=ls(), file=new_filename)
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
    }
  } else {

    print("L>0")

    for (j in jx:N){

      print(c("i",i,"j",j))

      if(save_seed==TRUE){

        saved_seed <- .Random.seed
        save(list=ls(), file=new_filename)
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

    survivors1[ix] <- sum(Gx1); survivors2[ix] <- sum(Gx2)

      if((survivors1[ix] == 0) | (survivors2[ix] == 0)){

        return(0)
      } else if(ix < m) {

        indices1 <- numeric(N); indices2 <- numeric(N)
        no_ind <- 0
        for(k in 1:N){
          new_indices <- coupled_resampling(Gx1, Gx2)
          indices1[k] <- new_indices$I1; indices2[k] <- new_indices$I2
          no_ind <- no_ind + new_indices$type
        }
        x1_vals <- x1_sample[indices1]; x2_vals <- x2_sample[indices2]
      }
      prop_ind[ix] <- no_ind/N

    if(ix < m){
      for (i in (ix+1):m){

      z_b <- levels[i]

      for (j in 1:N){

        print(c("i",i,"j",j))

        if(save_seed==TRUE){

          saved_seed <- .Random.seed
          save(list=ls(), file=new_filename)
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
        x1_vals <- x1_sample[indices1]; x2_vals <- x2_sample[indices2]
      }
      prop_ind[i] <- no_ind/N
      }
    }
    
    print(survivors1); print(survivors2)    

    proportion_independent <<- prop_ind
    out <- prod(survivors1/N) - prod(survivors2/N)
    s1 <<- survivors1; s2 <<- survivors2
  }

  v <- c(out,L)

  results_name <- paste(paste("survival-results/task",task,sep="_"),".RData",sep="")
  save(list=ls(), file=results_name)

  return(v)
}
