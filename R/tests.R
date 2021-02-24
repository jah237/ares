ex1_compare <- function(no_trials,no_particles,no_levels,save_seed,random_seed=1){

  n <- no_trials

  L_gen <- function(){
    rnbinom(1,3,0.5)
  }

  L_dens <- function(x){
    dnbinom(x,3,0.5)
  }

  lambda <- function(n){
    return(matrix(1.5,nrow=n,ncol=2))
  }

  xi <- function(x){return(min(x))}

  levels <- 2^(seq(1,length=no_levels,by=0.5))
  
  set.seed(random_seed)
  v <- matrix(0,nrow=no_trials,ncol=2)

  for(i in 1:no_trials){
    print(c("trial",i))
    trial <-  mlpf(L_gen, L_dens, N=50, d=2, lambda, xi=xi, z_A=0, levels, plot=FALSE,save_seed,random_seed)
    v[i,1] <- trial$p
    v[i,2] <- trial$L
  }

  name <- paste(paste("results/task",random_seed,sep="_"),".RData",sep="")
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
    L_samples[i] <- v[2]
  }

  return(list(p=p_estimates, L=L_samples))
}
