gathb <- function(Ns,L,l_start,samples){
  
  Ncol <- rep(as.character(Ns), each=L-(l_start-1))
  Lcol <- rep(l_start:L,length(Ns))
  vcol <- numeric((L-(l_start-1))*length(Ns))
  
  for(i in seq_along(Ns)){
    print(c("i",i))
    N <- Ns[i]
    
    for(l in l_start:L){
      print(c("l",l))
      record <- numeric(samples)
      
      for(j in seq_along(record)){
        N_str <- as.integer(substr(as.character(N),1,2))
        load(paste(paste(paste("beta-",N_str,"-results/L",sep=""),l,"task",j,sep="_"),".RData",sep=""))
        record[j] <- out
      }
      
      var <- var(record)
      vcol[(i-1)*(L-(l_start-1))+(l-(l_start-1))] <- log2(var)
    }
  }
  
  print(Ncol); print(Lcol); print(vcol)
  
  df <- data.frame("l" = Lcol, "v" = vcol, "N"=Ncol)
  p <- ggplot(df, aes(l,v, group=N, col=N)) + ylab("log2(variance)")
  p2 <- p + geom_point() + geom_smooth(method="lm",se=FALSE) 
  p2
}
