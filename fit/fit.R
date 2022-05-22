fit.rstan <- function( save.dir, file.name, model, dat, loadExisting=TRUE,
                       seed=1, chains=4, refresh=1600, iter=4000,
                       adapt_delta=0.8, max_treedepth=10, multi.C="" ){
  
  save.path.fit <- file.path(save.dir, paste(file.name, multi.C, "_stanfit.RData", sep=""))
  print(save.path.fit)
  start.time <- Sys.time()
  if ( !file.exists(save.path.fit) ){
    
    stanfit <- rstan::sampling(object=model, data=dat, 
                       iter=(iter+1500), warmup=1500, chain=chains, cores=chains,
                       refresh=refresh, seed=seed,
                       save_warmup=TRUE,
                       pars=c("V"), include=FALSE)
    # save model fit
    save( stanfit, file=save.path.fit )
  } else {
    if (loadExisting){
      print("already exists, loading stanfit...")
      load( file=save.path.fit )
      
    } else{
      print("already exists, return")
      stanfit<-NA
    }
  }
  end.time <- Sys.time()
  print(end.time - start.time)
  t.diff <- difftime( end.time, start.time, units=c("mins") )
  
  return( list(fit=stanfit, 
               t.diff=t.diff, path=save.path.fit, name=file.name, seed=seed ) )
}