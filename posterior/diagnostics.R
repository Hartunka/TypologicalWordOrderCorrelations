load.stanfit <- function(save.dir, file.name, cVar=""){
  load( file=file.path(save.dir, paste(file.name, cVar, "_stanfit.RData", sep="")) )
  return(stanfit)
}

diag.param.dens <- function( fit, corr ){
  
  logit.plot <- bayesplot::mcmc_areas(stanfit, pars=c("z[1]", "z[2]", corr)) + 
    labs(title="logit scale")
  
  prob.plot <- bayesplot::mcmc_areas(stanfit, pars=c("z[1]", "z[2]", corr), 
                                     transformations=rethinking::inv_logit, prob=0.89 ) + 
    labs(title="probability scale",
         subtitle="t = rethinking::inv_logit")
  return( list(logit=logit.plot, prob=prob.plot) )
}

get.bridge <- function( model, lang.name, trait.name.1, trait.name.2, dat, cores=1, multi.C="", trait.var="-0" ){
  
  # load fit model
  save.dir <- file.path( fit.dir, lang )
  data.name <- paste( lang.name, "_" ,trait.name.1, "-", trait.name.2, sep="" )
  file.name <- file.path( paste(model@model_name, data.name, trait.var, sep="") )
  
  if(file.exists(file.path(post.dir, "bridge", lang.name, paste(file.name, multi.C, "_bridge-warp3.RData", sep="")))){
    print( file.path(post.dir, "bridge", lang.name, paste(file.name, multi.C, "_bridge-warp3.RData", sep="")) )
    load( file=file.path(post.dir, "bridge", lang.name, paste(file.name, multi.C, "_bridge-warp3.RData", sep="")) )
  } else {
    print( file.path(post.dir, "bridge", lang.name, paste(file.name, multi.C, "_bridge-warp3.RData", sep="")) )
    stanfit <- load.stanfit( save.dir, file.name, cVar=multi.C  )
    
    # create empty stanfit.model
    cat("load empty model...\n")
    #stan.model <- stan_model( file=file.path( models.dir, 'categorical', m.name ),
    #                             model_name=dep )
    stan.model.fit <- sampling( object=model, data=dat, iter=0, chains=1, verbose=FALSE )
    
    # bridgesampling
    cat("bridgesampling...\n")
    bridge <- bridge_sampler( samples=stanfit, stanfit_model=stan.model.fit,
                              silent=TRUE, cores=cores, method="warp3" )
    cat("save...\n")
    save(bridge, file=file.path(post.dir, "bridge", lang.name, paste(file.name, multi.C, "_bridge-warp3.RData", sep="")) )
  }
  return(bridge)
}

#' compute leave-one-out cross validation approximation for a model
#' 
#' @param fit fit model to be evaluated (stanfit)
#' @return fit.loo : list containing evaluations for model returned from \code{loo::loo}
#' @export
#' 
eval.loo <- function( fit, cores=1 ){
  
  cat( "...loo...\n" )
  
  log_lik <- loo::extract_log_lik(fit, merge_chains=FALSE)
  r_eff <- loo::relative_eff( exp(log_lik), cores=cores )
  fit.loo <- loo::loo( log_lik, r_eff=r_eff, cores=cores )
  
  return( fit.loo )
}

get.LOOCs <- function( lang.name, model_name, trait.name.1, trait.name.2, cores=1, multi.C="" ){
  
  # load fit model
  save.dir <- file.path( fit.dir, lang )
  data.name <- paste( lang.name, "_" ,trait.name.1, "-", trait.name.2, sep="" )
  file.name <- file.path( paste( model_name, data.name, "-0", multi.C, sep="" ) )
  
  cat(file.name, "\n")
  if(file.exists(file.path(post.dir, "looc", lang.name, paste(file.name, "_looc.RData", sep="")))){
    
    load( file=file.path(post.dir, "looc", lang.name, paste(file.name, "_looc.RData", sep="")) )
    
  } else {
    
    stanfit <- load.stanfit(save.dir, file.name, cVar="" )
    # LOOC
    loo <- eval.loo( stanfit, cores=cores )
    save( loo, file=file.path( post.dir, "looc", lang.name, paste(file.name, multi.C, "_looc.RData", sep="") ) )
  }
  return(loo)
}

get.WAICs <- function( model.1.name, model.2.name, lang.name, trait.name.1, trait.name.2, multi.C="", dep="_", model.3.name=NA ){
  
  # load fit model
  save.dir <- file.path( fit.dir, lang )
  data.name <- paste( lang.name, "_" ,trait.name.1, "-", trait.name.2, sep="" )
  
  file.name.waic <- file.path( paste( model.1.name, data.name, "-0", multi.C, "_waic.RData", sep="" ) )
  
  cat(file.path(post.dir, "waic", lang.name, file.name.waic), "\n")
  if(file.exists(file.path(post.dir, "waic", lang.name, file.name.waic))){
    
    load( file=file.path(post.dir, "waic", lang.name, file.name.waic) )
    
  } else {
    file.name.1 <- file.path( paste( model.1.name, data.name, "-0", sep="" ) )
    file.name.2 <- file.path( paste( model.2.name, data.name, "-0", sep="" ) )
    
    cat("loading", file.name.1, "...\n")
    stanfit.1 <- load.stanfit(save.dir, file.name.1, cVar=multi.C)
    cat("loading", file.name.2, "...\n")
    stanfit.2 <- load.stanfit(save.dir, file.name.2, cVar=multi.C)
    
    if (is.na(model.3.name)){
      # WAIC
      print("waic...")
      waic <- rethinking::compare( stanfit.1, stanfit.2 )
    } else {
      file.name.3 <- file.path( paste( model.3.name, "_", data.name, "-0", sep="" ) )
      cat("loading", file.name.3, "...\n")
      stanfit.3 <- load.stanfit(save.dir, file.name.3, cVar=multi.C)
      
      print("waic...")
      waic <- rethinking::compare( stanfit.1, stanfit.2, stanfit.3 )
    }
    
    save( waic, file=file.path( post.dir, "waic", lang.name, file.name.waic ) )
  }
  return( waic )
}
