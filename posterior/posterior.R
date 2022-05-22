## ------------ plot posterior parameters ------------

plot.multiR <- function( fit, nMat ){
  
  iterations <- fit$metadata()$iter_sampling
  chains <- fit$num_chains()
  draws.df <- posterior::as_draws_df( fit$draws() )
  
  # read correlation and log likelihood for each iteration
  # and average over absolute value of chains
  log_lik <- matrix(ncol=nMat, nrow=iterations)
  R <- matrix(ncol=nMat, nrow=iterations)
  
  prog <- txtProgressBar(min=1, max=iterations, style=3)
  for (i in 1:iterations) {
    
    log.lik.i <- apply( as.matrix(
      posterior::subset_draws(draws.df, "log_lik", iteration=i)), 2, function(x) log(mean(exp(x))) )
    log_lik[i,] <- log.lik.i[1:nMat]
    
    rs <- apply( as.matrix(
      posterior::subset_draws(draws.df, "R", iteration=i)), 2, function(x) mean(abs(x)) )
    rs <- rs[1:(4*nMat)]
    R[i,] <- matrix(rs, ncol=nMat, byrow=TRUE)[2,]
    
    setTxtProgressBar(prog, i)
  }
  close(prog)
  
  log_lik <- apply( log_lik, 2, function(x) log(mean(exp(x))) )
  R <- apply( R, 2, mean )
  
  hpdi.r <- rethinking::HPDI(R)
  hpdi.log_lik <- rethinking::HPDI(log_lik)
  
  r.density <- plot( density(R), xlab="R[1,2]", main="Correlation density")
  rethinking::shade(density(R), c(hpdi.r[1], hpdi.r[2]))
  
  abline(v=mean(R), col="green")
  
  r.plot <- plot(R, log_lik, xlim=c(0,1))
  rect(xleft=hpdi.r[1], ybottom=hpdi.log_lik[1], xright=hpdi.r[2], ytop=hpdi.log_lik[2],
       lty="blank", col=rethinking::col.alpha("black",0.15))
  
  return( list(density=r.density, plot=r.plot) )
}


#' plot distribution of correlations of multiVar models
#' 
#' extract correlation parameter (R[,1,2]) and log likelihood
#'
plot.multiR.stanfit <- function( stanfit, nMat, chains=4, keep_log=TRUE ){
  
  # draw posterior samples from fit model
  samples <- rstan::extract( stanfit, permuted=FALSE )
  
  # read correlation and log likelihood for each iteration
  # and average over absolute value of chains
  log_lik <- vector(mode="numeric", length=nMat)
  R <- vector(mode="numeric", length=nMat)
  
  for (i in 1:nMat) {
    
    log_lik.i <- vector(mode="numeric", length=4)
    R.i <- vector(mode="numeric", length=4)

    for (j in 1:chains) {
      chain.i.j <- as.data.frame(samples[i,j,])
      
      par.name <- paste("log_lik[", i, "]", sep="")
      log_lik.i[j] <- chain.i.j[par.name,]
      par.name <- paste("R[", i, ",1,2]", sep="")
      R.i[j] <- abs(chain.i.j[par.name,])
    }
    log_lik[i] <- log( mean( exp(log_lik.i) ) )
    R[i] <- mean(R.i)
  }
  
  hpdi.r <- unlist(rethinking::HPDI(R))
  hpdi.log_lik <- rethinking::HPDI(log_lik)
  
  r.density <- plot( density(R), xlab="R[1,2]", main="Correlation density")
  rethinking::shade(density(R), c(hpdi.r[1], hpdi.r[2]))
  
  abline(v=mean(R), col="green")

  lik <- ifelse( rep(keep_log, nMat), log_lik, exp(log_lik) )
  r.plot <- plot(R, log_lik)
  rect(xleft=hpdi.r[1], ybottom=hpdi.log_lik[1], xright=hpdi.r[2], ytop=hpdi.log_lik[2],
       lty="blank", col=rethinking::col.alpha("black",0.15))

  return( list(density=r.density, plot=r.plot) )
}

plot.multiSigma <- function( fit, nMat, keep_log=TRUE ){
  
  iterations <- fit$metadata()$iter_sampling
  chains <- fit$num_chains()
  draws.df <- posterior::as_draws_df( fit$draws() )
  
  # read correlation and log likelihood for each iteration
  # and average over absolute value of chains
  log_lik <- matrix(ncol=nMat, nrow=iterations)
  sigma <- matrix(ncol=nMat, nrow=iterations)
  
  prog <- txtProgressBar(min=1, max=iterations, style=3)
  for (i in 1:iterations) {
    
    log.lik.i <- apply( as.matrix(
      posterior::subset_draws(draws.df, "log_lik", iteration=i)), 2, function(x) log(mean(exp(x))) )
    log_lik[i,] <- log.lik.i[1:nMat]
    
    sigmas <- apply( as.matrix(
      posterior::subset_draws(draws.df, "sigma", iteration=i)), 2, mean )
    sigma[i,] <- sigmas[1:nMat]
    
    setTxtProgressBar(prog, i)
  }
  close(prog)
  
  log_lik <- apply( log_lik, 2, function(x) log(mean(exp(x))) )
  sigma <- apply( sigma, 2, mean )
  
  hpdi.sigma <- rethinking::HPDI(sigma)
  hpdi.log_lik <- rethinking::HPDI(log_lik)
  
  sigma.density <- plot( density(sigma), xlab="R[1,2]", main="Correlation density")
  rethinking::shade(density(sigma), c(hpdi.sigma[1], hpdi.sigma[2]))
  
  abline(v=mean(sigma), col="green")
  
  #lik <- ifelse( rep(keep_log, nMat), log_lik, exp(log_lik) )
  sigma.plot <- plot(sigma, log_lik, xlim=c(0,1))
  rect(xleft=hpdi.sigma[1], ybottom=hpdi.log_lik[1], xright=hpdi.sigma[2], ytop=hpdi.log_lik[2],
       lty="blank", col=rethinking::col.alpha("black",0.15))
  
  return( list(density=sigma.density, plot=sigma.plot) )
}



## ------------ plot predictions ------------

plot.dens <- function( samples, x.breaks=seq(from=0, to=1, by=0.25) ){
  
  pi.89 <- rethinking::PI( samples )
  pi.50 <- rethinking::PI( samples, prob=0.5 )
  
  sample.dens <- density( samples )
  dens.df <- with( sample.dens, data.frame( x, y ) )
  # start plot
  dens.plot <- qplot( x, y, data=dens.df, geom="line" ) + 
    scale_x_continuous( breaks=x.breaks )
  
  # add 89% PI
  if ( pi.89[1]==pi.89[2] ){
    dens.plot <- dens.plot + geom_vline( xintercept=pi.89[1], color="red", size=1.5 )
  } else{
    dens.plot <- dens.plot + 
      geom_ribbon( data=subset(dens.df, x>pi.89[1] & x<pi.89[2]), 
                   aes(ymax=y), ymin=0, fill="red", colour=NA, alpha=0.5 )
  }
  # add 50% PI
  if ( pi.50[1]==pi.50[2] ){
    dens.plot <- dens.plot + geom_vline( xintercept=pi.50[1], color="blue", size=1.5 )
  } else{
    dens.plot <- dens.plot + 
      geom_ribbon( data=subset(dens.df, x>pi.50[1] & x<pi.50[2]), 
                   aes(ymax=y), ymin=0, fill="blue", colour=NA, alpha=0.5 )
  }
  # add mean
  dens.plot <- dens.plot + geom_vline( xintercept=mean(samples) )
  
  return( dens.plot )
}


plot.sim.scat <- function( fit, nMat, n, r, traits, lang.name, trait.name.1, trait.name.2,
                           ylim=c(-0.5, 1.5), trait.values=c(1/3, 2/3, 1) ){
  
  iterations <- fit$metadata()$iter_sampling
  draws.df <- posterior::as_draws_df( fit$draws() )
  
  x.list <- vector(mode="list", length=n*r)
  ll.mat <- matrix(nrow=nMat, ncol=iterations)
  for (i in 1:(n*r)) {
    x.list[[i]] <- matrix(nrow=nMat, ncol=iterations)
  }
  
  prog <- txtProgressBar(min=1, max=iterations, style=3)
  
  for (iteration in 1:iterations) {
    
    x_array <- apply( as.matrix( posterior::subset_draws(draws.df, "x_pred", iteration=iteration)), 
                      2, function(x) mean(x) ) # average over chains
    
    x_array <- x_array[1:(length(x_array)-3)] # avoid meta-data .chain, .iteration, .draw
    x.mat <- matrix(data=x_array, ncol=nMat, byrow=TRUE)
    
    for (i in 1:(n*r)) {
      x.list[[i]][,iteration] <- x.mat[i,]
    }
    
    ll_array <- apply( as.matrix( posterior::subset_draws(draws.df, "log_lik", iteration=iteration)), 
                       2, function(x) log(mean(exp(x))) ) # average over chains
    ll.mat[,iteration] <- ll_array[1:(length(ll_array)-3)] # avoid meta-data .chain, .iteration, .draw
    
    setTxtProgressBar(prog, iteration)
  }
  close(prog)
  
  ll.means <- apply( ll.mat, 1, function(x) log(mean(exp(x))) )
  ll.means <- sort( ll.means )
  
  x.seq <- seq(from=min(ll.means), to=max(ll.means), length.out=nMat )
  
  for (i in 1:n) {
    
    x.i.mat <- x.list[[i]]
    # x.seq <- seq(from=1, to=nMat, by=1)
    y.means <- apply(x.i.mat, 1, mean)
    hpdis <- apply(x.i.mat, 1, rethinking::HPDI)
    pis <- apply(x.i.mat, 1, rethinking::PI)
    
    plot(x=ll.means, y=y.means, ylim=ylim,
         main=paste(lang.name, ".", i, ": ", trait.name.1, sep=""), 
         xlab="log.lik", ylab="posterior predictions for x")
    
    lines( x.seq, rep(trait.values[1], nMat), lty=1, lwd=1 )
    lines( x.seq, rep(trait.values[2], nMat), lty=1, lwd=1 )
    lines( x.seq, rep(trait.values[3], nMat), lty=1, lwd=1 )
    lines( ll.means , y.means, col="green" )
    lines( x.seq , rep(traits[i], nMat), col="red", lty=5, lwd=1 )
    text( x=max(ll.means), y=traits[i], labels=trait.name.1, col="red" )
    
    rethinking::shade(pis, x.seq)
  }
  
  
  for (i in 1:n) {
    j <- i + n
    
    x.i.mat <- x.list[[j]]
    # x.seq <- seq(from=1, to=nMat, by=1)
    y.means <- apply(x.i.mat, 1, mean)
    hpdis <- apply(x.i.mat, 1, rethinking::HPDI)
    pis <- apply(x.i.mat, 1, rethinking::PI)
    
    plot(x=ll.means, y=y.means, ylim=ylim,
         main=paste(lang.name, ".", i, ": ", trait.name.1, sep=""), 
         xlab="log.lik", ylab="posterior predictions for x")
    
    lines( x.seq, rep(trait.values[1], nMat), lty=1, lwd=1 )
    lines( x.seq, rep(trait.values[2], nMat), lty=1, lwd=1 )
    lines( x.seq, rep(trait.values[3], nMat), lty=1, lwd=1 )
    lines( ll.means , y.means, col="green" )
    lines( x.seq , rep(traits[j], nMat), col="red", lty=5, lwd=1 )
    text( x=max(ll.means), y=traits[j], labels=trait.name.1, col="red" )
    
    rethinking::shade(pis, x.seq)
  }
}


plot.sim.scat.McElreath <- function( fit, nMat, n, traits, lang.name, trait.name,
                                     ylim=c(-0.5, 1.5), trait.values=c(1/3, 2/3, 1) ){
  
  iterations <- fit$metadata()$iter_sampling
  draws.df <- posterior::as_draws_df( fit$draws() )
  
  x.list <- vector(mode="list", length=n)
  ll.mat <- matrix(nrow=nMat, ncol=iterations)
  for (i in 1:n) {
    x.list[[i]] <- matrix(nrow=nMat, ncol=iterations)
  }
  
  prog <- txtProgressBar(min=1, max=iterations, style=3)
  
  for (iteration in 1:iterations) {
    
    x_array <- apply( as.matrix( posterior::subset_draws(draws.df, "x_pred", iteration=iteration)), 
                      2, function(x) mean(x) ) # average over chains
    
    x_array <- x_array[1:(length(x_array)-3)] # avoid meta-data .chain, .iteration, .draw
    x.mat <- matrix(data=x_array, ncol=nMat, byrow=TRUE)
    
    for (i in 1:n) {
      x.list[[i]][,iteration] <- x.mat[i,]
    }
    
    ll_array <- apply( as.matrix( posterior::subset_draws(draws.df, "log_lik", iteration=iteration)), 
                       2, function(x) log(mean(exp(x))) ) # average over chains
    ll.mat[,iteration] <- ll_array[1:(length(ll_array)-3)] # avoid meta-data .chain, .iteration, .draw
    
    setTxtProgressBar(prog, iteration)
  }
  close(prog)
  
  ll.means <- apply( ll.mat, 1, function(x) log(mean(exp(x))) )
  ll.means <- sort( ll.means )
  
  x.seq <- seq(from=min(ll.means), to=max(ll.means), length.out=nMat )
  
  for (i in 1:n) {
    
    x.i.mat <- x.list[[i]]
    y.means <- apply(x.i.mat, 1, mean)
    hpdis <- apply(x.i.mat, 1, rethinking::HPDI)
    pis <- apply(x.i.mat, 1, rethinking::PI)
    
    plot(x=ll.means, y=y.means, ylim=ylim,
         main=paste(lang.name, ".", i, ": ", trait.name, sep=""), 
         xlab="log.lik", ylab="posterior predictions for x")
    
    lines( x.seq, rep(trait.values[1], nMat), lty=1, lwd=1 )
    lines( x.seq, rep(trait.values[2], nMat), lty=1, lwd=1 )
    lines( x.seq, rep(trait.values[3], nMat), lty=1, lwd=1 )
    lines( ll.means , y.means, col="green" )
    lines( x.seq , rep(traits[i], nMat), col="red", lty=5, lwd=1 )
    text( x=max(ll.means), y=traits[i], labels=trait.name, col="red" )
    
    rethinking::shade(pis, x.seq)
  }
}

## ------------ simulate draws from posterior and plot from those draws ------------ 

#' returns a list of posterior predictions based on posterior draws from binomial model
#' 
#' @param fit fitted CmdStanMCMC object
#' @param iterations number of iterations to use for predictions
#' @param dat data list used for fitting the model
#' @return list, where each element is a matrix 
#'               with one column for each language and one row for each iteration;
#'               the list contains one such matrix for each cov-Matrix in the data
#'
sim.x.r <- function( fit, iterations, M, N, chains ){
  
  p <- posterior::as_draws_df( rethinking::inv_logit( rstan::extract(fit, pars="p", permuted=FALSE) ) )

  pp.List <- lapply( 1:(2*N), function(n) matrix(ncol=M, nrow=iterations) )
  
  prog <- txtProgressBar(min=1, max=(2*N), style=3)
  
  for (n in 1:(2*N)) {
    for (m in 1:(M)) {
      p.Mat <- matrix(nrow=iterations, ncol=chains)
      for (chain in 1:chains) {
        curr <- paste("p[", m, ",", n, "]", sep="")
        p.Mat[,chain] <- unlist(posterior::subset_draws( p, chain=chain )[curr])[1:iterations]
      }
      p_mn <- apply(p.Mat, 2, mean)
      pp_mn <- rbinom( n=iterations, size=1, prob=p_mn )
      pp.List[[n]][,m] <- pp_mn
    }
    setTxtProgressBar(prog, n)
  }
  close(prog)
  
  return(pp.List)
}

plot.sim.x <- function( pp.Mat ){
  for ( m in 1:ncol(pp.Mat) ) {
    
    m.mean <- mean(pp.Mat[,m])
    m.dens <- density(pp.Mat[,m])
    m.HPDI <- rethinking::PI(pp.Mat[,m])
  }
}