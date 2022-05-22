## run main script to fit pairwise trait models for given family

# load functions ----
source( 'paths.R', encoding="utf8" )
source( file.path(fit.dir, 'fit.R'), encoding="utf8" )
source( file.path(dataprep.dir, 'dataprep.R'), encoding="utf8" )
source( file.path(post.dir, 'diagnostics.R'), encoding="utf8" )

library(rstan)
library(bridgesampling)

# prepare models and metaparameters ----
chains <- 14
iter <- 1500

print("parse models...")
model.univ <- rstan::stan_model( 
  file.path( models.dir, 'Binom-Y-nonCentered-fullLang-multiC-univ.stan' ),
  model_name="Y-binom-univ-dep_" )
print("...universal done")
model.univ.indep <- rstan::stan_model( 
  file.path( models.dir, 'Binom-Y-nonCentered-fullLang-multiC-univ-indep.stan' ),
  model_name="Y-binom-univ-indep_" )
print("...universal-indep done")
model.lin <- rstan::stan_model( 
  file.path( models.dir, 'Binom-Y-nonCentered-fullLang-multiC-lin.stan' ),
  model_name="Y-binom-lin_" )
print("...lineage done")

lang <- "full"
trait.names <- c( "AN", "PN", "ND", "NG", "NNum", "VO", "NRc", "VS" )
t.pairs <- c()

# ---- main loop ----

for (trait.name.1 in trait.names){
  for (trait.name.2 in trait.names){
    pair = paste( trait.name.1, trait.name.2 )
    
    if ( (trait.name.1 != trait.name.2) & !(pair %in% t.pairs) ){
      print( paste("----", trait.name.1, "-", trait.name.2, "----") )
      
      t.pairs = c( t.pairs, pair, paste(trait.name.2, trait.name.1) )
        
      # ---- prepare data ----
        
      # setup save
      save.dir <- file.path(fit.dir, lang)
      if (!dir.exists(save.dir)) { dir.create(save.dir) }
      model.name <- "Y-binom"
      data.name <- paste( lang, "_" ,trait.name.1, "-", trait.name.2, sep="" )
      
      # get data:
      lang.data <- get.data.full( c(trait.name.1, trait.name.2), avg=TRUE, other=0, n_avg=5, getC=TRUE )
      dat <- list(
        N = lang.data$N,
        x = lang.data$traits,
        M = 5,
        C = lang.data$Cs,
        # Hyperparameters
        eta = 0.8,
        lambda = 1,
        mu_z = 0,
        sigma_z = 1.5
      )
      remove(lang.data)
      
      # ---- prepare metrics ----
      
      save.dir.looc <- file.path(post.dir, "looc", lang)
      save.dir.waic <- file.path(post.dir, "waic", lang)
      save.dir.bf <- file.path(post.dir, "bridge", lang)
      if (!dir.exists(save.dir.looc)) { dir.create(save.dir.looc) }
      if (!dir.exists(save.dir.waic)) { dir.create(save.dir.waic) }
      if (!dir.exists(save.dir.bf)) { dir.create(save.dir.bf) }
      
      # ============== universally dependent binomial model ==============
        
      print( "-- universally dependent --" )
      print(model.univ@model_name)
      file.name <- file.path( paste(model.name, "-univ-dep_", data.name, "-0", sep="") )
      
      fit.rstan( save.dir=save.dir, file.name=file.name, model=model.univ, dat=dat, 
                 loadExisting=FALSE, seed=1,
                 chains=chains, iter=iter, refresh=2000,
                 multi.C="_5" )
      print("...fit")

      # ---- eval ----
      get.bridge( model=model.univ, lang.name=lang, 
                  trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                  dat=dat, cores=chains,
                  multi.C="_5" )
      print("...bridge")
      
      get.LOOCs( lang.name=lang, model_name=model.univ@model_name,
                 trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                 cores=chains,
                 multi.C="_5" )
      print("...looc")

      # ============== universally independent binomial model ==============
      
      print( "-- universally independent --" )
      print(model.univ.indep@model_name)
      file.name <- file.path( paste(model.name, "-univ-indep_", data.name, "-0", sep="") )
      
      fit.rstan( save.dir=save.dir, file.name=file.name, model=model.univ.indep, dat=dat, 
                 loadExisting=FALSE, seed=1,
                 chains=chains, iter=iter, refresh=2000,
                 multi.C="_5" )
      print("...fit")

      # ---- eval ----
      get.bridge( model=model.univ.indep, lang.name=lang, 
                  trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                  dat=dat, cores=chains,
                  multi.C="_5" )
      print("...bridge")
      
      get.LOOCs( lang.name=lang, model_name=model.univ.indep@model_name,
                 trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                 cores=chains,
                 multi.C="_5" )
      print("...looc")

      # ============== lineage specific binomial model ==============
        
      print("")
      print( "- lineage specific" )
      print(model.lin@model_name)
      file.name <- file.path( paste(model.name, "-lin_", data.name, "-0", sep="") )
      
      fit.rstan( save.dir=save.dir, file.name=file.name, model=model.lin, dat=dat, 
                 loadExisting=FALSE, seed=1,
                 chains=chains, iter=iter, refresh=2000,
                 multi.C="_5" )
      print("...fit")

      # ---- eval ----
      get.bridge( model=model.lin, lang.name=lang, 
                  trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                  dat=dat, cores=chains,
                  multi.C="_5" )
      print("...bridge")
      
      get.LOOCs( lang.name=lang, model_name=model.lin@model_name,
                 trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                 cores=chains,
                 multi.C="_5" )
      print("...looc")
      
      # ============== WAIC ==============
      print(" -- WAIC --")
      
      if (!dir.exists(file.path(post.dir, "waic", lang))){
        dir.create(file.path(post.dir, "waic", lang))
      }
      
      model_name.univ="Y-binom-univ-dep_"
      model_name.indep="Y-binom-univ-indep_"
      model_name.lin="Y-binom-lin_"
      
      print("universal - lineage")
      get.WAICs( model.1.name=model_name.univ,
                 model.2.name=model_name.lin,
                 lang.name=lang,
                 trait.name.1=trait.name.1, trait.name.2=trait.name.2,
                 multi.C="_5", dep="")
      
      print("universal dependent - independent")
      get.WAICs( model.1.name=model_name.univ,
                 model.2.name=model_name.indep,
                 lang.name=lang,
                 trait.name.1=trait.name.1, trait.name.2=trait.name.2,
                 multi.C="_5", dep="")
      print("...done")
      
      remove( dat )
    }
  }
}
