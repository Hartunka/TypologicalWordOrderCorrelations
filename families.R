## run main script to fit pairwise trait models for given family

# load functions ----
source('paths.R', encoding="utf8" )
source( file.path(fit.dir, 'fit.R'), encoding="utf8" )
source( file.path(dataprep.dir, 'dataprep.R'), encoding="utf8" )
source( file.path(post.dir, 'diagnostics.R'), encoding="utf8" )

library(rstan)
library(bridgesampling)

# prepare models and metaparameters ----
chains <- 14
iter <- 1500

print("parse models...")
model.dep <- rstan::stan_model( 
  file.path( models.dir, 'Binom-multiC-Y-nonCentered-main.stan' ),
  model_name="Y-binom-dep_" )
print("...dep done")
model.indep <- rstan::stan_model( 
  file.path( models.dir, 'Binom-indep-multiC-Y-nonCentered-main.stan' ),
  model_name="Y-binom-indep_" )
print("...indep done")

lang.fam.names <- c( "Afro-Asiatic", "Atlantic-Congo", "Austronesian", "Indo-European",
                     "Nuclear_Trans_New_Guinea", "Pama-Nyungan", "Pano-Tacanan", "Sino-Tibetan",
                     "Algic", "Athapaskan-Eyak-Tlingit", "Austroasiatic", "Cariban", "Central_Sudanic",
                     "Chibchan", "Dravidian", "Mande", "Mongolic", "Nakh-Daghestanian",
                     "Nuclear_Torricelli", "Otomanguean", "Salishan", "Sepik", "Siouan", "Tupian",
                     "Turkic", "Uto-Aztecan", "Arawakan", "Chibchan", "Gunwinyguan", "Nilotic",
                     "Pano-Tacanan", "Sko", "Surmic", "Tai-Kadai", "Tucanoan",
                     "Uralic")
trait.names <- c( "AN", "PN", "ND", "NG", "NNum", "VO", "NRc", "VS" )
t.pairs <- c()

# ---- main loop ----

for (trait.name.1 in trait.names){
  for (trait.name.2 in trait.names){
    pair = paste( trait.name.1, trait.name.2 )
    
    if ( (trait.name.1 != trait.name.2) & !(pair %in% t.pairs) ){
      print( paste("----", trait.name.1, "-", trait.name.2, "----") )
      
      t.pairs = c( t.pairs, pair, paste(trait.name.2, trait.name.1) )
      
      for ( i in 1:length(lang.fam.names) ){
        
        # ---- prepare data ----
          
        lang <- lang.fam.names[i]
        
        print( paste("---", lang, "---") )
        
        # setup save
        save.dir <- file.path("fit", "fit", lang)
        if (!dir.exists(save.dir)) { dir.create(save.dir) }
        model.name <- "Y-binom"
        data.name <- paste( lang, "_" ,trait.name.1, "-", trait.name.2, sep="" )
        
        # get data:
        lang.data <- get.data( lang, c(trait.name.1, trait.name.2), bin=TRUE, other=0, avg=TRUE, n_avg=5 )
        dat <- list(
          N = lang.data$N,
          x = lang.data$traits,
          M = 5,
          C = lang.data$C,
          dMat = designMatrix(lang.data$N, 2),
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
        
        # ============== dependent binomial model ==============
          
        print( "-- dependent --" )
        print(model.dep@model_name)
        file.name <- file.path( paste(model.name, "-dep_", data.name, "-0", sep="") )
        
        fit.rstan( save.dir=save.dir, file.name=file.name, model=model.dep, dat=dat, 
                   loadExisting=FALSE, seed=1,
                   chains=chains, iter=iter, refresh=2000,
                   multi.C="_5" )
        print("...fit")
        gc(verbose=FALSE)
        
        # ---- eval ----
        get.bridge( model=model.dep, lang.name=lang, 
                    trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                    dat=dat, cores=chains,
                    multi.C="_5" )
        print("...bridge")
        gc(verbose=FALSE)
        get.LOOCs( lang.name=lang, model_name=model.dep@model_name,
                   trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                   cores=chains,
                   multi.C="_5" )
        print("...looc")
        gc(verbose=FALSE)
        
        # ============== independent binomial model ==============
        
        print("")
        print( "- independent" )
        print(model.indep@model_name)
        file.name <- file.path( paste(model.name, "-indep_", data.name, "-0", sep="") )
        
        fit.rstan( save.dir=save.dir, file.name=file.name, model=model.indep, dat=dat, 
                   loadExisting=FALSE, seed=1,
                   chains=chains, iter=iter, refresh=2000,
                   multi.C="_5" )
        print("...fit")
        gc(verbose=FALSE)
        
        # ---- eval ----
        get.bridge( model=model.indep, lang.name=lang, 
                    trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                    dat=dat, cores=chains,
                    multi.C="_5" )
        print("...bridge")
        gc(verbose=FALSE)
        get.LOOCs( lang.name=lang, model_name=model.indep@model_name,
                   trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                   cores=chains,
                   multi.C="_5" )
        print("...looc")
        gc(verbose=FALSE)

        # ============== WAIC ==============
        print(" -- WAIC --")

        model_name.dep <- "Y-binom-dep_"
        model_name.indep <- "Y-binom-indep_"

        print("dependent - independent")
        get.WAICs( model.1.name=model_name.dep,
                   model.2.name=model_name.indep,
                   lang.name=lang,
                   trait.name.1=trait.name.1, trait.name.2=trait.name.2,
                   multi.C="_5", dep="")

        remove( dat )
        gc(verbose=FALSE)
      }
    }
  }
}