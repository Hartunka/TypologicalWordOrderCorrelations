## run main script to fit pairwise trait models for given family

# load functions ----
source('paths.R', encoding="utf8" )
source( file.path(fit.dir, 'fit.R'), encoding="utf8" )
source( file.path(dataprep.dir, 'dataprep.R'), encoding="utf8" )
source( file.path(post.dir, 'diagnostics.R'), encoding="utf8" )

library(rstan)
library(bridgesampling)

model.name <- "Y-binom"
model_name.dep <- "Y-binom-dep_"
model_name.indep <- "Y-binom-indep_"

print("parse models...")
model.1 <- rstan::stan_model( 
  file.path( models.dir, 'Binom-multiC-Y-nonCentered-main.stan' ),
  model_name=model_name.dep )
print("...dep done")
model.2 <- rstan::stan_model( 
  file.path( models.dir, 'Binom-indep-multiC-Y-nonCentered-main.stan' ),
  model_name=model_name.indep )
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

multi.C <- "_5"

# ---- main loop ----

for ( i in 1:length(lang.fam.names) ){
  lang <- lang.fam.names[i]
  print( paste("---", lang, "---") )
  
  save.dir.tab <- file.path(post.dir, "rDats", lang)
  if (!dir.exists(save.dir.tab)) { dir.create(save.dir.tab) }
  
  trait.pairs <- vector(mode="character", length=28)
  values.looc <- vector(mode="numeric", length=28)
  values.waic <- vector(mode="numeric", length=28)
  wWaic.1 <- vector(mode="numeric", length=28)
  wWaic.2 <- vector(mode="numeric", length=28)
  values.bf.1 <- vector(mode="numeric", length=28)
  values.bf.2 <- vector(mode="numeric", length=28)
  j <- 1
  
  t.pairs <- c()
  for (trait.name.1 in trait.names){
    for (trait.name.2 in trait.names){
      pair <- paste( trait.name.1, trait.name.2 )
    
      if ( (trait.name.1 != trait.name.2) & !(pair %in% t.pairs) ){
        print( paste("----", trait.name.1, "-", trait.name.2, "----") )
      
        t.pairs <- c( t.pairs, pair, paste(trait.name.2, trait.name.1) )
        trait.pairs[j] <- pair
        
        ## ---- LOOC ----
        print(" -- LOOC --")
        
        loo.1 <- get.LOOCs( lang.name=lang, model_name=model_name.dep,
                            trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                            multi.C=multi.C )
        loo.2 <- get.LOOCs( lang.name=lang, model_name=model_name.indep,
                            trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                            multi.C=multi.C )
                            
        looc <- loo::loo_compare( loo.1, loo.2)
        looc.df <- as.data.frame( looc )
                            
        if ( looc.df["model1","elpd_diff"] == looc.df[2,1] ){
          value.looc <- looc.df[2,1]
        } else {
          value.looc <- - looc.df[2,1]
        }
        values.looc[j] <- value.looc
        
        ## ---- WAIC ----
        print(" -- WAIC --")
        
        if (!dir.exists(file.path(post.dir, "waic", lang))){
          dir.create(file.path(post.dir, "waic", lang))
        }
          
        waic <- get.WAICs( model.1.name=model_name.dep,
                           model.2.name=model_name.indep,
                           lang.name=lang,
                           trait.name.1=trait.name.1, trait.name.2=trait.name.2,
                           multi.C=multi.C, dep="")
        waic.df <- as.data.frame( waic )
          
        values.waic[j] <- waic.df["stanfit.1","WAIC"] - waic.df["stanfit.2","WAIC"]
          
        wWaic.1[j] <- waic.df["stanfit.1","weight"]
        wWaic.2[j] <- waic.df["stanfit.2","weight"]
          
        
        ## ---- Bayes Factors ----
        print(" -- Bayes Factors --")
        
        dependent <- get.bridge( model=model.1, lang.name=lang, 
                                 trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                                 dat=list(), multi.C=multi.C )
        independent <- get.bridge( model=model.2, lang.name=lang, 
                                   trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
                                   dat=list(), multi.C=multi.C )
        
        values.bf.1[j] <- bf(dependent, independent)$bf
        values.bf.2[j] <- bf(independent, dependent)$bf
        
        j <- j+1
        gc(verbose=FALSE)
      }
    }
  }
  
  print("saving...")
  # save looc
  looc.Y.dep.5.0.df <- data.frame(trait.pairs=as.factor(trait.pairs), loocs=values.looc)
  f_looc_df <- paste(lang, "_avgC-5_0_Y_loocdf.RData", sep="")
  save(looc.Y.dep.5.0.df, file=file.path(save.dir.tab, f_looc_df))
  
  # save waic
  waic.Y.dep.5.0.df <- data.frame(trait.pairs=as.factor(trait.pairs), dWAIC=values.waic)
  f_waic_df <- paste(lang, "_avgC-5_0_Y_waicdf.RData", sep="")
  save(waic.Y.dep.5.0.df, file=file.path(save.dir.tab, f_waic_df))
  
  # save waic weights
  wWAIC.df.a <- data.frame( trait.pairs=as.factor(trait.pairs), WAIC.weight=wWaic.1, 
                            model=rep("dependent", length(wWaic.1)) )
  wWAIC.df.b <- data.frame( trait.pairs=as.factor(trait.pairs), WAIC.weight=wWaic.2,
                            model=rep("independent", length(wWaic.2)) )
  waic.W.Y.dep.5.0.df <- rbind(wWAIC.df.a, wWAIC.df.b)
  f_wWaic_df <- paste(lang, "_avgC-5_0_Y_wWaicdf.RData", sep="")
  save(waic.W.Y.dep.5.0.df, file=file.path(save.dir.tab, f_wWaic_df))
  
  # save bf
  bf.Y.dep.5.0.df <- data.frame(trait.pairs=as.factor(trait.pairs), bayes.factors=values.bf.1)
  f_bf_df <- paste(lang, "_avgC-5_0_Y_bfdf.RData", sep="")
  save(bf.Y.dep.5.0.df, file=file.path(save.dir.tab, f_bf_df))
}

