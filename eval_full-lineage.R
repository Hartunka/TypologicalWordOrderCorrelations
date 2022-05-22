## run main script to fit pairwise trait models for given family

# load functions ----
source('paths.R', encoding="utf8" )
source( file.path(fit.dir, 'fit.R'), encoding="utf8" )
source( file.path(dataprep.dir, 'dataprep.R'), encoding="utf8" )
source( file.path(post.dir, 'diagnostics.R'), encoding="utf8" )

library(rstan)
library(bridgesampling)

model.name <- "Y-binom"
model_name.univ="Y-binom-univ-dep_"
model_name.lin="Y-binom-lin_"

print("parse models...")
model.univ <- rstan::stan_model( 
  file.path( models.dir, 'categorical', 'Binom-Y-nonCentered-fullLang-multiC-univ.stan' ),
  model_name="Y-binom-univ-dep_" )
print("...universal done")
model.lin <- rstan::stan_model( 
  file.path( models.dir, 'categorical', 'Binom-Y-nonCentered-fullLang-multiC-lin.stan' ),
  model_name="Y-binom-lin_" )
print("...lineage done")

lang <- "full"
trait.names <- c( "AN", "PN", "ND", "NG", "NNum", "VO", "NRc", "VS" )
t.pairs <- c()

multi.C <- "_5"

# ---- main loop ----

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
    pair = paste( trait.name.1, trait.name.2 )
    
    if ( (trait.name.1 != trait.name.2) & !(pair %in% t.pairs) ){
      print( paste("----", trait.name.1, "-", trait.name.2, "----") )
      
      t.pairs = c( t.pairs, pair, paste(trait.name.2, trait.name.1) )
      trait.pairs[j] <- pair
        
      ## ---- LOOC ----
      #print(" -- LOOC --")
      #
      #loo.1 <- get.LOOCs( lang.name=lang, model_name=model_name.univ,
      #                    trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
      #                    multi.C="_5_5" )
      #loo.2 <- get.LOOCs( lang.name=lang, model_name=model_name.lin,
      #                    trait.name.1=trait.name.1, trait.name.2=trait.name.2, 
      #                    multi.C="_5_5" )
      #                      
      #looc <- loo::loo_compare( loo.1, loo.2)
      #print("looc: univ, lin")
      #print(looc)
      #looc.df <- as.data.frame( looc )
      #                      
      #if ( looc.df["model1","elpd_diff"] == looc.df[2,1] ){
      #  value.looc <- looc.df[2,1]
      #} else {
      #  value.looc <- - looc.df[2,1]
      #}
      #values.looc[j] <- value.looc
      #  
      ## ---- WAIC ----
      print(" -- WAIC --")
      
      if (!dir.exists(file.path(post.dir, "waic", lang))){
        dir.create(file.path(post.dir, "waic", lang))
      }
          
      waic <- get.WAICs( model.1.name=model_name.univ,
                         model.2.name=model_name.lin,
                         lang.name=lang,
                         trait.name.1=trait.name.1, trait.name.2=trait.name.2,
                         multi.C=multi.C, dep="")
      print(waic)
      waic.df <- as.data.frame( waic )
          
      values.waic[j] <- waic.df["stanfit.1","WAIC"] - waic.df["stanfit.2","WAIC"]
          
      wWaic.1[j] <- waic.df["stanfit.1","weight"]
      wWaic.2[j] <- waic.df["stanfit.2","weight"]
          
        
      ## ---- Bayes Factors ----
      #print(" -- Bayes Factors --")
      #
      ## get data:
      #lang.data <- get.data.full( c(trait.name.1, trait.name.2), avg=TRUE, other=0, n_avg=5, getC=TRUE )
      #dat <- list(
      #  N = lang.data$N,
      #  x = lang.data$traits,
      #  M = 5,
      #  C = lang.data$Cs,
      #  # Hyperparameters
      #  eta = 0.8,
      #  lambda = 1,
      #  mu_z = 0,
      #  sigma_z = 1.5
      #)
      #remove(lang.data)
      #
      #universal <- get.bridge( model=model.univ, lang.name=lang, 
      #                         trait.name.1=trait.name.1, trait.name.2=trait.name.2, cores=14,
      #                         dat=dat, multi.C=multi.C )
      #lineage <- get.bridge( model=model.lin, lang.name=lang, 
      #                           trait.name.1=trait.name.1, trait.name.2=trait.name.2, cores=14,
      #                           dat=dat, multi.C=multi.C )
      #
      #print("universal:")
      #print(universal)
      #print("lineage:")
      #print(lineage)
      #print(bf(universal, lineage))
      #values.bf.1[j] <- bf(universal, lineage)$bf
      #print(bf(lineage, universal))
      #values.bf.2[j] <- bf(lineage, universal)$bf
      #
      j <- j+1
      gc(verbose=FALSE)
    }
  }
}
  
print("saving...")
# save looc
#looc.Y.full.5.0.df <- data.frame(trait.pairs=as.factor(trait.pairs), loocs=values.looc)
#f_looc_df <- paste(lang, "_avgC-5_0_v3_Y_loocdf.RData", sep="")
#save(looc.Y.full.5.0.df, file=file.path(save.dir.tab, f_looc_df))
  
# save waic
#waic.Y.full.5.0.df <- data.frame(trait.pairs=as.factor(trait.pairs), dWAIC=values.waic)
#f_waic_df <- paste(lang, "_avgC-5_0_v3_Y_waicdf.RData", sep="")
#save(waic.Y.full.5.0.df, file=file.path(save.dir.tab, f_waic_df))
  
# save waic weights
#wWAIC.df.a <- data.frame( trait.pairs=as.factor(trait.pairs), WAIC.weight=wWaic.1, 
#                          model=rep("universal", length(wWaic.1)) )
#wWAIC.df.b <- data.frame( trait.pairs=as.factor(trait.pairs), WAIC.weight=wWaic.2,
#                          model=rep("lineage", length(wWaic.2)) )
#waic.W.Y.full.5.0.df <- rbind(wWAIC.df.a, wWAIC.df.b)
#f_wWaic_df <- paste(lang, "_avgC-5_0_v3_Y_wWaicdf.RData", sep="")
#save(waic.W.Y.full.5.0.df, file=file.path(save.dir.tab, f_wWaic_df))
  
# save bf
#bf.Y.full.5.0.df <- data.frame(trait.pairs=as.factor(trait.pairs), bayes.factors=values.bf.1)
#f_bf_df <- paste(lang, "_avgC-5_0_v3_Y_bfdf.RData", sep="")
#save(bf.Y.full.5.0.df, file=file.path(save.dir.tab, f_bf_df))

#bf.Y.full.5.0.df.rev <- data.frame(trait.pairs=as.factor(trait.pairs), bayes.factors=values.bf.2)
#f_bf_df <- paste(lang, "_avgC-5_0_Y_v3_bfdf_reverse.RData", sep="")
#save(bf.Y.full.5.0.df.rev, file=file.path(save.dir.tab, f_bf_df))


