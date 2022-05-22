install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
install.packages("loo", "bridgesampling", 
                 dependencies=c("Depends", "Imports", "LinkingTo"))
install.packages(c("coda","mvtnorm","devtools","dagitty"))
devtools::install_github("rmcelreath/rethinking")









