#source('paths.R', encoding="utf8")
#source('libraries.R', encoding="utf8")
#load.lib(prep=TRUE)
#library(ape)

tree.glot <- list(
  "Afro-Asiatic", "Algic", "Arawakan", "Athapaskan-Eyak-Tlingit", "Atlantic-Congo", "Austroasiatic", "Austronesian",
  "Cariban", "Central_Sudanic", "Chibchan",
  "Dravidian", "Gunwinyguan", "Indo-European", 
  "Mande", "Mongolic",
  "Nakh-Daghestanian", "Nilotic", "Nuclear_Torricelli", "Nuclear_Trans_New_Guinea",
  "Otomanguean",
  "Pama-Nyungan", "Pano-Tacanan",
  "Salishan", "Sepik", "Sino-Tibetan", "Siouan", "Sko", "Surmic",
  "Tai-Kadai", "Tucanoan", "Tupian", "Turkic",
  "Uralic", "Uto-Aztecan"
  )

# ======== read features from csv ========

#' read language features from csv file "charMtx.csc"
#' 
#' filters for features of one specified language family 
#' or for all families for which relation information is given in list \code{tree.glot}
#' 
#' transforms value types to numeric - if bin=TRUE, binary - values, such that
#' | "-" | 1 | 0/1 |
#' | "0" | 2 | 0   |
#' | "1" | 3 | 1   |
#' 
#' @param lang.name the name of language family to get features for, or "all" for unfiltered (character)
#' @param bin if TRUE features are to be read as binary (boolean)
#' @param other if \code{bin=TRUE}, is the 'other' category "-" to be read as 0 or 1 (numeric)
#' @return feats.lang : a data.frame containing:
#' | X    | family.subfam.lang name    |                          |                          |             |
#' | AN   | noun-adjective order       | "0" : noun-adjective     | "1" : adjective-noun     | "-" : other |
#' | PN   | adposition-noun order      | "0" : prepositions       | "1" : adpositions        | "-" : other |
#' | ND   | noun-demonstrative order   | "0" : noun-demonstrative | "1" : demonstrative-noun | "-" : other |
#' | NG   | noun-genetive order        | "0" : noun-genetive      | "1" : genetive-noun      | "-" : other |
#' | NNum | noun-numeral order         | "0" : noun-numeral       | "1" : numeral-noun       | "-" : other |
#' | VO   | verb-object order          | "0" : verb-object        | "1" : object-verb        | "-" : other |
#' | NRc  | noun-relative clause order | "0" : noun-relative      | "1" : relative-noun      | "-" : other |
#' | VS   | verb-subject order         | "0" : verb-subject       | "1" : subject-verb       | "-" : other |
#' | glot | glottolog lang family name |                          |                          |             |
#' @export
#' 
read.feats <- function( lang.name, bin=FALSE, other=0 ){
  
  feats <- read.csv( file.path( 
    data.dir,     # path to directory containing data, defined in "paths.R"
    "charMtx.csv" # name of file containing the feature information
    ) )
  # renaming feature name to avoid naming conflict with NA
  names(feats)[names(feats) == 'NA.'] <- 'AN'
  
  if (lang.name=="all"){
    # filter all languages, for which no relation information is given
    feats.filtered <- feats[feats$glot %in% tree.glot,]
    # sort by language family
    order.glot <- order(feats.filtered$glot)
    feats.lang <- feats.filtered[order.glot, ]
  } else {
  # filter specified language
    feats.lang <- feats[feats$glot == lang.name,]
  }
  
  # features as numerical, "-" -> 1, "0" -> 2, "1" -> 3
  # whether they're read as strings or as factor seems to depend on the R version
  if (is.factor(feats.lang$AN)){  # =< R 3.5
    for (column in 2:( ncol(feats.lang)-1) ) {
      feats.lang[column] <- as.numeric( feats.lang[[column]] )
    }
  } else if (is.character(feats.lang$AN)){ # >= R 4.0
    for (column in 2:( ncol(feats.lang)-1) ) {
      # transform first to factor, to fit with how csv is read in older version
      feats.lang[column] <- factor( feats.lang[[column]], c("-", "0", "1") )
      feats.lang[column] <- as.numeric( feats.lang[[column]] )
    }
  }
  # features as binary, "-"/1 -> 0|1, "0"/2 -> 0, "1"/3 -> 1
  if (bin){
    for (column in 2:( ncol(feats.lang)-1) ) {
      feats.lang[[column]][ feats.lang[[column]] != 3-other] <- 0+other
      feats.lang[[column]][ feats.lang[[column]] == 3-other] <- 1-other
    }
  }
  
  return(feats.lang)
}


# ======== read trees ========

#' read phylogenetic tree data for specified language family 
#' from file in data directory specified in paths.R
#'
#' @param lang.name name of language family (character),
#'                  respective tree file is expecteded to be named '<lang.name>.posterior.tree'
#' TODO: do I actually need this?:
#' @param feats.lang.names array of names of languages of family specified with \code{lang.name} (character)
#'                         phylogenetic information for languages not in \code{feats.lang.names} will be filtered out
#' @param phyl if TRUE, return trees of class phylo, else covariance matrices (boolean)
#' @return if \code{phyl=FALSE} trees.C : list of covariance matrices 
#'                                       for specified language family derived from phylogenetic trees
#'         else trees.trimmed : list of phylo, containing phylogenetic trees for family
#' @export
#'
read.trees <- function(lang.name, feats.lang.names, 
                       avg=FALSE, n_avg=1 ){
                       
  trees <- ape::read.tree(
    file.path( data.dir, "trees", paste(lang.name, ".posterior.tree", sep="" ) )
    )
  # filter trees to only contain languages for which feature information is given
  trees.trimmed <- lapply( trees, 
                           function(i) { ape::keep.tip( i, feats.lang.names ) }
                           )
  if (!avg) {
    trees.C <- list()
    for ( i in 1:length( trees.trimmed ) ){
      suppressWarnings({
        tree <- trees.trimmed[[i]]
        tree.C <- ape::vcv( ape::corBrownian(1, tree), corr=FALSE )[feats.lang.names, feats.lang.names]
        trees.C[[i]] <- tree.C
      })
    }
  return( trees.C )
  } else{
    nC <- length(trees.trimmed)
    n.split <- ceiling(nC/n_avg)
    tree.splits <- split(trees.trimmed, rep(1:ceiling(nC/n.split), each=n.split, length.out=nC))
    avg.Cs <- list()
    
    for ( h in 1:length(tree.splits) ) {
      tree.split <- tree.splits[[h]]
      
      nrc <- length(tree.split[[1]]$tip.label)
      sum.C <- matrix(data=rep(0, nrc*nrc),  nrow=nrc, ncol=nrc)
      for ( i in 1:length( tree.split ) ) {
        suppressWarnings({
          tree <- tree.split[[i]]
          C.i <- ape::vcv( ape::corBrownian(1, tree), corr=FALSE )[feats.lang.names, feats.lang.names]
        })
        sum.C = sum.C + C.i
      }
      avg.Cs[[h]] = sum.C / length( tree.split )
    }
    return(avg.Cs)
  }
}

#' create covariance matrix for all language families (in data directory), 
#' or a list of languages from input
#' 
#' for each language family the values are the average over all trees for the family
#' 
#' @return if \code{phyl=FALSE} trees.C : covariance matrix for all or specified languages families,
#'         else trees: list of multiPhylo object per language
#' @export
#' 
read.trees.multi.Lang <- function(feats, avg=TRUE, n_avg=1){
  
  fam.names <- unique(feats$glot)
  trees <- list()
  for (i in 1:length(fam.names)) {
    lang.names <- as.character(feats[feats$glot==fam.names[i],]$X)
    if (is.factor(lang.names)){
      lang.names <- as.character(lang.names)
    }
    trees_i <- read.trees(fam.names[i], lang.names, avg=avg, n_avg=n_avg)
    if (!avg){
      set.seed(1)
      trees[[i]] <- sample(trees_i, max.C)
      set.seed(NULL)
    } else {
      trees[[i]] = trees_i
    }
  }
  
  n.Lang <- nrow(feats)
  full.Cs <- list()
  for (n in 1:n_avg) {
    # initiate empty matrix, to be filled with languages' phylogenetic distances
    full.C.n <- matrix( data=0, nrow=n.Lang, 
                                ncol=n.Lang )
    prev <- 0
    for ( i in 1:length(trees) ){
      # trees for i^th family
      trees.i <- trees[[i]]

      # insert values of current family into full matrix
      curr <- nrow(trees.i[[1]])
      for (j in 1:curr) {
        for (k in 1:curr) {
          full.C.n[j+prev,k+prev] <- trees.i[[n]][j,k]
        }
      }
      prev <- prev + curr
    }
    full.Cs[[n]] <- full.C.n
  }
  return(full.Cs)
}

get.data.full <- function( trait.names, other=0, avg=TRUE, n_avg=1, getCs=TRUE){
  
  max.C <- ifelse(avg, 1, max.C)
  
  feats <- read.feats("all", bin=TRUE, other=other)
  fam.names <- unique(feats$glot)
  
  lang.names <- c()
  traits <- c()
  Ns <- c()
  for (name in fam.names) {
    fam.feats <- feats[feats$glot==name,]
    lang.names <- c( lang.names, fam.feats$X )
    Ns <- c(Ns, nrow(fam.feats))
    traits <- c( traits, as.numeric(unlist(fam.feats[trait.names[1]])) )
    traits <- c( traits, as.numeric(unlist(fam.feats[trait.names[2]])) )
  }
  
  if (getCs){
    Cs <- read.trees.multi.Lang(feats, avg=avg, n_avg=n_avg )
  } else { Cs <- NA }
  
  return(list(Ns=Ns, traits=traits, Cs=Cs, nMat=n_avg, lang.names=lang.names))
}
