This is the code for the experiments described in the paper 
"Kai Hartung, Gerhard Jäger, Munir Georges, Sören Gröttrup, *Typological Word Order Correlations with Logistic Brownian Motion*", submitted to the SIGTYP 2022 Workshop (https://sigtyp.github.io/workshop.html)

The script families.R runs the stan models for the single language models. The script full.R runs the model variants over all language families.
The script eval_families.R runs the comparisons between dependent and independent variants for the fitted single language models. 
The scripts eval_full-dep.R, eval_full-lin_indep.R, eval_full-lineage.R run the comparisons between dependent - independent, lineage-specific - independent and lineage-specific - dependent respectively, for the models over all language families.
The comparisons include Bayes Factors, WAIC and LOOIC.

The necessary libraries are listed in libraries.R and can be installed by running the script. 
The main project path in the paths.R might has to be adapted locally.

The word-order data can be found in dataprep/dat/charMtx.csv. The phylogenetic trees are in the directory dataprep/dat/trees.
