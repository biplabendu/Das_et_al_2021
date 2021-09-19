## README: cflo_go_enrichment -------

# goal: using the cflo_go_enrichment() function
# author: Biplabendu Das

# CLEAN THE WORKING ENVIRONMENT
rm(list = ls())


# INSTALL/LOAD NECESSARY PACKAGES
require("pacman") # if this code does not install/load, try install.packages("pacman")
pacman::p_load(pheatmap, tidyverse, viridis, ggplot2)


# Step 00: provide the path to the "functions" folder 
#
# !! VERY IMPORTANT !! 
# !! CHANGE THIS TO YOUR path/to/folder !!
path_to_function = "/Users/biplabendudas/Documents/GitHub/Das_et_al_2021"


# Step 01: load the function into your R-Studio environment
source(paste0(path_to_function,"/functions/enrichment_analysis.R"))
#
## You should see four functions loaded onto your environment,
## we will need the cflo_go_enrichment() function to run GO enrichments


# Step 02: Using the function to run enrichment
#
### part 1: load your genes of interest ---- 
#
### here, I am loading all the data from Das et al. 2021 (bioRxiv))
load(paste0(path_to_function,"/functions/func_data/TC5_core_datasets.RData"))
### sample 1000 random genes to run GO enrichment on
mock.genes <- sample(cflo.annots.exp[[1]], 1000)
### check the first couple of the gene names
head(mock.genes)
#
### part 2: run enrichment analysis and save results ----
#
mock.enriched.gos <- cflo_go_enrichment(geneset = mock.genes, function.dir = path_to_function)
# View the first 5 rows of the results
mock.enriched.gos[1:5,]


# Step 03: Plot the results
go_enrichment_plot(data = mock.enriched.gos,
                   function.dir = path_to_function,
                   # change your False Discovery Rate here
                   fdr=10)

