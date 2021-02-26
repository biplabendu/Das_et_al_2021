
# Housekeeping ----------------------------------
rm(list = ls())

## Load packages ----------
pacman::p_load(pheatmap, dendextend, tidyverse, viridis)
pacman::p_load(RSQLite, tidyverse, dbplyr, DT, conflicted)

# set conflict preference (housekeeping to make sure functions work as expected)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("layout", "plotly")

conflict_prefer("setdiff", "BiocGenerics")
conflict_prefer("union", "BiocGenerics")
conflict_prefer("intersect", "BiocGenerics")


# Load functions ----------------------------------------------------------

# my functions
source("./functions/plot_zscores.R")
source("./functions/enrichment_analysis.R")
source("./functions/theme_publication.R")


## Load files --------------
# - database (db)
db <- dbConnect(RSQLite::SQLite(), "./data/TC5_data.db")
src_dbi(db)
# load the annotation file for Cflo
cflo.annots.exp <- tbl(db, "annot_fpkm") %>% collect()

# Get a list of all rhythmic genes ----------------------------------------

# Circadian genes (period = 24h)
## Foragers
for24 <-
  tbl(db, "ejtk_all") %>% 
  filter(caste == "for" & rhy == "yes") %>% 
  select(gene_name, GammaP) %>% collect() %>% arrange(GammaP) %>%
  select(gene_name) %>% pull()
## Nurses
nur24 <-
  tbl(db, "ejtk_all") %>% 
  filter(caste == "nur" & rhy == "yes") %>% 
  select(gene_name, GammaP) %>% collect() %>% arrange(GammaP) %>% 
  select(gene_name) %>% pull()

# Ultradian genes (period = 8h)
## Foragers
for8 <- 
  tbl(db, "ejtk_8h_all") %>% 
  filter(caste == "for" & rhy == "yes") %>% 
  select(gene_name, GammaP) %>% collect() %>% arrange(GammaP) %>% 
  select(gene_name) %>% pull()
## Nurses
nur8 <- 
  tbl(db, "ejtk_8h_all") %>% 
  filter(caste == "nur" & rhy == "yes") %>% 
  select(gene_name, GammaP) %>% collect() %>% arrange(GammaP) %>% 
  select(gene_name) %>% pull()

# Ultradian genes (period = 12h)
## Foragers
for12 <- 
  tbl(db, "ejtk_12h_all") %>% 
  filter(caste == "for" & rhy == "yes") %>% 
  select(gene_name, GammaP) %>% collect() %>% arrange(GammaP) %>% 
  select(gene_name) %>% pull()
## Nurses
nur12 <- 
  tbl(db, "ejtk_12h_all") %>% 
  filter(caste == "nur" & rhy == "yes") %>% 
  select(gene_name, GammaP) %>% collect() %>% arrange(GammaP) %>% 
  select(gene_name) %>% pull()



# Supp. Excel File 2 ----------------------------------------------

# Save the gene lists

# for24
cflo.annots.exp %>% 
  filter(gene_name %in% for24) %>% 
  write.csv(., file = "./results/gene_lists/for24_all.csv")
# nur24
cflo.annots.exp %>% 
  filter(gene_name %in% nur24) %>% 
  write.csv(., file = "./results/gene_lists/nur24_all.csv")

# for8
cflo.annots.exp %>% 
  filter(gene_name %in% for8) %>% 
  write.csv(., file = "./results/gene_lists/for8_all.csv")
# nur8
cflo.annots.exp %>% 
  filter(gene_name %in% nur8) %>% 
  write.csv(., "./results/gene_lists/nur8_all.csv")

# for12
cflo.annots.exp %>% 
  filter(gene_name %in% for12) %>% 
  write.csv(., "./results/gene_lists/for12_all.csv")
# nur12
cflo.annots.exp %>% 
  filter(gene_name %in% nur12) %>% 
  write.csv(., "./results/gene_lists/nur12_all.csv")



# Supp. Excel File 4 ------------------------------------------------------

# for24-nur8
cflo.annots.exp %>% 
  filter(gene_name %in% intersect(for24, nur8)) %>% 
  write.csv(., file = "./results/gene_lists/for24nur8_all.csv")

intersect(for24, nur8) %>% 
  cflo_go_enrichment(bg=bg.genes) %>% 
  write.csv(., file = "./results/enrichment_results/for24nur8_all.csv")

