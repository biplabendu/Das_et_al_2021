
# Supp. Excel 2 ---------------------------------------------------
  
rm(list = ls())
set.seed(420)
  
# load libraries
pacman::p_load(pheatmap, dendextend, tidyverse, viridis)
pacman::p_load(RSQLite, tidyverse, dbplyr, DT, conflicted)
# set conflict preference
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("layout", "plotly")

# load database
db <- dbConnect(RSQLite::SQLite(), "/Users/biplabendudas/Documents/GitHub/R-scripts_zombie_ant_lab/RSQLite/sql_dbs/TC5_data.db")

# Retrieve the gene annotations:
cflo.annots <-
  tbl(db, "annot_fpkm") %>% 
  # select(gene_name, blast_annotation = old_annotation, signalP, TMHMM, GOs, pfams) %>% 
  collect()

# Retrieve the eJTK data
ejtk.24 <- tbl(db, "ejtk_all") %>% collect()
ejtk.12 <- tbl(db, "ejtk_12h_all") %>% collect()
ejtk.8 <- tbl(db, "ejtk_8h_all") %>% collect()

## make gene lists:
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

# split the ejtk results into foragers and nurses
# Circadian genes - 24-hour
ejtk.for24 <- ejtk.24 %>% filter(caste == "for") %>% arrange(GammaP) 
ejtk.nur24 <- ejtk.24 %>% filter(caste == "nur") %>% arrange(GammaP) 
# Ultradian genes - 12-hour
ejtk.for12 <- ejtk.12 %>% filter(caste == "for") %>% arrange(GammaP) 
ejtk.nur12 <- ejtk.12 %>% filter(caste == "nur") %>% arrange(GammaP) 
# Ultradian genes - 8-hour
ejtk.for8 <- ejtk.8 %>% filter(caste == "for") %>% arrange(GammaP) 
ejtk.nur8 <- ejtk.8 %>% filter(caste == "nur") %>% arrange(GammaP) 



# Supp 2 - Worksheet 1 ----------------------------------------------------

# eTJK output for ALL tested forager genes, including their gene number and 
# normalized expression levels for each time point, sort based on significance
  
ejtk.for24 %>% 
    left_join((cflo.annots %>% 
                 select(gene_name, X2F:X24F)),
              by="gene_name") %>% 
    write.csv(., "./results/supp_files/ejtk24_for_all.csv")
  

# Supp 2 - Worksheet 2 ----------------------------------------------------

# eTJK output for ALL tested nurse genes, including their gene number and 
# normalized expression levels for each time point, sort based on significance
  
ejtk.nur24 %>% 
    left_join((cflo.annots %>% 
                 select(gene_name, X2N:X24N)),
              by="gene_name") %>% 
    write.csv(., "./results/supp_files/ejtk24_nur_all.csv")
  

# Cluster for24 and nur24 genes ------------------------------------------
  
## load zscore datasets
cflo.zscores.for <- tbl(db, "zscore_for") %>% collect()
cflo.zscores.nur <- tbl(db, "zscore_nur") %>% collect()

# Filter the zscores to keep only circadian genes
cflo.rhy.exp.for <- cflo.zscores.for %>% 
  filter(gene_name %in% for24) %>% as.data.frame()
cflo.rhy.exp.nur <- cflo.zscores.nur %>% 
  filter(gene_name %in% nur24) %>% as.data.frame()

# Set genes as rownames and convert it into a matrix
rownames(cflo.rhy.exp.for) = cflo.rhy.exp.for$gene_name
cflo.rhy.exp.for <- as.matrix(cflo.rhy.exp.for[-1])
rownames(cflo.rhy.exp.nur) = cflo.rhy.exp.nur$gene_name
cflo.rhy.exp.nur <- as.matrix(cflo.rhy.exp.nur[-1])

# Hierarchical clustering of the for24 and nur24 genesets
my_hclust_gene.for <- hclust(dist(cflo.rhy.exp.for), method = "complete")
my_hclust_gene.nur <- hclust(dist(cflo.rhy.exp.nur), method = "complete")

# Make annotations for the gene-clusters 
## Foragers
my_gene_col.for <- cutree(tree = as.dendrogram(my_hclust_gene.for), k = 4)
my_gene_col.for <- data.frame(for_cluster = my_gene_col.for)
my_gene_col.for$gene_name <- row.names(my_gene_col.for)

## Nurses
my_gene_col.nur <- cutree(tree = as.dendrogram(my_hclust_gene.nur), k = 4)
my_gene_col.nur <- data.frame(nur_cluster = my_gene_col.nur)
my_gene_col.nur$gene_name <- row.names(my_gene_col.nur)
  
  
# Supp 2 - Worksheet 3 ----------------------------------------------------

# Gene numbers with normalized expression levels of circadian forager genes, 
# divided into the 4 hierarchical clusters.  

my_gene_col.for %>% 
  left_join((cflo.annots %>% 
               select(gene_name, X2F:X24F)),
            by="gene_name") %>%
  select(gene_name, everything()) %>% 
  arrange(for_cluster) %>% 
  write.csv(., "./results/supp_files/for24_clusters.csv")


# Supp 2 - Worksheet 4 ----------------------------------------------------

# Gene numbers with normalized expression levels of circadian nurse genes, 
# divided into the 4 hierarchical clusters.

my_gene_col.nur %>% 
  left_join((cflo.annots %>% 
               select(gene_name, X2N:X24N)),
            by="gene_name") %>%
  select(gene_name, everything()) %>% 
  arrange(nur_cluster) %>% 
  write.csv(., "./results/supp_files/nur24_clusters.csv")

# Supp 2 - Worksheet 5 ----------------------------------------------------

# GO enrichment results for circadian genes in foragers (for-24h) and nurses (nur-24h) 
# that peak during the day (day-peaking clusters) and night (night-peaking clusters). 
# Also includes the enrichment results for day-peaking cluster of overlapping for-24h and 
# nur-24h genes (for-24h-nur-24h).

# Note: enrichments of for-24h and nur-24h clusters have been performed in "plotting_pheatmaps.R"

# for-24h-nur-24h, day-peaking cluster: Cluster1
# Load the enrichment function
source(file = "./functions/enrichment_analysis.R")

# Clustering for24-nur24 --------------------------------------------------
# Subset the intersecting geneset
for24.nur24 <- cflo.zscores.for %>% 
  filter(gene_name %in% intersect(for24, nur24)) %>% 
  left_join(cflo.zscores.nur) %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "gene_name")

# Set genes as rownames and convert it into a matrix
for24.nur24 <- as.matrix(for24.nur24[-1])  
# cluster
my_hclust_gene.for24.nur24 <- hclust(dist(for24.nur24), 
                                     method = "complete")
# Make annotations for the heatmaps
my_gene_col.for24.nur24 <- cutree(tree = as.dendrogram(my_hclust_gene.for24.nur24), k = 4)
my_gene_col.for24.nur24 <- data.frame(cluster = my_gene_col.for24.nur24)
my_gene_col.for24.nur24$gene_name <- row.names(my_gene_col.for24.nur24)

# for-24h-nur-24h genes and their expression
my_gene_col.for24.nur24 %>% 
  left_join((cflo.annots %>% 
               select(gene_name, X2F:X24N)),
            by="gene_name") %>%
  select(gene_name, everything()) %>% 
  arrange(cluster) %>% 
  write.csv(., "./results/supp_files/for24nur24_clusters.csv")

# Enrichment for all the day-peaking for-24h-nur-24h genes (Cluster 1)
my_gene_col.for24.nur24 %>% 
  filter(cluster==1) %>% 
  pull(gene_name) %>%
  cflo_go_enrichment() %>% 
  write.csv(., "./results/supp_files/for24nur24_daypeaking_cluster1.csv")


# Supp. Excel 3 -----------------------------------------------------------


# Supp 3 - Worksheet 1 ----------------------------------------
# eTJK output for 12h periodicity for ALL tested forager genes, 
# including their gene number and normalized expression levels for each time point, 
# sort based on significance

ejtk.for12 %>% 
  arrange(GammaP) %>% 
  left_join((cflo.annots %>% 
               select(gene_name, X2F:X24F)),
            by="gene_name") %>% 
  write.csv(., "./results/supp_files/ejtk12_for_all.csv")


# Supp 3 - Worksheet 2 ----------------------------------------
# eTJK output for 12h periodicity ALL tested nurse genes, 
# including their gene number and normalized expression levels for each time point,
# sort based on significance

ejtk.nur12 %>% 
  arrange(GammaP) %>% 
  left_join((cflo.annots %>% 
               select(gene_name, X2N:X24N)),
            by="gene_name") %>% 
  write.csv(., "./results/supp_files/ejtk12_nur_all.csv")


# Supp 3 - Worksheet 3 ----------------------------------------------
# eTJK output for 8h periodicity for ALL tested forager genes, 
# including their gene number and normalized expression levels for each time point, 
# sort based on significance

ejtk.for8 %>% 
  arrange(GammaP) %>% 
  left_join((cflo.annots %>% 
               select(gene_name, X2F:X24F)),
            by="gene_name") %>% 
  write.csv(., "./results/supp_files/ejtk8_for_all.csv")


# Supp 3 - Worksheet 4 ----------------------------------------------
# eTJK output for 8h periodicity ALL tested nurse genes, 
# including their gene number and normalized expression levels for each time point

ejtk.nur8 %>%
  arrange(GammaP) %>% 
  left_join((cflo.annots %>% 
               select(gene_name, X2N:X24N)),
            by="gene_name") %>% 
  write.csv(., "./results/supp_files/ejtk8_nur_all.csv")



# Supp. Excel 5 -----------------------------------------------------------


# Supp 5 - Worksheet 1 ----------------------------------------------------
# Gene numbers with normalized expression levels of circadian forager genes, 
# and ultradian 8h nurse genes that clustered with per.

# What I am doing:
# - show all for24h-nur8h genes
# - indicate cluster identity using a cluster column

# Subset the for24h-nur8h genes
for24.nur8 <- cflo.zscores.for %>% 
  filter(gene_name %in% intersect(for24, nur8)) %>% 
  left_join(cflo.zscores.nur) %>% 
  as.data.frame()

# Set genes as rownames and convert it into a matrix
rownames(for24.nur8) <- for24.nur8[[1]]
for24.nur8 <- as.matrix(for24.nur8[-1])  

# cluster
my_hclust_gene.for24.nur8 <- hclust(dist(for24.nur8), method = "complete")
my_gene_col.for24.nur8 <- cutree(tree = as.dendrogram(my_hclust_gene.for24.nur8), k = 4)
my_gene_col.for24.nur8 <- data.frame(cluster = my_gene_col.for24.nur8)
my_gene_col.for24.nur8$gene_name <- row.names(my_gene_col.for24.nur8)



# Identify cluster identity of Period gene --------------------------------

# Plot for24h-nur8h as heatmap and infer cluster identity 

# Let's make the column annotation data frame
my_sample_col <- data.frame(caste = rep(c("foragers","nurses"), c(12,12)),
                            phase = rep(rep(c("light", "dark", "light"), c(5,6,1)),2))
row.names(my_sample_col) <- colnames(for24.nur8)

## For combined plot, use the following palatte
my_colour = list(
  caste = c(foragers = "#F23030", nurses = "#1A80D9"),
  phase = c(light = "white", dark = "#403F3D"),
  cluster = c("#A6A056", "#260101", "#F2A25C", "#732F29")
  # overlap = c(no = "white", yes = "#D92818")
)


# Color scale
my.breaks = seq(min(for24.nur8), max(for24.nur8), by=0.06)

png("./results/figures/figure_4/for24nur8_all_v2.png", width = 800, height = 1400, res = 300)
for24.nur8 %>%
  pheatmap(show_rownames = F,
           show_colnames = F,
           # I want to see which FORAGER cluster is the gene from
           annotation_row = my_gene_col.for24.nur8 %>% select(cluster),
           # both light phase and ant caste identity
           annotation_col = my_sample_col,
           cutree_rows = 4,
           # cutree_cols = 2,
           annotation_colors = my_colour,
           border_color=F,
           cluster_cols = F,
           cluster_rows = T,
           breaks = my.breaks,
           ## color scheme borrowed from:
           color = inferno(length(my.breaks) - 1),
           # treeheight_row = 0,
           # treeheight_col = 0,
           # remove the color scale or not
           # main = paste0("Foragers - circadian genes \n (n=", nrow(cflo.rhy.exp.for), " genes)"),
           ## annotation legend
           annotation_legend = T,
           ## Color scale
           legend = T)
dev.off()
# Done.
# Cluster 1 represents the for24h-nur8h genes that cluster with Period

# Let's save the gene_names, cluster identity, annotation and expression to a csv
my_gene_col.for24.nur8 %>% 
  remove_rownames() %>% 
  arrange(for24nur8_cluster) %>% 
  left_join(cflo.annots, by="gene_name") %>% 
  select(gene_name, for24nur8_cluster, everything()) %>% 
  write.csv(., "./results/supp_files/supp_file_5/for24_nur8_clusters_exp.csv")


# Supp 5 - Worksheet 2 ----------------------------------------------------

# Results of GO enrichment analysis for Per-like DRG cluster, 
# oscillating every 24h in foragers and every 8h in nurses 
# (corresponding to Supplementary Figure 6) 
# note: supp. figure 6 will not be included!

# As we saw above, the Cluster 1 contains the for24h-nur8h genes clustered with Period

# run enrichment and save the results to a csv file

my_gene_col.for24.nur8 %>% 
  filter(for24nur8_cluster == 1) %>% # cluster 1 has 215 genes
  pull(gene_name) %>% 
  cflo_go_enrichment() %>% 
  write.csv(., "./results/supp_files/supp_file_5/for24_nur8_Cluster1_enrichedGos.csv")




# Supp. Excel 6 -----------------------------------------------------------


# Supp 6 - Worksheet 1 ----------------------------------------------------
# Please note, the genes involved in mammalian circadian entrainment were 
# obtained from KEGG pathways (KEGG pathway: hsa04713)

# x <- read.delim(pipe("pbpaste"))
# x %>% 
#   mutate(gene_name = cflo_gene_name) %>% 
#   select(-cflo_gene_name) %>% 
#   left_join(cflo.annots, by="gene_name") %>% 
#   rename(Cflo_ortholog = gene_name) %>% 
#   arrange(Cflo_ortholog) %>%
#   write.csv(., "./results/supp_files/supp_file_6/Cflo_ortho_mammal_circa_entrain.csv")



# Supp. Excel 7 -----------------------------------------------------------


# Supp 6 - Worksheet 1 ----------------------------------------------------
# Results of DEG analysis using LimoRhyde, for all genes tested

# The following shows all the results of the DEG analysis,
# sorted by abs(logFC)

tbl(db, "TC5_degs_all") %>% 
  collect() %>% 
  arrange(desc(abs(logFC))) %>% 
  write.csv(., "./results/supp_files/supp_file_7/limorhyde_DEG_all_tested.csv")


# Supp 7 - Worksheet 2 ----------------------------------------------------
# All DEGs “Of these DEGs, 34 were significantly higher expressed in forager brains, 
# and the remaining 47 were higher expressed in nurses (Fig. 8A-B; Supp. Excel 6).”
tbl(db, "TC5_degs_all") %>% 
  filter(significant=="yes") %>% 
  collect() %>% 
  arrange(upregulation, desc(abs(logFC))) %>% 
  left_join((cflo.annots %>% 
               select(gene_name, old_annotation, X2F:X2N)),
            by = "gene_name") %>% 
  write.csv(., "./results/supp_files/supp_file_7/limorhyde_DEG_sig.csv")


# Supp 7 - Worksheet 3 ----------------------------------------------------
# Enriched GOs for nurse-DEGs (there were no enriched GOs in forager-DEGs)

tbl(db, "TC5_degs_all") %>% 
  filter(significant=="yes") %>% 
  collect() %>% 
  filter(upregulation == "nur") %>% 
  pull(gene_name) %>% 
  cflo_go_enrichment() %>% 
  write.csv(., "./results/supp_files/supp_file_7/nurseDEG_enrichedGOs.csv")

  

