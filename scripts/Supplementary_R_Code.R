
set.seed(420)
rm(list = ls())

## Load packages ----------
pacman::p_load(pheatmap, dendextend, tidyverse, viridis)
pacman::p_load(RSQLite, tidyverse, dbplyr, DT, conflicted)

# set conflict preference (housekeeping to make sure functions work as expected)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("layout", "plotly")


## Load files --------------
# - database (db)
db <- dbConnect(RSQLite::SQLite(), "./data/TC5_data.db")
src_dbi(db)
# - db::cflo.annots
cflo.annots.exp <- tbl(db, "annot_fpkm") %>% collect()

## Load functions ----------
# - cflo_go_enrichment
# - go_enrichment_plot
source(file = "./functions/enrichment_analysis.R")
source(file = "./functions/plot_zscores.R")


#########################################################################################
## Supplementary Figure 2 ## "TC5_expressed_genes::line_103"
############################
# Run enrichment for unexpressed and lowly-expressed genes
# and save the results as csvs & plot the results.
############################

# list of all protein coding genes from Cflo_v7
all.Cflo.genes <- 
	cflo.annots.exp %>% 
	pull(gene_name) %>% 
	unique()

# list of all Cflo genes that have non-zero expression in either forager or nurse brains
some.expression <-
  cflo.annots.exp %>% 
  dplyr::select(gene_name, X2F:X24N) %>% 
  na.omit() %>% 
  filter_at(vars(starts_with("X")), any_vars(. != 0)) %>% 
  pull(gene_name)

# list of all Cflo "Expressed" genes (>1 FPKM for either foragers or nurses)
expressed <- 
  cflo.annots.exp %>% 
  dplyr::select(gene_name, X2F:X24N) %>% 
  na.omit() %>% 
  filter_at(vars(starts_with("X")), any_vars(. > 1)) %>% 
  pull(gene_name)


# no_expression -----------------------------------------------------------

# list of genes that show "No expression"
no.expression <- setdiff(all.Cflo.genes, some.expression)

# Save list of genes with "no_expression"
cflo.annots.exp %>% 
  filter(gene_name %in% no.expression) %>% 
  write.csv(., file = "./results/gene_lists/no_exp.csv")
# run enrichment analysis for no.expression geneset against all Cflo genes
no.expression %>% 
  cflo_go_enrichment(bg=all.Cflo.genes) %>% 
  write.csv(., file = "./results/enrichment_results/no_exp.csv")
# plot the enrichment results
no.expression %>%
	cflo_go_enrichment(bg=all.Cflo.genes) %>% 
	go_enrichment_plot(clean = "no")


# low_expression ----------------------------------------------------------

# list of genes that show "Low expression"
low.expression <- setdiff(some.expression, expressed)

# Save list to a csv
cflo.annots.exp %>% 
  filter(gene_name %in% low.expression) %>% 
  write.csv(., file = "./results/gene_lists/low_exp.csv")
# Save enrichment results to a csv
low.expression %>% 
  cflo_go_enrichment(bg=all.Cflo.genes) %>%
  write.csv(., file = "./results/enrichment_results/low_exp.csv")
# plot the enrichment results
low.expression %>% 
  cflo_go_enrichment(bg=all.Cflo.genes) %>%
  go_enrichment_plot(clean="no")
  

#########################################################################################
## Supplementary Figure 3 ## "formatting_feature_table_27Apr20::line_1162"
############################
# Run enrichment for genes expressed only in forager or nurse brains
# and save the results as csvs & plot the results.
############################

## List of expressed genes
# Foragers
dat.f <- cflo.annots.exp %>%
  dplyr::select(gene_name, X2F:X24F) %>%
  # omit the genes that have NAs
  na.omit() %>%
  # remove genes that have at least 1 fpkm for at least one time point
  filter_at(vars(starts_with("X")), any_vars(. >= 1))
# Nurses
dat.n <- cflo.annots.exp %>%
  dplyr::select(gene_name, X2N:X24N) %>%
  # omit the genes that have NAs
  na.omit() %>%
  # remove genes that have at least 1 fpkm for at least one time point
  filter_at(vars(starts_with("X")), any_vars(. >= 1))

# "expressed" genes - foragers
	bg.genes.f <- dat.f %>% pull(gene_name)
# "expressed" genes - foragers
	bg.genes.n <- dat.n %>% pull(gene_name) 
# background "expressed" genes - all
	bg.genes <- unique(c(bg.genes.f, bg.genes.n))
  ## "bg.genes" will be used as the background geneset for 
	## further functional enrichment analyses
	
## What are the uniquely expressed genes in the forager ant brain?
## FORAGERS
# save the list to a csv
cflo.annots.exp %>% 
  filter( gene_name %in% setdiff(bg.genes.f, bg.genes.n)) %>% 
  write.csv(., "./results/gene_lists/uniquely_exp_for.csv")
# run enrichment 
for.only.enriched <- setdiff(bg.genes.f, bg.genes.n) %>% 
    cflo_go_enrichment(bg = bg.genes) 
# plot enrichment
for.only.enriched %>% 
  go_enrichment_plot(fdr = 5, clean = "no") # set clean = "yes" to remove text annotations

## NURSES
# save the list to a csv
cflo.annots.exp %>% 
  filter( gene_name %in% setdiff(bg.genes.n, bg.genes.f)) %>% 
  write.csv(., "./results/gene_lists/uniquely_exp_nur.csv")
# run enrichment 
nur.only.enriched <- 
   setdiff(bg.genes.n, bg.genes.f) %>%
   cflo_go_enrichment(bg = bg.genes)
# plot enrichment
nur.only.enriched %>% 
  go_enrichment_plot(fdr = 5, clean = "no")


#########################################################################################
## Figure 2 ## "TC5_heatmaps2"
############################
# perform hierarchical clustering of for-24h and nur-24h into four clusters;
# plot time-course heatmaps for the clustered for-24h and nur-24h genesets
# Identify the day-peaking and night-peaking clusters visually.
############################

## Load all the rhythmic genesets 
## Note, ordered according to their p-value; highly rhythmic at the top.
# Circadian genes (period = 24h)
tbl(db, "ejtk_all") %>% head()
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
tbl(db, "ejtk_8h_all") %>% head()
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
tbl(db, "ejtk_12h_all") %>% head()
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

# Make annotations for the heatmaps
my_gene_col.for <- cutree(tree = as.dendrogram(my_hclust_gene.for), k = 4) # four clusters
my_gene_col.for <- data.frame(overlap = ifelse(rownames(cflo.rhy.exp.for) %in% nur24, "yes", "no"),
                              cluster = my_gene_col.for)

my_gene_col.nur <- cutree(tree = as.dendrogram(my_hclust_gene.nur), k = 4)
my_gene_col.nur <- data.frame(overlap = ifelse(rownames(cflo.rhy.exp.nur) %in% for24, "yes", "no"),
                              cluster = my_gene_col.nur)


# I’ll add some column annotations and create the heatmap.
# Annotations for:
# 1. Is the sample collected during the light or dark phase? 
my_sample_col.for <- data.frame(phase = rep(c("light", "dark", "light"), c(5,6,1)))
row.names(my_sample_col.for) <- colnames(cflo.rhy.exp.for)
my_sample_col.nur <- data.frame(phase = rep(c("light", "dark", "light"), c(5,6,1)))
row.names(my_sample_col.nur) <- colnames(cflo.rhy.exp.nur)

# Manual color palette
my_colour.for = list(
  phase = c(light = "#F2E205", dark = "#010440"),
  cluster = viridis::cividis(100)[c(10,90,60,30)],
  overlap = c(no = "white", yes = "#D92818"))
my_colour.nur = list(
  phase = c(light = "#F2E205", dark = "#010440"),
  # cluster = viridis::cividis(100)[c(90,60,30,10)],
  overlap = c(no = "white", yes = "#D92818"))

# Color scale
my.breaks = seq(-3, max(cflo.rhy.exp.for), by=0.06)

# Let's plot!
for.rhy.heat <- pheatmap(cflo.rhy.exp.for, show_rownames = F, show_colnames = F,
                         annotation_row = my_gene_col.for[,c("cluster","overlap")], 
                         annotation_col = my_sample_col.for,
                         cutree_rows = 4,
                         cutree_cols = 2,
                         annotation_colors = my_colour.for,
                         border_color=FALSE,
                         cluster_cols = F,
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

## The following code will sort the clusters so we have the night-peaking genes at the top
# install.packages("dendsort")
library(dendsort)
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
mat_cluster_cols <- sort_hclust(my_hclust_gene.nur)
# plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")


nur.rhy.heat <- pheatmap(cflo.rhy.exp.nur, show_rownames = F, show_colnames = F,
                         annotation_row = my_gene_col.nur %>% select(cluster, overlap), 
                         annotation_col = my_sample_col.nur,
                         cutree_rows = 2,
                         cutree_cols = 2,
                         annotation_colors = my_colour.nur,
                         border_color=FALSE,
                         cluster_cols = F,
                         # cluster_rows = mat_cluster_cols,
                         breaks = my.breaks,
                         ## color scheme borrowed from: 
                         color = inferno(length(my.breaks) - 1),
                         # treeheight_row = 0, 
                         # treeheight_col = 0,
                         # remove the color scale or not
                         # main = paste0("Nurses - circadian genes \n (n=", nrow(cflo.rhy.exp.nur), " genes)"),
                         ## annotation legend
                         annotation_legend = T,
                         ## Color scale
                         legend = T)

# Note, the day-peaking cluster in foragers is cluster 2, but in nurses it is cluster 1
# The night-peaking cluster in foragers is cluster 1, but in nurses it's cluster 2 & 3

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
my_gene_col.for24.nur24 <- cutree(tree = as.dendrogram(my_hclust_gene.for24.nur24), 
                                 k = 4)
my_gene_col.for24.nur24 <- data.frame(cluster = my_gene_col.for24.nur24)

## For combined plot, use the following palatte
my_colour = list(
  caste = c(foragers = "#F23030", nurses = "#1A80D9"),
  phase = c(light = "#F2E205", dark = "#010440")
  # cluster = viridis::cividis(100)[c(10,90,60,30)]
  # overlap = c(no = "white", yes = "#D92818")
)

for24.nur24.rhy.heat <- 
  for24.nur24 %>% 
  pheatmap(show_rownames = F, show_colnames = F,
           # I want to see which FORAGER cluster is the gene from
           annotation_row = my_gene_col.for24.nur24 %>% select(cluster),
           # both light phase and ant caste identity
           # annotation_col = my_sample_col,
           # cutree_rows = 4,
           # cutree_cols = 2,
           annotation_colors = my_colour,
           border_color=FALSE,
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

# Note, Cluster 1 is the day-peaking gene cluster in for-24h-nur-24h.
# Perform enrichment of the day-peaking shared circadian genes:
for24.nur24.daypeaking.cluster1 <- 
  my_gene_col.for24.nur24 %>%
  rownames_to_column(var = "gene") %>%
  filter(cluster %in% c(1)) %>% 
  pull(gene) %>% 
  cflo_go_enrichment(bg=bg.genes) 

# Plot the enrichment results
for24.nur24.daypeaking.cluster1 %>% 
  go_enrichment_plot()

#########################################################################################

#########################################################################################
## Supplementary Figure 4 ## "TC5_heatmaps2::line_329"
############################
# Run enrichment for for circadian genes in foragers (for-24h) and nurses (nur-24h) 
# that peak during the day (day-peaking clusters) and night (night-peaking clusters);
# save the results as csvs & plot the enrichment results.
# Also, plot the z-scores for the different gene clusters.
############################

### FORAGERS ###

## day-peaking | cluster 2 ##
for24.daypeaking.cluster2 <- 
  my_gene_col.for %>%
  rownames_to_column(var = "gene") %>%
  filter(cluster == 2) %>%
  pull(gene) %>%
  cflo_go_enrichment(bg = bg.genes)
  
# save results #
write.csv(for24.daypeaking.cluster2, file = "./results/enrichment_results/rhythmic_genes/for24_daypeaking_cluster2.csv")
# save plot with enrichment results
png("./results/enrichment_results/rhythmic_genes/for24_daypeaking_cluster2.png",
     width = 18, height = 30, units = "cm",
     res = 300)
  # plot enrichment results   
  for24.daypeaking.cluster2 %>%
   go_enrichment_plot()
dev.off()


## night-peaking | cluster 1 ##
for24.nightpeaking.cluster1 <- 
  my_gene_col.for %>%
  rownames_to_column(var = "gene") %>%
  filter(cluster == 1) %>%
  pull(gene) %>%
  cflo_go_enrichment(bg = bg.genes)

# save results #
write.csv(for24.nightpeaking.cluster1, file = "./results/enrichment_results/rhythmic_genes/for24_nightpeaking_cluster1.csv")
png("./results/enrichment_results/rhythmic_genes/for24_nightpeaking_cluster1.png",
    width = 18, height = 30, units = "cm",
    res = 300)
	for24.nightpeaking.cluster1 %>%
	  go_enrichment_plot()
dev.off()


### NURSES ###

## day-peaking | cluster 1 ##
nur24.daypeaking.cluster1 <- 
  my_gene_col.nur %>% 
  rownames_to_column(var = "gene") %>% 
  filter(cluster %in% c(1)) %>% 
  pull(gene) %>% 
  cflo_go_enrichment(bg = bg.genes)

# save results #
write.csv(nur24.daypeaking.cluster1, file = "./results/enrichment_results/rhythmic_genes/nur24_daypeaking_cluster1.csv")
png("./results/enrichment_results/rhythmic_genes/nur24_daypeaking_cluster1.png",
    width = 18, height = 8, units = "cm",
    res = 300)
 nur24.daypeaking.cluster1 %>%
  go_enrichment_plot()
dev.off()


## night-peaking | cluster 2 ##
nur24.nightpeaking.cluster2 <- 
  my_gene_col.nur %>% 
  rownames_to_column(var = "gene") %>% 
  filter(cluster %in% c(2)) %>%
  pull(gene) %>% 
  cflo_go_enrichment(bg = bg.genes)

## night-peaking | cluster 3 ##
nur24.nightpeaking.cluster3 <- 
  my_gene_col.nur %>% 
  rownames_to_column(var = "gene") %>% 
  filter(cluster %in% c(3)) %>%
  pull(gene) %>% 
  cflo_go_enrichment(bg = bg.genes)

# save results #
write.csv(nur24.nightpeaking.cluster2, file = "./results/enrichment_results/rhythmic_genes/nur24_nightpeaking_cluster2.csv")
write.csv(nur24.nightpeaking.cluster3, file = "./results/enrichment_results/rhythmic_genes/nur24_nightpeaking_cluster3.csv")

png("./results/enrichment_results/rhythmic_genes/nur24_nightpeaking_cluster3.png",
    width = 18, height = 8, units = "cm",
    res = 300)
 nur24.nightpeaking.cluster3 %>%
  go_enrichment_plot()
dev.off()

#########################################################################################



