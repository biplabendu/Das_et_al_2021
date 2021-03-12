rm(list = ls())


# Housekeeping ------------------------------------------------------------

pacman::p_load(pheatmap, dendextend, tidyverse, viridis)
pacman::p_load(RSQLite, tidyverse, dbplyr, DT, conflicted)

# set conflict preference
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("layout", "plotly")

db <- dbConnect(RSQLite::SQLite(), 
                "./data/TC5_data.db")
src_dbi(db)

# Load all the rhythmic genesets 
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

# Get a glimpse
cflo.zscores.for %>% head()

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
set.seed(420)
my_hclust_gene.for <- hclust(dist(cflo.rhy.exp.for), method = "complete")
my_hclust_gene.nur <- hclust(dist(cflo.rhy.exp.nur), method = "complete")


# Make annotations for the heatmaps
my_gene_col.for <- cutree(tree = as.dendrogram(my_hclust_gene.for), k = 4) # four clusters
my_gene_col.for <- data.frame(cluster = my_gene_col.for)

my_gene_col.nur <- cutree(tree = as.dendrogram(my_hclust_gene.nur), k = 4)
my_gene_col.nur <- data.frame(cluster = my_gene_col.nur)


# I’ll add some column annotations and create the heatmap.

# Annotations for:
# 1. Is the sample collected during the light or dark phase? 
my_sample_col.for <- data.frame(phase = rep(c("light", "dark", "light"), c(5,6,1)))
row.names(my_sample_col.for) <- colnames(cflo.rhy.exp.for)
my_sample_col.nur <- data.frame(phase = rep(c("light", "dark", "light"), c(5,6,1)))
row.names(my_sample_col.nur) <- colnames(cflo.rhy.exp.nur)


# Specify your color palette -----
my_colour.for = list(
  phase = c(light = "white", dark = "#403F3D"),
  cluster = viridis::cividis(100)[c(10,35,65,80)])
# overlap = c(no = "white", yes = "#D92818"))
my_colour.nur = list(
  phase = c(light = "white", dark = "#403F3D"),
  cluster = viridis::viridis(100)[c(10,35,65,80)])
# overlap = c(no = "white", yes = "#D92818"))

# Color scale
my.breaks = seq(-3, max(cflo.rhy.exp.for), by=0.06)

# Let's plot!
png("./results/for24_all.png", width = 1000, height = 1600, res = 300)
for.rhy.heat <- pheatmap(cflo.rhy.exp.for, show_rownames = F, show_colnames = F,
                         annotation_row = my_gene_col.for, 
                         annotation_col = my_sample_col.for,
                         cutree_rows = 4,
                         cutree_cols = 2,
                         annotation_colors = my_colour.for,
                         border_color=T,
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
dev.off()

## The following code will add border to your Phase column annotation
# library(grid)
# grid.ls(grid.force()) # "col_annotation" looks like it's the one to edit
# grid.gedit("col_annotation", gp = gpar(col="grey60"))




## The following code will sort the clusters so we have the night-peaking genes at the top
# install.packages("dendsort")
library(dendsort)
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
mat_cluster_cols <- sort_hclust(my_hclust_gene.nur)
# plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")

png("./results/nur24_all.png", width = 1000, height = 1600, res = 300)
nur.rhy.heat <- pheatmap(cflo.rhy.exp.nur, show_rownames = F, show_colnames = F,
                         annotation_row = my_gene_col.nur, 
                         annotation_col = my_sample_col.nur,
                         cutree_rows = 2,
                         cutree_cols = 2,
                         annotation_colors = my_colour.nur,
                         border_color=FALSE,
                         cluster_cols = F,
                         cluster_rows = mat_cluster_cols,
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
dev.off()
# Note, the day-peaking cluster in foragers is 2, but in nurses it is cluster 1
# Note, cluster 2 and 3 might be a single cluster, but not as strong as forager night-peaking cluster 

## for-24h-nur-24h ----------------------------------------------
intersect(for24, nur24) %>% 
  cflo_heatmap(annotation = F, cluster.r = T)

# Subset the intersecting geneset
for24.nur24 <- cflo.zscores.for %>% 
  filter(gene_name %in% intersect(for24, nur24)) %>% 
  left_join(cflo.zscores.nur) %>% 
  as.data.frame()

# Set genes as rownames and convert it into a matrix
for24.nur24 <- as.matrix(for24.nur24[-1])  
# cluster
my_hclust_gene.for24.nur24 <- hclust(dist(for24.nur24), 
                                     method = "complete")

# Let's make the column annotation data frame
my_sample_col <- data.frame(caste = rep(c("foragers","nurses"), c(12,12)),
                            phase = rep(rep(c("light", "dark", "light"), c(5,6,1)),2))
row.names(my_sample_col) <- colnames(for24.nur24)

## For combined plot, use the following palatte
my_colour = list(
  caste = c(foragers = "#F23030", nurses = "#1A80D9"),
  phase = c(light = "white", dark = "#403F3D")
  # cluster = viridis::cividis(100)[c(10,90,60,30)]
  # overlap = c(no = "white", yes = "#D92818")
)

# Color scale
my.breaks = seq(-3, max(for24.nur24), by=0.06)

png("./results/for24nur24_all.png", width = 1600, height = 1200, res = 300)
for.nur.rhy.heat <- 
  for24.nur24 %>% 
  pheatmap(show_rownames = F, show_colnames = F,
           # I want to see which FORAGER cluster is the gene from
           # annotation_row = my_gene_col.for %>% select(cluster),
           # both light phase and ant caste identity
           annotation_col = my_sample_col,
           # cutree_rows = 4,
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


# Enrichment - for24nur24 --------------------------------------------------------------

# Get the enrichment for all for-24h-nur-24h
# cluster 1(+maybe 3) - day-peaking; cluster 2(+maybe 4) - night-peaking; 
for24.nur24.daypeaking.cluster1 <- 
  my_gene_col.for24.nur24 %>%
  rownames_to_column(var = "gene") %>%
  filter(cluster %in% c(1)) %>% 
  pull(gene) %>% 
  cflo_go_enrichment(bg=bg.genes)

## Supplementary Excel 1
write.csv(for24.nur24.daypeaking.cluster1,
          file = "~/OneDrive - University of Central Florida/BD-TC5/TC5_Paper/results_csv/enrichment_results/go_enrichments/rhythmic_genes/For24_Nur24/for24nur24_daypeaking_cluster1.csv")

# png("~/OneDrive - University of Central Florida/BD-TC5/TC5_Paper/results_figs/enrichments/for24_nur24_day_cluster1.png",
#     width = 14, height = 10, units = "cm",
#     res = 300)
# for24.nur24.daypeaking.cluster1 %>% 
#   go_enrichment_plot()
# dev.off()


# Plot for24nur8 ----------------------------------------------------------

# Subset the intersecting geneset
for24.nur8 <- cflo.zscores.for %>% 
  filter(gene_name %in% intersect(for24, nur8)) %>% 
  left_join(cflo.zscores.nur) %>% 
  as.data.frame()

# Set genes as rownames and convert it into a matrix
rownames(for24.nur8) <- for24.nur8[[1]]
for24.nur8 <- as.matrix(for24.nur8[-1])  
# cluster
my_hclust_gene.for24.nur8 <- hclust(dist(for24.nur8), method = "complete")

# Let's make the column annotation data frame
my_sample_col <- data.frame(caste = rep(c("foragers","nurses"), c(12,12)),
                            phase = rep(rep(c("light", "dark", "light"), c(5,6,1)),2))
row.names(my_sample_col) <- colnames(for24.nur8)

## For combined plot, use the following palatte
my_colour = list(
  caste = c(foragers = "#F23030", nurses = "#1A80D9"),
  phase = c(light = "white", dark = "#403F3D"),
  cluster = viridis::cividis(100)[c(10,35,65,80)]
  # overlap = c(no = "white", yes = "#D92818")
)

# Color scale
my.breaks = seq(min(for24.nur8), max(for24.nur8), by=0.06)

png("./results/for24nur8_all.png", width = 800, height = 1400, res = 300)
for24nur8.heat <- 
  for24.nur8 %>% 
  pheatmap(show_rownames = F, 
           show_colnames = F,
           # I want to see which FORAGER cluster is the gene from
           annotation_row = my_gene_col.for %>% select(cluster),
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
           annotation_legend = F,
           ## Color scale
           legend = T)
dev.off()


# Plot daily oscillations of period ---------------------------------------

# Load the zplot function
source(file="./functions/plot_zscores.R")

png("./results/period.png", width = 1600, height = 800, res = 300)
c("LOC105256454") %>% # this is the gene_name for Period
  zplot(lwd=2, alpha=0.8)
dev.off()


# Additional File XX ------------------------------------------------------

# plot daily expressions for the following genes:

# 1. Dbt
png("./results/Dbt.png", width = 1600, height = 800, res = 300)
c("LOC105255207") %>% # this is the gene_name for Period
  zplot(lwd=2, alpha=0.8)
dev.off()

# 2. Nemo
png("./results/Nemo.png", width = 1600, height = 800, res = 300)
c("LOC105248529") %>% # this is the gene_name for Period
  zplot(lwd=2, alpha=0.8)
dev.off()

# 3. Foraging
png("./results/Foraging.png", width = 1600, height = 800, res = 300)
c("LOC105255628") %>% # this is the gene_name for Period
  zplot(lwd=2, alpha=0.8)
dev.off()

# 4. Slowpoke
png("./results/Slowpoke.png", width = 1600, height = 800, res = 300)
c("LOC105258647") %>% # this is the gene_name for Period
  zplot(lwd=2, alpha=0.8)
dev.off()
 