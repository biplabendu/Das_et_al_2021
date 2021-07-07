
# NOTE: DID NOT use for rebuttal.

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

# - db::cflo.annots
cflo.annots.exp <- tbl(db, "annot_fpkm") %>% collect()

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
bg.genes <- intersect(bg.genes.f, bg.genes.n)


## load zscore datasets
cflo.zscores.for <- tbl(db, "zscore_for") %>% collect()
cflo.zscores.nur <- tbl(db, "zscore_nur") %>% collect()

# Get a glimpse
cflo.zscores.for %>% dim()

# Filter the zscores to keep only circadian genes
cflo.rhy.exp.for <- cflo.zscores.for %>% 
  filter(gene_name %in% bg.genes) %>% as.data.frame()
  # filter(gene_name %in% intersect(for24,nur24)) %>% as.data.frame()
cflo.rhy.exp.nur <- cflo.zscores.nur %>% 
  filter(gene_name %in% bg.genes) %>% as.data.frame()
  # filter(gene_name %in% intersect(for24,nur24)) %>% as.data.frame()


cflo.dat <- cflo.rhy.exp.for %>% 
  left_join(cflo.rhy.exp.nur, by="gene_name")


# Set genes as rownames and convert it into a matrix
# rownames(cflo.rhy.exp.for) = cflo.rhy.exp.for$gene_name
# cflo.rhy.exp.for <- as.matrix(cflo.rhy.exp.for[-1])
# rownames(cflo.rhy.exp.nur) = cflo.rhy.exp.nur$gene_name
# cflo.rhy.exp.nur <- as.matrix(cflo.rhy.exp.nur[-1])

rownames(cflo.dat) = cflo.dat$gene_name
cflo.dat <- as.matrix(cflo.dat[-1])

# # Hierarchical clustering of the for24 and nur24 genesets
# set.seed(420)
# my_hclust_gene.for <- hclust(dist(cflo.rhy.exp.for), method = "complete")
# my_hclust_gene.nur <- hclust(dist(cflo.rhy.exp.nur), method = "complete")
# 
# 
# 
# # Make annotations for the heatmaps
# my_gene_col.for <- cutree(tree = as.dendrogram(my_hclust_gene.for), k = 4) # four clusters
# my_gene_col.for <- data.frame(cluster = my_gene_col.for)
# 
# my_gene_col.nur <- cutree(tree = as.dendrogram(my_hclust_gene.nur), k = 4)
# my_gene_col.nur <- data.frame(cluster = my_gene_col.nur)
# 
# 
# # I’ll add some column annotations and create the heatmap.
# 
# # Annotations for:
# # 1. Is the sample collected during the light or dark phase? 
# my_sample_col.for <- data.frame(phase = rep(c("light", "dark", "light"), c(5,6,1)))
# row.names(my_sample_col.for) <- colnames(cflo.rhy.exp.for)
# my_sample_col.nur <- data.frame(phase = rep(c("light", "dark", "light"), c(5,6,1)))
# row.names(my_sample_col.nur) <- colnames(cflo.rhy.exp.nur)
# 
# 
# # Specify your color palette -----
# my_colour.for = list(
#   phase = c(light = "white", dark = "#403F3D"),
#   cluster = viridis::cividis(100)[c(10,35,65,80)])
# # overlap = c(no = "white", yes = "#D92818"))
# my_colour.nur = list(
#   phase = c(light = "white", dark = "#403F3D"),
#   cluster = viridis::viridis(100)[c(10,35,65,80)])
# # overlap = c(no = "white", yes = "#D92818"))




# Color scale
my.breaks = seq(min(cflo.dat), 
                max(cflo.dat), 
                by=0.06)

# Let's plot!
# png("./results/for24_all.png", width = 1000, height = 1600, res = 300)
pheatmap(cflo.dat, 
         show_rownames = F, 
         show_colnames = T,
         # annotation_row = my_gene_col.for, 
         # annotation_col = my_sample_col.for,
         # cutree_rows = 4,
         # cutree_cols = 3,
         # annotation_colors = my_colour.for,
         border_color=T,
         cluster_cols = T,
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
# dev.off()

## The following code will add border to your Phase column annotation
# library(grid)
# grid.ls(grid.force()) # "col_annotation" looks like it's the one to edit
# grid.gedit("col_annotation", gp = gpar(col="grey60"))



# Plot results of hierarchical clustering (method="complete") only
plot(hclust(dist(t(cflo.dat)), method = "complete"))
     