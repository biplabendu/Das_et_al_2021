# 1 Detecting DEGs -------------------------------------------------------------

# DEGs = Differentially Expressed Genes (while accounting for the time course fluctuations)

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


## Load files --------------
# - database (db)
db <- dbConnect(RSQLite::SQLite(), "./data/TC5_data.db")
src_dbi(db)
# - db::cflo.annots
cflo.annots.exp <- tbl(db, "annot_fpkm") %>% collect()
# - db::expressed_genes (expressed for at least half time points)
gs1 <- tbl(db, "expressed_genes") %>% filter(exp_half_for=="yes") %>% collect() %>% pull(gene_name)
gs2 <- tbl(db, "expressed_genes") %>% filter(exp_half_nur=="yes") %>% collect() %>% pull(gene_name)


# Let's load all functions:
# # limorhyde functions
# source(system.file('extdata', 'vignette_functions.R', package = 'limorhyde'))

# my functions
source("./functions/plot_zscores.R")
# source("~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/vennDia.R")
source("./functions/enrichment_analysis.R")
source("./functions/theme_publication.R")


# ## 1.1 LimoRhyde --------------------------------------
# 
# # Code borrowed from: http://127.0.0.1:20841/library/limorhyde/doc/introduction.html
# 
# # Let's load the libraries required for running Limorhyde
# # library('annotate')
# library('data.table')
# # library('foreach')
# # library('GEOquery')
# 
# # library('knitr')
# library('limma')
# library('limorhyde')
# 
# ### 1.1.1 Format the meta-data ----------------
# # load the meta-data
# sm.tc5 <- read.csv("./data/TC5_meta_data.csv",
#                    stringsAsFactors = F)
# 
# # Let's format the columns in the right data-type
# sm.tc5$genotype <- as.factor(sm.tc5$genotype)
# sm.tc5$cond <- as.factor(sm.tc5$cond)
# sm.tc5$time <- as.numeric(sm.tc5$time)
# sm.tc5$batch <- as.factor(sm.tc5$batch)
# sm.tc5$LD <- as.factor(sm.tc5$LD)
# sm.tc5$location <- as.factor(sm.tc5$location)
# 
# # Let's get a glimpse of the metadata
# sm.tc5 %>%
#   glimpse()
# 
# # Next we use limorhyde to calculate time_cos and time_sin, which are based on the first
# #harmonic of a Fourier decomposition of the time column, and append them to the sm data frame.
# sm.tc5 = cbind(sm.tc5, limorhyde(sm.tc5$time, 'time_'))
# # convert the dataframe into a data.table
# sm.tc5 <- data.table(sm.tc5)
# # check that it worked
# sm.tc5[1:5, ]
# 
# ### 1.1.2 Format the expression data ----------------
# 
# c.for <- cflo.annots.exp %>%
#   dplyr::select(gene_name,X2F:X24F) %>%
#   na.omit() %>%
#   # remove genes that have no expression for any sample (for and nur)
#   filter_at(vars(starts_with("X")), any_vars(. > 0)) %>%
#   # keep only the genes that have at least 1 fpkm for at least half the timepoints
#   filter(gene_name %in% unique(c(gs1,gs2))) # gs1 - foragers; gs2 - nurses
# 
# c.nur <- cflo.annots.exp %>%
#   dplyr::select(gene_name,X2N:X24N) %>%
#   na.omit() %>%
#   # remove genes that have no expression for any sample (for and nur)
#   filter_at(vars(starts_with("X")), any_vars(. > 0)) %>%
#   # keep only the genes that have at least 1 fpkm for at least half the timepoints
#   filter(gene_name %in% unique(c(gs1,gs2)))
# 
# 
# # rename the columns to match the names in the meta-data
# #Foragers
# oldnames.for <- names(c.for[,-1])
# newnames.for <- c("2F","4F","6F","8F","10F","12F","14F","16F","18F","20F","22F","24F")
# emat.for <- c.for %>%
#   rename_at(vars(oldnames.for), ~ newnames.for)
# # take a look if it worked
# emat.for[1:5,]
# 
# #Nurses
# oldnames.nur <- names(c.nur[,-1])
# newnames.nur <- c("2N","4N","6N","8N","10N","12N","14N","16N","18N","20N","22N","24N")
# emat.nur <- c.nur %>%
#   rename_at(vars(oldnames.nur), ~ newnames.nur)
# # take alook if it worked
# emat.nur[1:5,]
# 
# # Merge the for and nur files into one
# #Merge emat.for and emat.nur by gene_name
# emat.tc5 <- emat.for %>%
#   # keep only the genes that are expressed in both foragers and nurses
#   # (Note: there are 4 genes that are expressed in foragers but not in nurses)
#   filter(gene_name %in% emat.nur$gene_name) %>%
#   left_join(emat.nur, by = "gene_name")
# # Did it work?
# head(emat.tc5)
# 
# 
# 
# # Create the log2 transformed input matrix for identifying DEGs ----------------------------
# emat.tc5 <- as.data.frame(emat.tc5)
# rownames(emat.tc5) <- emat.tc5[,1]
# emat.tc5 <- emat.tc5[,-1]
# summary(emat.tc5)
# # Let's remove NAs
# emat.tc5 <- na.omit(emat.tc5) # n(genes) = 9,361
# 
# # Need to make the emat.tc5 into a matrix.
# emat.tc5 <- data.matrix(emat.tc5)
# emat.tc5[1:5,1:5]
# class(emat.tc5) # is a matrix with samples=columns and genes=rows
# 
# ## log2 transform the data
# emat.tc5 <- log2(emat.tc5 + 1)
# 
# 
# # 1.2 Identify DEGs ----------------------------------------------
# # Identify DEGs:
# ### Differential expression is based on the coefficient for cond in a linear model with no 
# ### interaction terms. We pass all genes to limma, but keep results only for non-differentially 
# ### rhythmic genes, and adjust for multiple testing accordingly.
# 
# Set threshold for q-value (BH adjusted p-value)
q.threshold <- 0.05 # currently, using 5% FDR
log2.foldchange <- 1 # thus, any gene with a 2^(log2.foldchange) fold change in it's expression
# 
# # 1.2.1 all rhy genes -----------------------------------------------
# 
# # Use the subsetted emat.tc5.all.rhy to find DEGs
# design.deg = model.matrix(~ cond + time_cos + time_sin, data = sm.tc5)
# 
# fit = lmFit(emat.tc5, design.deg)
# fit = eBayes(fit, trend = TRUE)
# # Take a look at the coefficients table
# fit$coefficients %>% head()
# 
# deLimma.deg = data.table(topTable(fit, coef = 2, number = Inf), keep.rownames = TRUE)
# setnames(deLimma.deg, 'rn', 'gene_name')
# deLimma.deg[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
# setorderv(deLimma.deg, 'adj.P.Val')
# 
# deLimma.deg %>%
#   arrange(adj.P.Val) %>%
#   head()


# Save the results of the DEG analysis 
# save(deLimma.deg, 
#      file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/DEGs/TC5_DEG_results.RData")
##
## Done.

# Load the results of the DEG analysis (dataframe: deLimma.deg)
load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/DEGs/TC5_DEG_results.RData")

# Filter the results to keep only the significant genes
tc5.all.DEGs <- 
  deLimma.deg %>% 
  arrange(adj.P.Val) %>% 
  mutate(sig = as.factor(ifelse(adj.P.Val < q.threshold & abs(logFC) >= log2.foldchange, "yes", "no"))) %>% 
  mutate(for_direction = as.factor(ifelse(sig=="yes", ifelse( logFC > 0, "down", "up" ), "NA")))

head(tc5.all.DEGs)

## How many DEGs - 5% FDR and ≥ 1 fold change in gene expression
tc5.all.DEGs %>% 
  filter(adj.P.Val < q.threshold) %>%
  filter(abs(logFC) >= 1) %>% # change the criteria here for top DEG or all DEG (logFC≥1)
  filter(sig == "yes") %>% 
  group_by(for_direction) %>% 
  summarise(total = n())
# n = 72 DEGs (1% FDR; 1 log2 fold change)
## n = 81 DEGs (5% FDR; 1 log2 fold change) ## THIS IS WHAT WE ARE GOING WITH
# n = 86 DEGs (10% FDR; 1 log2 fold change)



# Volcano plot -------
library(viridis)
# png("~/OneDrive - University of Central Florida/BD-TC5/TC5_Paper/results_figs/volcano_for_nur_DEGs.png",
#     width = 14, height = 10, units = "cm",
#     res = 300)
ggplot(tc5.all.DEGs) +
  #geom_hline(yintercept = -log10(0.05), col="red", alpha=0.6) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), color=sig), size = 1.5, alpha = 0.5) +
  labs(x = expression(log[2]*' fold-change'), y = expression(-log[10]*' '*q[DE]),
       title = "") +
  # xlim(c(-50,50)) +
  theme_Publication() +
  scale_color_viridis(discrete = T, direction = -1, option = "viridis") 
# annotation_custom(grid::textGrob(paste0( nrow(tc5.all.DEGs[tc5.all.DEGs$sig == "yes", ]),
#                                          " DEGs ",
#                                         "(tested: ", nrow(tc5.all.DEGs),"genes)")),
#                   xmin = -4, xmax=2, ymin = 7, ymax = 8)
# dev.off()

# Enrichment analysis | all the genes --------
tc5.all.DEGs %>% 
  filter(sig == "yes") %>% 
  pull(gene_name) %>% 
  cflo_go_enrichment() %>% 
  go_enrichment_plot(fdr = 5, clean = "no")

sig.DEGs.tc5.all <-
  tc5.all.DEGs %>% 
  filter(sig=="yes") %>% 
  arrange(desc(abs(logFC)))

sig.DEGs.tc5.all %>%
  pull(gene_name) %>% 
  head(50) %>%
  cflo_heatmap(show_rownames = T)

###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-
#-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-##-###-##-###-##-###-##
###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-


# Plot the fold-change with column annotations ----------------------------

### How many of the DEGs are also rhythmic genes

## Load all the rhythmic genesets 

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


for.rhy.genes <- unique(c(for24,for12,for8))
nur.rhy.genes <- unique(c(nur24,nur12,nur8))


# Fig 13: Upset - DEGs, rhythmic genes ------------------------------------
# example of list input (list of named vectors)
listInput <- list(
  for_DEG = tc5.all.DEGs %>% filter(for_direction == "up") %>% pull(gene_name),
  nur_DEG = tc5.all.DEGs %>% filter(for_direction == "down") %>% pull(gene_name),
  for_Rhy = for.rhy.genes,
  nur_Rhy = nur.rhy.genes
)

library(UpSetR)
library(viridis)
caste.col <- c("#F23030","#1A80D9")
upset(fromList(listInput), order.by = "freq")



# Plot heatmap with gene name and annotations -----------------------------

## make data.frame with logfold data
sig.DEGs.tc5.all %>% glimpse()
sig.DEGs.tc5.all %>% tail() # ordered in desc by the abs(logFC)

## Make a custom annot that pastes LOC# and blast annotation
# load gene to annotation file
load(file = "./functions/func_data/gene_to_annot.RData")
gene_to_annot %>% head()

# # Let's write the DEG gene_names and annotations to a file, so we can simplify blast annotations
# sig.DEGs.tc5.all %>%
#   select(gene_name) %>% 
#   left_join(gene_to_annot, by="gene_name") %>% 
#   write.csv(., "./results/gene_lists/for_nur_DEG_all.csv")

# let's read the simplified annotation file 
for_nur_DEG <- read.csv("./results/gene_lists/for_nur_DEG_all_2.csv", header = T, stringsAsFactors = F)
for_nur_DEG <- for_nur_DEG[-1]

deg.dat <- 
  sig.DEGs.tc5.all %>% 
  select(gene_name, logFC, for_direction) %>% 
  mutate(abslogFC = (abs(logFC))) %>%
  left_join(for_nur_DEG, by="gene_name") %>% 
  # mutate(annot1 = ifelse(gene_name==annot, NA,annot)) %>% 
  mutate(annot2 = ifelse(is.na(annot),gene_name,
                         stringr::str_trunc(paste0(gene_name, ", ", annot), 
                                            60)))  # how many characters do you want to keep


## Specify how many genes you want to do
n.degs <- nrow(deg.dat)
# n.degs <- 50
dat <- deg.dat %>% 
  head(n.degs)

# Load genes that code for the most abundant TF in Cflo
tf.genes <- read.csv("./data/most_abundant_TF_Cflo.csv", header = F, stringsAsFactors = F, na.strings = c("NA"))
tf.genes <- na.omit(tf.genes)
tf.genes <- tf.genes[,1, drop=T]

# Make annotations for the heatmaps
my_gene_col.for <- data.frame(
   Rhy = ifelse(dat$gene_name %in% c(for.rhy.genes,nur.rhy.genes), "yes", "no"),
   # Rhy = ifelse(dat$gene_name %in% c(nur24,nur12,nur8), "yes", "no"),
   DRG = ifelse(dat$gene_name %in% c(intersect(for24,nur12),
                                      intersect(for24,nur8),
                                      intersect(for12,nur24),
                                      intersect(for12,nur8),
                                      intersect(for8,nur24),
                                      intersect(for8,nur12)), "yes","no"),
   Troph = ifelse(dat$gene_name %in% tf.genes, "yes", "no"),
   DEG = ifelse(dat$for_direction=="down","nurse","forager")) %>% 
  select(Troph, DRG, Rhy, DEG)

row.names(my_gene_col.for) <- dat$annot2

# Change rownames to unique annot2 values
dat <- dat %>% select(annot2, abslogFC)
rownames(dat) <- dat$annot2
dat <- as.matrix(dat[-1])

# Manual color palette
my_colour.for = list(
  Rhy = c(yes = "#D9BF73", no = "white"),
  # Rhy_Nur = c(yes = "#F2E205", no = "#010440"),
  DRG = c(yes = "#A66249", no = "white"),
  Troph = c(yes = "#6C8C26", no = "white"),
  DEG = c(nurse = "#1A80D9", forager = "#F23030"))

# Color scale
my.breaks = seq(0, max(dat), by=0.06)

# Let's plot!
deg.all <- pheatmap(dat, 
               show_rownames = T, 
               show_colnames = F,
               annotation_row = my_gene_col.for, 
               # annotation_col = my_sample_col.for,
               # cutree_rows = 4,
               # cutree_cols = 2,
               annotation_colors = my_colour.for,
               border_color=T,
               cluster_cols = F,
               cluster_rows = F,
               breaks = my.breaks,
               ## color scheme borrowed from: 
               color = inferno(length(my.breaks) - 1),
               # treeheight_row = 0, 
               # treeheight_col = 0,
               # remove the color scale or not
               # main = paste0("Top 50", " DEGs"),
               ## annotation legend
               annotation_legend = T,
               ## cell width
               cellwidth = 20,
               cellheight = 10,
               ## Color scale
               legend = T)

save_pheatmap_png(deg.all, 
                  filename = "./results/figures/figure_5/all_degs_5.png",
                  width = 2000,
                  height = 3600,
                  res = 300)


# Plot daily rhythms of DEGs ----------------------------------------------

## All DRGs and Vitellogenin

deg.dat %>% glimpse()

drgs <- c(intersect(for24,nur12),
               intersect(for24,nur8),
               intersect(for12,nur24),
               intersect(for12,nur8),
               intersect(for8,nur24),
               intersect(for8,nur12))

degs.drgs <- deg.dat %>% 
  filter(gene_name %in% drgs) %>% 
  pull(gene_name)
# png("./results/degs_drgs.png", width = 2000, height = 1800, res = 300)
# degs.drgs %>% 
#   zplot() %>% 
#   pluck(1) %>% 
#   multi.plot(rows = 3, cols = 2)
# dev.off()

degs.drgs.vg <- c(degs.drgs,"LOC105254427")
# png("./results/degs_drgs_vg.png", width = 2000, height = 1800, res = 300)
# degs.drgs.vg %>% 
#   zplot() %>% 
#   pluck(1) %>% 
#   multi.plot(rows = 3, cols = 2)
# dev.off()

## all the degs that have been discussed in text

degs.in.text <- c("LOC105250772", # arylphorin subunit alpha
                  "LOC105256665") # med12

degs.drgs.vg.intext <- c(degs.drgs.vg, degs.in.text)                  
# png("./results/degs_drgs_vg_intext.png", width = 2000, height = 2400, res = 300)
degs.drgs.vg.intext %>% 
  zplot() %>% 
  pluck(1) %>% 
  multi.plot(rows = 4, cols = 2)
# dev.off()


# Volcano plot with gene annotations --------------------------------------
png("./results/volcano_DEGs.png",
    width = 14, height = 10, units = "cm",
    res = 300)
  
  ggplot(tc5.all.DEGs) +
  
    # plot all points and indicate significant ones
    geom_point(aes(x = logFC, y = -log10(adj.P.Val), color=sig), size = 1.5, alpha = 0.7) +
  
    # highlight the genes that have been plotted separately
    geom_point(data = (tc5.all.DEGs %>% filter(gene_name %in% degs.drgs.vg.intext)),
               aes(x=logFC, y=-log10(adj.P.Val)),
               shape=24,
               # fill="#8850BF",
               color="#304019",
               size=4,
               alpha=1) +
  
    # add axis-labels  
    labs(x = expression(log[2]*' fold-change'), y = expression(-log[10]*' '*q[DE]),
         title = "") +
    # xlim(c(-50,50)) +
    theme_Publication() +
    scale_color_viridis(discrete = T, direction = -1, option = "viridis") 
    # # add gene annotation
    # annotate(geom = "text", 
    #          x = -4, 
    #          y = 5, 
    #          label = "italic(Cdk4)", parse=T,
    #          hjust = "left") +
    # ggrepel::geom_text_repel(data = . %>% filter(gene_name %in% degs.drgs.vg.intext),
    #                         aes(x=logFC, y=-log10(adj.P.Val), 
    #                             label=gene_name),
    #                             # label= (data.frame(gene_name=degs.drgs.vg.intext) %>% 
    #                             #           left_join(gene_to_annot, by="gene_name") %>% 
    #                             #           pull(annot) %>% stringr::str_trunc(40))),
    #                         # arrow = grid::arrow(length = unit(0.02, "npc")),
    #                         box.padding = 1,
    #                         # direction = "x",
    #                         show.legend = F,
    #                         alpha=0.9) 
dev.off()  
    