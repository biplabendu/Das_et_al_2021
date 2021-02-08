# 1 Detecting DEGs -------------------------------------------------------------

# DEGs = Differentially Expressed Genes (while accounting for the time course fluctuations)

rm(list = ls())

# Load libraries
pacman::p_load(conflicted, 
               tidyverse)

# Set preference for functions
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("setdiff", "BiocGenerics")
conflict_prefer("union", "BiocGenerics")
conflict_prefer("intersect", "BiocGenerics")


# Let's load the datasets:
load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/TC5_core_datasets.RData")
# load genes that are "expressed" in the Cflo brain (≥ 1 fpkm)
load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/expressed_genes/cflo_expressed_genes.RData")
# load the genes that are expressed (≥ 1 fpkm) for at least half time points
load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/expressed_genes/cflo_expressed_halftimepoints.RData")
# load the gene to annotation file
load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/gene_to_annot.RData")
# load the blastp_summary with rhythmic data
load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/reciprocal_blast/TC5_cflo_blastp_rhy_summ.RData")
## load the list of rhy genes
load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/rhy_genes/TC5_rhy_genes.RData")


# Let's load all functions:
# limorhyde functions
source(system.file('extdata', 'vignette_functions.R', package = 'limorhyde'))
# my functions
source("~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/plot_zscores.R")
source("~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/vennDia.R")
source("~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/enrichment_analysis.R")
source("~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/theme_publication.R")


## 1.1 LimoRhyde --------------------------------------

# # Code borrowed from: http://127.0.0.1:20841/library/limorhyde/doc/introduction.html
# 
# # Let's load the libraries required for running Limorhyde
# # library('annotate')
# library('data.table')
# library('foreach')
# # library('GEOquery')
# library('ggplot2')
# library('knitr')
# library('limma')
# library('limorhyde')
# 
# ### 1.1.1 Format the meta-data ----------------
# # load the meta-data
# sm.tc5 <- read.csv("/Users/biplabendudas/OneDrive - University of Central Florida/BD-TC5/1212LD_24h/LimoRhyde/TC5_meta_data.csv",
#                    stringsAsFactors = F)
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
#   filter(gene_name %in% unique(c(gs1,gs2)))
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
# oldnames.for <- names(cflo.zscores.for[,-1])
# newnames.for <- c("2F","4F","6F","8F","10F","12F","14F","16F","18F","20F","22F","24F")
# emat.for <- c.for %>%
#   rename_at(vars(oldnames.for), ~ newnames.for)
# # take a look if it worked
# emat.for[1:5,]
# 
# #Nurses
# oldnames.nur <- names(cflo.zscores.nur[,-1])
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
# 
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


# 1.2 Identify DEGs ----------------------------------------------
# Identify DEGs:
### Differential expression is based on the coefficient for cond in a linear model with no 
### interaction terms. We pass all genes to limma, but keep results only for non-differentially 
### rhythmic genes, and adjust for multiple testing accordingly.

# Set threshold for q-value (BH adjusted p-value)
q.threshold <- 0.05 # currently, using 5% FDR
log2.foldchange <- 1 # thus, any gene with a 2^(log2.foldchange) fold change in it's expression

# 1.2.1 all rhy genes -----------------------------------------------

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
  # filter(abs(logFC) >= 2) %>% # change the criteria here for top DEG or all DEG (logFC≥1)
  filter(sig == "yes") %>% 
  # pull(gene_name) %>% 
  # exp.plot(log = T) %>% 
  # pluck(1) %>% 
  # multi.plot(rows = 5, cols = 5)
  # the foraging gene
  # filter(geneId == "LOC105255628") # is not sig DE
  group_by(for_direction) %>% 
  summarise(total = n())
# n = 25 DEGs (5% FDR; 2 log2 fold change)
# n = 42 DEGs (5% FDR; 1.5 log2 fold change)
## n = 81 DEGs (5% FDR; 1 log2 fold change) ## THIS IS WHAT WE ARE GOING WITH
# n = 72 DEGs (1% FDR; 1 log2 fold change)
# n = 86 DEGs (10% FDR; 1 log2 fold change)


# Plot top DEGs -----------------------------------------------------------

# Forager-biased
png("~/OneDrive - University of Central Florida/BD-TC5/TC5_Paper/results_figs/top_DEG_for.png",
    width = 18, height = 5, units = "cm",
    res = 300)
exp.plot(c("LOC105249672", #cyclin-dependent kinase 4
           "LOC105256576"), # glycine N-methyltransferase
         log = T) %>% pluck(1) %>% 
  multi.plot(rows = 1, cols = 2)  
dev.off()  

# Nurse-biased
png("~/OneDrive - University of Central Florida/BD-TC5/TC5_Paper/results_figs/top_DEG_nur2.png",
    width = 18, height = 5, units = "cm",
    res = 300)
exp.plot(
  # c("LOC105249045", # protein CREG1
  #          "LOC105252207"), # regucalcin
  c("LOC105254427", # vitellogenin
    "LOC105253142"), # venom carboxylesterase-6
  log = T) %>% pluck(1) %>% 
  multi.plot(rows = 1, cols = 2)  
dev.off()  


# Histogram of adj p-values
ggplot(tc5.all.DEGs, aes(as.numeric(adj.P.Val))) +
  # geom_histogram(binwidth = 0.01, col="white", alpha=0.8) +
  geom_density() +
  geom_vline(xintercept = 0.05, col="red", size=1, alpha=0.8, lty=2) +
  theme_bw() +
  xlim(0,1) +
  #ylim(0,3000) +
  ggtitle("") +
  xlab("BH-corrected p-values") +
  annotation_custom(grid::textGrob(paste0("number of DEGs: ", 
                                          nrow(tc5.all.DEGs[tc5.all.DEGs$sig == "yes", ]),
                                          "/", nrow(tc5.all.DEGs))),
                    xmin = 0.1, xmax=0.9, ymin = 400, ymax = 500) +
  theme_Publication()


# Volcano plot 
library(viridis)
png("~/OneDrive - University of Central Florida/BD-TC5/TC5_Paper/results_figs/volcano_for_nur_DEGs.png",
    width = 14, height = 10, units = "cm",
    res = 300)
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
dev.off()

# Enrichment analysis | all the genes --------
tc5.all.DEGs %>% 
  filter(sig == "yes") %>% 
  pull(gene_name) %>% 
  cflo_go_enrichment() %>% 
  go_enrichment_plot(fdr = 5)
# No enriched terms at 5% FDR and 4 fold change
# no enriched terms at 5% FDR and 3 fold change


# heatmap of the top 20 DEGs | all rhy genes ------
tc5.all.DEGs %>%
  filter(sig == "yes") %>% 
  arrange(adj.P.Val) %>% 
  # head(20) %>%
  pull(gene_name) %>% 
  cflo_heatmap(cluster.r = T, cluster.c = F,
               show_rownames = F,
               annotation = F,
               title = "")

# # zplot (plotting z-scores)
# tc5.all.DEGs %>% 
#   filter(sig=="yes") %>% 
#   head(20) %>% 
#   pull(gene_name) %>% 
#   zplot() %>% 
#   pluck(1) %>% 
#   multi.plot(rows = 5, cols = 4)

# exp.plot (plotting raw or log-transformed expression)
tc5.all.DEGs %>% 
  filter(sig=="yes") %>% 
  head(20) %>% 
  pull(gene_name) %>% 
  exp.plot(log = T) %>% 
  pluck(1) %>% 
  multi.plot(rows = 5, cols = 4)

sig.DEGs.tc5.all <- tc5.all.DEGs %>% filter(sig=="yes") %>% pull(gene_name)


## Compare DEGs with Rhythmic genes and expressed genes in Cflo ----------

## For+Nur
venndiagram(x= sig.DEGs.tc5.all,
            y= all.rhy.genes,
            z= bg.genes,
            unique=T, title="DEGs", 
            labels=c("DEG", "Rhythmic", "Expressed"), 
            lines=1:3, lcol=1:3, 
            diacol=1, plot=T, type="3")


# Manipulation - DEGs -----------------------------------------------

load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/cflo_ophio_degs.RData")

# set threshold for log2FC to identify DEG
### Will et al. 2020 used log2FC ≥ 1 (i.e. 2 fold change) and ≥ 4 fpkm
### Here, I am using the log2FC ≥ 1 and ≥ 1 fpkm (comparable geneset)
ophio.deg.logfc <- 1

sig.DEGs.ophio.alive.control <-
  cflo_Ian_cuffdiff %>% 
  # only the alive vs control 
  filter(sample_1 == "Alive" & sample_2 == "Control") %>% 
  # only the significant DEGs
  filter(significant == "yes") %>%
  # double check that at least one of the values was ≥ 1 fpkm
  filter(value_1 >= 1 | value_2 >= 1) %>%
  mutate(logFC = as.numeric(log2.fold_change.)) %>% 
  filter(abs(logFC) >= ophio.deg.logfc) %>%
  pull(gene)

sig.DEGs.ophio.alive.dead <-
  cflo_Ian_cuffdiff %>% 
  filter(sample_1 == "Alive" & sample_2 == "Dead") %>% 
  filter(significant == "yes") %>%
  filter(value_1 >= 1 | value_2 >= 1) %>%
  mutate(logFC = as.numeric(log2.fold_change.)) %>% 
  filter(abs(logFC) >= ophio.deg.logfc) %>%
  pull(gene)

sig.DEGs.ophio.all <- union(sig.DEGs.ophio.alive.control, sig.DEGs.ophio.alive.dead)

# Compare the DEGs (for/nur, alive/dead, alive/control)
venndiagram(x= sig.DEGs.tc5.all,
            y= sig.DEGs.ophio.alive.control,
            unique=T, title="DEGs", 
            labels=c("ForVsNur", "AliveVsControl"), 
            lines=1:2, lcol=1:2, 
            diacol=1, plot=T, type="2")


# Fig 9: Upset diagram - DEGs - DOL vs Manipulation -----------------------

## how many are up or down regulated? How does the direction compare?
# reformat the ophio data
ophio.data <- cflo_Ian_cuffdiff %>% 
  select(gene, sample_1, sample_2, value_1:annotation) %>% 
  # only the alive vs control 
  filter(sample_1 == "Alive" & sample_2 == "Control") %>% 
  # only the significant DEGs
  filter(significant == "yes") %>%
  # double check that at least one of the values was ≥ 1 fpkm
  filter(value_1 >= 1 | value_2 >= 1) %>%
  mutate(logFC = as.numeric(log2.fold_change.)) %>% 
  # DEGs ≥ 1 log2FC
  filter(abs(logFC) >= ophio.deg.logfc) %>%
  # make a column indicating if the gene is up- or down-regulated at the biting phase
  mutate(ophio_direction = ifelse(significant == "yes", 
                                  ifelse(logFC > 0, "down", "up"), 
                                  "NA"))

## make an upset diagram (think VennDiagram)
# example of list input (list of named vectors)
listInput <- list(for_up = tc5.all.DEGs %>% filter(for_direction == "up") %>% pull(gene_name),
                  nur_up = tc5.all.DEGs %>% filter(for_direction == "down") %>% pull(gene_name),
                  ophio_up = ophio.data %>% filter(ophio_direction == "up") %>% pull(gene),
                  ophio_down = ophio.data %>% filter(ophio_direction == "down") %>% pull(gene))

library(UpSetR)
library(viridis)
caste.col <- c("#F23030","#1A80D9")
upset(fromList(listInput), order.by = "freq")
upset(fromList(listInput), 
      number.angles = 0, point.size = 3, line.size = 1.5, 
      mainbar.y.label = "# Genes", sets.x.label = "DEGs", 
      text.scale = c(1.5, 1.3, 1.5, 1, 1.5, 1.5),
      keep.order = T,
      sets.bar.color = viridis(1),
      # adding queries
      query.legend = "bottom",
      queries = list( # all forager only sets
        list(query = intersects, 
             params = list("for_up", "ophio_up"),
             active = T, color=caste.col[1],
             query.name="upregulated in foragers")))


## All DEGs involved in (behavioral) DOL
sig.DEGs.tc5.all %>%
  cflo_go_enrichment(bg = bg.genes) %>%
  go_enrichment_plot(fdr = 5)

# 
# # let's plot the first 20
# intersect(sig.DEGs.tc5.all, sig.DEGs.ophio.alive.control) %>% 
#   exp.plot(log=T) %>% 
#   pluck(1) %>% 
#   multi.plot(rows = 5, cols = )

# Let's do a heatmap
intersect(sig.DEGs.tc5.all, sig.DEGs.ophio.alive.control) %>% 
  cflo_heatmap(cluster.r = T, cluster.c = F,
               show_rownames = F,
               annotation = F,
               title = "")

# Enrichment for the overlapped genes
intersect(sig.DEGs.tc5.all, sig.DEGs.ophio.alive.control) %>%
  cflo_go_enrichment(bg = bg.genes) %>%
  go_enrichment_plot(fdr = 5)
#
## Done. No enriched terms.


# Compare the DEGs (for/nur, alive/dead, alive/control)
venndiagram(x= sig.DEGs.tc5.all,
            y= sig.DEGs.ophio.alive.control,
            z= sig.DEGs.ophio.alive.dead,
            unique=T, title="DEGs", 
            labels=c("ForVsNur", "AliveVsControl", "AliveVsDead"), 
            lines=1:3, lcol=1:3, 
            diacol=1, plot=T, type="3")

## Note: There are 25 DEGs common between the three conditions
sig.DEGs.ophio.tc5 <- intersect(sig.DEGs.tc5.all, intersect(sig.DEGs.ophio.alive.control, sig.DEGs.ophio.alive.dead))

gene_to_annot %>% 
  filter(gene_name %in% sig.DEGs.ophio.tc5)

# Let's do a heatmap
sig.DEGs.ophio.tc5 %>% 
  cflo_heatmap(cluster.r = T, cluster.c = F,
               show_rownames = F,
               annotation = F,
               title = "")

### Let's look at any evidence of manipulation, either alive-control or alive-dead
venndiagram(x= sig.DEGs.tc5.all,
            y= sig.DEGs.ophio.all,
            unique=T, title="DEGs", 
            labels=c("ForVsNur", "Manipulation"), 
            lines=1:2, lcol=1:2, 
            diacol=1, plot=T, type="2")
## Note: There are 57 DEGs common between for-nur (behavioral DOL) and during manipulation

sig.DEGs.ophio.all.tc5.all <- 
  intersect(sig.DEGs.tc5.all, sig.DEGs.ophio.all)

## Heatmaps
sig.DEGs.ophio.all.tc5.all %>% 
  cflo_heatmap(cluster.r = T, cluster.c = F,
               show_rownames = F,
               annotation = F,
               title = "")
## Enriched GO terms?
sig.DEGs.ophio.all.tc5.all %>% 
  cflo_go_enrichment(bg = bg.genes) %>% 
  go_enrichment_plot(fdr=5)


gene_to_annot %>% 
  filter(gene_name %in% sig.DEGs.ophio.all.tc5.all) %>% view()

# cflo.blastp.rhy.summ %>% 
# filter(cflo_gene %in% sig.DEGs.ophio.tc5) %>%
#   filter(!duplicated(cflo_gene)) %>% 
#   view()

###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-
#-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-##-###-##-###-##-###-##
###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-


### How many of the DEGs are also rhythmic genes
## load the list of rhy genes
load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/rhy_genes/TC5_rhy_genes.RData")

# Compare the DEGs and Rhy genes
venndiagram(x= sig.DEGs.tc5.all,
            y= sig.DEGs.ophio.alive.control,
            z= all.rhy.genes,
            unique=T, title="DEGs", 
            labels=c("ForVsNur", "AliveVsControl", "Rhythmic"), 
            lines=1:3, lcol=1:3, 
            diacol=1, plot=T, type="3")
# Note:
## 11 rhythmic genes that have potential role in behavior manipulation (alive-control) as well as behavioral DOL.
## 16/81 for-nur DEGs are rhythmic in at least one caste
## 135/806 alive-control DEGs are rhythmic in at least one caste
intersect(sig.DEGs.tc5.all, intersect(sig.DEGs.ophio.alive.control, all.rhy.genes)) %>% 
  cflo_heatmap(cluster.r = T, cluster.c = F,
               show_rownames = T,
               annotation = T,
               title = "")


# What about genes with any signal of manipulation?
venndiagram(x= sig.DEGs.tc5.all,
            y= sig.DEGs.ophio.all,
            z= all.rhy.genes,
            unique=T, title="DEGs", 
            labels=c("DOL", "Manipulation", "Rhythmic"), 
            lines=1:3, lcol=1:3, 
            diacol=1, plot=T, type="3")
# Note:
## 12 rhythmic genes that have potential role in behavior manipulation as well as behavioral DOL.

intersect(sig.DEGs.tc5.all, intersect(sig.DEGs.ophio.all, all.rhy.genes)) %>% 
  cflo_heatmap(cluster.r = T, cluster.c = F,
               show_rownames = T,
               annotation = T,
               title = "")

# Keep only the genes that are rhythmic in both foragers and nurses
venndiagram(x= sig.DEGs.tc5.all,
            y= sig.DEGs.ophio.all,
            z= intersect(for.rhy.genes, nur.rhy.genes), # this is what changes
            unique=T, title="DEGs", 
            labels=c("DOL", "Manipulation", "Rhythmic_in_both"), 
            lines=1:3, lcol=1:3, 
            diacol=1, plot=T, type="3")
venndiagram(x= sig.DEGs.ophio.alive.control,
            y= sig.DEGs.ophio.alive.dead,
            z= intersect(for.rhy.genes, nur.rhy.genes), # this is what changes
            unique=T, title="DEGs", 
            labels=c("Control-Biting", "Biting-Dead", "Rhythmic_in_both"), 
            lines=1:3, lcol=1:3, 
            diacol=1, plot=T, type="3")
## Note: 
# 30 genes commonly rhythmic in for-nur are also involved in manipulation.
# 20 of those are involved in manipulated biting
# 2 are differentially expressed between control-biting as well as biting-dead

# what are these 2 genes?
intersect(intersect(sig.DEGs.ophio.alive.control,sig.DEGs.ophio.alive.dead),
          intersect(for.rhy.genes, nur.rhy.genes)) %>% 
  cflo_heatmap(cluster.r = T, cluster.c = F,
               show_rownames = T,
               annotation = T,
               title = "")
# what are these 20 genes?
intersect(sig.DEGs.ophio.alive.control,
          intersect(for.rhy.genes, nur.rhy.genes)) %>% 
  cflo_heatmap(cluster.r = T, cluster.c = F,
               show_rownames = T,
               annotation = T,
               title = "")


gene_to_annot %>% 
  filter(gene_name %in% intersect(sig.DEGs.ophio.all,intersect(for.rhy.genes, nur.rhy.genes)))
cflo.blastp.rhy.summ %>% 
  filter(cflo_gene %in% intersect(sig.DEGs.ophio.all,intersect(for.rhy.genes, nur.rhy.genes))) %>% 
  view()



# Fig 13: Upset - DEGs, rhythmic genes ------------------------------------
# example of list input (list of named vectors)
listInput <- list(
  # for_up = tc5.all.DEGs %>% filter(for_direction == "up") %>% pull(gene_name),
  # nur_up = tc5.all.DEGs %>% filter(for_direction == "down") %>% pull(gene_name),
  # ophio_up = ophio.data %>% filter(ophio_direction == "up") %>% pull(gene),
  # ophio_down = ophio.data %>% filter(ophio_direction == "down") %>% pull(gene),
  ophio_DEG = sig.DEGs.ophio.alive.control,
  for_rhy = setdiff(for.rhy.genes, intersect(for.rhy.genes, nur.rhy.genes)),
  nur_rhy = setdiff(nur.rhy.genes, intersect(for.rhy.genes, nur.rhy.genes)),
  both_rhy = intersect(for.rhy.genes, nur.rhy.genes)
)

library(UpSetR)
library(viridis)
caste.col <- c("#F23030","#1A80D9")
# upset(fromList(listInput), order.by = "freq")
upset(fromList(listInput), 
      number.angles = 0, point.size = 3, line.size = 1.5, 
      mainbar.y.label = "# genes", sets.x.label = "total", 
      text.scale = c(2, 1.5, 1.5, 1.5, 2, 2),
      # keep.order = F,
      sets.bar.color = viridis(1),
      # adding queries
      query.legend = "bottom",
      queries = list( # all forager only sets
        list(query = intersects, 
             params = list("ophio_DEG","for_rhy"),
             active = T, color=caste.col[1],
             query.name="")),
      keep.order=T,
      main.bar.color = viridis(100)[50])


# Random queries ----------------------------------------------------------

# all Cflo genes (n=13808)
cflo.annots.exp %>% select(gene_name, old_annotation) %>% view()

# all expressed Cflo genes (n = )
cflo.annots.exp %>% select(gene_name, old_annotation) %>% 
  filter(gene_name %in% bg.genes) %>% 
  view()


# all DOL-DEGs
tc5.all.DEGs %>% head()
# filter(sig == "yes") %>% 
select(gene_name) %>% 
  left_join((cflo.annots.exp %>%
               select(gene_name, old_annotation, enzyme_names)), 
            by = "gene_name") %>% 
  view()

## Get all the DEGs
tc5.all.DEGs %>%
  filter(sig == "yes") %>% 
  mutate(DE = ifelse(for_direction == "down", 
                     "nurse-biased", "forager-biased"),
         fold_change = 2^abs(logFC)) %>%
  select(DE, gene_name, fold_change, adj.P.Val) %>% 
  left_join((cflo.annots.exp %>%
               select(gene_name, old_annotation, signalP, TMHMM, GOs, pfams)), 
            by = "gene_name") %>% 
  arrange(DE, desc(fold_change)) %>% 
  view()

## Get the dmel orthologs of Cflo DOL-related DE genes
cflo.pogo.summ %>% view()
filter(cflo_gene %in% (tc5.all.DEGs %>% 
                         filter(sig == "yes") %>% 
                         pull(gene_name))) %>% view()
pull(cflo_gene) %>% 
  unique() %>% 
  length()

