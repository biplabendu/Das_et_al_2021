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
# limorhyde functions
source(system.file('extdata', 'vignette_functions.R', package = 'limorhyde'))
# my functions
source("./functions/plot_zscores.R")
# source("~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/vennDia.R")
source("./functions/enrichment_analysis.R")
source("./functions/theme_publication.R")


## 1.1 LimoRhyde --------------------------------------

# Code borrowed from: http://127.0.0.1:20841/library/limorhyde/doc/introduction.html

# Let's load the libraries required for running Limorhyde
# library('annotate')
library('data.table')
# library('foreach')
# library('GEOquery')

# library('knitr')
library('limma')
library('limorhyde')

### 1.1.1 Format the meta-data ----------------
# load the meta-data
sm.tc5 <- read.csv("./data/TC5_meta_data.csv",
                   stringsAsFactors = F)

# Let's format the columns in the right data-type
sm.tc5$genotype <- as.factor(sm.tc5$genotype)
sm.tc5$cond <- as.factor(sm.tc5$cond)
sm.tc5$time <- as.numeric(sm.tc5$time)
sm.tc5$batch <- as.factor(sm.tc5$batch)
sm.tc5$LD <- as.factor(sm.tc5$LD)
sm.tc5$location <- as.factor(sm.tc5$location)

# Let's get a glimpse of the metadata
sm.tc5 %>%
  glimpse()

# Next we use limorhyde to calculate time_cos and time_sin, which are based on the first
#harmonic of a Fourier decomposition of the time column, and append them to the sm data frame.
sm.tc5 = cbind(sm.tc5, limorhyde(sm.tc5$time, 'time_'))
# convert the dataframe into a data.table
sm.tc5 <- data.table(sm.tc5)
# check that it worked
sm.tc5[1:5, ]

### 1.1.2 Format the expression data ----------------

c.for <- cflo.annots.exp %>%
  dplyr::select(gene_name,X2F:X24F) %>%
  na.omit() %>%
  # remove genes that have no expression for any sample (for and nur)
  filter_at(vars(starts_with("X")), any_vars(. > 0)) %>%
  # keep only the genes that have at least 1 fpkm for at least half the timepoints
  filter(gene_name %in% unique(c(gs1,gs2))) # gs1 - foragers; gs2 - nurses

c.nur <- cflo.annots.exp %>%
  dplyr::select(gene_name,X2N:X24N) %>%
  na.omit() %>%
  # remove genes that have no expression for any sample (for and nur)
  filter_at(vars(starts_with("X")), any_vars(. > 0)) %>%
  # keep only the genes that have at least 1 fpkm for at least half the timepoints
  filter(gene_name %in% unique(c(gs1,gs2)))


# rename the columns to match the names in the meta-data
#Foragers
oldnames.for <- names(c.for[,-1])
newnames.for <- c("2F","4F","6F","8F","10F","12F","14F","16F","18F","20F","22F","24F")
emat.for <- c.for %>%
  rename_at(vars(oldnames.for), ~ newnames.for)
# take a look if it worked
emat.for[1:5,]

#Nurses
oldnames.nur <- names(c.nur[,-1])
newnames.nur <- c("2N","4N","6N","8N","10N","12N","14N","16N","18N","20N","22N","24N")
emat.nur <- c.nur %>%
  rename_at(vars(oldnames.nur), ~ newnames.nur)
# take alook if it worked
emat.nur[1:5,]

# Merge the for and nur files into one
#Merge emat.for and emat.nur by gene_name
emat.tc5 <- emat.for %>%
  # keep only the genes that are expressed in both foragers and nurses
  # (Note: there are 4 genes that are expressed in foragers but not in nurses)
  filter(gene_name %in% emat.nur$gene_name) %>%
  left_join(emat.nur, by = "gene_name")
# Did it work?
head(emat.tc5)



# Create the log2 transformed input matrix for identifying DEGs ----------------------------
emat.tc5 <- as.data.frame(emat.tc5)
rownames(emat.tc5) <- emat.tc5[,1]
emat.tc5 <- emat.tc5[,-1]
summary(emat.tc5)
# Let's remove NAs
emat.tc5 <- na.omit(emat.tc5) # n(genes) = 9,361

# Need to make the emat.tc5 into a matrix.
emat.tc5 <- data.matrix(emat.tc5)
emat.tc5[1:5,1:5]
class(emat.tc5) # is a matrix with samples=columns and genes=rows

## log2 transform the data
emat.tc5 <- log2(emat.tc5 + 1)


# 1.2 Identify DEGs ----------------------------------------------
# Identify DEGs:
### Differential expression is based on the coefficient for cond in a linear model with no 
### interaction terms. We pass all genes to limma, but keep results only for non-differentially 
### rhythmic genes, and adjust for multiple testing accordingly.

# Set threshold for q-value (BH adjusted p-value)
q.threshold <- 0.05 # currently, using 5% FDR
log2.foldchange <- 1 # thus, any gene with a 2^(log2.foldchange) fold change in it's expression

# 1.2.1 all rhy genes -----------------------------------------------

# Use the subsetted emat.tc5.all.rhy to find DEGs
design.deg = model.matrix(~ cond + time_cos + time_sin, data = sm.tc5)

fit = lmFit(emat.tc5, design.deg)
fit = eBayes(fit, trend = TRUE)
# Take a look at the coefficients table
fit$coefficients %>% head()

deLimma.deg = data.table(topTable(fit, coef = 2, number = Inf), keep.rownames = TRUE)
setnames(deLimma.deg, 'rn', 'gene_name')
deLimma.deg[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
setorderv(deLimma.deg, 'adj.P.Val')

deLimma.deg %>%
  arrange(adj.P.Val) %>%
  head()


# Save the results of the DEG analysis 
# save(deLimma.deg, 
#      file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/DEGs/TC5_DEG_results.RData")
##
## Done.

# # Load the results of the DEG analysis (dataframe: deLimma.deg)
# load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/DEGs/TC5_DEG_results.RData")

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
  go_enrichment_plot(fdr = 5)

sig.DEGs.tc5.all <- tc5.all.DEGs %>% filter(sig=="yes") %>% pull(gene_name)



###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-
#-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-##-###-##-###-##-###-##
###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-


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
                     "nurse-DEG", "forager-DEG"),
         fold_change = 2^abs(logFC)) %>%
  select(DE, gene_name, fold_change, adj.P.Val) %>% 
  left_join((cflo.annots.exp %>%
               select(gene_name, old_annotation, signalP, TMHMM, GOs, pfams)), 
            by = "gene_name") %>% 
  arrange(DE, desc(fold_change)) %>% 
  view()


