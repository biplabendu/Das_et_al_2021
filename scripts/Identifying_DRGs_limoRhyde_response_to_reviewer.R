# 1 Detecting DRGs -------------------------------------------------------------

# DRGs = Differentially Rhythmic Genes (differences in phase or amplitude)

rm(list = ls())

# ## Load packages ----------
# pacman::p_load(pheatmap, dendextend, tidyverse, viridis)
# pacman::p_load(RSQLite, tidyverse, dbplyr, DT, conflicted)
# 
# # set conflict preference (housekeeping to make sure functions work as expected)
# conflict_prefer("select", "dplyr")
# conflict_prefer("filter", "dplyr")
# conflict_prefer("layout", "plotly")
# 
# conflict_prefer("setdiff", "BiocGenerics")
# conflict_prefer("union", "BiocGenerics")
# conflict_prefer("intersect", "BiocGenerics")
# conflict_prefer("desc", "dplyr")
# 
# ## Load files --------------
# # - database (db)
# db <- dbConnect(RSQLite::SQLite(), "./data/TC5_data.db")
# src_dbi(db)
# # - db::cflo.annots
# cflo.annots.exp <- tbl(db, "annot_fpkm") %>% collect()
# # - db::expressed_genes (expressed for at least half time points)
# gs1 <- tbl(db, "expressed_genes") %>% filter(exp_half_for=="yes") %>% collect() %>% pull(gene_name)
# gs2 <- tbl(db, "expressed_genes") %>% filter(exp_half_nur=="yes") %>% collect() %>% pull(gene_name)
# 
# 
# # Let's load all functions:
# # # limorhyde functions
# # source(system.file('extdata', 'vignette_functions.R', package = 'limorhyde'))
# 
# # my functions
# source("./functions/plot_zscores.R")
# # source("~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/vennDia.R")
# source("./functions/enrichment_analysis.R")
# source("./functions/theme_publication.R")
# 
# 
# # ## 1.1 LimoRhyde --------------------------------------
# # 
# # # Code borrowed from: http://127.0.0.1:20841/library/limorhyde/doc/introduction.html
# # 
# # # Let's load the libraries required for running Limorhyde
# library('annotate')
# library('data.table')
# library('foreach')
# # # library('GEOquery')
# # 
# library('knitr')
# library('limma')
# library('limorhyde')
# # 
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
#   arrange(batch) %>% 
#   head()
# # 
# # Next we use limorhyde to calculate time_cos and time_sin, which are based on the first
# # harmonic of a Fourier decomposition of the time column, and append them to the sm data frame.
# sm.tc5 = cbind(sm.tc5, limorhyde(sm.tc5$time, 'time_'))
# # convert the dataframe into a data.table
# sm.tc5 <- data.table(sm.tc5)
# # check that it worked
# sm.tc5[1:5, ]
# # 
# ### 1.1.2 Format the expression data ----------------
# #
# c.for <- cflo.annots.exp %>%
#   dplyr::select(gene_name,X2F:X24F) %>%
#   na.omit() %>%
#   # remove genes that have no expression for any sample (for and nur)
#   filter_at(vars(starts_with("X")), any_vars(. > 0)) %>%
#   # keep only the genes that have at least 1 fpkm for at least half the timepoints
#   filter(gene_name %in% unique(c(gs1,gs2))) # gs1 - foragers; gs2 - nurses
# #
# c.nur <- cflo.annots.exp %>%
#   dplyr::select(gene_name,X2N:X24N) %>%
#   na.omit() %>%
#   # remove genes that have no expression for any sample (for and nur)
#   filter_at(vars(starts_with("X")), any_vars(. > 0)) %>%
#   # keep only the genes that have at least 1 fpkm for at least half the timepoints
#   filter(gene_name %in% unique(c(gs1,gs2)))
# # 
# # 
# # rename the columns to match the names in the meta-data
# #Foragers
# oldnames.for <- names(c.for[,-1])
# newnames.for <- c("2F","4F","6F","8F","10F","12F","14F","16F","18F","20F","22F","24F")
# emat.for <- c.for %>%
#   rename_at(vars(oldnames.for), ~ newnames.for)
# # take a look if it worked
# emat.for[1:5,]
# # 
# #Nurses
# oldnames.nur <- names(c.nur[,-1])
# newnames.nur <- c("2N","4N","6N","8N","10N","12N","14N","16N","18N","20N","22N","24N")
# emat.nur <- c.nur %>%
#   rename_at(vars(oldnames.nur), ~ newnames.nur)
# # take alook if it worked
# emat.nur[1:5,]
# # 
# # Merge the for and nur files into one
# # (Merge emat.for and emat.nur by gene_name)
# emat.tc5 <- emat.for %>%
#   # keep only the genes that are expressed in both foragers and nurses
#   # (Note: there are 4 genes that are expressed in foragers but not in nurses)
#   filter(gene_name %in% emat.nur$gene_name) %>%
#   left_join(emat.nur, by = "gene_name")
# # Did it work?
# head(emat.tc5)
# # 
# # 
# # 
# # # Create the log2 transformed input matrix for identifying DEGs ----------------------------
# emat.tc5 <- as.data.frame(emat.tc5)
# rownames(emat.tc5) <- emat.tc5[,1]
# emat.tc5 <- emat.tc5[,-1]
# summary(emat.tc5)
# # Let's remove NAs
# emat.tc5 <- na.omit(emat.tc5) # n(genes) = 9,361
# # 
# # Need to make the emat.tc5 into a matrix.
# emat.tc5 <- data.matrix(emat.tc5)
# emat.tc5[1:5,1:5]
# class(emat.tc5) # is a matrix with samples=columns and genes=rows
# # 
# ## log2 transform the data
# emat.tc5 <- log2(emat.tc5 + 1)
# # 
# #
# ## Load 24h rhythmic genesets 
# #
# ## Foragers
# for24 <-
#   tbl(db, "ejtk_all") %>% 
#   filter(caste == "for" & rhy == "yes") %>% 
#   select(gene_name, GammaP) %>% collect() %>% arrange(GammaP) %>%
#   select(gene_name) %>% pull()
# ## Nurses
# nur24 <-
#   tbl(db, "ejtk_all") %>% 
#   filter(caste == "nur" & rhy == "yes") %>% 
#   select(gene_name, GammaP) %>% collect() %>% arrange(GammaP) %>% 
#   select(gene_name) %>% pull()
# #
# # # 1.2 Identify DRGs ----------------------------------------------
# # # Identify DRGs:
# # ### Differential rhythmicity is based on statistical interactions between cond and time
# # ### components. We pass all genes in emat.tc5 to limma (whose Empirical Bayes does best 
# # ### with many genes), but keep results only for rhythmic genes, and adjust for multiple
# # ### testing accordingly.
# # 
# # Set threshold for q-value (BH adjusted p-value)
# q.threshold <- 0.05 # currently, using 5% FDR
# # log2.foldchange <- 1 # thus, any gene with a 2^(log2.foldchange) fold change in it's expression
# # 
# # # 1.2.1 all rhy genes -----------------------------------------------
# # 
# # Use the subsetted emat.tc5.all.rhy to find DEGs
# design.drg = model.matrix(~ cond * (time_cos + time_sin), data = sm.tc5)
# # 
# fit = lmFit(emat.tc5, design.drg)
# fit = eBayes(fit, trend = TRUE)
# # Take a look at the coefficients table
# fit$coefficients %>% head()
# # 
# deLimma.drg = data.table(topTable(fit, coef = 5:6, number = Inf), keep.rownames = TRUE)
# setnames(deLimma.drg, 'rn', 'gene_name')
# # Subset to keep only the overlapping 24h-rhythmic genes in foragers and nurses
# deLimma.drg = deLimma.drg %>% filter(gene_name %in% intersect(for24, nur24))
# deLimma.drg[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
# setorderv(deLimma.drg, 'adj.P.Val')
# #
# deLimma.drg %>%
#   head()
#
# Save the results of the DRG analysis
# save(deLimma.drg,
#      file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/DEGs/TC5_DRG_limorhyde_results.RData")
##
## Done.

# Load the results of the DEG analysis (dataframe: deLimma.deg)
load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/DEGs/TC5_DRG_limorhyde_results.RData")


# Save results of limorhyde_DRG analysis to a csv file --------------------

deLimma.drg %>% write.csv("./results/gene_lists/for_nur_DRG_limorhyde_output.csv", row.names = F)
