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


# 01_Check point ----------------------------------------------------------


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

emat.tc5 %>% class()
emat.tc5 %>% head()

# 02_Check point ----------------------------------------------------------

# function to extract "p" consecutive numbers from a vector "x"
sample_consec_p <- function(x, p){
  firstIndex = sample(1:(length(x)-p+1), 1)
  x[firstIndex:(firstIndex + p -1)]
}

### Simualtion Begins ------

# summary.results <- list()
sim.results.wo.FC <- data.frame()
sim.results.with.FC <- data.frame()

time.points <- 2:12
iterations <- 100

# Set threshold for q-value (BH adjusted p-value)
q.threshold <- 0.05 # currently, using 5% FDR
log2.foldchange <- 1 # thus, any gene with a 2^(log2.foldchange) fold change in it's expression


for (i in 1:length(time.points)){
  
  
  for (j in 1:iterations) {
  
    # print(c(j,i))
    
  # list of all sampling time points
  all.time.points <- 2*(1:12)
  
  # draw random time-points (n = time.points)
  # zeitgeber.times <- sort(sample(all.time.points, time.points[i]))
  zeitgeber.times <- sample_consec_p(all.time.points, time.points[i])
  zeitgeber.times.for <- paste0(zeitgeber.times, "F")
  zeitgeber.times.nur <- paste0(zeitgeber.times, "N")
  
  # extract forager and nurse data for the random zeitgeber times
  emat.data <- emat.tc5 %>% select(c(zeitgeber.times.for, zeitgeber.times.nur))
  
  # adjust metadata accordingly
  sm.data <- sm.tc5 %>% filter(sample %in% c(zeitgeber.times.for, zeitgeber.times.nur) )
  
  ### Perform DGE analysis ###
  
  
  # Need to make the emat.tc5 into a matrix.
  emat.data <- data.matrix(emat.data)
  
  # log2 transform the data
  emat.data <- log2(emat.data + 1)
  
  
  # create the design matrix
  # design.deg = model.matrix(~ cond, data = sm.data)
  design.deg = model.matrix(~ cond + time_cos + time_sin, data = sm.data)
  
  # run the DGE analysis and save the data
  fit = lmFit(emat.data, design.deg)
  fit = eBayes(fit, trend = TRUE)
  # # Take a look at the coefficients table
  # fit$coefficients %>% head()
  
  deLimma.deg = data.table(topTable(fit, coef = 2, number = Inf), keep.rownames = TRUE)
  setnames(deLimma.deg, 'rn', 'gene_name')
  deLimma.deg[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
  setorderv(deLimma.deg, 'adj.P.Val')
  
  # deLimma.deg %>%
  #   arrange(adj.P.Val) %>%
  #   head()
  
  # Filter the results to keep only the significant genes
  all.DEGs <- 
    deLimma.deg %>% 
    arrange(adj.P.Val) %>% 
    mutate(sig = as.factor(ifelse(adj.P.Val < q.threshold & abs(logFC) >= log2.foldchange, "yes", "no"))) %>% 
    mutate(for_direction = as.factor(ifelse(sig=="yes", ifelse( logFC > 0, "down", "up" ), "NA")))
  
  # head(all.DEGs)
  
  ## How many DEGs at 5% FDR and no cut-off for fold-change
  sig.DEGs.wo.FC <- 
    all.DEGs %>% 
    filter(adj.P.Val < q.threshold)
    
  ## How many DEGs at 5% FDR and ≥ 2 fold-change in gene expression
  sig.DEGs.with.FC <- 
    all.DEGs %>% 
    filter(adj.P.Val < q.threshold) %>%
    filter(abs(logFC) >= log2.foldchange)
  
  # all.DEGs.summ
  
  # # Save the summary file
  # summary.results[[i]] <- all.DEGs.summ$total
  
  # ## What are these DEGs (gene names)?
  # sig.DEGs <-
  #   all.DEGs %>%
  #   filter(sig == "yes") %>%
  #   pull(gene_name)
  
  sim.results.wo.FC[j,i] <- nrow(sig.DEGs.wo.FC)
  sim.results.with.FC[j,i] <- nrow(sig.DEGs.with.FC)
  
  }
}

# plot(sim.results[,1])

# 
names(sim.results.with.FC) <- paste0("timepoints_",time.points)
names(sim.results.wo.FC) <- paste0("timepoints_",time.points)
# rownames(sim.results) <- paste0("random_sampling_", 1:iterations)
# 
# sim.results %>% view()

## Format the results for plotting
sim.results.with.FC %>% 
  pivot_longer(cols = starts_with("timepoints_"),
               names_to = "time_points", 
               values_to = "numDEGs") %>%
  mutate(time_points = factor(time_points, levels=c(names(sim.results.with.FC)))) %>% 

  ## Plot the results as boxplots      
  ggplot() + 
    geom_boxplot(aes(x=time_points, y=numDEGs)) +
    coord_flip() +
    theme_Publication() +
    ggtitle(paste0("#DEGs with abs(fold-change) ≥ ", 2^log2.foldchange)) +
    xlab("number of consecutive timepoints sampled")


## Format the results for plotting
sim.results.wo.FC %>% 
  pivot_longer(cols = starts_with("timepoints_"),
               names_to = "time_points", 
               values_to = "numDEGs") %>%
  mutate(time_points = factor(time_points, levels=c(names(sim.results.wo.FC)))) %>% 
  
  ## Plot the results as boxplots      
  ggplot() + 
  geom_boxplot(aes(x=time_points, y=numDEGs)) +
  coord_flip() +
  theme_Publication() +
  ggtitle("#DEGs with NO cut-off for fold-change ") +
  xlab("number of consecutive timepoints sampled") 
  
