# setwd("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/reciprocal_blast/feature_tables")

rm(list = ls())

# Load libraries
pacman::p_load(conflicted, tidyverse)

# Set preference for functions
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("setdiff", "BiocGenerics")
conflict_prefer("union", "BiocGenerics")
conflict_prefer("intersect", "BiocGenerics")


# RefSeq to Gene ----------------------------------------------------------
## Make the refseq to gene name (symbol) file for each of the organisms
## 
# 1. Cflo --------------------------------------------------------------------

# # Read the feature_table
# cflo.ft <- read.csv("cflo.csv",
#                       stringsAsFactors = F,
#                       header = T,
#                       na.strings = c(""," ","NA"))
# 
# cflo.ft %>% glimpse()
# 
# cflo.ft.formatted <- 
#   cflo.ft %>% 
#   select(feature,
#          class,
#          product_accession,
#          related_accession,
#          cflo_gene = symbol,
#          geneID = GeneID,
#          cflo_name = name)
# 
# cflo.ft.formatted %>% glimpse()
# 
# cflo.ft.formatted$cflo_gene <- as.factor(cflo.ft.formatted$cflo_gene)
# 
# ## Number of unique genes?
#   cflo.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature == "gene" & class == "protein_coding") %>% 
#   # how many unique genes are there?
#   pull(cflo_gene) %>% 
#   unique() 
#   # 13,861 unique protein coding genes
# 
# ## number of unique proteins?
# cflo.proteins <-
#   cflo.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature != "gene") %>%
#   # how many unique genes are there?
#   mutate(protein = as.factor(product_accession)) %>% 
#   pull(protein) %>% 
#   unique() %>% 
#   as.character()
# 
# head(cflo.proteins)
# 
# # Let's make a RefSeq accession number to gene_name mapping
# cflo.refseq_to_gene <- 
#   cflo.ft.formatted %>% 
#   select(product_accession, 
#          cflo_gene, 
#          cflo_name) %>% 
#   filter(product_accession %in% cflo.proteins) %>% 
#   filter(!is.na(product_accession))
# 
# cflo.refseq_to_gene %>% head()


# 2. Sinv --------------------------------------------------------------------

# # Read the feature_table
# sinv.ft <- read.csv("sinv.csv",
#                       stringsAsFactors = F,
#                       header = T,
#                       na.strings = c(""," ","NA"))
# 
# sinv.ft %>% glimpse()
# 
# sinv.ft.formatted <- 
#   sinv.ft %>% 
#   select(feature,
#          class,
#          product_accession,
#          related_accession,
#          sinv_gene = symbol,
#          geneID = GeneID,
#          sinv_name = name)
# 
# sinv.ft.formatted %>% glimpse()
# sinv.ft.formatted$sinv_gene <- as.factor(sinv.ft.formatted$sinv_gene)
# 
# ## Number of unique genes?
# sinv.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature == "gene" & class == "protein_coding") %>% 
#   # how many unique genes are there?
#   pull(sinv_gene) %>% 
#   unique() %>% 
#   length()
# 
# ## number of unique proteins?
# sinv.proteins <-
#   sinv.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature != "gene") %>%
#   # how many unique genes are there?
#   mutate(protein = as.factor(product_accession)) %>% 
#   pull(protein) %>% 
#   unique() %>% 
#   as.character() 
# 
# head(sinv.proteins)
# 
# # Let's make a RefSeq accession number to gene_name mapping
# sinv.refseq_to_gene <- 
#   sinv.ft.formatted %>% 
#   select(product_accession, 
#          sinv_gene, 
#          sinv_name) %>% 
#   filter(product_accession %in% sinv.proteins) %>% 
#   filter(!is.na(product_accession))
# 
# sinv.refseq_to_gene %>% head()

# 3. amel --------------------------------------------------------------------

# # Read the feature_table
# amel.ft <- read.csv("amel.csv",
#                     stringsAsFactors = F,
#                     header = T,
#                     na.strings = c(""," ","NA"))
# 
# amel.ft %>% glimpse()
# 
# amel.ft.formatted <- 
#   amel.ft %>% 
#   select(feature,
#          class,
#          product_accession,
#          related_accession,
#          amel_gene = symbol,
#          geneID = GeneID,
#          amel_name = name)
# 
# amel.ft.formatted %>% glimpse()
# amel.ft.formatted$amel_gene <- as.factor(amel.ft.formatted$amel_gene)
# 
# ## Number of unique genes?
# amel.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature == "gene" & class == "protein_coding") %>% 
#   # how many unique genes are there?
#   pull(amel_gene) %>% 
#   unique() %>% 
#   length()
# 
# ## number of unique proteins?
# amel.proteins <-
#   amel.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature != "gene") %>%
#   # how many unique genes are there?
#   mutate(protein = as.factor(product_accession)) %>% 
#   pull(protein) %>% 
#   unique() %>% 
#   as.character() 
# 
# head(amel.proteins)
# 
# # Let's make a RefSeq accession number to gene_name mapping
# amel.refseq_to_gene <- 
#   amel.ft.formatted %>% 
#   select(product_accession, 
#          amel_gene, 
#          amel_name) %>% 
#   filter(product_accession %in% amel.proteins) %>% 
#   filter(!is.na(product_accession))
# 
# amel.refseq_to_gene %>% head()


# 4. dmel --------------------------------------------------------------------

# # Read the feature_table
# dmel.ft <- read.csv("dmel.csv",
#                     stringsAsFactors = F,
#                     header = T,
#                     na.strings = c(""," ","NA"))
# 
# dmel.ft %>% glimpse()
# 
# dmel.ft.formatted <- 
#   dmel.ft %>% 
#   select(feature,
#          class,
#          product_accession,
#          related_accession,
#          dmel_gene = symbol,
#          geneID = GeneID,
#          dmel_name = name)
# 
# dmel.ft.formatted %>% glimpse()
# dmel.ft.formatted$dmel_gene <- as.factor(dmel.ft.formatted$dmel_gene)
# 
# ## Number of unique genes?
# dmel.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature == "gene" & class == "protein_coding") %>% 
#   # how many unique genes are there?
#   pull(dmel_gene) %>% 
#   unique() %>% 
#   length()
# 
# ## number of unique proteins?
# dmel.proteins <-
#   dmel.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature != "gene") %>%
#   # how many unique genes are there?
#   mutate(protein = as.factor(product_accession)) %>% 
#   pull(protein) %>% 
#   unique() %>% 
#   as.character() 
# 
# head(dmel.proteins)
# 
# # Let's make a RefSeq accession number to gene_name mapping
# dmel.refseq_to_gene <- 
#   dmel.ft.formatted %>% 
#   select(product_accession, 
#          dmel_gene, 
#          dmel_name) %>% 
#   filter(product_accession %in% dmel.proteins) %>% 
#   filter(!is.na(product_accession))
# 
# dmel.refseq_to_gene %>% head()


# 5. mmus --------------------------------------------------------------------

# # Read the feature_table
# mmus.ft <- read.csv("mmus.csv",
#                     stringsAsFactors = F,
#                     header = T,
#                     na.strings = c(""," ","NA"))
# 
# mmus.ft %>% glimpse()
# 
# mmus.ft.formatted <- 
#   mmus.ft %>% 
#   select(feature,
#          class,
#          product_accession,
#          related_accession,
#          mmus_gene = symbol,
#          geneID = GeneID,
#          mmus_name = name)
# 
# mmus.ft.formatted %>% glimpse()
# mmus.ft.formatted$mmus_gene <- as.factor(mmus.ft.formatted$mmus_gene)
# 
# ## Number of unique genes?
# mmus.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature == "gene" & class == "protein_coding") %>% 
#   # how many unique genes are there?
#   pull(mmus_gene) %>% 
#   unique() %>% 
#   length()
# # n = 22,075
# 
# ## number of unique proteins?
# mmus.proteins <-
#   mmus.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature != "gene") %>%
#   # how many unique genes are there?
#   mutate(protein = as.factor(product_accession)) %>% 
#   pull(protein) %>% 
#   unique() %>% 
#   as.character() 
# 
# head(mmus.proteins)
# 
# # Let's make a RefSeq accession number to gene_name mapping
# mmus.refseq_to_gene <- 
#   mmus.ft.formatted %>% 
#   select(product_accession, 
#          mmus_gene, 
#          mmus_name) %>% 
#   filter(product_accession %in% mmus.proteins) %>% 
#   filter(!is.na(product_accession))
# 
# mmus.refseq_to_gene %>% head()


# 6. human --------------------------------------------------------------------

# # Read the feature_table
# human.ft <- read.csv("human.csv",
#                     stringsAsFactors = F,
#                     header = T,
#                     na.strings = c(""," ","NA"))
# 
# human.ft %>% glimpse()
# 
# human.ft.formatted <- 
#   human.ft %>% 
#   select(feature,
#          class,
#          product_accession,
#          related_accession,
#          human_gene = symbol,
#          geneID = GeneID,
#          human_name = name)
# 
# human.ft.formatted %>% glimpse()
# human.ft.formatted$human_gene <- as.factor(human.ft.formatted$human_gene)
# 
# ## Number of unique genes?
# human.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature == "gene" & class == "protein_coding") %>% 
#   # how many unique genes are there?
#   pull(human_gene) %>% 
#   unique() %>% 
#   length()
# # n = 19,826
# 
# ## number of unique proteins?
# human.proteins <-
#   human.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature != "gene") %>%
#   # how many unique genes are there?
#   mutate(protein = as.factor(product_accession)) %>% 
#   pull(protein) %>% 
#   unique() %>% 
#   as.character() 
# 
# head(human.proteins)
# 
# # Let's make a RefSeq accession number to gene_name mapping
# human.refseq_to_gene <- 
#   human.ft.formatted %>% 
#   select(product_accession, 
#          human_gene, 
#          human_name) %>% 
#   filter(product_accession %in% human.proteins) %>% 
#   filter(!is.na(product_accession))
# 
# human.refseq_to_gene %>% head()


# 7. agam --------------------------------------------------------------------

# # Read the feature_table
# agam.ft <- read.csv("agam.csv",
#                      stringsAsFactors = F,
#                      header = T,
#                      na.strings = c(""," ","NA"))
# 
# agam.ft %>% glimpse()
# 
# agam.ft.formatted <- 
#   agam.ft %>% 
#   select(feature,
#          class,
#          product_accession,
#          related_accession,
#          agam_gene = symbol,
#          geneID = GeneID,
#          agam_name = name)
# 
# agam.ft.formatted %>% glimpse()
# agam.ft.formatted$agam_gene <- as.factor(agam.ft.formatted$agam_gene)
# 
# ## Number of unique genes?
# agam.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature == "gene" & class == "protein_coding") %>% 
#   # how many unique genes are there?
#   pull(agam_gene) %>% 
#   unique() %>% 
#   length()
# # n = 1053 (this is a bit odd!)
# 
# ## number of unique proteins?
# agam.proteins <-
#   agam.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature != "gene") %>%
#   # how many unique genes are there?
#   mutate(protein = as.factor(product_accession)) %>% 
#   pull(protein) %>% 
#   unique() %>% 
#   as.character() 
# 
# head(agam.proteins)
# 
# # Let's make a RefSeq accession number to gene_name mapping
# agam.refseq_to_gene <- 
#   agam.ft.formatted %>% 
#   select(product_accession, 
#          agam_gene, 
#          agam_name) %>% 
#   filter(product_accession %in% agam.proteins) %>% 
#   filter(!is.na(product_accession))
# 
# agam.refseq_to_gene %>% head()


# 8. mono --------------------------------------------------------------------

# # Read the feature_table
# mono.ft <- read.csv("mono.csv",
#                     stringsAsFactors = F,
#                     header = T,
#                     na.strings = c(""," ","NA"))
# 
# mono.ft %>% glimpse()
# 
# mono.ft.formatted <- 
#   mono.ft %>% 
#   select(feature,
#          class,
#          product_accession,
#          related_accession,
#          mono_gene = symbol,
#          geneID = GeneID,
#          mono_name = name)
# 
# mono.ft.formatted %>% glimpse()
# mono.ft.formatted$mono_gene <- as.factor(mono.ft.formatted$mono_gene)
# 
# ## Number of unique genes?
# mono.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature == "gene" & class == "protein_coding") %>% 
#   # how many unique genes are there?
#   pull(mono_gene) %>% 
#   unique() %>% 
#   length()
# # n = 12,651 
# 
# ## number of unique proteins?
# mono.proteins <-
#   mono.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature != "gene") %>%
#   # how many unique genes are there?
#   mutate(protein = as.factor(product_accession)) %>% 
#   pull(protein) %>% 
#   unique() %>% 
#   as.character() 
# 
# head(mono.proteins)
# 
# # Let's make a RefSeq accession number to gene_name mapping
# mono.refseq_to_gene <- 
#   mono.ft.formatted %>% 
#   select(product_accession, 
#          mono_gene, 
#          mono_name) %>% 
#   filter(product_accession %in% mono.proteins) %>% 
#   filter(!is.na(product_accession))
# 
# mono.refseq_to_gene %>% head()


# 9. pogo --------------------------------------------------------------------

# # Read the feature_table
# pogo.ft <- read.csv("pogo.csv",
#                     stringsAsFactors = F,
#                     header = T,
#                     na.strings = c(""," ","NA"))
# 
# pogo.ft %>% glimpse()
# 
# pogo.ft.formatted <- 
#   pogo.ft %>% 
#   select(feature,
#          class,
#          product_accession,
#          related_accession,
#          pogo_gene = symbol,
#          geneID = GeneID,
#          pogo_name = name)
# 
# pogo.ft.formatted %>% glimpse()
# pogo.ft.formatted$pogo_gene <- as.factor(pogo.ft.formatted$pogo_gene)
# 
# ## Number of unique genes?
# pogo.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature == "gene" & class == "protein_coding") %>% 
#   # how many unique genes are there?
#   pull(pogo_gene) %>% 
#   unique() %>% 
#   length()
# # n = 11,348
# 
# ## number of unique proteins?
# pogo.proteins <-
#   pogo.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature != "gene") %>%
#   # how many unique genes are there?
#   mutate(protein = as.factor(product_accession)) %>% 
#   pull(protein) %>% 
#   unique() %>% 
#   as.character() 
# 
# head(pogo.proteins)
# 
# # Let's make a RefSeq accession number to gene_name mapping
# pogo.refseq_to_gene <- 
#   pogo.ft.formatted %>% 
#   select(product_accession, 
#          pogo_gene, 
#          pogo_name) %>% 
#   filter(product_accession %in% pogo.proteins) %>% 
#   filter(!is.na(product_accession))
# 
# pogo.refseq_to_gene %>% head()


# 10. temno --------------------------------------------------------------------

# ## Note: It does not have a feature table, instead I am using the protein table available on NCBI
# 
# # Read the protein_table
# temno.ft <- read.csv("./protein_tables/temno.csv",
#                     stringsAsFactors = F,
#                     header = T,
#                     na.strings = c(""," ","NA"))
# 
# temno.ft %>% glimpse()
# 
# temno.ft.formatted <- 
#   temno.ft %>% 
#   select(
#          # feature,
#          # class,
#          # product_accession,
#           product_accession = Protein.product,
#          # related_accession,
#          # temno_gene = symbol,
#           temno_gene = Locus.tag,
#          # geneID = GeneID,
#          # temno_name = name
#           temno_name = Protein.Name)
# 
# temno.ft.formatted %>% glimpse()
# temno.ft.formatted$temno_gene <- as.factor(temno.ft.formatted$temno_gene)
# 
# ## Number of unique genes?
# temno.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   # filter(feature == "gene" & class == "protein_coding") %>% 
#   # how many unique genes are there?
#   pull(temno_gene) %>% 
#   unique() %>% 
#   length()
# # n = 13,028 
# 
# ## number of unique proteins?
# temno.proteins <-
#   temno.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   # filter(feature != "gene") %>%
#   # how many unique genes are there?
#   mutate(protein = as.factor(product_accession)) %>% 
#   pull(protein) %>% 
#   unique() %>% 
#   as.character() 
# 
# head(temno.proteins)
# 
# # Let's make a RefSeq accession number to gene_name mapping
# temno.refseq_to_gene <- 
#   temno.ft.formatted %>% 
#   select(product_accession, 
#          temno_gene, 
#          temno_name) %>% 
#   filter(product_accession %in% temno.proteins) %>% 
#   filter(!is.na(product_accession))
# 
# temno.refseq_to_gene %>% head()


# 11. clonal --------------------------------------------------------------------

# # Read the feature_table
# clonal.ft <- read.csv("clonal.csv",
#                     stringsAsFactors = F,
#                     header = T,
#                     na.strings = c(""," ","NA"))
# 
# clonal.ft %>% glimpse()
# 
# clonal.ft.formatted <- 
#   clonal.ft %>% 
#   select(feature,
#          class,
#          product_accession,
#          related_accession,
#          clonal_gene = symbol,
#          geneID = GeneID,
#          clonal_name = name)
# 
# clonal.ft.formatted %>% glimpse()
# clonal.ft.formatted$clonal_gene <- as.factor(clonal.ft.formatted$clonal_gene)
# 
# ## Number of unique genes?
# clonal.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature == "gene" & class == "protein_coding") %>% 
#   # how many unique genes are there?
#   pull(clonal_gene) %>% 
#   unique() %>% 
#   length()
# # n = 11,927
# 
# ## number of unique proteins?
# clonal.proteins <-
#   clonal.ft.formatted %>% 
#   # let's filter to keep only protein coding genes
#   filter(feature != "gene") %>%
#   # how many unique genes are there?
#   mutate(protein = as.factor(product_accession)) %>% 
#   pull(protein) %>% 
#   unique() %>% 
#   as.character() 
# 
# head(clonal.proteins)
# 
# # Let's make a RefSeq accession number to gene_name mapping
# clonal.refseq_to_gene <- 
#   clonal.ft.formatted %>% 
#   select(product_accession, 
#          clonal_gene, 
#          clonal_name) %>% 
#   filter(product_accession %in% clonal.proteins) %>% 
#   filter(!is.na(product_accession))
# 
# clonal.refseq_to_gene %>% head()


###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-
#-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-##-###-##-###-##-###-##
###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-

### save the protein to gene data for each organism
# save(cflo.refseq_to_gene, 
#      sinv.refseq_to_gene, amel.refseq_to_gene, dmel.refseq_to_gene, mmus.refseq_to_gene, human.refseq_to_gene,
#      agam.refseq_to_gene, mono.refseq_to_gene, pogo.refseq_to_gene, temno.refseq_to_gene, clonal.refseq_to_gene,
#      file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/reciprocal_blast/TC5_protein_to_gene.RData")
##
## Done.

# To load the protein to gene data
# load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/reciprocal_blast/TC5_protein_to_gene.RData")

###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-
#-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-##
###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-





# Read results of proteinortho5 ---------------------------------------

# 1. cflo v sinv -------------------------------------------------------------

# cflo.sinv.ortho <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/reciprocal_blast/proteinortho_results/cflo_v_sinv_1.csv",
#                             header = T,
#                             stringsAsFactors = F,
#                             na.strings = c(""," ","NA"))
# 
# cflo.sinv.ortho %>% glimpse()
# 
# cflo.sinv.ortho <- 
#   cflo.sinv.ortho %>% 
#   left_join(cflo.refseq_to_gene, by=c("cflo" = "product_accession")) %>% 
#   left_join(sinv.refseq_to_gene, by=c("sinv" = "product_accession"))
# 
# cflo.sinv.ortho %>% head()
# 
# cflo.sinv.ortho %>% 
#   pull(cflo_gene) %>% 
#   unique() %>% 
#   length() # 9,907 unique Cflo genes mapped


# 2. cflo v amel -------------------------------------------------------------

# cflo.amel.ortho <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/reciprocal_blast/proteinortho_results/cflo_v_amel_1.csv",
#                             header = T,
#                             stringsAsFactors = F,
#                             na.strings = c(""," ","NA"))
# 
# cflo.amel.ortho %>% glimpse()
# 
# cflo.amel.ortho <- 
#   cflo.amel.ortho %>% 
#   left_join(cflo.refseq_to_gene, by=c("cflo" = "product_accession")) %>% 
#   left_join(amel.refseq_to_gene, by=c("amel" = "product_accession"))
# 
# cflo.amel.ortho %>% head()
# 
# cflo.amel.ortho %>% 
#   pull(cflo_gene) %>% 
#   unique() %>% 
#   length() # 8,702 unique Cflo genes mapped


# 3. cflo v dmel -------------------------------------------------------------

# cflo.dmel.ortho <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/reciprocal_blast/proteinortho_results/cflo_v_dmel_1.csv",
#                             header = T,
#                             stringsAsFactors = F,
#                             na.strings = c(""," ","NA"))
# 
# cflo.dmel.ortho %>% glimpse()
# 
# cflo.dmel.ortho <- 
#   cflo.dmel.ortho %>% 
#   left_join(cflo.refseq_to_gene, by=c("cflo" = "product_accession")) %>% 
#   left_join(dmel.refseq_to_gene, by=c("dmel" = "product_accession"))
# 
# cflo.dmel.ortho %>% head()
# 
# cflo.dmel.ortho %>%
#   pull(cflo_gene) %>% 
#   unique() %>% 
#   length() # 6,032 unique Cflo genes mapped


# 4. cflo v mmus -------------------------------------------------------------

# cflo.mmus.ortho <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/reciprocal_blast/proteinortho_results/cflo_v_mmus_1.csv",
#                             header = T,
#                             stringsAsFactors = F,
#                             na.strings = c(""," ","NA"))
# 
# cflo.mmus.ortho %>% glimpse()
# 
# cflo.mmus.ortho <-
#   cflo.mmus.ortho %>%
#   left_join(cflo.refseq_to_gene, by=c("cflo" = "product_accession")) %>%
#   left_join(mmus.refseq_to_gene, by=c("mmus" = "product_accession"))
# 
# cflo.mmus.ortho %>% head()
# cflo.mmus.ortho %>% 
#   pull(cflo_gene) %>% 
#   unique() %>% 
#   length() # 5,430 unique Cflo genes mapped



# 5. cflo v human ------------------------------------------------------------

# cflo.human.ortho <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/reciprocal_blast/proteinortho_results/cflo_v_human_1.csv",
#                             header = T,
#                             stringsAsFactors = F,
#                             na.strings = c(""," ","NA"))
# 
# cflo.human.ortho %>% glimpse()
# 
# cflo.human.ortho <-
#   cflo.human.ortho %>%
#   left_join(cflo.refseq_to_gene, by=c("cflo" = "product_accession")) %>%
#   left_join(human.refseq_to_gene, by=c("human" = "product_accession"))
# 
# cflo.human.ortho %>% head()
# 
# cflo.human.ortho %>% 
#   pull(cflo_gene) %>% 
#   unique() %>% 
#   length() # 5,479 unique Cflo genes mapped



# 6. cflo v agam ------------------------------------------------------------

# cflo.agam.ortho <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/reciprocal_blast/proteinortho_results/cflo_v_agam_1.csv",
#                             header = T,
#                             stringsAsFactors = F,
#                             na.strings = c(""," ","NA"))
# 
# cflo.agam.ortho %>% glimpse()
# 
# cflo.agam.ortho <-
#   cflo.agam.ortho %>%
#   left_join(cflo.refseq_to_gene, by=c("cflo" = "product_accession")) %>%
#   left_join(agam.refseq_to_gene, by=c("agam" = "product_accession"))
# 
# cflo.agam.ortho %>% head()
# 
# cflo.agam.ortho %>% 
#   pull(cflo_gene) %>% 
#   unique() %>% 
#   length() # 5,993 unique Cflo genes mapped



# 7. cflo v mono ------------------------------------------------------------

# cflo.mono.ortho <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/reciprocal_blast/proteinortho_results/cflo_v_mono_1.csv",
#                             header = T,
#                             stringsAsFactors = F,
#                             na.strings = c(""," ","NA"))
# 
# cflo.mono.ortho %>% glimpse()
# 
# cflo.mono.ortho <-
#   cflo.mono.ortho %>%
#   left_join(cflo.refseq_to_gene, by=c("cflo" = "product_accession")) %>%
#   left_join(mono.refseq_to_gene, by=c("mono" = "product_accession"))
# 
# cflo.mono.ortho %>% head()
# 
# cflo.mono.ortho %>% 
#   pull(cflo_gene) %>% 
#   unique() %>% 
#   length() # 9,712 unique Cflo genes mapped


# 8. cflo v pogo ------------------------------------------------------------
# 
# cflo.pogo.ortho <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/reciprocal_blast/proteinortho_results/cflo_v_pogo_1.csv",
#                             header = T,
#                             stringsAsFactors = F,
#                             na.strings = c(""," ","NA"))
# 
# cflo.pogo.ortho %>% glimpse()
# 
# cflo.pogo.ortho <-
#   cflo.pogo.ortho %>%
#   left_join(cflo.refseq_to_gene, by=c("cflo" = "product_accession")) %>%
#   left_join(pogo.refseq_to_gene, by=c("pogo" = "product_accession"))
# 
# cflo.pogo.ortho %>% head()
# 
# cflo.pogo.ortho %>% 
#   pull(cflo_gene) %>% 
#   unique() %>% 
#   length() # 9,520 unique Cflo genes mapped

# 9. cflo v temno ------------------------------------------------------------

# cflo.temno.ortho <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/reciprocal_blast/proteinortho_results/cflo_v_temno_1.csv",
#                             header = T,
#                             stringsAsFactors = F,
#                             na.strings = c(""," ","NA"))
# 
# cflo.temno.ortho %>% glimpse()
# 
# cflo.temno.ortho <-
#   cflo.temno.ortho %>%
#   left_join(cflo.refseq_to_gene, by=c("cflo" = "product_accession")) %>%
#   left_join(temno.refseq_to_gene, by=c("temno" = "product_accession"))
# 
# cflo.temno.ortho %>% head()
# 
# cflo.temno.ortho %>% 
#   pull(cflo_gene) %>% 
#   unique() %>% 
#   length() # 7,202 unique Cflo genes mapped

# 10. cflo v clonal ------------------------------------------------------------

# cflo.clonal.ortho <- read.csv("~/Documents/GitHub/R-scripts_zombie_ant_lab/TC5/reciprocal_blast/proteinortho_results/cflo_v_clonal_1.csv",
#                              header = T,
#                              stringsAsFactors = F,
#                              na.strings = c(""," ","NA"))
# 
# cflo.clonal.ortho %>% glimpse()
# 
# cflo.clonal.ortho <-
#   cflo.clonal.ortho %>%
#   left_join(cflo.refseq_to_gene, by=c("cflo" = "product_accession")) %>%
#   left_join(clonal.refseq_to_gene, by=c("clonal" = "product_accession"))
# 
# cflo.clonal.ortho %>% head()
# 
# cflo.clonal.ortho %>% 
#   pull(cflo_gene) %>% 
#   unique() %>% 
#   length() # 9,907 unique Cflo genes mapped

# How many unique Cflo genes mapped to each of the organisms?
cflo.sinv.ortho %>% distinct(cflo_gene, .keep_all = T) %>% nrow()  #unique Cflo genes = 9,907
cflo.amel.ortho %>% distinct(cflo_gene, .keep_all = T) %>% nrow()  #unique Cflo genes = 8,702
cflo.dmel.ortho %>% distinct(cflo_gene, .keep_all = T) %>% nrow()  #unique Cflo genes = 6,032
cflo.mmus.ortho %>% distinct(cflo_gene, .keep_all = T) %>% nrow()  #unique Cflo genes = 5,430
cflo.human.ortho %>% distinct(cflo_gene, .keep_all = T) %>% nrow() #unique Cflo genes = 5,479
cflo.agam.ortho %>% distinct(cflo_gene, .keep_all = T) %>% nrow() #unique Cflo genes = 5,993
cflo.mono.ortho %>% distinct(cflo_gene, .keep_all = T) %>% nrow() #unique Cflo genes = 9,712
cflo.pogo.ortho %>% distinct(cflo_gene, .keep_all = T) %>% nrow() #unique Cflo genes = 9,520
cflo.temno.ortho %>% distinct(cflo_gene, .keep_all = T) %>% nrow() #unique Cflo genes = 7,202
cflo.clonal.ortho %>% distinct(cflo_gene, .keep_all = T) %>% nrow() #unique Cflo genes = 9,481

###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-
#-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-##-###-##-###-##-###-##
###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-

### save the blastp results with protein-gene mapping
# save(cflo.sinv.ortho, cflo.amel.ortho, cflo.dmel.ortho, cflo.mmus.ortho, cflo.human.ortho,
#      cflo.agam.ortho, cflo.mono.ortho, cflo.pogo.ortho, cflo.temno.ortho, cflo.clonal.ortho,
#      file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/reciprocal_blast/TC5_cflo_blastp_orthologs.RData")
##
## Done.

###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-
#-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-##-###-##-###-##-###-##
###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-




### Gene to gene summary tables ----

# 1. cflo v sinv -------
# cflo.sinv.summ <-
#   cflo.sinv.ortho %>% 
#     # keep only the gene names and gene descriptions
#   select(cflo_gene,sinv_gene, cflo_name, sinv_name) %>% 
#     # keep only the unique gene-to-gene mappings (note: isoforms will be lost)
#   distinct(cflo_gene, sinv_gene, .keep_all = T)

# 2. cflo v amel -------
# cflo.amel.summ <-
#   cflo.amel.ortho %>% 
#   # keep only the gene names and gene descriptions
#   select(cflo_gene,amel_gene, cflo_name, amel_name) %>% 
#   # keep only the unique gene-to-gene mappings (note: isoforms will be lost)
#   distinct(cflo_gene, amel_gene, .keep_all = T)

# 3. cflo v dmel -------
# cflo.dmel.summ <-
#   cflo.dmel.ortho %>% 
#   # keep only the gene names and gene descriptions
#   select(cflo_gene,dmel_gene, cflo_name, dmel_name) %>% 
#   # keep only the unique gene-to-gene mappings (note: isoforms will be lost)
#   distinct(cflo_gene, dmel_gene, .keep_all = T)

# 4. cflo v mmus -------
# cflo.mmus.summ <-
#   cflo.mmus.ortho %>%
#   # keep only the gene names and gene descriptions
#   select(cflo_gene,mmus_gene, cflo_name, mmus_name) %>%
#   # keep only the unique gene-to-gene mappings (note: isoforms will be lost)
#   distinct(cflo_gene, mmus_gene, .keep_all = T)

# 5. cflo v human -------
# cflo.human.summ <-
#   cflo.human.ortho %>%
#   # keep only the gene names and gene descriptions
#   select(cflo_gene,human_gene, cflo_name, human_name) %>%
#   # keep only the unique gene-to-gene mappings (note: isoforms will be lost)
#   distinct(cflo_gene, human_gene, .keep_all = T)

# 6. cflo v agam -------
# cflo.agam.summ <-
#   cflo.agam.ortho %>%
#   # keep only the gene names and gene descriptions
#   select(cflo_gene,agam_gene, cflo_name, agam_name) %>%
#   # keep only the unique gene-to-gene mappings (note: isoforms will be lost)
#   distinct(cflo_gene, agam_gene, .keep_all = T)

# 7. cflo v mono -------
# cflo.mono.summ <-
#   cflo.mono.ortho %>%
#   # keep only the gene names and gene descriptions
#   select(cflo_gene,mono_gene, cflo_name, mono_name) %>%
#   # keep only the unique gene-to-gene mappings (note: isoforms will be lost)
#   distinct(cflo_gene, mono_gene, .keep_all = T)

# 8. cflo v pogo -------
# cflo.pogo.summ <-
#   cflo.pogo.ortho %>%
#   # keep only the gene names and gene descriptions
#   select(cflo_gene,pogo_gene, cflo_name, pogo_name) %>%
#   # keep only the unique gene-to-gene mappings (note: isoforms will be lost)
#   distinct(cflo_gene, pogo_gene, .keep_all = T)

# 9. cflo v temno -------
# cflo.temno.summ <-
#   cflo.temno.ortho %>%
#   # keep only the gene names and gene descriptions
#   select(cflo_gene,temno_gene, cflo_name, temno_name) %>%
#   # keep only the unique gene-to-gene mappings (note: isoforms will be lost)
#   distinct(cflo_gene, temno_gene, .keep_all = T)

# 10. cflo v clonal -------
# cflo.clonal.summ <-
#   cflo.clonal.ortho %>%
#   # keep only the gene names and gene descriptions
#   select(cflo_gene,clonal_gene, cflo_name, clonal_name) %>%
#   # keep only the unique gene-to-gene mappings (note: isoforms will be lost)
#   distinct(cflo_gene, clonal_gene, .keep_all = T)

###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-
#-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-##-###-##-###-##-###-##
###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-

### save the summary tables for gene-to-gene mapping
# save(cflo.sinv.summ, cflo.amel.summ, cflo.dmel.summ, cflo.mmus.summ, cflo.human.summ,
#      cflo.agam.summ, cflo.mono.summ, cflo.pogo.summ, cflo.temno.summ, cflo.clonal.summ,
#      file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/reciprocal_blast/TC5_cflo_blastp_summ.RData")
##
## Done.

# load(file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/reciprocal_blast/TC5_cflo_blastp_summ.RData")

###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-
#-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-##-###-##-###-##-###-##
###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-





# # Combine summary tables -------------
# # for Cflo, sinv, amel, dmel, mmus and human
# cflo.blastp.summ <- 
#   cflo.sinv.summ %>% 
#   # combine sinv summary table to amel summ. table
#   left_join(cflo.amel.summ %>% select(-cflo_name), by="cflo_gene") %>% 
#   left_join(cflo.dmel.summ %>% select(-cflo_name), by="cflo_gene") %>% 
#   left_join(cflo.mmus.summ %>% select(-cflo_name), by="cflo_gene") %>%
#   left_join(cflo.human.summ %>% select(-cflo_name), by="cflo_gene") %>%
#   # reorder the columns
#   select(cflo_gene, cflo_name, sinv_name, 
#          amel_name, 
#          dmel_name, 
#          mmus_name, human_name, 
#          everything()) 
# 
# # How many Cflo genes have an ortholog in at least one of the other organisms (sinv, amel, dmel, mmus or human)?
# cflo.blastp.summ %>% distinct(cflo_gene, .keep_all = T) %>% nrow() # n = 9,907 genes have an ortholog
# 
# ## Get the list of Cflo genes that have an ortholog in at least one of the organisms
# cflo.blastp.genes <- cflo.blastp.summ %>% pull(cflo_gene) %>% as.character() %>% unique()

## Save this file as a csv and .RData 
# write.csv(cflo.blastp.summ,
#           row.names = F,
#           file = "~/OneDrive - University of Central Florida/BD-TC5/TC5_Paper/results_csv/reciprocal_blast/cflo_blastp_combined_summary.csv")
## 
## Done.

### save the summary tables for gene-to-gene mapping
# save(cflo.blastp.summ,
#      file = "~/Documents/GitHub/R-scripts_zombie_ant_lab/Functions/data/reciprocal_blast/cflo_blastp_combined_summary.RData")
##
## Done.

### To load the gene-to-gene mapping for Cflo, sinv, amel, dmel, mmus, and human
load(file = "./data/cflo_blastp_combined_summary.RData")

###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-
#-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-##-###-##-###-##-###-##
###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-

