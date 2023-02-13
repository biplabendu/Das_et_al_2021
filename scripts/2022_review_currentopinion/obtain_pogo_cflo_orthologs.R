pacman::p_load(dbplyr, tidyverse)
blast.db <- dbConnect(RSQLite::SQLite(), "./data/blast_data.db")
src_dbi(blast.db)

pogo.cflo <- 
  tbl(blast.db, "pogo") %>% 
  collect()

clonal.cflo <- 
  tbl(blast.db, "clonal") %>% 
  collect()

mono.cflo <- 
  tbl(blast.db, "mono") %>% 
  collect()

sinv.cflo <- 
  tbl(blast.db, "sinv") %>% 
  collect()

temno.cflo <- 
  tbl(blast.db, "temno") %>% 
  collect

n_orthologs <- function(sp1, sp2) {
  foo <- 
    pogo.cflo %>% 
    select({{sp1}}, {{sp2}}) %>% 
    na.omit() %>% 
    group_by({{sp1}}) %>% 
    summarize(ortholog=n()) 
  
  
  foo %>%
    select(n_orthologs = ortholog) %>% 
    group_by(n_orthologs) %>% 
    summarize(n_genes=n()) %>% 
    print()
  
  bar <- 
    foo %>% 
    ggplot() +
    geom_density(aes(ortholog)) +
    theme_Publication(20) +
    scale_x_continuous(limits = c(0,5))
    # scale_y_continuous(limits = c(0,10))
  # labs(title=paste0("number of ", as.character(sp1), " orthologs in ", as.character(sp2)))
  return(bar)
}

# n_orthologs(cflo_gene, pogo_gene)
n_orthologs(pogo_gene, cflo_gene)

## save the file
# pogo.cflo %>% 
#   write.csv(file="./data/gordon_data/pogo_cflo_orthologs.csv", row.names = F)
# temno.cflo %>%
#   write.csv(file="./data/gordon_data/temno_cflo_orthologs.csv", row.names = F)
