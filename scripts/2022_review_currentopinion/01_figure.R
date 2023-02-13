cflo.sinv.amel <- 
  tbl(db, "blastp_rhy_summary") %>% 
  select(cflo_gene, rhy_for, rhy_nur, cflo_name, sinv_gene, amel_gene) %>% 
  collect() %>% 
  distinct()

n_orthologs <- function(sp1, sp2) {
  foo <- 
    cflo.sinv.amel %>% 
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
    scale_x_continuous(limits = c(0,5)) +
    scale_y_continuous(limits = c(0,10))
    # labs(title=paste0("number of ", as.character(sp1), " orthologs in ", as.character(sp2)))
  return(bar)
}

n_orthologs(cflo_gene, sinv_gene)
n_orthologs(sinv_gene, cflo_gene)

n_orthologs(cflo_gene, amel_gene)
n_orthologs(amel_gene, cflo_gene)


  
