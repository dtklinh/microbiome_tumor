## get all taxID

Lst_taxID <- NCTs %>% taxa_names()
write.table(Lst_taxID, "./script/AllTaxID.txt", quote = F, row.names = F, col.names = F)

## process All_taxonomy file

Lines <- read.csv("./script/All_taxonomy.txt", sep = "\t", header = F)
##phylum;class;order;family;genus;species
Taxa_tab <- data.frame(TaxID = character(0),
                       superkingdom = character(0),
                       phylum = character(0),
                       class = character(0),
                       order = character(0),
                       family = character(0),
                       genus = character(0),
                       species = character(0))

## a  iterated each line
tmp <- NULL
for(i in 1:nrow(Lines)) {
  row <- Lines[i,]
  taxID <- row[[1]]
  content <- strsplit(row[[2]], ";")[[1]]
  colnam <- strsplit(row[[3]], ";")[[1]]
  sub_content <- lapply(c("superkingdom", "phylum", "class", "order", "family", "genus", "species"), function(x){
    idx <- which(colnam == x)
    if(length(idx)==0){return("unclassified")}
    return(content[idx])
  }) %>% unlist()
  Taxa_tab[nrow(Taxa_tab)+1, ] <- c(taxID, sub_content)
}

## put back to Phyloseq
NCTs <- readRDS("./physeqs/NCTs_16S_Bac.rds")
library(tidyverse)
df_taxa <- Taxa_tab %>% remove_rownames() %>% column_to_rownames(var = "TaxID")
tax_table(NCTs) <-  tax_table(as.matrix(df_taxa))
saveRDS(NCTs, "./physeqs/NCTs_CorrectedTaxa.rds")

## rename taxa rank
### class
df_taxa_correct <- df_taxa %>% 
  mutate(class_co = ifelse(class=="unclassified", paste0(phylum, "_", "unclassified"), class)) %>% 
  mutate(order_co = ifelse(order=="unclassified", 
                           ifelse(grepl("unclassified", class_co, fixed = TRUE), class_co, paste0(class_co, "_", "unclassified")), order)) %>% 
  mutate(family_co = ifelse(family=="unclassified", 
                           ifelse(grepl("unclassified", order_co, fixed = TRUE), order_co, paste0(order_co, "_", "unclassified")), family)) %>% 
  mutate(genus_co = ifelse(genus=="unclassified", 
                            ifelse(grepl("unclassified", family_co, fixed = TRUE), family_co, paste0(family_co, "_", "unclassified")), genus)) %>% 
  select(superkingdom, phylum, class = class_co, order = order_co, family = family_co, genus = genus_co, species)
  
## put back
tax_table(NCTs) <-  tax_table(as.matrix(df_taxa_correct))
saveRDS(NCTs, "./physeqs/NCTs_CorrectedTaxa_AdjustName.rds")
