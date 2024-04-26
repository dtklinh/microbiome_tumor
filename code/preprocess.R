## get all taxID
library(phyloseq)


FF <- readRDS("./raw_data/Bulk_FFPE/combined_bulk_physeq_object.rds")
FFPE <- readRDS("./raw_data/Bulk_FFPE/combined_ffpe_physeq_object.rds")
NCT <- readRDS("./raw_data/Bulk_FFPE/combined_ntc_physeq_object.rds")

Lst_taxID <-union( NCT %>% taxa_names(), FFPE %>% taxa_names())
Lst_taxID <- union(Lst_taxID, FF %>% taxa_names())
write.table(Lst_taxID, "./script/FFPE_FF_TaxID.txt", quote = F, row.names = F, col.names = F)

## process All_taxonomy file

Lines <- read.csv("./script/FFPEvsFF_taxonomy.txt", sep = "\t", header = F)
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

##df_taxa[NCT %>% taxa_names(),] %>% head()
tax_table(NCT) <-  tax_table(as.matrix(df_taxa[NCT %>% taxa_names(),]))
tax_table(FFPE) <-  tax_table(as.matrix(df_taxa[FFPE %>% taxa_names(),]))
tax_table(FF) <-  tax_table(as.matrix(df_taxa[FF %>% taxa_names(),]))

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
tax_table(NCT) <-  tax_table(as.matrix(df_taxa_correct[NCT %>% taxa_names(),]))
tax_table(FFPE) <-  tax_table(as.matrix(df_taxa_correct[FFPE %>% taxa_names(),]))
tax_table(FF) <-  tax_table(as.matrix(df_taxa_correct[FF %>% taxa_names(),]))

saveRDS(NCT, "./raw_data/Bulk_FFPE/NCT_TaxaAdj.rds")
saveRDS(FFPE, "./raw_data/Bulk_FFPE/FFPE_TaxaAdj.rds")
saveRDS(FF, "./raw_data/Bulk_FFPE/FF_TaxaAdj.rds")
