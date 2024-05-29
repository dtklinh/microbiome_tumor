## analysis 
library(mbtools)
bracken_rel <- read.table("./raw_data/Bracken/LIB-CP1309a-03-A-1_CP05665.bracken", sep = "\t", header = T)
br <- bracken_to_phyloseq("./raw_data/Bracken/LIB-CP1309a-03-A-1_CP05665.bracken")

TaxIDs <- bracken_rel %>% 
  dplyr::select(taxonomy_id)

write.table(TaxIDs, "./script/Sample03_TaxIDs.txt", quote = F, col.names = F, row.names = F, sep = "\t")
