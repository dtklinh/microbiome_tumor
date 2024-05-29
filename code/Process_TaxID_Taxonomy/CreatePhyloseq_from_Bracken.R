## get all taxID
library(phyloseq)
library(tidyverse)
library(microViz)

df_bracken <- read.table("../../raw_data/Bracken/LIB-CP1309a-03-A-1_CP05665.bracken", sep = "\t", header = T)

df_bracken <- df_bracken %>% 
  filter(taxonomy_id !=9606)
## process All_taxonomy file

Lines <- read.csv("../../script/Sample03_taxonomy.txt", sep = "\t", header = F)
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
##tmp <- NULL
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

df_taxa <- Taxa_tab %>% remove_rownames() %>% column_to_rownames(var = "TaxID")

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
  dplyr::select(superkingdom, phylum, class = class_co, order = order_co, family = family_co, genus = genus_co, species)
## put back to Phyloseq
#NCTs <- readRDS("./physeqs/NCTs_16S_Bac.rds")

sm_data <- data.frame(SampleID = c("FF_NT_1"),
                          SampleType = "Normal",
                          Storage="FF")
sm_data <- sm_data %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "SampleID")

otu_tab <- df_bracken %>% 
  dplyr::select(taxonomy_id, FF_NT_1 = new_est_reads) %>% 
  remove_rownames() %>% 
  column_to_rownames("taxonomy_id") 

pseq <- phyloseq(otu_table(otu_tab, taxa_are_rows = T),
                 tax_table(as.matrix(df_taxa_correct)),
                 sample_data(sm_data))



saveRDS(pseq, "FF_NT_1.rds")

my_cols <- c("#9d547c","#56ca63","#a357d6","cornflowerblue","#419d2a","sandybrown","red3","peachpuff","cyan","paleturquoise3","mistyrose","mediumpurple","mediumseagreen","mediumorchid","moccasin","orange4","olivedrab","midnightblue","papayawhip","palevioletred4","brown1","greenyellow","orchid","navy","darkred","navajowhite1","mistyrose1","grey85","#525fd6","red2","#8cbe3a","#c944aa","indianred3","#5ba557","#9e66cb","#c1b735","#6d82ec","grey25","#e69728","#6654b0","lightsalmon3","lightcyan1","khaki1","seagreen1","plum1","lightsteelblue1","palevioletred3","mintcream","magenta3","#799330","#da7fdf","#3c782c","#e44586","blue4","#63c996","#dc3f53","#49cbc8","#cf3f29","#4fabda","#da6c2b","#598bd1","#b78c24","#8d4191","#a0b971","slategray1","sienna","plum1","lightyellow1","lightskyblue3","linen","limegreen","cornsilk1","mediumaquamarine","gray14","gold3","darkviolet","#b2386a","#479d71","#ae4341","#2ba198","#e07557","#5361a3","#dda353","#aa98df","#5b6114","#dc89bf","#327243","slateblue1","#e57b94","#277257","#9b62a0","#bbab59","#98495a","#526229","#d8827d","#857624","gray40","#9a4a22","#7c7d46","mediumslateblue","lemonchiffon1","#e3a073","#9e6b33", "gray74","slateblue1","rosybrown3", "lawngreen","gainsboro","dodgerblue3","deeppink3","firebrick3", "orchid2", "olivedrab1", "ivory3", "darkseagreen", "bisque2", "darkgoldenrod2", "blue2", "skyblue", "seashell2", "turquoise", "tan1", "seagreen2", "palevioletred3", "linen", "steelblue4","ghostwhite","dodgerblue1","deeppink1","firebrick1", "limegreen", "purple3", "khaki3", "snow3", "darkslategray","darkorchid","lavender", "magenta2", "palegreen", "salmon", "maroon", "cyan2","#671408","#FAEBD7","#7FFFD4","#F0FFFF","#A52A2A","burlywood","cadetblue","#7FFF00","chocolate","cornsilk","slateblue1","#FF7F50","red1","#008B8B","darkgoldenrod1","darkolivegreen","darkorange4","white","hotpink","honeydew1","goldenrod2","darkgreen","oldlace","darkslategray3","navajowhite3","orchid4","gray25","#F0924D")


pseq %>% 
  comp_barplot(
    tax_level = "species", n_taxa = 15,
    bar_outline_colour = NA,
    sample_order = "bray",
    bar_width = 0.7,
    taxon_renamer = toupper
  ) + coord_flip()

