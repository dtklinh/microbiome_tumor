## Maasalin2

library(Maaslin2)

input_data = system.file("extdata", "HMP2_taxonomy.tsv", package="Maaslin2") # The abundance table file input_data
#InData <- read.csv(input_data, sep="\t", row.names = 1)

input_metadata = system.file("extdata", "HMP2_metadata.tsv", package="Maaslin2") # The metadata table file input_metadata
#InMetadata <- read.csv(input_metadata, sep="\t")

#get the pathway (functional) data - place holder
download.file("https://raw.githubusercontent.com/biobakery/biobakery_workflows/master/examples/tutorial/stats_vis/input/pathabundance_relab.tsv", "./pathabundance_relab.tsv")

df_input_data = read.table(file             = input_data,
                           header           = TRUE,
                           sep              = "\t", 
                           row.names        = 1,
                           stringsAsFactors = FALSE)
df_input_data[1:5, 1:5]

df_input_metadata = read.table(file             = input_metadata, 
                               header           = TRUE, 
                               sep              = "\t", 
                               row.names        = 1,
                               stringsAsFactors = FALSE)
df_input_metadata[1:5, ]

df_input_path = read.csv("./pathabundance_relab.tsv", 
                         sep              = "\t", 
                         stringsAsFactors = FALSE, 
                         row.names        = 1)
df_input_path[1:5, 1:5]

####
fit_data = Maaslin2(input_data     = df_input_data, 
                    input_metadata = df_input_metadata, 
                    min_prevalence = 0,
                    normalization  = "NONE",
                    output         = "OutDir/demo_output_df", 
                    fixed_effects  = c("diagnosis", "dysbiosis"),
                    reference      = c("diagnosis,nonIBD"))

## functional
#This can also be done with with the HUMAnN 3 untiliy `humann_split_stratified_table`
unstrat_pathways <-function(dat_path){
  temp = dat_path[!grepl("\\|",rownames(dat_path)),]
  return(temp)
}

df_input_path_2 = unstrat_pathways(df_input_path)

fit_func = Maaslin2(input_data     = df_input_path_2, 
                    input_metadata = df_input_metadata, 
                    output         = "demo_functional", 
                    fixed_effects  = c("diagnosis", "dysbiosis"),
                    reference      = c("diagnosis,nonIBD"),
                    min_abundance  = 0.0001,
                    min_prevalence = 0.1)
