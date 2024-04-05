### is.valid.batch
is.valid.batch <- function(batch_num, batch_type, meta_data ){
  tmp <- NULL
  if(grepl("dna", tolower(batch_type), fixed = TRUE)){
    tmp <- meta_data[meta_data$DNAex_round==batch_num, "TRUE_control"] %>% table() %>% length()
  }else if(grepl("pcr", tolower(batch_type), fixed = TRUE)){
    tmp <- meta_data[meta_data$PCR_round==batch_num, "TRUE_control"] %>% table() %>% length()
  }else{
    stop("Batch type must contain either dna or pcr.")
  }
  return(tmp==2)
}

Extract_SampleIdx <- function(batch_num, batch_type, meta_data){
  sub_meta_data <- NULL
  if(grepl("dna", tolower(batch_type), fixed = TRUE)){
    sub_meta_data <- meta_data[meta_data$DNAex_round == batch_num,]
  }else if(grepl("pcr", tolower(batch_type), fixed = TRUE)){
    sub_meta_data <- meta_data[meta_data$PCR_round == batch_num,]
  }else{
    stop("Batch type must contain either dna or pcr.")
  }
  Sample.TRUE.idx <- sub_meta_data[sub_meta_data$TRUE_control=="TRUE","uniqueID"]
  Sample.NCT.idx <- sub_meta_data[sub_meta_data$TRUE_control=="control","uniqueID"]
  return(list("TrueSam"=Sample.TRUE.idx, "NCTSam"=Sample.NCT.idx))
}


is.valid.nct.type <- function(batch_num, batch_type, meta_data){
  sub_meta_data <- NULL
  if(grepl("dna", tolower(batch_type), fixed = TRUE)){
    sub_meta_data <- meta_data[meta_data$DNAex_round == batch_num,]
  }else if(grepl("pcr", tolower(batch_type), fixed = TRUE)){
    sub_meta_data <- meta_data[meta_data$PCR_round == batch_num,]
  }else{
    stop("Batch type must contain either dna or pcr.")
  }
  n.nct.type <- c()
  #for buffer
  n.buffer <- length(which(sub_meta_data$sample_type == "buffer"))
  n.buffer <- !(n.buffer == 0)
  n.nct.type["buffer"] <- n.buffer
  #for paraffin
  n.paraffin <- length(which(sub_meta_data$sample_type == "paraffin"))
  n.paraffin <- !(n.paraffin == 0)
  n.nct.type["paraffin"] <- n.paraffin
  #for pcr_control
  n.pcr <- length(which(sub_meta_data$sample_type == "pcr_ctrl"))
  n.pcr <- !(n.pcr == 0)
  n.nct.type["pcr_ctrl"] <- n.pcr
  return(n.nct.type)
}

create.sub.NCT.Sample.type <- function(batch_num,batch_type,meta_data,typeVector,NCT.sample,i){
  sub_meta_data <- NULL
  if(grepl("dna", tolower(batch_type), fixed = TRUE)){
    sub_meta_data <- meta_data[meta_data$DNAex_round == batch_num,]
  }else if(grepl("pcr", tolower(batch_type), fixed = TRUE)){
    sub_meta_data <- meta_data[meta_data$PCR_round == batch_num,]
  }else{
    stop("Batch type must contain either dna or pcr.")
  }
  sub.NCT.Sample.type.idx <- sub_meta_data[sub_meta_data$sample_type == names(typeVector)[i],"uniqueID"]
  sub.NCT.Sample.type <- prune_samples(sample_data(NCT.sample) %>% pull(uniqueID) %in% sub.NCT.Sample.type.idx, NCT.sample)
  sub.NCT.Sample.type <- prune_taxa(taxa_sums(sub.NCT.Sample.type) >0, sub.NCT.Sample.type)
  return(sub.NCT.Sample.type)
} 

  
testTaxa <- function(True.Sample,NCT.sample){
  prev.NCT <- tibble(
    my_p = prevalence(NCT.sample, detection  = 0, sort = TRUE, count = FALSE),
    OTU =names(prevalence(NCT.sample , detection = 0, sort = TRUE, count = FALSE))
  )
  
  prev.TRUE <- tibble(
    my_x = prevalence(True.Sample, detection  = 0, sort = TRUE, count = TRUE),
    OTU =names(prevalence(True.Sample , detection = 0, sort = TRUE, count = TRUE))
  )
    
  Tab <- merge(x=prev.NCT, y=prev.TRUE, by.x = "OTU", by.y = "OTU",all.x = TRUE, all.y = FALSE) %>% tibble()
  Tab[is.na(Tab)] <- 0
  Tab$my_x <- as.integer(Tab$my_x)
  Tab$pval <- 0
  n <- nsamples(True.Sample)
  for (j in seq(1,nrow(Tab))) {
    Tab[j,"pval"] <- binom.test(x = unlist(Tab[j,"my_x"]),n = n,p = unlist(Tab[j,"my_p"]),alternative = "greater")$p.value
  }
  return(Tab)
}


Binom.Test.complex <- function(True.Sample, NCT.sample, batch_num, batch_type, meta_data){
  cont <- list()
  A <- is.valid.nct.type(batch_num, batch_type, meta_data)
  for (i in seq(1:length(A))) {
    if (A[i] == FALSE) {
      cont[names(A)[i]] <- NA
    }
    if (A[i] == TRUE) {
      sub.NCT.Sample.type <- create.sub.NCT.Sample.type(batch_num,batch_type,meta_data,typeVector = A,NCT.sample,i)
      df <- testTaxa(NCT.sample = sub.NCT.Sample.type,True.Sample)
      cont[names(A)[i]] <- list(as.character(unlist(df[df$pval > 0.05, "OTU"]$OTU)))
    }
  }
  cont <- unique(unlist(cont))
  cont <- cont[!is.na(cont)]
  return(list("Contaminants" = cont))
}


BinomTest_Wrapper <- function(batch_num, batch_type, meta_data, True.Sample, NCT.sample){
  #check if batch for given batch_num and batch_type is valid
  if(!is.valid.batch(batch_num, batch_type, meta_data)){
    return(NULL)}
  #create subset of samples that are element of batch
  Lst_idx <- Extract_SampleIdx(batch_num, batch_type, meta_data)
  tmp_idx_true <- Lst_idx$TrueSam
  sub.True.Sample <- prune_samples(sample_data(True.Sample) %>% pull(uniqueID) %in% tmp_idx_true ,True.Sample)
  sub.True.Sample <- prune_taxa(taxa_sums(sub.True.Sample)>0, sub.True.Sample)
  
  return(Binom.Test.complex(sub.True.Sample, NCT.sample, batch_num, batch_type, meta_data))
  }


Wrapper_Nejman <- function(True.Sample, NCT.Sample, metadata, x){
    #get batch numbers for DNAex
    tmp <- metadata$DNAex_round %>% as.numeric()
    DNAex_lst <- tmp[!is.na(tmp)] %>% unique()
    #get batch numbers for PCR
    tmp <- metadata$PCR_round %>% as.numeric()
    PCR_lst <- tmp[!is.na(tmp)] %>% unique()
    
    Lst_contaminant = c() # list of contaminant taxa, empty at the beginning
    tmp_lst_DNA <- lapply(DNAex_lst, BinomTest_Wrapper, batch_type="DNA", meta_data=metadata, 
                      True.Sample=True.Sample, NCT.sample=NCT.Sample)
    
    tmp_lst_DNA <- unique(unlist(tmp_lst_DNA))
    ## Second, we fix the condition PCR
    tmp_lst_PCR <- lapply(PCR_lst, BinomTest_Wrapper, batch_type="PCR", meta_data=metadata, 
                      True.Sample=True.Sample, NCT.sample=NCT.Sample)
    
    tmp_lst_PCR <- unique(unlist(tmp_lst_PCR))
    
    Lst_contaminant <- union(tmp_lst_DNA,tmp_lst_PCR)
    
    keep_taxa <- setdiff(taxa_names(x), Lst_contaminant)
    return(prune_taxa(keep_taxa, x))
  }

## filter by low abundance
filter_by_low_abundance <- function(phyloseqObj, threshold, A){
  physeq.abun.filt <- phyloseq::transform_sample_counts(phyloseqObj, function(x){x/sum(x)})
  noCont <- phyloseq::genefilter_sample(physeq.abun.filt, filterfun_sample(function(x) x > threshold), A= A)
  physeq.abun.filt <- phyloseq::prune_taxa(taxa = noCont,x = phyloseqObj)
  return(physeq.abun.filt)
}

filter_by_low_readcount <- function(phyloseqObj, threshold, A){
  noCont <- phyloseq::genefilter_sample(phyloseqObj, filterfun_sample(function(x) x > threshold), A= A)
  physeq.read.filt <- phyloseq::prune_taxa(taxa = noCont,x = phyloseqObj)
  return(physeq.read.filt)
}

## filtering sample from low read count
visualize_samples_low_reads <- function(PhyObj, type = "NTC"){
  # use to remove sample from NTC whose number of reads are less than a  certain threshold
  bac.count.ntc <- tibble(bac.count = sample_sums(PhyObj), 
                          sample = sample_data(PhyObj)$id) %>% 
    arrange (bac.count)
  tibble(
    ">150 reads" = length(which(bac.count.ntc$bac.count > 150)),
    ">250 reads" = length(which(bac.count.ntc$bac.count > 250)),
    ">500 reads" = length(which(bac.count.ntc$bac.count > 500)),
    ">750 reads" = length(which(bac.count.ntc$bac.count > 750)),
    ">1000 reads" = length(which(bac.count.ntc$bac.count > 1000)),
  )  %>% 
    kable (caption=paste0(type, "read counts"), booktabs =TRUE) %>%
    kable_classic(full_width = F, html_font = "Cambria") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
}

filter_samples_low_reads <- function(PhyObj, threshold){
  bac.count.ntc <- dplyr::tibble(bac.count = phyloseq::sample_sums(PhyObj), 
                          sample = phyloseq::sample_data(PhyObj)$id)
  id.ex <- bac.count.ntc %>% dplyr::filter(bac.count < threshold) %>% dplyr::pull (sample)
  physeq.read.filt <- phyloseq::subset_samples(PhyObj, !(id %in% id.ex))
  physeq.read.filt <- phyloseq::prune_taxa(phyloseq::taxa_sums(physeq.read.filt)>0, physeq.read.filt)
  return(physeq.read.filt)
}


## wrapper of low abundance filtering by PERfect package
Wrapper_FERfect <- function(PhyObj){
  otu = as(otu_table(PhyObj), "matrix")
  if(taxa_are_rows(PhyObj)){otu <- t(otu)}
  otu = as_tibble(otu)
  res_sim <- PERFect_sim(X = otu)
  # dim(res_sim$filtX) 
  ids.sim<- colnames(res_sim$filtX) 
  return(prune_taxa(ids.sim, PhyObj))
}

## determine a threshold for high prevalence filtering.
## return an array, how many percent of taxa in true sample I remove if I choose that threshold 
HighPrevalence_Data <- function(True.Sample, NCT.Sample){
  num_taxa <- ntaxa(True.Sample)
  # subset of taxa which in in overlap
  idx_overlap <- intersect(taxa_names(True.Sample) %>% unique(),
                           taxa_names(NCT.Sample) %>% unique())
  NCT.Sample.Subset <- prune_taxa(idx_overlap, NCT.Sample)
  prev.nct <- tibble(
    prev = prevalence(NCT.Sample.Subset, detection  = 0, sort = TRUE, count = FALSE),
    OTU =as.factor(names(prevalence(NCT.Sample.Subset , detection = 0, sort = TRUE, count = FALSE)))
  )
  x <- prev.nct$prev
  factorx <- factor(cut(x, breaks=nclass.Sturges(x)))
  xout <- as.data.frame(table(factorx)) %>% map_df(rev)
  xout <- mutate(xout, cumFreq = cumsum(Freq), relative = prop.table(Freq))
  xout <- xout %>% mutate(rel_tax=cumFreq/num_taxa)
  return(xout[,c(1,5)])
}

ConQur_apply <- function(True.Sample, metadata, batchName, batch_covar, batch_ref){
  otu = otu_table(True.Sample) %>% as.data.frame()
  if(taxa_are_rows(True.Sample)){
    otu <- t(otu)
  }
  metadata <- as.data.frame(metadata)
  keep_names <- c(batchName,batch_covar)
  keep_names <- keep_names[!duplicated(keep_names)]
  otu.merge <- merge(x=otu, y=metadata[, keep_names, drop = FALSE],  
                     by=0, all.x = TRUE, all.y = FALSE)
  otu.merge[,batchName] <- as.factor(otu.merge[,batchName])
  otu.merge[,batch_covar] <- as.factor(otu.merge[,batch_covar])
  rownames(otu.merge) <- otu.merge$Row.names
  otu.merge <- otu.merge[,-1]
  
  batchid <- otu.merge[,batchName]
  covar <- otu.merge[, batch_covar]
  taxa_tab <- otu.merge[,c(1:(taxa_names(True.Sample) %>% length()))]
  
  taxa_corrected = ConQuR(tax_tab=taxa_tab, batchid=batchid, covariates=covar, batch_ref=batch_ref)
  tmp_Sample <- True.Sample
  otu_table(tmp_Sample) <- otu_table(taxa_corrected %>% t(), taxa_are_rows = TRUE)
  return(tmp_Sample)
}

Wrapper_BatchEffectCorrection_ConQuR <- function(True.Sample, metadata){
  
}

Subtract_Species <- function(Sample1, Sample2){
  set_diff <- setdiff(Sample1 %>% otu_table() %>% rownames() %>% unique(),
                      Sample2 %>% otu_table() %>% rownames() %>% unique()
  )
  return(prune_taxa(set_diff, Sample1))
}

## Inspect sequencing depth
Inspect_SequencingDepth <- function(PhyloObj){
  bac.count <- as_tibble(colSums(as.data.frame(otu_table(PhyloObj))))%>% mutate (id=sample_data(PhyloObj)$id)
  bac.count %>% 
    arrange(value) %>% 
    kable (caption="Bacterial reads", booktabs =TRUE) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>% 
    kable_classic(full_width = F, html_font = "Cambria") %>%
    scroll_box(width = "100%", height = "600px")
}

ConQuR_applyCP <- function(phyloseq,meta_data,batch_type,batch_ref){
  otu <- as(otu_table(phyloseq), "matrix")
  if(taxa_are_rows(phyloseq)){otu <- t(otu)}
  if(grepl("dna", tolower(batch_type), fixed = TRUE)){
    sub.meta.data <- meta_data[meta_data$TRUE_control == "TRUE",]
    sub.meta.data <- sub.meta.data %>% mutate(sample_sidev2 = as.factor(sample_side),
                                              DNAex_roundv2 = as.factor(DNAex_round))
    batchid = sub.meta.data[,c("DNAex_roundv2")]
    covar = sub.meta.data[,c("sample_sidev2")]
  }
  if(grepl("pcr", tolower(batch_type), fixed = TRUE)){
    sub.meta.data <- meta_data[meta_data$TRUE_control == "TRUE",]
    sub.meta.data <- sub.meta.data %>% mutate(sample_sidev2 = as.factor(sample_side),
                                              PCR_roundv2 = as.factor(PCR_round))
    batchid = sub.meta.data[,c("PCR_roundv2")]
    covar = sub.meta.data[,c("sample_sidev2")]
  }
  options(warn = -2)
  taxa_corrected <- ConQuR(tax_tab = otu,batchid = batchid,covariates = covar,batch_ref = batch_ref)
  tmp_phy <- phyloseq
  otu_table(tmp_phy) <- otu_table(taxa_corrected %>% t(),taxa_are_rows = TRUE)
  Plot_PCoA(TAX = otu,factor = batchid,main="Pre-Correction, Bray Curtis")
  Plot_PCoA(TAX = taxa_corrected,factor = batchid ,main = "Post-Correction, Bray Curtis")
  return(tmp_phy)
}