AlphaPlot <- function(PhyloObj, index = "Observed", strata = "Treatment", y_label = "Observed species", add_legend = FALSE){
  #group.colors <- c(ctrl = "dodgerblue3", Gem = "firebrick2")
  condition_names <- sample_data(PhyloObj)[[strata]] %>% table() %>% names()
  group.colors <- my_cols[1:length(condition_names)]
  names(group.colors) <- condition_names
  #my.labels <- c("Ctrl", "Gem")
  my.labels <- condition_names
  rich_meta <- merge(PhyloObj %>% sample_data(), PhyloObj %>% estimate_richness(), by = "row.names")
  ##-- only using statistics test
  ob <- rich_meta %>% t_test(as.formula(paste0(index, " ~ ", strata))) %>% adjust_pvalue(method = "BH") %>%  add_significance("p.adj") %>% add_xy_position()
  
  ##-----
  p1 <- ggplot (rich_meta, aes_string (x=strata, y=index, fill=strata))+ 
    geom_boxplot()+
    geom_point (position=position_jitterdodge( jitter.width = 0.05))+
    theme_gray() + 
    scale_x_discrete(labels= my.labels)+
    xlab(strata)+
    ylab(y_label)+
    # facet_grid(.~tp)+
    #scale_x_discrete(labels= my.labels)+
    scale_fill_manual(values=group.colors, labels = c(""))+
    # ggtitle("KPC tumor vs. Healthy pancreas - Observed species") +
    theme(axis.text.y = element_text (size=11),
          axis.title = element_text(size=12, face="bold"))
  
  if(!add_legend){    
    p1 <- p1 + theme(legend.position = "none",
                     # legend.background = element_rect(size=0.25, linetype="solid", colour ="black"),
                     legend.key.size = unit(4,"mm"),
                     #axis.text.x = element_blank(),
                     plot.title = element_text(size = 12))
    
  } else {
    p1 <- p1 + theme(legend.text = element_text(size = 12),
                     #legend.title = element_blank(),
                     # legend.background = element_rect(size=0.25, linetype="solid", colour ="black"),
                     legend.key.size = unit(4,"mm"),
                     axis.text.x = element_blank(),
                     plot.title = element_text(size = 12))
  }
  p1 <- p1 + stat_pvalue_manual(ob, label = "p.adj", inherit.aes = FALSE, tip.length = 0.01)
  return(p1)
}
###--------------------
AlphaPlotWrapper <- function(PhyloObj, roundUp = TRUE){
  ## round up/down otu_table
  if(roundUp){
    otu_table(PhyloObj) <- PhyloObj %>% otu_table() %>% round()
  } else{
    otu_table(PhyloObj) <- PhyloObj %>% otu_table() %>% ceiling()
  }
  plt.1 <- AlphaPlot(PhyloObj, index = "Observed", y_label = "Observed species", add_legend = F)
  plt.2 <- AlphaPlot(PhyloObj, index = "Shannon", y_label = "Shannon index", add_legend = F)
  plt.3 <- AlphaPlot(PhyloObj, index = "InvSimpson", y_label = "Inv Simpson index", add_legend = F)
  return(list("Observed" = plt.1, "Shannon" = plt.2, "InvSimpson" = plt.3))
}
###----- Beta --------------
BetaPlot <- function(PhyloObjct, strata_f , n_per = 1000, type = NULL, title_method = "Original", dis_method = "bray", ordination_method = "PCoA"){
  bray_dist = phyloseq::distance(PhyloObjct, method=dis_method)
  ordination = ordinate(PhyloObjct, method=ordination_method, distance=bray_dist)
  
  pcoa1 <- paste("PCoA 1 [", round(ordination[[3]]$Relative_eig[1], digits = 3)*100, "%]", sep = "")
  pcoa2 <- paste("PCoA 2 [", round(ordination[[3]]$Relative_eig[2], digits = 3)*100, "%]", sep = "")
  
  ##p.adonis <- adonis2(bray_dist ~ sample_data(PhyloObjct)$Treatment)
  p.adonis <- adonis2(as.formula(paste0("bray_dist ~ ", strata_f)), data = PhyloObjct %>% sample_data() %>% as_tibble(), permutations = n_per)
  
  p <- case_when(
    p.adonis$`Pr(>F)`[1] > 0.05 ~ paste("p =", p.adonis$`Pr(>F)`[1], "n.s.", sep = " "),
    p.adonis$`Pr(>F)`[1] < 0.05 &  p.adonis$`Pr(>F)`[1] > 0.01 ~ paste("p =", p.adonis$`Pr(>F)`[1], "*", sep = " "),
    p.adonis$`Pr(>F)`[1] <= 0.01 & p.adonis$`Pr(>F)`[1] > 0.001  ~ paste("p =", p.adonis$`Pr(>F)`[1], "**", sep = " "),
    p.adonis$`Pr(>F)`[1] <= 0.001 ~ paste("p =",p.adonis$`Pr(>F)`[1], "***", sep = " "),
  )
  
  annotations <- data.frame(
    xpos = c(-Inf),
    ypos =  c(Inf),
    annotateText = p,
    hjustvar = c(-0.2) ,
    vjustvar = c(1.5))
  
  p1 <- plot_ordination(PhyloObjct, ordination, color = strata_f) +
    geom_point(aes_string(colour = strata_f, shape = type), size = 3) +
    #geom_point(aes(colour = .data[[strata_f]], shape = .data[[type]]), size = 3) +
    #geom_point(aes(colour = sym(strata_f), shape = sym(type)), size = 3)
    ##geom_text_repel(aes(label = id), size = 4) +
    theme(aspect.ratio=1) +
    theme_bw()+
    scale_color_brewer(palette = "Set1")+
    stat_ellipse() +
    xlab(pcoa1)+
    ylab(pcoa2)+
    theme(panel.grid =  element_blank())+
    ggtitle(paste0(title_method, "")) +
    theme (axis.text=element_text(size=14),
           axis.title=element_text(size=16,face="bold"),
           legend.text = element_text(size = 12),
           legend.title = element_blank())+
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size = 4.5, inherit.aes = FALSE)
  return(p1)
}
##### ---- Beta simple
BetaPlot_Simple <- function(PhyloObjct, strata_f = NULL, type = NULL, title_method = "Original", dis_method = "bray", ordination_method = "PCoA"){
  bray_dist = phyloseq::distance(PhyloObjct, method=dis_method)
  ordination = ordinate(PhyloObjct, method=ordination_method, distance=bray_dist)
  
  pcoa1 <- paste("PCoA 1 [", round(ordination[[3]]$Relative_eig[1], digits = 3)*100, "%]", sep = "")
  pcoa2 <- paste("PCoA 2 [", round(ordination[[3]]$Relative_eig[2], digits = 3)*100, "%]", sep = "")
  
  ##p.adonis <- adonis2(bray_dist ~ sample_data(PhyloObjct)$Treatment)
  
  
  
  
  p1 <- plot_ordination(PhyloObjct, ordination, color = strata_f, shape = type) +
    geom_point(aes_string(colour = strata_f, shape = type), size = 3) +
    #geom_text_repel(aes(label = id), size = 4) +
    theme(aspect.ratio=1) +
    theme_bw()+
    scale_color_brewer(palette = "Set1")+
    #stat_ellipse() +
    xlab(pcoa1)+
    ylab(pcoa2)+
    theme(panel.grid =  element_blank())+
    ggtitle(paste0(title_method)) +
    theme (axis.text=element_text(size=14),
           axis.title=element_text(size=16,face="bold"),
           legend.text = element_text(size = 12),
           legend.title = element_blank())
  # geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size = 4.5, inherit.aes = FALSE)
  return(p1)
}
##-------- alternative for plot_bar by phyloseq
plot_bar_2 <- function(physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, 
                       facet_grid = NULL, border_color = NA)
{
  mdf = psmelt(physeq)
  p = if(!is.null(fill) & fill %in% colnames(mdf)) { 
    ggplot(mdf, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]]))
  } else {
    ggplot(mdf, aes(x = .data[[x]], y = .data[[y]]))
  }
  p = p + geom_bar(stat = "identity", position = "stack", color = border_color)
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  p
}
## ---- wrench normalization wrapper
WrenchWrapper <- function(PhyloObjct, grp, roundUp = F){
  cnt_table <- PhyloObjct %>% otu_table()
  group <- PhyloObjct %>% sample_data() %>% pull(grp)
  w <- wrench(cnt_table, condition = group)
  
  # deseq.obj <- DESeqDataSetFromMatrix(cnt_table %>% as.data.frame(), DataFrame(group), ~group)
  # DESeq2::sizeFactors(deseq.obj) <- w$nf
  # cnt_table_normalized <- DESeq2::counts(deseq.obj, normalized=TRUE)
  
  norm_factors <- w$nf
  norm_counts <- sweep(cnt_table, 2, norm_factors, FUN = '/')
  if(roundUp){norm_counts <- norm_counts %>% round()}
  return(phyloseq(otu_table(norm_counts, taxa_are_rows = T), tax_table(PhyloObjct %>% tax_table()), sample_data(PhyloObjct %>% sample_data())))
}
## MAAsLin2 Wrapper
## Apply MaAslin2 for phyloseq object after decontam and normalized
MaAsLin2_Wrapper <- function(PhyLoObj, barcode, Treatment, OutDir){
  ## Normalization: rarefying (rar), Wrench (wrench)
  
  df_data <- PhyLoObj %>% otu_table()  %>% t() %>% as.data.frame()
  df_metadata <- PhyLoObj %>% sample_data() %>% as_tibble() %>% dplyr::select(all_of(barcode), all_of(Treatment)) %>% as.data.frame()
  row.names(df_metadata) <- df_metadata[[barcode]]
  fit.RmLowAbun.Rar <- Maaslin2(input_data = df_data,
                                input_metadata = df_metadata,
                                #output = "./MaAsLin2_OutDir/RmLowabun.Rar",
                                output = OutDir,
                                min_abundance = 0,
                                min_prevalence = 0,
                                normalization = "NONE",
                                fixed_effects = c("Treatment"),
                                reference = c("Treatment,ctrl"))
}
## rarefying data
rarefy_even_depth_wrapper <- function(PhyLoObj, seed=911){
  return(rarefy_even_depth(PhyLoObj, sample.size = PhyLoObj %>% sample_sums() %>% min(), rngseed = seed))
}
## V2
MaAsLin2_Wrapper2 <- function(df_data, df_metadata, Treatment, TreatMent_ref, OutDir){
  ## Normalization: rarefying (rar), Wrench (wrench)
  
  # df_data <- PhyLoObj %>% otu_table()  %>% t() %>% as.data.frame()
  # df_metadata <- PhyLoObj %>% sample_data() %>% as_tibble() %>% dplyr::select(all_of(barcode), all_of(Treatment)) %>% as.data.frame()
  # row.names(df_metadata) <- df_metadata[[barcode]]
  fit.RmLowAbun.Rar <- Maaslin2(input_data = df_data,
                                input_metadata = df_metadata,
                                #output = "./MaAsLin2_OutDir/RmLowabun.Rar",
                                output = OutDir,
                                min_abundance = 0,
                                min_prevalence = 0,
                                normalization = "NONE",
                                fixed_effects = c(Treatment),
                                #reference = c("Treatment,ctrl")
                                reference = c(paste0(Treatment, ",", TreatMent_ref))
                                )
}
## rarefying data
rarefy_even_depth_wrapper <- function(PhyLoObj, seed=911){
  return(rarefy_even_depth(PhyLoObj, sample.size = PhyLoObj %>% sample_sums() %>% min(), rngseed = seed))
}
###---------------------------------------------------------------
## plot MaAsLin2 results
MaAslin2_plot <- function(Path2Tab, Tax_tab, title_method=""){
  ## plot significant results
  ## Read table
  Tab <- read.table(Path2Tab, header = TRUE)
  if(dim(Tab)[1]==0){return(NULL)}
  Tab <- Tab %>% as_tibble() %>% 
    mutate(TaxaID = substr(feature, 2, nchar(feature))) %>% 
    left_join(., Tax_tab[, c("TaxaID", "species")], by = "TaxaID") %>% 
    mutate(legend = paste0(TaxaID, ":", species))
  
  Tab$legend = factor(Tab$legend,levels=Tab$legend[order(Tab$coef)])
  plt <- ggplot(Tab, aes(x=legend,y=coef,fill=coef>0))+
    geom_col() + coord_flip()+
    scale_fill_manual(values=c("blue","red"),
                      labels=c("negative","positive"))
  plt <- plt + xlab("") +
    ggtitle(paste0(title_method, " - MaAsLin2")) +
    theme (axis.title=element_text(size=12,face="bold"),
           legend.title = element_blank())
  return(plt)
}
## MaAslin_plot version 2
MaAslin2_plot_2 <- function(Path2Tab, coef_thres = 2.0, qval_thres = 0.05, Tax_tab, legend_rank = "family", title_method=""){
  Tab <- read.table(Path2Tab, header = TRUE)
  if(dim(Tab)[1]==0){return(NULL)}
  Tab <- Tab %>% as_tibble() %>% 
    dplyr::filter(coef >=coef_thres | coef <= -coef_thres) %>% 
    dplyr::filter(qval <= qval_thres) %>% 
    mutate(TaxaID = substr(feature, 2, nchar(feature))) %>% 
    left_join(., Tax_tab[, c("TaxaID", "species", legend_rank)], by = "TaxaID") %>% 
    mutate(legend = paste0(TaxaID, ":", species)) %>% 
    mutate(fam_col = Kolors[factor(.data[[legend_rank]]) %>% as.integer()])
  
  Tab$legend = factor(Tab$legend,levels=Tab$legend[order(Tab$coef)])
  Tab[[legend_rank]] = factor(Tab[[legend_rank]])
  KKK <- Tab[, c(legend_rank, "fam_col")] %>% unique()
  
  plt_data <- Tab %>% 
    #mutate(fam_col = Kolors[factor(.data[["family"]]) %>% as.integer()]) %>% 
    ggplot()+
    geom_col(aes(x=legend,y=coef,fill=coef>0))+
    coord_flip()+
    scale_fill_manual(values=c("blue","red"),
                      labels=c("negative","positive"))+
    theme(axis.text.y=element_blank(),
          axis.title = element_blank()
    )
  #coord_cartesian(expand = F)
  
  plt_anno <- Tab %>% 
    #mutate(fam_col = Kolors[factor(.data[["family"]]) %>% as.integer()]) %>% 
    ggplot()+
    geom_col(aes(x=1, y=legend, fill = .data[[legend_rank]]))+
    #scale_fill_brewer(type = "qual", )+
    scale_fill_manual(values = KKK$fam_col
    )+
    #scale_fill_identity(name = waiver(), labels = LETTERS[1:18])+
    theme(axis.title = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.text = element_blank()
    )
  
  plt_anno2 <- Tab %>% 
    #mutate(fam_col = Kolors[factor(.data[["family"]]) %>% as.integer()]) %>% 
    ggplot()+
    geom_col(aes(x=1, y=legend, fill = .data[[legend_rank]]))+
    #scale_fill_brewer(type = "qual", )+
    scale_fill_manual(values = KKK$fam_col
    )+
    #scale_fill_identity(name = waiver(), labels = LETTERS[1:18])+
    theme(axis.title = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks = element_blank(),
          legend.position="bottom"
          #legend.position = "none",
          #legend.text = element_blank()
    )
  m_legend <- get_legend(plt_anno2)
  return(list(gg_data = plt_data, gg_anno = plt_anno, gg_legend=m_legend))
}
###---------------------------------------------------------------
## plot ALDEx2 results
ALDex2_plot <- function(x.all, Tax_tab, p_val_thres = 0.05, title_method=""){
  Tab <- x.all %>% filter(we.ep <= p_val_thres) %>% 
    rownames_to_column(var="TaxaID") %>% 
    left_join(., Tax_tab[, c("TaxaID", "species")], by = "TaxaID") %>% 
    mutate(legend = paste0(TaxaID, ":", species))
  if(dim(Tab)[1]==0){return(NULL)}
  Tab$legend = factor(Tab$legend,levels=Tab$legend[order(Tab$diff.btw)])
  plt <- ggplot(Tab, aes(x=legend,y=diff.btw,fill=diff.btw>0))+
    geom_col() + coord_flip()+
    scale_fill_manual(values=c("blue","red"),
                      labels=c("negative","positive"))
  plt <- plt + xlab("") +
    ggtitle(paste0(title_method, " - ALDEx2")) +
    theme (axis.title=element_text(size=12,face="bold"),
           axis.text = element_text(size = 8),
           legend.title = element_blank())
  return(plt)
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
## filter taxa by low abundance
filter_by_low_abundance <- function(phyloseqObj, threshold, A){
  physeq.abun.filt <- phyloseq::transform_sample_counts(phyloseqObj, function(x){x/sum(x)})
  noCont <- phyloseq::genefilter_sample(physeq.abun.filt, filterfun_sample(function(x) x > threshold), A= A)
  physeq.abun.filt <- phyloseq::prune_taxa(taxa = noCont,x = phyloseqObj)
  return(physeq.abun.filt)
}

## Wrapper of Decontam
Wrapper_Decontam <- function(phylo, neg_f, thres = 0.1){
  is_contam <- phylo %>% 
    #ps_mutate(is_neg = ifelse(true.control == "true", FALSE, TRUE)) %>% 
    #filter_by_low_abundance(threshold = 1e-4, A = 1) %>% 
    #ntaxa() %>% 
    #Inspect_SequencingDepth() %>% 
    isContaminant(method="prevalence", neg=neg_f, threshold = thres, normalize = FALSE, detailed = FALSE) 
  return(prune_taxa(!is_contam, phylo))
}
#Plot Venn diagram for two sets
PlotVenn2Sets <- function(set1, set2, nameset1, nameset2){
  plt <- venn.diagram(
    x = list(set1, set2),
    category.names = c(nameset1 , nameset2),
    disable.logging = TRUE,
    filename = NULL
  )
  grid::grid.newpage()
  grid::grid.draw(plt)
}
## pairwise adonis for phyloseq object
ps_pwAdonis2 <- function(PhyloObjct, feature, m_strata = NULL, dis_method = "bray"){
  bray_dist = phyloseq::distance(PhyloObjct, method=dis_method)
  res_adonis <- pairwise.adonis2(x = as.formula(paste0("bray_dist ~ ", feature)),
                                 data = PhyloObjct %>% sample_data() %>% as_tibble() %>% as.data.frame(),
                                 strata = m_strata)
  return(res_adonis)
}
## format result from pairwise.adonis2 
reshape_adonis2 <- function(res_adonis){
  pw_name <- lapply(names(res_adonis)[2:length(res_adonis)], function(x){strsplit(x,"_vs_")[[1]]})
  dat <- data.frame(f1 = character(0), f2 = character(0), pVal = numeric(0))
  for (i in 1:length(pw_name)) {
    v1 <- pw_name[[i]][1]
    v2 <- pw_name[[i]][2]
    v3 <- res_adonis[[i+1]]$`Pr(>F)`[1]
    dat[nrow(dat)+1, ] <- c(v1, v2, v3)
  }
  dat <- dat %>% reshape(idvar = "f1", timevar = "f2", direction = "wide") %>% as_tibble() %>% 
    remove_rownames() %>% column_to_rownames(var = "f1")
  #dat %>% remove_rownames() %>%  %>% rename()
  return(dat)
}
## randomly subset samples
sample_ps <- function(ps, FUN = sample, ...){
  ids <- sample_names(ps)
  sampled_ids <- FUN(ids, ...)
  ps <- prune_samples(sampled_ids, ps)
  return(ps)
}
