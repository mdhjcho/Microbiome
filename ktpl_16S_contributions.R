### Compare taxonomic contributions to pathways between the 16S (biopsy) and MGS (stool) datarm(list=ls(all=TRUE)).

rm(list=ls(all.names=TRUE))

getwd()
setwd("/Users/Chohyunjeong/Dropbox/KTPL_2021/ktpl_working/ktpl_tables")
source("~/Dropbox/KTPL_2021/ktpl_working/hmp2_util_functions.R")
source("~/Dropbox/KTPL_2021/ktpl_working/picrust2_ms_functions.R")

# Read in stratified HMP2 pathway abundances
ktpl_16S_pathabun_strat <- data.frame(t(readRDS("count_tables/ktpl_pathabun_strat_filt_count_t.rds")), check.names=FALSE)

metadata_0m <- read.table("~/Dropbox/KTPL_2021/ktpl_working/metadata_0m.csv",
                            header=TRUE, sep=",", comment.char="", stringsAsFactors = FALSE)
rownames(metadata_0m) <- metadata_0m$SampleID

raw_R_NR_samples <- c(metadata_0m[which(metadata_0m$Rejection_1YA == "Rejection"), "SampleID"],
                      metadata_0m[which(metadata_0m$Rejection_1YA == "No Rejection"), "SampleID"])

R_NR_samples<-raw_R_NR_samples[which(raw_R_NR_samples  %in% colnames(ktpl_16S_pathabun_strat))]

ktpl_16S_pathabun_strat <- ktpl_16S_pathabun_strat[, R_NR_samples]
ktpl_16S_pathabun_strat_1 <- ktpl_16S_pathabun_strat[-which(rowSums(ktpl_16S_pathabun_strat) == 0), ] ## Error 


# Rearrange both tables to be contributions by genera.
asv_tax_in <- read.table("~/Dropbox/KTPL_2021/ktpl_working/taxonomy.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
asv_tax_in_levels <- add_tax_cols(asv_tax_in)
rownames(asv_tax_in_levels) <- asv_tax_in_levels$Feature.ID


picrust2_orig_col <- colnames(ktpl_16S_pathabun_strat)

ktpl_16S_pathabun_strat$pathway <- as.character(sapply(rownames(ktpl_16S_pathabun_strat), function(x) { stri_split(str=x, regex="\\|")[[1]][1] } ))
ktpl_16S_pathabun_strat$sequence <- as.character(sapply(rownames(ktpl_16S_pathabun_strat), function(x) { stri_split(str=x, regex="\\|")[[1]][2] } ))
ktpl_16S_pathabun_strat$genus <- gsub("^.*g__", "g__", asv_tax_in_levels[ktpl_16S_pathabun_strat$sequence, "Genus"])
ktpl_16S_pathabun_strat$genus[which(ktpl_16S_pathabun_strat$genus == "g__")] <- "unclassified"
ktpl_16S_pathabun_strat <- ktpl_16S_pathabun_strat[, -which(colnames(ktpl_16S_pathabun_strat) == "sequence")]

# Re-convert both to relative abundance.
ktpl_16S_pathabun_strat[, 1:(ncol(ktpl_16S_pathabun_strat) - 2)] <- data.frame(sweep(ktpl_16S_pathabun_strat[, 1:(ncol(ktpl_16S_pathabun_strat) - 2)],
                                                                                     2, colSums(ktpl_16S_pathabun_strat[, 1:(ncol(ktpl_16S_pathabun_strat) - 2)]), '/'), check.names = FALSE) * 100

ktpl_16S_pathabun_strat_genus_sum <- aggregate(. ~ genus + pathway, FUN=sum, data=ktpl_16S_pathabun_strat)


##################

# Also get collapsed genera without requiring that final genus level be classified and save for other analyses.
ktpl_16S_pathabun_strat_full <- data.frame(t(readRDS("count_tables/ktpl_pathabun_strat_filt_count_t.rds")), check.names=FALSE)
ktpl_16S_pathabun_strat_full <- ktpl_16S_pathabun_strat_full[, R_NR_samples]
ktpl_16S_pathabun_strat_full <- ktpl_16S_pathabun_strat_full[-which(rowSums(ktpl_16S_pathabun_strat_full) == 0), ] # Error, Not done

ktpl_16S_pathabun_strat_full$pathway <- as.character(sapply(rownames(ktpl_16S_pathabun_strat_full), function(x) { stri_split(str=x, regex="\\|")[[1]][1] } ))
ktpl_16S_pathabun_strat_full$sequence <- as.character(sapply(rownames(ktpl_16S_pathabun_strat_full), function(x) { stri_split(str=x, regex="\\|")[[1]][2] } ))
ktpl_16S_pathabun_strat_full$genus <- asv_tax_in_levels[ktpl_16S_pathabun_strat_full$sequence, "Genus"]
ktpl_16S_pathabun_strat_full <- ktpl_16S_pathabun_strat_full[, -which(colnames(ktpl_16S_pathabun_strat_full) == "sequence")]
ktpl_16S_pathabun_strat_full[, 1:(ncol(ktpl_16S_pathabun_strat_full) - 2)] <- data.frame(sweep(ktpl_16S_pathabun_strat_full[, 1:(ncol(ktpl_16S_pathabun_strat_full) - 2)],
                                                                                               2, colSums(ktpl_16S_pathabun_strat_full[, 1:(ncol(ktpl_16S_pathabun_strat_full) - 2)]), '/'), check.names = FALSE) * 100
ktpl_16S_pathabun_strat_full_genus_sum <- aggregate(. ~ genus + pathway, FUN=sum, data=ktpl_16S_pathabun_strat_full)

#saveRDS(object = hmp2_16S_pathabun_strat_full_genus_sum, file = "results_out/hmp2_16S_pathabun_strat_full_genus_sum.rds")


library(ALDEx2)
library(ggplot2)
library(stringi)
library(reshape2)
library(sinaplot)
library(cowplot)
install.packages("ggbeeswarm")
library(ggbeeswarm)

# Read in stratified pathways (--per-seq-contrib option used).
pathabun_in_strat<-readRDS("count_tables/ktpl_pathabun_strat_filt_count_t.rds")

# Transpose stratified table and re-separate pathways and sequences as separate columns.
pathabun_in_strat <- data.frame(t(pathabun_in_strat), check.names = FALSE)
orig_col <- colnames(pathabun_in_strat)
pathabun_in_strat$pathway <- as.character(sapply(rownames(pathabun_in_strat), function(x) { stri_split(str=x, regex="\\|")[[1]][1] } ))
pathabun_in_strat$sequence <- as.character(sapply(rownames(pathabun_in_strat), function(x) { stri_split(str=x, regex="\\|")[[1]][2] } ))
pathabun_in_strat <- pathabun_in_strat[, c("pathway", "sequence", orig_col)]


# Unstrat pathways and phenotypes.
pathabun_in_unstrat <- readRDS("count_tables/ktpl_pathabun_filt_count_t.rds")
pathabun_in_unstrat <- data.frame(t(pathabun_in_unstrat), check.names = FALSE)
pheno_in_unstrat <- readRDS("count_tables/ktpl_pheno_count.rds")


# Read in ASVs:
asv_abun <- readRDS("count_tables/ktpl_biom_count.rds")


# Only keep pathways in PICRUSt2 mapfile (i.e. prokaryotic pathways only) that aren't superpathways or engineered.
descrip_gzfile <- gzfile('/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/description_mapfiles/metacyc_pathways_info.txt.gz', 'rt')

pathway_descrip <- read.table(descrip_gzfile, header=FALSE, sep="\t", row.names=1, comment.char="", quote="", stringsAsFactors = FALSE)

close(descrip_gzfile)

### Test for pathways and phenotypes that differ between nonIBD and CD individuals (at lenient cut-off).
R_NR_pathabun_in <- pathabun_in_unstrat[, raw_R_NR_samples[which(raw_R_NR_samples %in% colnames(pathabun_in_unstrat))]]
R_NR_pathabun_in_filt <- R_NR_pathabun_in[which(rowSums(R_NR_pathabun_in > 0) >= 0.33 * ncol(R_NR_pathabun_in)), ]
R_NR_pathabun_in_filt_aldex <- aldex(round(R_NR_pathabun_in_filt), metadata_0m[colnames(R_NR_pathabun_in_filt), "Rejection_1YA"], effect=TRUE)
rownames(R_NR_pathabun_in_filt_aldex)[which(R_NR_pathabun_in_filt_aldex$wi.eBH < 0.1)]
# 45 pathway (q <0.1)

R_NR_pheno_in <- pheno_in_unstrat[, raw_R_NR_samples[which(raw_R_NR_samples %in% colnames(pheno_in_unstrat))]]
R_NR_pheno_in_filt <- R_NR_pheno_in[which(rowSums(R_NR_pheno_in > 0) >= 0.33 * ncol(R_NR_pheno_in )), ]
R_NR_pheno_in_filt_aldex <- aldex(round(R_NR_pheno_in_filt), metadata_0m[colnames(R_NR_pheno_in_filt), "Rejection_1YA"], effect=TRUE)
rownames(R_NR_pheno_in_filt_aldex)[which(R_NR_pheno_in_filt_aldex$wi.eBH < 0.18)]
# q <0.2, many pathway difference

### Identify ASVs and taxa that differ between control and CD individuals (at lenient cut-off).
R_NR_asv_abun <- asv_abun[, raw_R_NR_samples[which(raw_R_NR_samples %in% colnames(pathabun_in_unstrat))]]
R_NR_asv_abun_taxa <- abun_by_taxa_level(asv_tax_in_levels, R_NR_asv_abun)

tail(R_NR_asv_abun_taxa)

R_NR_asv_all_levels <- rbind(R_NR_asv_abun, R_NR_asv_abun_taxa$Species,
                             R_NR_asv_abun_taxa$Genus, R_NR_asv_abun_taxa$Family,
                             R_NR_asv_abun_taxa$Order, R_NR_asv_abun_taxa$Class, R_NR_asv_abun_taxa$Phylum) # abundance table that contain all (asv, phylum~species)


R_NR_asv_all_levels_filt <-  R_NR_asv_all_levels[which(rowSums(R_NR_asv_all_levels > 0) >= 0.33 * ncol(R_NR_asv_all_levels)), ]

R_NR_asv_all_levels_filt_aldex <- aldex(R_NR_asv_all_levels_filt, metadata_0m[colnames(R_NR_pheno_in_filt), "Rejection_1YA"], effect=TRUE)
sig_taxa <- rownames(R_NR_asv_all_levels_filt_aldex)[which(R_NR_asv_all_levels_filt_aldex$wi.eBH < 0.2)]
sig_taxa

# [1] "a7e39fb3e455440bb2ee624c04a2fbc4" ->d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia-Shigella                                                                                  
# [2] "981d989dd8301b017834b0d7fb1eb3e1" -> d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia-Shigella                                                                                   
# [3] "43629854914aac4a4e21fcf1a4720365"  -> d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia-Shigella                                                                                 
# [4] "d__Bacteria; p__Firmicutes; c__Negativicutes; o__Acidaminococcales; f__Acidaminococcaceae; g__Phascolarctobacterium"
# [5] "d__Bacteria; p__Firmicutes; c__Negativicutes; o__Acidaminococcales; f__Acidaminococcaceae"                          
# [6] "d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae"                 
# [7] "d__Bacteria; p__Firmicutes; c__Negativicutes; o__Acidaminococcales"                                                 
# [8] "d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales"                                        
# [9] "d__Bacteria; p__Firmicutes; c__Negativicutes"   


# Plot these significant hits in separate script.
# ktpl_aldex2_hits
R_NR_pathabun_in_filt_relab <- data.frame(sweep(R_NR_pathabun_in, 2, colSums(R_NR_pathabun_in), FUN="/")) * 100
R_NR_asv_abun_relab <- data.frame(sweep(R_NR_asv_abun, 2, colSums(R_NR_asv_abun), FUN="/")) * 100
R_NR_asv_abun_Species_relab <- data.frame(sweep(R_NR_asv_abun_taxa$Species, 2, colSums(R_NR_asv_abun_taxa$Species), FUN="/")) * 100
R_NR_asv_abun_Genus_relab <- data.frame(sweep(R_NR_asv_abun_taxa$Genus, 2, colSums(R_NR_asv_abun_taxa$Genus), FUN="/")) * 100
R_NR_asv_abun_Phylum_relab <- data.frame(sweep(R_NR_asv_abun_taxa$Phylum, 2, colSums(R_NR_asv_abun_taxa$Phylum), FUN="/")) * 100
R_NR_asv_abun_Order_relab <- data.frame(sweep(R_NR_asv_abun_taxa$Order, 2, colSums(R_NR_asv_abun_taxa$Order), FUN="/")) * 100


ASV_samples_meta <- metadata_0m[colnames(R_NR_asv_abun), ]
ASV_R_samples <- ASV_samples_meta[which(ASV_samples_meta$Rejection_1YA == "Rejection"), "SampleID"]
ASV_NR_samples <- ASV_samples_meta[which(ASV_samples_meta$Rejection_1YA == "No Rejection"), "SampleID"]

sig_taxa_df <- data.frame(asv_sig=as.numeric(R_NR_asv_abun_relab["a7e39fb3e455440bb2ee624c04a2fbc4",]), # which ASV should be selected?, combine 3 ASV? 
                          s_sig=as.numeric(R_NR_asv_abun_Species_relab["d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia-Shigella; s__",]),
                          g_sig=as.numeric(R_NR_asv_abun_Genus_relab["d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia-Shigella",]),
                          p_sig=as.numeric(R_NR_asv_abun_Phylum_relab["d__Bacteria; p__Proteobacteria",]),
                          o_sig=as.numeric(R_NR_asv_abun_Order_relab["d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales",]),
                          rejection=metadata_0m[colnames(R_NR_asv_abun_relab), "Rejection_1YA"])


asv_boxplot <- ggplot(sig_taxa_df, aes(x=rejection, y=asv_sig, fill=rejection)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=1, dodge.width=0.8) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  ggtitle("ASV for ") +
  scale_fill_manual(values=c("white", "grey")) +
  theme(legend.position = "none") 

asv_boxplot

species_boxplot <- ggplot(sig_taxa_df, aes(x=rejection, y=s_sig, fill=rejection)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=1, dodge.width=0.8) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  ggtitle("g__Escherichia-Shigella; s__") +
  scale_fill_manual(values=c("white", "grey")) +
  theme(legend.position = "none") 

species_boxplot 

genus_boxplot <- ggplot(sig_taxa_df, aes(x=rejection, y=g_sig, fill=rejection)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=1, dodge.width=0.8) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  ggtitle("g__Escherichia-Shigella") +
  scale_fill_manual(values=c("white", "grey")) +
  theme(legend.position = "none") 

genus_boxplot

phylum_boxplot <- ggplot(sig_taxa_df, aes(x=rejection, y=p_sig, fill=rejection)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=1, dodge.width=0.8) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  ggtitle("p__Proteobacteria") +
  scale_fill_manual(values=c("white", "grey")) +
  theme(legend.position = "none") 

phylum_boxplot

order_boxplot <- ggplot(sig_taxa_df, aes(x=rejection, y=o_sig, fill=rejection)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=1, dodge.width=0.8) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  ggtitle("o__Enterobacterales") +
  scale_fill_manual(values=c("white", "grey")) +
  theme(legend.position = "none") 

order_boxplot



pdf(file = "~/Dropbox/KTPL_2021/ktpl_working/ktpl_figures/ktpl_aldex2_entero.pdf", width=10, height=8)

plot_grid(phylum_boxplot, order_boxplot, genus_boxplot, labels=c("a", "b", "c"))

dev.off()


## phascolarctobacterium

sig_taxa_df1 <- data.frame(s_sig=as.numeric(R_NR_asv_abun_Species_relab["d__Bacteria; p__Firmicutes; c__Negativicutes; o__Acidaminococcales; f__Acidaminococcaceae; g__Phascolarctobacterium; s__",]),
                          g_sig=as.numeric(R_NR_asv_abun_Genus_relab["d__Bacteria; p__Firmicutes; c__Negativicutes; o__Acidaminococcales; f__Acidaminococcaceae; g__Phascolarctobacterium",]),
                          o_sig=as.numeric(R_NR_asv_abun_Order_relab["d__Bacteria; p__Firmicutes; c__Negativicutes; o__Acidaminococcales",]),
                          p_sig=as.numeric(R_NR_asv_abun_Phylum_relab["d__Bacteria; p__Firmicutes",]),
                          rejection=metadata_0m[colnames(R_NR_asv_abun_relab), "Rejection_1YA"])


order_boxplot1 <- ggplot(sig_taxa_df1, aes(x=rejection, y=o_sig, fill=rejection)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=1, dodge.width=0.8) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  ggtitle("o__Acidaminococcales") +
  scale_fill_manual(values=c("white", "grey")) +
  theme(legend.position = "none") 


species_boxplot1 <- ggplot(sig_taxa_df1, aes(x=rejection, y=s_sig, fill=rejection)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=1, dodge.width=0.8) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  ggtitle("g__Phascolarctobacterium; s__") +
  scale_fill_manual(values=c("white", "grey")) +
  theme(legend.position = "none") 

genus_boxplot1 <- ggplot(sig_taxa_df1, aes(x=rejection, y=g_sig, fill=rejection)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=1, dodge.width=0.8) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  ggtitle("g__Phascolarctobacterium") +
  scale_fill_manual(values=c("white", "grey")) +
  theme(legend.position = "none") 


phylum_boxplot1 <- ggplot(sig_taxa_df1, aes(x=rejection, y=p_sig, fill=rejection)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=1, dodge.width=0.8) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  ggtitle("p__Firmicutes") +
  scale_fill_manual(values=c("white", "grey")) +
  theme(legend.position = "none") 

pdf(file = "~/Dropbox/KTPL_2021/ktpl_working/ktpl_figures/ktpl_aldex2_hits_phasco.pdf", width=10, height=8)

plot_grid( phylum_boxplot1, order_boxplot1, genus_boxplot1, labels=c("a", "b", "c"))

dev.off()


# Based on significant taxa instead limit testing to functions contributed by taxa in the 2 groups:
# (1) Those that are at higher RA in NR and (2) those are higher RA in R

# Determine ASVs which fall within the significant taxa categories.
R_sig_higher_ASVs <- c("a7e39fb3e455440bb2ee624c04a2fbc4", "981d989dd8301b017834b0d7fb1eb3e1", "43629854914aac4a4e21fcf1a4720365")
R_sig_higher_ASVs <- c(R_sig_higher_ASVs, rownames(asv_tax_in_levels)[which(asv_tax_in_levels$Species == "d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia-Shigella; s__")])
R_sig_higher_ASVs <- c(R_sig_higher_ASVs, rownames(asv_tax_in_levels)[which(asv_tax_in_levels$Genus == "d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia-Shigella")])
# 315 ASVs

NR_sig_higher_ASVs <- rownames(asv_tax_in_levels)[which(asv_tax_in_levels$Genus == "d__Bacteria; p__Firmicutes; c__Negativicutes; o__Acidaminococcales; f__Acidaminococcaceae; g__Phascolarctobacterium")]
# 97 ASVs

R_sig_higher_out <- compare_func_ratios_by_asv_groups(strat_table = pathabun_in_strat,
                                                           ASVs_of_interest = R_sig_higher_ASVs,
                                                           sample_group1 = ASV_NR_samples[which(ASV_NR_samples %in% colnames(pathabun_in_strat))],
                                                           sample_group2 = ASV_R_samples[which(ASV_R_samples %in% colnames(pathabun_in_strat))])


NR_sig_higher_out <- compare_func_ratios_by_asv_groups(strat_table = pathabun_in_strat,
                                                       ASVs_of_interest = NR_sig_higher_ASVs,
                                                       sample_group1 = ASV_NR_samples[which(ASV_NR_samples %in% colnames(pathabun_in_strat))],
                                                       sample_group2 = ASV_R_samples[which(ASV_R_samples %in% colnames(pathabun_in_strat))])


R_sig_higher_out_wil <- R_sig_higher_out$wilcox
R_sig_higher_out_wil$feature <- paste("Rejection-higher", R_sig_higher_out_wil$feature, sep="_")

NR_sig_higher_out_wil <- NR_sig_higher_out$wilcox
NR_sig_higher_out_wil$feature <- paste("NoRejection-higher", NR_sig_higher_out_wil$feature, sep="_")

combined_out_wil <- rbind(R_sig_higher_out_wil, NR_sig_higher_out_wil)
combined_out_wil$fdr <- p.adjust(combined_out_wil$wilcox_p, "fdr")

combined_out_wil_sig_0.1 <- combined_out_wil[which(combined_out_wil$fdr < 0.1),]
combined_out_wil_sig_0.06 <- combined_out_wil[which(combined_out_wil$fdr < 0.06),]

combined_out_wil_sig_NR_higher <- gsub("NoRejection-higher_", "", combined_out_wil_sig_0.06[grep("NoRejection-higher", combined_out_wil_sig_0.06$feature), "feature"])
combined_out_wil_sig_R_higher <- gsub("Rejection-higher_" , "", combined_out_wil_sig_0.1[grep("Rejection-higher", combined_out_wil_sig_0.1$feature), "feature"])

# Make plots of Proteobacteria/other sig pathways and Clostridia/other sig pathways
NR_sig_higher_ratio_prep <- data.frame(t(NR_sig_higher_out$ratio[combined_out_wil_sig_NR_higher,]), check.names=FALSE)
NR_sig_higher_ratio_prep$sample <- rownames(NR_sig_higher_ratio_prep)
NR_sig_higher_ratio_prep$Rejection_1YA <- metadata_0m[rownames(NR_sig_higher_ratio_prep), "Rejection_1YA"]
NR_sig_higher_ratio_prep_melt <- melt(NR_sig_higher_ratio_prep)
NR_sig_higher_ratio_prep_melt$descrip <- paste(NR_sig_higher_ratio_prep_melt$variable, pathway_descrip[as.character(NR_sig_higher_ratio_prep_melt$variable), "V2"], sep=": ")
NR_sig_higher_ratio_prep_melt$log2ratio <- log2(NR_sig_higher_ratio_prep_melt$value)

NR_sig_higher_ratio_prep_melt$descrip <- factor(NR_sig_higher_ratio_prep_melt$descrip,
                                                levels=c("ANAEROFRUCAT-PWY: homolactic fermentation", "PWY-6897: thiamin salvage II",                                                          "PWY-7315: dTDP-N-acetylthomosamine biosynthesis"))

NR_sig_higher_ratio_plot <- ggplot(NR_sig_higher_ratio_prep_melt, aes(x=descrip, y=log2ratio, fill=Rejection_1YA)) +
  geom_boxplot(width=0.75, outlier.shape = NA) +
  coord_flip() +
  scale_fill_manual(values=c("black", "grey")) +
  xlab("Pathway") +
  ylab("log2((Contributed by Phascolarctobacterium + 1)/(Contributed by other + 1))") +
  labs(fill="Rejection_1YA") +
  theme(legend.position = c(0.8, 0.2),
        legend.background = element_rect(color = "black", 
                                         fill = "white", size = 0.2, linetype = "solid"))

NR_sig_higher_ratio_plot


# Figure out ranking of significant features based on absolute mean differences.
R_sig_higher_out_wil$renamed_feat <- gsub("Rejection-higher_", "", R_sig_higher_out_wil$feature)
R_sig_higher_out_wil_sig <- R_sig_higher_out_wil[which(R_sig_higher_out_wil$renamed_feat %in% combined_out_wil_sig_R_higher),]
R_sig_higher_out_wil_sig_mean_diff_order <- order(R_sig_higher_out_wil_sig$mean_diff)
R_sig_higher_out_wil_sig$full_descrip <- paste(R_sig_higher_out_wil_sig$renamed_feat, pathway_descrip[as.character(R_sig_higher_out_wil_sig$renamed_feat), "V2"], sep=": ")

R_sig_higher_ratio_prep <- data.frame(t(R_sig_higher_out$ratio[combined_out_wil_sig_R_higher,]), check.names=FALSE)
R_sig_higher_ratio_prep$sample <- rownames(R_sig_higher_ratio_prep)
R_sig_higher_ratio_prep$rejection <- ASV_samples_meta[rownames(R_sig_higher_ratio_prep), "Rejection_1YA"]
R_sig_higher_ratio_prep_melt <- melt(R_sig_higher_ratio_prep)
R_sig_higher_ratio_prep_melt$descrip <- paste(R_sig_higher_ratio_prep_melt$variable, pathway_descrip[as.character(R_sig_higher_ratio_prep_melt$variable), "V2"], sep=": ")
R_sig_higher_ratio_prep_melt$descrip <- factor(R_sig_higher_ratio_prep_melt$descrip, levels=c(R_sig_higher_out_wil_sig$full_descrip[R_sig_higher_out_wil_sig_mean_diff_order]))
R_sig_higher_ratio_prep_melt$log2ratio <- log2(R_sig_higher_ratio_prep_melt$value)

R_sig_higher_ratio_plot <- ggplot(R_sig_higher_ratio_prep_melt, aes(x=descrip, y=log2ratio, fill=rejection)) + geom_boxplot(outlier.shape=NA) +
  coord_flip() + scale_fill_manual(values=c("black", "grey")) + xlab("Pathway") + ylab("log2((Contributed by significant Escherichia-Shigella)/(Contributed by other))")

R_sig_higher_ratio_plot 



# Plot stacked barcharts of main contributing to PWY-7315 and PWY0-1533

ktpl_16S_pathabun_strat_genus_sum_full<-ktpl_16S_pathabun_strat_full_genus_sum

# Remove "Bacteria" from genus string:
ktpl_16S_pathabun_strat_genus_sum_full$genus <- gsub("d__Bacteria; ", "", ktpl_16S_pathabun_strat_genus_sum_full$genus)

# Identify genera contributing most abundance across both pathways.
ktpl_16S_pathabun_strat_genus_sum_PWY_7315 <- ktpl_16S_pathabun_strat_genus_sum_full[which(ktpl_16S_pathabun_strat_genus_sum_full$pathway == "PWY-7315"), ]
ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt <- melt(ktpl_16S_pathabun_strat_genus_sum_PWY_7315)
ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt_tmp <- ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt
ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt_tmp <- ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt_tmp[, -3]
ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt_by_genera <- aggregate(value ~ genus + pathway, data=ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt_tmp, FUN=sum)

ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt_by_genera <- ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt_by_genera[with(ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt_by_genera, order(value, decreasing = TRUE)),]

PWY_7315_top_genera <- head(ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt_by_genera$genus, 10)

ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY <- ktpl_16S_pathabun_strat_genus_sum_full[which(ktpl_16S_pathabun_strat_genus_sum_full$pathway == "ANAEROFRUCAT-PWY"), ]
ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY_melt <- melt(ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY)
ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY_melt_tmp <- ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY_melt
ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY_melt_tmp <- ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY_melt_tmp[, -3]
ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY_melt_by_genera <- aggregate(value ~ genus + pathway, data=ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY_melt_tmp, FUN=sum)

ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY_melt_by_genera <- ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY_melt_by_genera[with(ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY_melt_by_genera, order(value, decreasing = TRUE)),]

ANAEROFRUCAT_PWY_top_genera <- head(ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY_melt_by_genera$genus, 10)

overall_top_genera <- sort(unique(c(ANAEROFRUCAT_PWY_top_genera, PWY_7315_top_genera)))

#qual_col <- c('#e6194b', '#3cb44b', 'yellow', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', 'greenyellow', '#fabebe', '#008080', '#e6beff',
#              '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', 'royalblue1', 'grey')

qual_col <- c('#e6194b', '#3cb44b', 'yellow', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', 'greenyellow', '#fabebe', '#008080', '#e6beff','#9a6324', '#fffac8','grey')

ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt$genus_clean <- ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt$genus
ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY_melt$genus_clean <- ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY_melt$genus

ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt[which(! ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt$genus_clean %in%  PWY_7315_top_genera), "genus_clean"] <- "Other"
ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY_melt[which(! ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY_melt$genus_clean %in%  ANAEROFRUCAT_PWY_top_genera), "genus_clean"] <- "Other"

ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt$genus_clean <- factor(ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt$genus_clean, levels=c(overall_top_genera, "Other"))
ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY_melt$genus_clean <- factor(ktpl_16S_pathabun_strat_genus_sum_ANAEROFRUCAT_PWY_melt$genus_clean, levels=c(overall_top_genera, "Other"))


# First get plot of only a shared legend.
ktpl_16S_pathabun_strat_genus_sum_TMP <- melt(ktpl_16S_pathabun_strat_genus_sum_full)
ktpl_16S_pathabun_strat_genus_sum_TMP <- ktpl_16S_pathabun_strat_genus_sum_TMP[which(ktpl_16S_pathabun_strat_genus_sum_TMP$pathway %in% c("PWY-7315", "ANAEROFRUCAT-PWY")),]

ktpl_16S_pathabun_strat_genus_sum_TMP[which(! ktpl_16S_pathabun_strat_genus_sum_TMP$genus %in%  overall_top_genera), "genus"] <- "Other"

ktpl_16S_pathabun_strat_genus_sum_TMP$genus <- factor(ktpl_16S_pathabun_strat_genus_sum_TMP$genus, levels=c(overall_top_genera, "Other"))

tmp_plot <- ggplot(ktpl_16S_pathabun_strat_genus_sum_TMP, aes(x=variable, y=value, fill=genus)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=qual_col) + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        legend.text=element_text(size=8)) +
  labs(fill="Genus")

tmp_plot
grobs <- ggplotGrob(tmp_plot)$grobs
stacked_legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

ktpl_16S_pathabun_strat_genus_sum_TMP$rejection <- metadata_0m[match(ktpl_16S_pathabun_strat_genus_sum_TMP$variable,metadata_0m$SampleID), "Rejection_1YA"]


tmp_plot1 <- ggplot(ktpl_16S_pathabun_strat_genus_sum_TMP, aes(x=rejection, y=value, fill=genus)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=qual_col) + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
            legend.text=element_text(size=8)) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  ggtitle("PWY-7315: dTDP-N-acetylthomosamine biosynthesis") +
  labs(fill="Genus")

tmp_plot1



# Get samples ordered by total relative abundance.
ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt_by_sample <- aggregate(value ~ variable, data=ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt, FUN=sum)
ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt_by_sample <- ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt_by_sample[with(ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt_by_sample, order(value, decreasing = FALSE)),]

ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt$variable <- factor(ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt$variable,
                                                                   levels=ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt_by_sample$variable)

PWY_7315_col <- qual_col[which(levels(ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt$genus_clean) %in% ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt$genus_clean)]

PWY_7315_stacked <- ggplot(ktpl_16S_pathabun_strat_genus_sum_PWY_7315_melt, aes(x=variable, y=value, fill=genus_clean)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=PWY_7315_col) + 
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text.x=element_blank(),
                     axis.text=element_text(size=12),
                     axis.title=element_text(size=14),
                     plot.title = element_text(hjust=0.2, vjust=-10, face="bold")) +
  ylab("Relative Abundance (%)") +
  xlab("Sample") +
  ggtitle("PWY-7315: Chondroitin sulfate degradation I (bacterial)") +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill=FALSE)



#  additional analysis
#1. Supp Clostridiales contributed functions.R
#2. Supp_CD_validation.R

getwd()
write.csv(R_NR_asv_all_levels_filt_aldex,"R_NR_asv_all_levels_filt_aldex.csv")
write.csv(R_NR_pathabun_in_filt_aldex,"R_NR_pathabun_in_filt_aldex.csv")
