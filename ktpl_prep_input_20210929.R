### Commands to prep all files for ktpl analyses.
### i.e. transform, filter, and subset tables and write out to be used for subsequent analyses.

rm(list=ls(all=TRUE))

source("~/Dropbox/KTPL_2021/ktpl_working/hmp2_util_functions.R")

# Read in metadata_0m table.

metadata_0m <- read.table("~/Dropbox/KTPL_2021/ktpl_working/metadata_0m.csv",
                            header=TRUE, sep=",", comment.char="", stringsAsFactors = FALSE, check.names = TRUE)

### Prep PICRUSt2 pathway predictions ###
# Read in stratified per-seq-contrib pathway abundances.
ktpl_pathabun_strat <- read.table("~/Dropbox/KTPL_2021/ktpl_working/path_abun_strat.tsv",
                                  header=TRUE, sep="\t", stringsAsFactors = FALSE)

# Get unstratified pathway abundances based on the per-seq-contrib table (so that contributors will sum to unstratified abundance).
ktpl_pathabun_strat_tmp <- ktpl_pathabun_strat[, -which(colnames(ktpl_pathabun_strat) == "sequence")] # remove sequence
ktpl_pathabun <- aggregate(. ~ pathway, data=ktpl_pathabun_strat_tmp, FUN=sum) #sum pathway
rm(ktpl_pathabun_strat_tmp)
rownames(ktpl_pathabun) <- ktpl_pathabun$pathway
ktpl_pathabun <- ktpl_pathabun[, -which(colnames(ktpl_pathabun) == "pathway")]


# Remove any pathways classified as "engineered" or "superpathways".
descrip_gzfile <- gzfile('~/Dropbox/KTPL_2021/ktpl_working/metacyc_pathways_info.txt.gz', 'rt')

pathway_descrip <- read.table(descrip_gzfile, header=FALSE, sep="\t", row.names=1, comment.char="", quote="", stringsAsFactors = FALSE)

close(descrip_gzfile)

pathway_descrip_subset <- pathway_descrip[rownames(ktpl_pathabun),, drop=FALSE]

path2remove_i <- c(grep("superpathway", pathway_descrip_subset$V2), grep("engineered", pathway_descrip_subset$V2))

path2remove <- rownames(pathway_descrip_subset)[path2remove_i]

ktpl_pathabun_filt <- ktpl_pathabun[-which(rownames(ktpl_pathabun) %in% path2remove), ]

ktpl_pathabun_filt_count <- ktpl_pathabun_filt # file that removed 'superpathway' 'engineered'


# Convert to relative abundance and perform asin transformation.
ktpl_pathabun_filt <- data.frame(sweep(ktpl_pathabun_filt, 2, colSums(ktpl_pathabun_filt), '/'), check.names=FALSE)
ktpl_pathabun_filt_asin <- asinTransform(ktpl_pathabun_filt)
ktpl_pathabun_filt_asin_t<-data.frame(t(ktpl_pathabun_filt_asin),check.names=FALSE)


# Split into ileum and rectum samples and re-name columns to be participant ids.(Not done)
# Sort mapfiles so that 1st week should be first (in case of duplicates from same patient) and the deduplicate.
# Also transpose final abundance tables.
biopsy_16S_meta <- hmp2_metadata[which(hmp2_metadata$data_type == "biopsy_16S"), ]

biopsy_16S_meta_Ileum <- biopsy_16S_meta[which(biopsy_16S_meta$biopsy_location == "Ileum"), ]
biopsy_16S_meta_Ileum <- biopsy_16S_meta_Ileum[with(biopsy_16S_meta_Ileum, order(Participant.ID, week_num)),]
biopsy_16S_meta_Ileum <- biopsy_16S_meta_Ileum[-which(duplicated(biopsy_16S_meta_Ileum$Participant.ID)), ]
biopsy_16S_meta_Ileum <- biopsy_16S_meta_Ileum[which(biopsy_16S_meta_Ileum$External.ID %in% colnames(hmp2_pathabun_filt_asin)), ]

hmp2_pathabun_filt_asin_Ileum <- data.frame(t(hmp2_pathabun_filt_asin[, biopsy_16S_meta_Ileum$External.ID]), check.names=FALSE)
rownames(hmp2_pathabun_filt_asin_Ileum) <- biopsy_16S_meta_Ileum$Participant.

biopsy_16S_meta_Rectum <- biopsy_16S_meta[which(biopsy_16S_meta$biopsy_location == "Rectum"), ]
biopsy_16S_meta_Rectum <- biopsy_16S_meta_Rectum[with(biopsy_16S_meta_Rectum, order(Participant.ID, week_num)),]
biopsy_16S_meta_Rectum <- biopsy_16S_meta_Rectum[-which(duplicated(biopsy_16S_meta_Rectum$Participant.ID)), ]
biopsy_16S_meta_Rectum <- biopsy_16S_meta_Rectum[which(biopsy_16S_meta_Rectum$External.ID %in% colnames(hmp2_pathabun_filt_asin)), ]

hmp2_pathabun_filt_asin_Rectum <- data.frame(t(hmp2_pathabun_filt_asin[, biopsy_16S_meta_Rectum$External.ID]), check.names=FALSE)
rownames(hmp2_pathabun_filt_asin_Rectum) <- biopsy_16S_meta_Rectum$Participant.ID


##### Also prep stratified pathway predictions. ##### 
# Note that this data wont be transformed since we're interested in the relative abundance of contributors.

ktpl_pathabun_strat_filt <- ktpl_pathabun_strat[-which(ktpl_pathabun_strat$pathway %in% path2remove), ]
rownames(ktpl_pathabun_strat_filt) <- paste(ktpl_pathabun_strat_filt$pathway, ktpl_pathabun_strat_filt$sequence, sep="|")
ktpl_pathabun_strat_filt <- ktpl_pathabun_strat_filt[, -which(colnames(ktpl_pathabun_strat_filt) %in% c("pathway", "sequence"))]

ktpl_pathabun_strat_filt_count <- ktpl_pathabun_strat_filt

# Convert to relative abundance
ktpl_pathabun_strat_filt <- data.frame(sweep(ktpl_pathabun_strat_filt, 2, colSums(ktpl_pathabun_strat_filt), '/'), check.names=FALSE)

ktpl_pathabun_strat_filt_t <- data.frame(t(ktpl_pathabun_strat_filt), check.names=FALSE)

# Subset to Ileum and Rectum samples. (Not done)
hmp2_pathabun_strat_filt_Ileum <- data.frame(t(hmp2_pathabun_strat_filt[, biopsy_16S_meta_Ileum$External.ID]), check.names=FALSE)
rownames(hmp2_pathabun_strat_filt_Ileum) <- biopsy_16S_meta_Ileum$Participant.ID

hmp2_pathabun_strat_filt_Rectum <- data.frame(t(hmp2_pathabun_strat_filt[, biopsy_16S_meta_Rectum$External.ID]), check.names=FALSE)
rownames(hmp2_pathabun_strat_filt_Rectum) <- biopsy_16S_meta_Rectum$Participant.ID

####### Prep ASV relative abundances ####### 

ktpl_biom <- read.table("~/Dropbox/KTPL_2021/ktpl_working/ktpl_biom.tsv",
                        header=TRUE, comment.char="", sep="\t", row.names=1)

ktpl_biom_count <- ktpl_biom

ktpl_biom <- data.frame(sweep(ktpl_biom, 2, colSums(ktpl_biom), '/'), check.names = FALSE)

# Also read in taxa classifications for this subset of ASVs.
ktpl_asv_taxa <- read.table("~/Dropbox/KTPL_2021/ktpl_working/taxonomy.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE)

ktpl_asv_taxa_subset <- ktpl_asv_taxa[which(ktpl_asv_taxa$Feature.ID %in% rownames(ktpl_biom)), ]

####### Prep IMG phenotype predictions ####### 
ktpl_pheno <- read.table("~/Dropbox/KTPL_2021/ktpl_working//pred_metagenome_unstrat.tsv",
                         header=TRUE,row.names=1, sep="\t")

ktpl_pheno_count <- ktpl_pheno

ktpl_pheno <- data.frame(sweep(ktpl_pheno, 2, colSums(ktpl_pheno), '/'), check.names = FALSE)

ktpl_pheno_asin <- asinTransform(ktpl_pheno)


##### Also get count tables formatted for pathways, ASVs, and phenos.
ktpl_pathabun_filt_count_t <- data.frame(t(ktpl_pathabun_filt_count), check.names=FALSE)

ktpl_pathabun_strat_filt_count_t <- data.frame(t(ktpl_pathabun_strat_filt_count), check.names=FALSE)

ktpl_biom_count_t <- ktpl_biom_count[, biopsy_16S_meta_Ileum$External.ID]
colnames(hmp2_biom_count_Ileum) <- biopsy_16S_meta_Ileum$Participant.ID


# Write out full tables.
saveRDS(object = ktpl_pathabun_strat_filt_t,
        file = "~/Dropbox/KTPL_2021/ktpl_working/prepped_tables/ktpl_pathabun_strat_filt_t.rds")

saveRDS(object = ktpl_pathabun_filt_asin_t,
        file = "~/Dropbox/KTPL_2021/ktpl_working/prepped_tables/ktpl_pathabun_filt_asin_t.rds")

saveRDS(object = ktpl_biom,
        file = "~/Dropbox/KTPL_2021/ktpl_working/prepped_tables/ktpl_biom.rds")

saveRDS(object = ktpl_pheno_asin,
        file = "~/Dropbox/KTPL_2021/ktpl_working/prepped_tables/ktpl_pheno_asin.rds")


# Also save tables with raw counts.
saveRDS(object = ktpl_pathabun_filt_count_t,
        file = "~/Dropbox/KTPL_2021/ktpl_working/prepped_tables/ktpl_pathabun_filt_count_t.rds")

saveRDS(object = ktpl_pathabun_strat_filt_count_t,
        file = "~/Dropbox/KTPL_2021/ktpl_working/prepped_tables/ktpl_pathabun_strat_filt_count_t.rds")

saveRDS(object = ktpl_biom_count,
        file = "~/Dropbox/KTPL_2021/ktpl_working/count_tables/ktpl_biom_count.rds")

saveRDS(object = ktpl_pheno_count,
        file = "~/Dropbox/KTPL_2021/ktpl_working/count_tables/ktpl_pheno_count.rds")
