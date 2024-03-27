#############################################################
# 
# Author: Julia Corell
# Date: January 25, 2024
# Description: miRNA target identification and KEGG 
# enrichment analysis
#
#############################################################

# Load the necessary libraries
suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(ggpubr))

# Paths
MTI_file <- "../../Additional_data/hsa_MTI.xlsx" # Get the file from mirTarBase (release 9.0)
DE_results <- "/path/to/results/07-DE"
output_dir <- "/path/to/results/09-miRNA_and_tiny_figures"
anno_file <- "/path/to/annotation_results/unitas_simplified_annnotation/curated_unitas_withtiny_annotation.tsv"

#############################################################
# MTI Filtering
#############################################################

MTI_df <- readxl::read_excel(MTI_file)

# We know that the stronger experimental evidences are: reporter assay,
# western blot and qPCR
#MTI_df_exp <- MTI_df[, c("miRTarBase ID", "miRNA", "Target Gene", "Experiments")]
MTI_df_exp <- MTI_df

# Separate each experiment type
MTI_df_exp$Reporter <- ifelse(grepl("Luciferase reporter assay", MTI_df_exp$Experiments), 1, 0)
MTI_df_exp$Western <- ifelse(grepl("Western blot", MTI_df_exp$Experiments), 1, 0)
MTI_df_exp$qPCR <- ifelse(grepl("qRT-PCR", MTI_df_exp$Experiments), 1, 0)

# Sum the strong validation experiments
MTI_df_exp$Sum_exp <- rowSums(MTI_df_exp[, c("Reporter", "Western", "qPCR")])

# Get only the MTIs with at least 2 exp. evidences
MTI_df_filt <- MTI_df_exp %>%
  filter(Sum_exp >= 2)

# Some entries are duplicated. Filter duplicates.
MTI_df_final <- MTI_df_filt[!duplicated(MTI_df_filt$`miRTarBase ID`),]


#############################################################
# Consistent miRNA
#############################################################

# Read files
anno_df <- read.table(anno_file, sep="\t", header=TRUE)
DE_T1_acute <- read.table(paste0(DE_results, "/results_DE_T1acute.tsv"), sep="\t", header=TRUE)
DE_T1_mild <- read.table(paste0(DE_results, "/results_DE_T1mild.tsv"), sep="\t", header=TRUE)
DE_T2_acute <- read.table(paste0(DE_results, "/results_DE_T2acute.tsv"), sep="\t", header=TRUE)
DE_T2_mild <- read.table(paste0(DE_results, "/results_DE_T2mild.tsv"), sep="\t", header=TRUE)
miRNA_binary_highconf <- read.table("/path/to/results/09-miRNA_and_tiny_figures/miRNA_binary_table_highconf.tsv",
                                    sep="\t", header = TRUE)


# Filter DE sequences
padj.cutoff <- 0.05
DE_T1_acute_sig <- DE_T1_acute %>%
  dplyr::filter(padj < padj.cutoff)
DE_T1_mild_sig <- DE_T1_mild %>%
  dplyr::filter(padj < padj.cutoff)
DE_T2_acute_sig <- DE_T2_acute %>%
  dplyr::filter(padj < padj.cutoff)
DE_T2_mild_sig <- DE_T2_mild %>%
  dplyr::filter(padj < padj.cutoff)

# Annotate the seqs
DE_T1_acute_sig_ann <- merge(DE_T1_acute_sig, anno_df, by.x = "seq", by.y = "Seq", all.x = TRUE)
DE_T1_mild_sig_ann <-  merge(DE_T1_mild_sig, anno_df, by.x = "seq", by.y = "Seq", all.x = TRUE)
DE_T2_acute_sig_ann <- merge(DE_T2_acute_sig, anno_df, by.x = "seq", by.y = "Seq", all.x = TRUE)
DE_T2_mild_sig_ann <-  merge(DE_T2_mild_sig, anno_df, by.x = "seq", by.y = "Seq", all.x = TRUE)

# Create the family column
# T1_acute
DE_T1_acute_miR <- DE_T1_acute_sig_ann %>%
  filter(grepl("miRNA$", RNA_class))
DE_T1_acute_miR$Family <- DE_T1_acute_miR$RNA_ID

DE_T1_acute_miR$Family <- str_remove(DE_T1_acute_miR$Family, "-[35]p")
DE_T1_acute_miR$Family <- str_remove(DE_T1_acute_miR$Family, "\\(.+\\)")
DE_T1_acute_miR$Family <- str_remove(DE_T1_acute_miR$Family, "\\-[1-9]$")
DE_T1_acute_miR$Family <- str_remove(DE_T1_acute_miR$Family, "[a-z]{1,}$")
DE_T1_acute_miR$Family <- paste0("hsa-", DE_T1_acute_miR$Family)

# T1_mild
DE_T1_mild_miR <- DE_T1_mild_sig_ann %>%
  filter(grepl("miRNA$", RNA_class))
DE_T1_mild_miR$Family <- DE_T1_mild_miR$RNA_ID

DE_T1_mild_miR$Family <- str_remove(DE_T1_mild_miR$Family, "-[35]p")
DE_T1_mild_miR$Family <- str_remove(DE_T1_mild_miR$Family, "\\(.+\\)")
DE_T1_mild_miR$Family <- str_remove(DE_T1_mild_miR$Family, "\\-[1-9]$")
DE_T1_mild_miR$Family <- str_remove(DE_T1_mild_miR$Family, "[a-z]{1,}$")
DE_T1_mild_miR$Family <- paste0("hsa-", DE_T1_mild_miR$Family)

# T2 acute
DE_T2_acute_miR <- DE_T2_acute_sig_ann %>%
  filter(grepl("miRNA$", RNA_class))
DE_T2_acute_miR$Family <- DE_T2_acute_miR$RNA_ID

DE_T2_acute_miR$Family <- str_remove(DE_T2_acute_miR$Family, "-[35]p")
DE_T2_acute_miR$Family <- str_remove(DE_T2_acute_miR$Family, "\\(.+\\)")
DE_T2_acute_miR$Family <- str_remove(DE_T2_acute_miR$Family, "\\-[1-9]$")
DE_T2_acute_miR$Family <- str_remove(DE_T2_acute_miR$Family, "[a-z]{1,}$")
DE_T2_acute_miR$Family <- paste0("hsa-", DE_T2_acute_miR$Family)

# T2 mild
DE_T2_mild_miR <- DE_T2_mild_sig_ann %>%
  filter(grepl("miRNA$", RNA_class))
DE_T2_mild_miR$Family <- DE_T2_mild_miR$RNA_ID

DE_T2_mild_miR$Family <- str_remove(DE_T2_mild_miR$Family, "-[35]p")
DE_T2_mild_miR$Family <- str_remove(DE_T2_mild_miR$Family, "\\(.+\\)")
DE_T2_mild_miR$Family <- str_remove(DE_T2_mild_miR$Family, "\\-[1-9]$")
DE_T2_mild_miR$Family <- str_remove(DE_T2_mild_miR$Family, "[a-z]{1,}$")
DE_T2_mild_miR$Family <- paste0("hsa-", DE_T2_mild_miR$Family)

# Get only high-confidence miRNAs
DE_T1_acute_miR_highconf <- DE_T1_acute_miR[!grepl("\\(", DE_T1_acute_miR$RNA_ID),]
DE_T2_acute_miR_highconf <- DE_T2_acute_miR[!grepl("\\(", DE_T2_acute_miR$RNA_ID),]
DE_T1_mild_miR_highconf <- DE_T1_mild_miR[!grepl("\\(", DE_T1_mild_miR$RNA_ID),]
DE_T2_mild_miR_highconf <- DE_T2_mild_miR[!grepl("\\(", DE_T2_mild_miR$RNA_ID),]

#############################################################
# Get the table with only the consistent miRNAs

miRNA_binary_highconf_filt <- miRNA_binary_highconf[apply(miRNA_binary_highconf[, -1], 1, function(x) all(x == 1)),]
miRNA_consist_fam <- miRNA_binary_highconf_filt$Family
#Remove hsa-miR-203
miRNA_consist_fam <- miRNA_consist_fam[-which(miRNA_consist_fam == "hsa-miR-203")]

for (time in c("T1", "T2")){
  for (sev in c("acute", "mild")){
    df_cont <- get(paste0("DE_", time,"_", sev, "_miR_highconf"))
    #df_cont_filt <- df_cont[, c("seq", "lfcShrunk", "Family")]
    df_cont_filt <- df_cont %>% filter(Family %in% miRNA_consist_fam)
    df_cont_filt$Contrast <- paste0(time, "_", sev)
    df_cont_filt$Severity <- sev
    df_cont_filt$Time <- time
    name_df <- paste0("DE_", time, "_", sev, "_filt")
    assign(name_df, df_cont_filt)
  }
}

#############################################################
# Intersect MTI and miRNAs
#############################################################

# We will study acute and mild patients miRNAs separatedly

# Intersect the common IDs in T1 an T2 
# Acute patients
miRNA_acute <- data.frame("miRNA_ID" = NA)

for (fam in miRNA_consist_fam){
  temp_T1 <- sort(unique(DE_T1_acute_filt[DE_T1_acute_filt$Family == fam, "RNA_ID"]))
  temp_T2 <- sort(unique(DE_T2_acute_filt[DE_T2_acute_filt$Family == fam, "RNA_ID"]))
  temp_int <- intersect(temp_T1, temp_T2)
  temp_df <- as.data.frame(temp_int)
  colnames(temp_df) <- colnames(miRNA_acute)
  miRNA_acute <- rbind(miRNA_acute, temp_df)
}

# Remove first row with NAs
miRNA_acute <- miRNA_acute[!is.na(miRNA_acute$miRNA_ID),, drop = FALSE]
# Add hsa- to the name
miRNA_acute$miRNA_ID <- paste0("hsa-", miRNA_acute$miRNA_ID)

# Mild patients
miRNA_mild <- data.frame("miRNA_ID" = NA)

for (fam in miRNA_consist_fam){
  temp_T1 <- sort(unique(DE_T1_mild_filt[DE_T1_mild_filt$Family == fam, "RNA_ID"]))
  temp_T2 <- sort(unique(DE_T2_mild_filt[DE_T2_mild_filt$Family == fam, "RNA_ID"]))
  temp_int <- intersect(temp_T1, temp_T2)
  temp_df <- as.data.frame(temp_int)
  colnames(temp_df) <- colnames(miRNA_mild)
  miRNA_mild <- rbind(miRNA_mild, temp_df)
}

# Remove first row with NAs
miRNA_mild <- miRNA_mild[!is.na(miRNA_mild$miRNA_ID),, drop = FALSE]
# Add hsa- to the name
miRNA_mild$miRNA_ID <- paste0("hsa-", miRNA_mild$miRNA_ID)


#############################################################
# Intersect miRNAs and MTIs

# There are some mismatches between the annotations
miRNA_acute[miRNA_acute$miRNA_ID == "hsa-let-7f-1-5p", "miRNA_ID"] <- "hsa-let-7f-5p"
miRNA_mild[miRNA_mild$miRNA_ID == "hsa-miR-16-1-5p", "miRNA_ID"] <- "hsa-miR-16-5p"

# Merge
miRNA_acute_MTI <- merge(miRNA_acute, MTI_df_final, by.x = "miRNA_ID", by.y = "miRNA", all.x = TRUE)
miRNA_mild_MTI <- merge(miRNA_mild, MTI_df_final, by.x = "miRNA_ID", by.y = "miRNA", all.x = TRUE)

# Get the list of genes
acute_genes <- unique(miRNA_acute_MTI[!is.na(miRNA_acute_MTI$`Target Gene`), "Target Gene"])
mild_genes <- unique(miRNA_mild_MTI[!is.na(miRNA_mild_MTI$`Target Gene`), "Target Gene"])
mild_genes <- str_to_upper(mild_genes)

#############################################################
# Supplemental table

# Severe
# Add the condition
miRNA_acute_MTI$Condition <- "Severe"
# Remove the columns used for the selection of valid MTIs
miRNA_acute_merged <- miRNA_acute_MTI
miRNA_acute_merged[, c("Reporter", "Western", "qPCR", "Sum_exp")] <- list(NULL)
# Reorder the columns
miRNA_acute_merged <- miRNA_acute_merged[, c(10,1,2,3,4,5,6,7,8,9)]

# Moderate
# Add the condition
miRNA_mild_MTI$Condition <- "Moderate"
# Remove the columns used for the selection of valid MTIs
miRNA_mild_merged <- miRNA_mild_MTI
miRNA_mild_merged[, c("Reporter", "Western", "qPCR", "Sum_exp")] <- list(NULL)
# Reorder the columns
miRNA_mild_merged <- miRNA_mild_merged[, c(10,1,2,3,4,5,6,7,8,9)]

# Concatenate both df
miRNA_targets_supp <- rbind(miRNA_acute_merged, miRNA_mild_merged)
# Write the supplemental table
write.table(miRNA_targets_supp, paste0(output_dir, "/miRNA_targets_supp.tsv"),
            sep="\t", row.names = FALSE)


#############################################################
# KEGG enrichment
#############################################################

# Acute
acute_kegg <- enrichKEGG(gene = acute_genes_entrez$ENTREZID,
                         organism = 'hsa',
                         pvalueCutoff = 0.05)

acute_kegg_dotplot <- dotplot(acute_kegg, title = "Severe", font.size = 10)
ggsave(paste0(output_dir, "/acute_kegg_plot.png"))

# Mild
mild_kegg <- enrichKEGG(gene = mild_genes_entrez$ENTREZID,
                        organism = 'hsa',
                        pvalueCutoff = 0.05)

mild_kegg_dotplot <- dotplot(mild_kegg, title = "Moderate", font.size = 10)
ggsave(paste0(output_dir, "/mild_kegg_plot.png"))


# Arrange individual plots together
both_keggs <- ggarrange(acute_kegg_dotplot, mild_kegg_dotplot)
ggsave(paste0(output_dir, "/acute_and_mild_individual_kegg_plot.png"), height = 5.5, width = 12.5)

