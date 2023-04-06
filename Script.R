############################################
### Import libraries
############################################
library(DESeq2)
library(pheatmap)
library(biomaRt)
library(EnhancedVolcano)
library(Gviz)
library(rtracklayer)
library(GenomicFeatures)

############################################
### Import and process data
############################################
metadata <- read.csv("samples.csv")
metadata$type <- relevel(as.factor(metadata$type), ref = 'stemcell')
direct <- paste(metadata$experiment, "/", metadata$filename, sep = "")
data <- lapply(direct, read.table, header = TRUE)
names(data) <- metadata$sampleID

counts <- cbind(data$stemcell1$expected_count,
                data$stemcell2$expected_count,
                data$endodermal1$expected_count,
                data$endodermal2$expected_count,
                data$mesodermal1$expected_count,
                data$mesodermal2$expected_count,
                data$ectodermal1$expected_count,
                data$ectodermal2$expected_count)

counts.int <- apply(counts, 2, as.integer)
rownames(counts.int) <- data$endodermal1$gene_id
counts.int.filter <- counts.int[rowSums(counts.int) > 80,]

############################################
### Run DESeq2 and store results
############################################
dds <- DESeqDataSetFromMatrix(countData = counts.int.filter,
                              colData = metadata,
                              design = ~ type)

dds <- DESeq(dds)
res <- results(dds)
resultsNames(dds)

res.ecto <- results(dds, name = "type_ectodermal_vs_stemcell")
df <- as.data.frame(res.ecto)
df$gene_id <- rownames(df)
write.csv(df, file = paste('Results/ecto.csv', sep = ''), row.names = F)

res.endo <- results(dds, name = "type_endodermal_vs_stemcell")
df <- as.data.frame(res.endo)
df$gene_id <- rownames(df)
write.csv(df, file = paste('Results/endo.csv', sep = ''), row.names = F)

res.meso <- results(dds, name = "type_mesodermal_vs_stemcell")
df <- as.data.frame(res.meso)
df$gene_id <- rownames(df)
write.csv(df, file = paste('Results/meso.csv', sep = ''), row.names = F)
  
############################################
### Query genes using biomaRt
############################################
mart <- useEnsembl("ensembl","hsapiens_gene_ensembl")

gene_id_ecto <- sapply(rownames(res.ecto), function(x) strsplit(x, "[.]")[[1]][1])
info_ecto <- getBM(c("ensembl_gene_id","hgnc_symbol", 'chromosome_name', 'start_position', 'end_position', 'strand'),
           "ensembl_gene_id",
           gene_id_ecto,
           mart)
write.csv(info_ecto, file = paste('Gene_info/info_ecto.csv', sep = ''), row.names = F)

gene_id_endo <- sapply(rownames(res.endo), function(x) strsplit(x, "[.]")[[1]][1])
info_endo <- getBM(c("ensembl_gene_id","hgnc_symbol", 'chromosome_name', 'start_position', 'end_position', 'strand'),
                   "ensembl_gene_id",
                   gene_id_endo,
                   mart)
write.csv(info_endo, file = paste('Gene_info/info_endo.csv', sep = ''), row.names = F)

gene_id_meso <- sapply(rownames(res.meso), function(x) strsplit(x, "[.]")[[1]][1])
info_meso <- getBM(c("ensembl_gene_id","hgnc_symbol", 'chromosome_name', 'start_position', 'end_position', 'strand'),
                   "ensembl_gene_id",
                   gene_id_meso,
                   mart)
write.csv(info_meso, file = paste('Gene_info/info_meso.csv', sep = ''), row.names = F)

############################################
### Save counts for heatmap
############################################
colnames(counts.int.filter) <- metadata$sampleID
write.csv(counts.int.filter, file = paste('counts.csv', sep = ''), row.names = rownames(counts.int.filter))
