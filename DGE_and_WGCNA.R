setwd('D:\\Internship Draft\\Maki')
load(".RData")

if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  install.packages("AnnotationDbi")
}
if (!requireNamespace("hgu133plus2.db", quietly = TRUE)) {
  BiocManager::install("hgu133plus2.db")
}

install.packages("devtools")
devtools::install_github("kevinblighe/CorLevelPlot",force = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
BiocManager::install("KEGGREST")
BiocManager::install("clusterProfiler",force = TRUE)
BiocManager::install("arrayQualityMetrics")
BiocManager::install("VennDiagram")

# Install additional required packages
install.packages("doParallel")
install.packages("foreach")
install.packages("splitstackshape")
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  BiocManager::install("WGCNA")
}
library(WGCNA)
library(doParallel)
library(tibble)
library(foreach)
library(CorLevelPlot)
library(AnnotationDbi)
library(hgu133plus2.db)
library(KEGGREST)
library(KEGG.db)
library(affy)
library(oligo)
library(Biobase)
library(GEOquery)
library(dplyr)
library(clusterProfiler)
library(arrayQualityMetrics)
library(splitstackshape)
library(tidyr)
library(limma)
library('org.Hs.eg.db')
library(VennDiagram)
library(WGCNA)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(hgu133plus2.db) 

# Data preprocessing
untar("GSE4183_RAW.tar")
celFiles <- list.celfiles()
affyRAW <- read.celfiles(celFiles)
hist(affyRAW)
eset <- oligo::rma(affyRAW)
hist(eset)
write.exprs(eset,file="data1.txt")

# Annotation
mydata <- read.delim("GSE4183_family.soft", check.names = FALSE)
abc <- data.frame(mydata)
a <- cSplit(abc,"Gene.Symbol", "//")
write.table(a, file="imp1.txt",sep = "\t", row.names=FALSE,quote = FALSE)

data1 <- read.delim("imp1.txt",check.names = FALSE)
data2 <- read.delim("data1.txt", check.names=FALSE)

combined <- left_join(data1, data2, by ="ID")
write.csv(combined,"annotated.csv")

# Differential Expression Analysis
pData(eset)
targets <- readTargets("GSE4183_target_file.csv", sep = ",")
condition <- factor(targets$Condition, levels = c("colon_normal", "colon_adenoma", "colon_CRC", "colon_IBD"))
design <- model.matrix(~ 0 + condition)
colnames(design) <- c("normal", "adenoma", "CRC", "IBD")
print(design)

# Fit linear model
fit <- lmFit(eset, design)

# Differential expression analysis for adenoma vs normal
contrast.matrix_adenoma <- makeContrasts(adenoma - normal, levels = design)
fit_adenoma <- contrasts.fit(fit, contrast.matrix_adenoma)
fit_adenoma <- eBayes(fit_adenoma)
res_adenoma_vs_normal <- topTable(fit_adenoma, number=Inf)
write.table(res_adenoma_vs_normal, "adenoma_vs_normal.txt", sep="\t")

# Differential expression analysis for CRC vs normal
contrast.matrix_CRC <- makeContrasts(CRC - normal, levels = design)
fit_CRC <- contrasts.fit(fit, contrast.matrix_CRC)
fit_CRC <- eBayes(fit_CRC)
res_CRC_vs_normal <- topTable(fit_CRC, number=Inf)
write.table(res_CRC_vs_normal, "CRC_vs_normal.txt", sep="\t")

# Differential expression analysis for IBD vs normal
contrast.matrix_IBD <- makeContrasts(IBD - normal, levels = design)
fit_IBD <- contrasts.fit(fit, contrast.matrix_IBD)
fit_IBD <- eBayes(fit_IBD)
res_IBD_vs_normal <- topTable(fit_IBD, number=Inf)
write.table(res_IBD_vs_normal, "IBD_vs_normal.txt", sep="\t")

# Find common DE genes
# Load DE results
res1 <- read.delim("IBD_vs_normal.txt", check.names = FALSE)
res2 <- read.delim("adenoma_vs_normal.txt", check.names = FALSE)
res3 <- read.delim("CRC_vs_normal.txt", check.names = FALSE)

# Define significance threshold
significance_threshold <- 0.05
logFC_threshold <- 1

# Filter DE genes for each dataset
de_genes1 <- rownames(res1[res1$P.Value < significance_threshold & abs(res1$logFC) > logFC_threshold, ])
de_genes2 <- rownames(res2[res2$P.Value < significance_threshold & abs(res2$logFC) > logFC_threshold, ])
de_genes3 <- rownames(res3[res3$P.Value < significance_threshold & abs(res3$logFC) > logFC_threshold, ])

# Convert to character vectors
de_genes1 <- as.character(de_genes1)
de_genes2 <- as.character(de_genes2)
de_genes3 <- as.character(de_genes3)

# Find common DE genes
common_de_genes <- Reduce(intersect, list(de_genes1, de_genes2, de_genes3))

# Print common DE genes
print("Common Differentially Expressed Genes:")
print(common_de_genes)

# Pathway Enrichment Analysis
# Annotate DEGs
res_adenoma_vs_normal <- read.delim("adenoma_vs_normal.txt", check.names = FALSE)
res_CRC_vs_normal <- read.delim("CRC_vs_normal.txt", check.names = FALSE)
res_IBD_vs_normal <- read.delim("IBD_vs_normal.txt", check.names = FALSE)

# Annotate DEGs according to logFC and adjusted P-value
res_adenoma_vs_normal <- res_adenoma_vs_normal %>%
  mutate(diffexpressed = case_when(
    logFC > 0 & adj.P.Val < 0.05 ~ 'UP',
    logFC < 0 & adj.P.Val < 0.05 ~ 'DOWN',
    adj.P.Val > 0.05 ~ 'NO'
  ))

res_CRC_vs_normal <- res_CRC_vs_normal %>%
  mutate(diffexpressed = case_when(
    logFC > 0 & adj.P.Val < 0.05 ~ 'UP',
    logFC < 0 & adj.P.Val < 0.05 ~ 'DOWN',
    adj.P.Val > 0.05 ~ 'NO'
  ))

res_IBD_vs_normal <- res_IBD_vs_normal %>%
  mutate(diffexpressed = case_when(
    logFC > 0 & adj.P.Val < 0.05 ~ 'UP',
    logFC < 0 & adj.P.Val < 0.05 ~ 'DOWN',
    adj.P.Val > 0.05 ~ 'NO'
  ))

# Save annotated files
write.csv(res_adenoma_vs_normal, "adenoma_vs_normal_annotated.csv")
write.csv(res_CRC_vs_normal, "CRC_vs_normal_annotated.csv")
write.csv(res_IBD_vs_normal, "IBD_vs_normal_annotated.csv")

# Map probe IDs to gene symbols
probe_ids_adenoma <- rownames(res_adenoma_vs_normal)
gene_symbols_adenoma <- mapIds(hgu133plus2.db, keys = probe_ids_adenoma, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
res_adenoma_vs_normal$gene_symbol <- gene_symbols_adenoma

probe_ids_crc <- rownames(res_CRC_vs_normal)
gene_symbols_crc <- mapIds(hgu133plus2.db, keys = probe_ids_crc, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
res_CRC_vs_normal$gene_symbol <- gene_symbols_crc

probe_ids <- rownames(res_IBD_vs_normal)
gene_symbols_ibd <- mapIds(hgu133plus2.db, keys = probe_ids, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
res_IBD_vs_normal$gene_symbol <- gene_symbols_ibd


# Convert probe IDs to gene symbols for common DE genes
gene_symbols <- mapIds(hgu133plus2.db, keys = common_de_genes, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
######################################################################
# PLOTS FOR DE ANALYSIS


# Function to create volcano plot
create_volcano_plot <- function(res_data, title) {
  # Add a column for significant + fold change status
  res_data <- res_data %>%
    mutate(
      sig = case_when(
        adj.P.Val < 0.05 & abs(logFC) > 1 ~ "Significant",
        TRUE ~ "Not Significant"
      ),
      label = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, gene_symbol, "")
    )
  
  # Create volcano plot
  ggplot(res_data, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = sig), alpha = 0.6) +
    geom_text(aes(label = label), size = 2, vjust = -0.5) +
    scale_color_manual(values = c("Not Significant" = "grey", 
                                  "Significant" = "red")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    theme_minimal() +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10 Adjusted P-value"
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    )
}

# Function to create heatmap with gene symbols
create_heatmap <- function(eset, condition, top_n_genes = 50, title) {
  # Get expression data
  expr_data <- exprs(eset)
  
  # Calculate variance for each gene
  var_genes <- apply(expr_data, 1, var)
  
  # Select top variable genes
  top_var_genes <- names(sort(var_genes, decreasing = TRUE)[1:top_n_genes])
  expr_data_subset <- expr_data[top_var_genes, ]
  
  # Convert probe IDs to gene symbols
  gene_symbols <- mapIds(hgu133plus2.db, 
                         keys = rownames(expr_data_subset), 
                         column = "SYMBOL", 
                         keytype = "PROBEID", 
                         multiVals = "first")
  
  # Replace any NA symbols with original probe IDs
  gene_symbols[is.na(gene_symbols)] <- names(gene_symbols)[is.na(gene_symbols)]
  
  # Scale the data
  expr_data_scaled <- t(scale(t(expr_data_subset)))
  
  # Assign gene symbols as row names
  rownames(expr_data_scaled) <- gene_symbols
  
  # Create annotation data frame
  annotation_col <- data.frame(
    Condition = condition,
    row.names = colnames(expr_data)
  )
  
  # Create color schemes
  ann_colors <- list(
    Condition = c(
      colon_normal = "lightblue",
      colon_adenoma = "orange",
      colon_CRC = "red",
      colon_IBD = "purple"
    )
  )
  
  # Generate heatmap
  pheatmap(
    expr_data_scaled,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_rownames = TRUE,
    show_colnames = FALSE,
    main = title,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    scale = "none",
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
    fontsize_row = 8  # Adjust this value if needed for better readability
  )
}

#volcano plots
pdf("volcano_plots.pdf", width = 10, height = 12)
par(mfrow = c(3,1))
print(create_volcano_plot(res_adenoma_vs_normal, "Adenoma vs Normal"))
print(create_volcano_plot(res_CRC_vs_normal, "CRC vs Normal"))
print(create_volcano_plot(res_IBD_vs_normal, "IBD vs Normal"))
dev.off()
#png vocalno plots
png("volcano_plots_full.png", width = 1000, height = 1500, res = 300)  # Adjust height for better spacing
par(mfrow = c(3, 1))
print(create_volcano_plot(res_adenoma_vs_normal, "Adenoma vs Normal"))
print(create_volcano_plot(res_CRC_vs_normal, "CRC vs Normal"))
print(create_volcano_plot(res_IBD_vs_normal, "IBD vs Normal"))
dev.off()
# Function to create heatmap with unique gene symbols
# Function to create heatmap with unique gene symbols
create_heatmap <- function(eset, condition, top_n_genes = 50, title) {
  # Get expression data
  expr_data <- exprs(eset)
  
  # Remove AFFX control probes
  expr_data <- expr_data[!grepl("^AFFX-", rownames(expr_data)), ]
  
  # Convert probe IDs to gene symbols
  gene_symbols <- mapIds(hgu133plus2.db, 
                         keys = rownames(expr_data), 
                         column = "SYMBOL", 
                         keytype = "PROBEID", 
                         multiVals = "first")
  
  # Create a data frame with probe IDs, gene symbols and expression data
  expr_df <- data.frame(
    probe_id = rownames(expr_data),
    gene_symbol = gene_symbols,
    expr_data,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Remove rows with NA gene symbols
  expr_df <- expr_df[!is.na(expr_df$gene_symbol), ]
  
  # Calculate row variances
  expr_df$variance <- apply(expr_df[, -c(1:2)], 1, var)
  
  # For duplicate gene symbols, keep probe with highest variance
  expr_df <- expr_df %>%
    group_by(gene_symbol) %>%
    arrange(desc(variance)) %>%
    slice(1) %>%
    ungroup()
  
  # Sort by variance to get top variable genes
  expr_df <- expr_df[order(expr_df$variance, decreasing = TRUE), ]
  
  # Select top n genes
  expr_df <- expr_df[1:min(top_n_genes, nrow(expr_df)), ]
  
  # Extract expression data for heatmap
  expr_mat <- as.matrix(expr_df[, -c(1:2, ncol(expr_df))])
  rownames(expr_mat) <- expr_df$gene_symbol
  
  # Scale the data
  expr_mat_scaled <- t(scale(t(expr_mat)))
  
  # Create annotation data frame
  annotation_col <- data.frame(
    Condition = condition,
    row.names = colnames(expr_mat)
  )
  
  # Create color schemes
  ann_colors <- list(
    Condition = c(
      colon_normal = "lightblue",
      colon_adenoma = "orange",
      colon_CRC = "red",
      colon_IBD = "purple"
    )
  )
  
  # Generate heatmap
  pheatmap(
    expr_mat_scaled,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_rownames = TRUE,
    show_colnames = FALSE,
    main = title,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    scale = "none",
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
    fontsize_row = 8
  )
}

# Generate heatmaps
pdf("heatmaps_unique_genes.pdf", width = 10, height = 12)
create_heatmap(eset, condition, top_n_genes = 50, 
               "Top 50 Variable Genes Across Conditions")

# Create condition-specific heatmaps
create_heatmap(eset[, condition %in% c("colon_normal", "colon_adenoma")],
               condition[condition %in% c("colon_normal", "colon_adenoma")],
               title = "Normal vs Adenoma")

create_heatmap(eset[, condition %in% c("colon_normal", "colon_CRC")],
               condition[condition %in% c("colon_normal", "colon_CRC")],
               title = "Normal vs CRC")

create_heatmap(eset[, condition %in% c("colon_normal", "colon_IBD")],
               condition[condition %in% c("colon_normal", "colon_IBD")],
               title = "Normal vs IBD")
dev.off()
#######################################################################
# KEGG Pathway Analysis
kegg_enrich <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = "hsa")
summary(kegg_enrich)

# GO Enrichment Analysis
go_enrich <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", 
                      ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
summary(go_enrich)

# Clean up gene symbols and prepare for visualization
gene_symbols <- na.omit(gene_symbols)
go_enrich <- enrichGO(gene = gene_symbols,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "BP",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

# Visualizations
barplot(go_enrich, showCategory = 10, title = "Top 10 GO Terms Enriched")
dotplot(go_enrich, showCategory = 10, title = "Top 10 GO Terms Enriched")
cnetplot(go_enrich, showCategory = 5, foldChange = NULL, layout = "kk", title = "GO Enrichment Map")

# KEGG visualization
entrez_ids <- na.omit(entrez_ids)
kk <- enrichKEGG(gene = entrez_ids$ENTREZID,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)

barplot(kk, showCategory = 10, title = "Top 10 KEGG Pathways Enriched")
dotplot(kk, showCategory = 10, title = "Top 10 KEGG Pathways Enriched")

# GSEA Analysis
geneList_adenoma <- res_adenoma_vs_normal$logFC
names(geneList_adenoma) <- res_adenoma_vs_normal$gene_symbol
geneList_adenoma <- sort(geneList_adenoma, decreasing = TRUE)

geneList_CRC <- res_CRC_vs_normal$logFC
names(geneList_CRC) <- res_CRC_vs_normal$gene_symbol
geneList_CRC <- sort(geneList_CRC, decreasing = TRUE)

geneList_IBD <- res_IBD_vs_normal$logFC
names(geneList_IBD) <- res_IBD_vs_normal$gene_symbol
geneList_IBD <- sort(geneList_IBD, decreasing = TRUE)

# GSEA for each comparison
gsea_GO_adenoma <- gseGO(geneList = geneList_adenoma, 
                         OrgDb = org.Hs.eg.db, 
                         keyType = "SYMBOL", 
                         ont = "BP", 
                         pvalueCutoff = 0.05, 
                         verbose = FALSE)

gsea_GO_CRC <- gseGO(geneList = geneList_CRC, 
                     OrgDb = org.Hs.eg.db, 
                     keyType = "SYMBOL", 
                     ont = "BP", 
                     pvalueCutoff = 0.05, 
                     verbose = FALSE)

gsea_GO_IBD <- gseGO(geneList = geneList_IBD, 
                     OrgDb = org.Hs.eg.db, 
                     keyType = "SYMBOL", 
                     ont = "BP", 
                     pvalueCutoff = 0.05, 
                     verbose = FALSE)

# Visualize GSEA results
dotplot(gsea_GO_adenoma, showCategory = 10, title = "GSEA GO Adenoma vs Normal")
dotplot(gsea_GO_CRC, showCategory = 10, title = "GSEA GO CRC vs Normal")
dotplot(gsea_GO_IBD, showCategory = 10, title = "GSEA GO IBD vs Normal")

# WGCNA Analysis
options(stringsAsFactors = FALSE)
expression_data <- as.data.frame(exprs(eset))
expression_data <- t(expression_data)
expression_data

gsg = goodSamplesGenes(expression_data)
summary(gsg)
gsg$allOK
table(gsg$goodGenes)

# Calculate variance for each gene
gene_variances <- apply(expression_data, 2, var)

# Filter genes based on variance
threshold <- quantile(gene_variances, 0.5)
filtered_expression_data <- expression_data[, gene_variances > threshold]

dim(filtered_expression_data)



############################################################################
#WGCNA Analysis

############################################################################

powers = c(1:20)
sft = pickSoftThreshold(
  filtered_expression_data,
  powerVector = powers,
  verbose = 5,
  networkType = "signed"
)

# 3. Then run your plotting code
# Plot Scale-Free Topology Fit
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit",
     type = "n", main = "Scale Independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.8, col = "blue", lty=2)

# Plot Mean Connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "b", pch = 19, col = "blue", main = "Mean Connectivity vs Power")


# Convert expression_data to norm.counts
norm.counts <- filtered_expression_data
norm.counts[] = sapply(norm.counts, as.numeric)  # Convert to numeric

soft_power <- 13
temp_cor <- cor
cor <- WGCNA::cor  # Temporarily reassign cor to WGCNA’s version

# Execute blockwiseModules
n_threads <- 16  # Use the full 16 threads

module_colors <- blockwiseModules(
  norm.counts, maxBlockSize = 14000, TOMType = "signed",
  power = soft_power, mergeCutHeight = 0.25,
  numericLabels = FALSE, randomSeed = 1234, verbose = 3,
  nThreads = n_threads  # Enable multithreading
)


# Look at sample relationships
sampleTree = hclust(dist(filtered_expression_data), method = "average")
plot(sampleTree)

# Check data distribution
boxplot(filtered_expression_data)


cor <- temp_cor


# Get the block assignment for each gene
blocks <- module_colors$blocks
print("Number of blocks:")
print(table(blocks))

# Plot dendrogram for each block separately
par(mfrow = c(length(unique(blocks)), 1))  # Set up plotting layout
for(block in unique(blocks)) {
  # Get genes in this block
  block_genes <- which(blocks == block)
  
  # Get corresponding colors for this block
  block_colors <- data.frame(
    unmerged = module_colors$unmergedColors[block_genes],
    merged = module_colors$colors[block_genes]
  )
  
  # Create plot title
  plot_title <- paste("Block", block, "Dendrogram and Module Colors")
  
  # Plot dendrogram for this block
  plotDendroAndColors(
    module_colors$dendrograms[[block]], 
    block_colors,
    c("unmerged", "merged"),
    main = plot_title,
    dendroLabels = FALSE,
    addGuide = TRUE,
    hang = 0.03,
    guideHang = 0.05
  )
}
par(mfrow = c(1, 1))  # Reset plotting layout

# Module Eigengenes analysis
module_eigengene <- module_colors$MEs
head(module_eigengene)

# Print gene counts for each module
print("Module sizes:")
print(table(module_colors$colors))

# Create diagnostic checks output
cat("Number of genes in dendrogram:", length(module_colors$dendrograms[[1]]$order), "\n")
cat("Number of genes in unmergedColors:", length(module_colors$unmergedColors), "\n")
cat("Number of genes in colors:", length(module_colors$colors), "\n")

# Plot final dendrogram with color matching
colors_to_plot <- data.frame(
  unmerged = module_colors$unmergedColors[module_colors$dendrograms[[1]]$order],
  merged = module_colors$colors[module_colors$dendrograms[[1]]$order]
)

plotDendroAndColors(
  module_colors$dendrograms[[1]], 
  colors_to_plot,
  c("unmerged", "merged"),
  dendroLabels = FALSE,
  addGuide = TRUE,
  hang = 0.03,
  guideHang = 0.05
)


##############################################################################3
blocks <- module_colors$blocks
print("Number of blocks:")
print(table(blocks))

# Plot dendrogram for each block in R graphics window
for (block in unique(blocks)) {
  block_genes <- which(blocks == block)
  block_colors <- data.frame(
    unmerged = module_colors$unmergedColors[block_genes],
    merged = module_colors$colors[block_genes]
  )
  
  plot_title <- paste("Block", block, "Dendrogram and Module Colors")
  plotDendroAndColors(
    module_colors$dendrograms[[block]], 
    block_colors,
    c("unmerged", "merged"),
    main = plot_title,
    dendroLabels = FALSE,
    addGuide = TRUE,
    hang = 0.03,
    guideHang = 0.05
  )
}



#############################################################################
#### Module Eigens 
module_eigengene <- module_colors$MEs


head(module_eigengene)

targets$Condition <- factor(targets$Condition,levels = c("colon_normal",  "colon_adenoma", "colon_CRC", "colon_IBD" ))
levels(targets$Condition)
disease_type = binarizeCategoricalColumns(targets$Condition,
                                          includePairwise=FALSE,
                                          includeLevelVsAll = TRUE,
                                          minCount = 1)
traits = targets
traits

traits = cbind(traits,disease_type)
traits
colnames(traits)
traits <- traits %>% select(-Condition)

# defining no. of genes and samples
nSamples = nrow(norm.counts)


nGenes = ncol(norm.counts)
module_trait_cor = cor(module_eigengene,traits,use="p")
module_trait_cor_pvals=corPvalueStudent(module_trait_cor,nSamples)

rownames(module_eigengene) <- sub("\\.CEL$", "", rownames(module_eigengene))

rownames(traits) <- rownames(module_eigengene)
traits
heatmap_data = merge(module_eigengene,traits,by="row.names")
head(heatmap_data)
dim(heatmap_data)
heatmap_data = heatmap_data %>%
  column_to_rownames(var="Row.names")
############################################################
###########################################################

CorLevelPlot(heatmap_data,
             x = names(heatmap_data)[31:33],  # Traits columns
             y = names(heatmap_data)[1:29],  # Module eigengene columns
             col = c("blue1", "skyblue", "white", "pink", "red"))



module_gene_map = as.data.frame(module_colors$colors)
blue_genes = module_gene_map %>%
  filter(`module_colors$colors`=="blue") %>%
  rownames()

head(blue_genes)
gene_symbols_blue <- mapIds(
  hgu133plus2.db, 
  keys = blue_genes, 
  column = "SYMBOL", 
  keytype = "PROBEID", 
  multiVals = "first"
)

common_de_genes_genesymbol <- mapIds(
  hgu133plus2.db, 
  keys = common_de_genes, 
  column = "SYMBOL", 
  keytype = "PROBEID", 
  multiVals = "first"
)



head(gene_symbols_blue)
head(common_de_genes_genesymbol)

# searching for common genes between DE and modules
# Convert to unique gene symbols (to remove duplicates)
unique_gene_symbols_blue <- unique(gene_symbols_blue)
unique_common_de_genes_genesymbol <- unique(common_de_genes_genesymbol)

# Find common gene symbols
common_genes_de_module <- intersect(unique_gene_symbols_blue, unique_common_de_genes_genesymbol)

# Display the common genes
common_genes_de_module

save.image("fixed.RData")



###############################################END##################
###############################################END##################
###############################################END##################
###############################################END##################
###############################################END##################
