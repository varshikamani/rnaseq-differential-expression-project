# 01_load_data.R: load GSE68086 TEP counts, select breast cancer vs healthy donors, filter, build DESeq2 object, VST + PCA
# Put GSE68086_TEP_data_matrix.txt (from GEO) in your working directory first.

counts <- read.table("GSE68086_TEP_data_matrix.txt", header = TRUE, row.names = 1, sep = "\t")
dim(counts)
head(colnames(counts))
unique(colnames(counts))
sample_types <- sub("^X?[0-9]*\\.?", "", colnames(counts))
sample_types <- sub("\\..*", "", sample_types)
table(sample_types)

# Find columns containing "Breast" anywhere in the name
breast_cols <- grep("Breast", colnames(counts), value = TRUE)
length(breast_cols)

# Find columns containing "HD" anywhere in the name  
healthy_cols <- grep("HD", colnames(counts), value = TRUE)
length(healthy_cols)

# Combine both groups into one list of column names
selected_cols <- c(breast_cols, healthy_cols)

# Subset the count matrix to only these samples
counts_subset <- counts[, selected_cols]
dim(counts_subset)

# Build the metadata: a label for each sample matching the column order
sample_group <- c(rep("breast_cancer", length(breast_cols)), 
                  rep("healthy", length(healthy_cols)))

metadata <- data.frame(
  sample = selected_cols,
  group = sample_group
)

metadata

all(colnames(counts_subset) == metadata$sample)

# Keep only genes with at least 10 total counts across all samples
keep <- rowSums(counts_subset) >= 10
counts_filtered <- counts_subset[keep, ]

dim(counts_subset)    # before filtering
dim(counts_filtered)  # after filtering

metadata$group <- factor(metadata$group)
levels(metadata$group)

library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ group
)

vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "group")