available_idents <- Idents(query)
print(available_idents)

Idents(query) <- query$scpred_prediction
de_markers <- FindMarkers(query, ident.1 = "DS3D", ident.2 = "EP3D", thresh.use=0, min.pct=0)
de_markers_sole <- de_markers

Idents(query) <- query$scpred_prediction
de_markers_sole <- FindMarkers(query, ident.1 = "DS3D", ident.2 = "EP3D")
de_markers_tabib <- de_markers_sole

filtered_de_markers_sole <- de_markers_sole %>%
  filter(p_val_adj <= 0.05)

filtered_de_markers_tabib <- de_markers_tabib %>%
  filter(p_val_adj <= 0.05)

write.csv(de_markers, file = "top_DEGs_predictedEPvsDS_sole.csv")

# Filter rows where p_val_adj is less than or equal to 0.05
filtered_de_markers <- de_markers %>%
  filter(p_val_adj <= 0.05)

# Optional: Print the filtered data frame
print(filtered_de_markers)

ggplot(de_markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = ifelse(p_val_adj > 0.05, "grey", ifelse(avg_log2FC > 0, "red", "black"))), alpha = 0.7) +
  scale_color_identity() +  # Use identity scale to remove legend
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "LogFC",
    y = "-log10(p_val_adj)"
  ) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "red")






  
# Check if APOE is present in the gene names of the Seurat object
if ("CDKN1A" %in% rownames(query@assays$RNA)) {
  
  # Extract APOE expression data
  apo_expression <- query@assays$RNA["CDKN1A", ]
  
  # Combine APOE expression data with metadata
  plot_data <- data.frame(Group = query@meta.data$scpred_prediction, APOE = apo_expression)
  
  # Create overlaid density plot for APOE expression with ggplot2
  ggplot(plot_data, aes(x = APOE, fill = as.factor(Group))) +
    geom_density(alpha = 0.5) +
    labs(title = "APOE Expression Density Comparison",
         x = "APOE Expression",
         y = "Density",
         fill = "Group") +
    scale_fill_manual(values = c("blue", "red", "green")) +  # Adjust colors based on unique values
    theme_minimal()
  
} else {
  print("APOE not found in the gene names of the Seurat object.")
}




# Check if APOE is present in the gene names of the Seurat object
if ("CDKN1A" %in% rownames(query@assays$RNA)) {
  
  # Extract APOE expression data
  apo_expression <- query@assays$RNA["CDKN1A", ]
  
  # Combine APOE expression data with metadata
  plot_data <- data.frame(Group = query@meta.data$scpred_prediction, APOE = apo_expression)
  
  # Function to estimate PDF using kernel density estimation
  pdf_estimate <- function(x) {
    density(x, adjust = 2)$y
  }
  
  # Create a line plot for estimated PDF of APOE expression with ggplot2
  ggplot(plot_data, aes(x = APOE, color = as.factor(Group))) +
    stat_function(fun = pdf_estimate, geom = "line", linewidth = 1) +
    labs(title = "APOE Expression PDF Comparison",
         x = "APOE Expression",
         y = "Probability Density",
         color = "Group") +
    scale_color_manual(values = c("blue", "red", "green")) +  # Adjust colors based on unique values
    theme_minimal()
  
} else {
  print("APOE not found in the gene names of the Seurat object.")
}


# Subset data for the gene of interest
cdkn1a_data <- query@assays$RNA_Seq["CDKN1A", , drop = FALSE]

# Create a data frame for plotting
plot_data <- data.frame(
  Expression = cdkn1a_data,
  Group = query$scpred_prediction
)

# Create a violin plot
ggplot(plot_data, aes(x = Group, y = Expression)) +
  geom_violin(fill = "lightblue", color = "blue", alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", color = "blue") +
  labs(
    title = "Violin Plot of CDKN1A Expression",
    x = "scpred_prediction",
    y = "CDKN1A Expression"
  ) +
  theme_minimal()


combined <- read.csv("/Users/bethanyhughes/Documents/SenPred_HDF/combined.csv")
combined_unique <- distinct(combined)


result <- de_markers_sole %>%
  filter(row.names(de_markers_sole) %in% combined_unique$GENE)



# Assuming 'query' is your Seurat object
# Assuming 'combined_unique' is your dataframe with a column 'GENE'

# Extract the genes of interest
genes_of_interest <- combined_unique$GENE

# Filter genes in Seurat object based on 'combined_unique' genes
query_subset <- subset(query, features = genes_of_interest)

# Calculate average expression for each gene across clusters (RNA assay)
avg_expr_all_genes <- AverageExpression(query_subset, assay = "RNA")

# Access the RNA data from the list
rna_data <- avg_expr_all_genes$RNA

# Transpose the matrix
rna_data <- t(rna_data)

# Create a SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(RNA = rna_data), colData = query_subset@meta.data)

# Calculate row means for each gene
avg_expr_df <- rowData(sce)
avg_expr_df$AverageExpression <- rowMeans(assay(sce, "RNA"))

# Add scpred_prediction information
avg_expr_df$scpred_prediction <- colData(sce)$scpred_prediction

# Subset the results for genes in 'combined_unique$GENE'
avg_expr_filtered <- avg_expr_df[avg_expr_df$gene_id %in% genes_of_interest, ]

# Print the resulting average expression values for the filtered genes
print(avg_expr_filtered)
