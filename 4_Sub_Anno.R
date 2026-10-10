library(RColorBrewer)
library(ComplexHeatmap)
library(ggthemes)
library(SingleCellExperiment)
library(SpatialExperiment)
library(tidyverse)

set.seed(100)

ZScore <- function(x) {
  value_sd <- sd(x, na.rm = TRUE)
  if (!is.finite(value_sd) || value_sd == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / value_sd
}

# Run in the following order
# Lymphocyte -> Myeloid -> Epithelial

cluster_label <- "Epithelial"

subcluster_spe_file <- file.path("./Output", cluster_label, paste0(tolower(cluster_label), "_cluster_spe.rds"))
annotation_file <- file.path("./Output", cluster_label, "annotation.csv")
celltype_spe_file <- "./Output/2_CellType/celltype_spe.rds"
corrected_spe_file <- "./Output/2_CellType/celltype_corrected_spe.rds"
merged_spe_file <- "./Output/2_CellType/celltype_subtype_spe.rds"
output_dir <- file.path("./Output", cluster_label)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

spe_sub <- readRDS(subcluster_spe_file)
annotation_table <- read.csv(annotation_file, check.names = FALSE, fileEncoding = "UTF-8")

required_columns <- c("pg_cluster", "celltype")
missing_columns <- setdiff(required_columns, colnames(annotation_table))

if (length(missing_columns) > 0) {
  stop("annotation.csv is missing required columns: ", paste(missing_columns, collapse = ", "))
}

duplicated_clusters <- annotation_table$pg_cluster[duplicated(annotation_table$pg_cluster)]

if (length(duplicated_clusters) > 0) {
  stop("Duplicated pg_cluster values were found in annotation.csv: ", paste(unique(duplicated_clusters), collapse = ", "))
}

if (cluster_label == "Myeloid") annotation_table$celltype <- gsub("MF", "MΦ", annotation_table$celltype)

cluster_index <- match(as.character(spe_sub$pg_cluster), as.character(annotation_table$pg_cluster))

if (any(is.na(cluster_index))) {
  missing_clusters <- unique(as.character(spe_sub$pg_cluster)[is.na(cluster_index)])
  stop("The following clusters are missing from annotation.csv: ", paste(missing_clusters, collapse = ", "))
}

spe_sub$celltype <- annotation_table$celltype[cluster_index]

if (cluster_label == "Lymphocyte") correction_labels <- c("Epithelial", "Mesenchymal")
if (cluster_label == "Myeloid") correction_labels <- character(0)
if (cluster_label == "Epithelial") correction_labels <- "Mesenchymal"


if (length(correction_labels) > 0) {
  if (cluster_label == "Lymphocyte") spe_corrected <- readRDS(celltype_spe_file)
  if (cluster_label == "Epithelial") spe_corrected <- readRDS(corrected_spe_file)
  
  correction_index <- which(spe_sub$celltype %in% correction_labels)
  spe_correction <- spe_sub[, correction_index]
  
  cell_index <- match(colnames(spe_correction), colnames(spe_corrected))
  
  if (any(is.na(cell_index))) {
    stop("Some correction cells were not found in the full SPE object.")
  }
  
  spe_corrected$celltype <- as.character(spe_corrected$celltype)
  spe_corrected$celltype[cell_index] <- as.character(spe_correction$celltype)
  
  lineage_index <- which(spe_corrected$celltype == "Lineage-")
  if (length(lineage_index) > 0) spe_corrected$celltype[lineage_index] <- "Mesenchymal"
  
  spe_corrected$celltype <- factor(spe_corrected$celltype)
  
  saveRDS(spe_corrected, corrected_spe_file)
}


remove_labels <- c("rm", correction_labels)
spe_sub_filtered <- spe_sub[, !spe_sub$celltype %in% remove_labels]

celltype_levels <- annotation_table %>% filter(!celltype %in% remove_labels) %>% pull(celltype) %>% unique()
spe_sub_filtered$celltype <- factor(spe_sub_filtered$celltype, levels = celltype_levels)

tableau_colors <- ggthemes_data$tableau$`color-palettes`$regular$`Tableau 20`$value
miller_stone_colors <- ggthemes_data$tableau$`color-palettes`$regular$`Miller Stone`$value

celltype_colors <- colorRampPalette(c(tableau_colors, rev(miller_stone_colors)))(length(celltype_levels))
names(celltype_colors) <- celltype_levels


celltype_median_expression <- as.data.frame(
  t(assay(spe_sub_filtered, "exprs"))[, rowData(spe_sub_filtered)$use_channel, drop = FALSE],
  check.names = FALSE
)

celltype_median_expression$celltype <- spe_sub_filtered$celltype

celltype_median_expression <- celltype_median_expression %>%
  group_by(celltype) %>%
  summarise(across(everything(), median), .groups = "drop") %>%
  column_to_rownames("celltype") %>%
  as.matrix()

celltype_median_expression <- celltype_median_expression[celltype_levels, , drop = FALSE]

heatmap_colors <- colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(101)

celltype_counts <- table(factor(spe_sub_filtered$celltype, levels = celltype_levels))
celltype_row_annotation <- rowAnnotation(counts = anno_barplot(as.numeric(celltype_counts), gp = grid::gpar(fill = celltype_colors)))

celltype_heatmap_unscaled <- Heatmap(celltype_median_expression,
                                     col = heatmap_colors,
                                     name = "Expression value",
                                     cluster_columns = FALSE,
                                     cluster_rows = FALSE,
                                     right_annotation = celltype_row_annotation,
                                     heatmap_legend_param = list(legend_height = grid::unit(4, "cm"),
                                                                 title_position = "lefttop-rot"))

png(file.path(output_dir, paste0(tolower(cluster_label), "_celltype_heatmap_unscaled.png")), width = 12000, height = 5000, res = 72 * 15)
draw(celltype_heatmap_unscaled)
dev.off()


scaled_celltype_expression <- apply(celltype_median_expression, 2, ZScore)

celltype_heatmap_scaled <- Heatmap(scaled_celltype_expression,
                                   col = heatmap_colors,
                                   name = "Z score",
                                   cluster_columns = FALSE,
                                   cluster_rows = FALSE,
                                   right_annotation = celltype_row_annotation,
                                   heatmap_legend_param = list(legend_height = grid::unit(4, "cm"),
                                                               title_position = "lefttop-rot"))

png(file.path(output_dir, paste0(tolower(cluster_label), "_celltype_heatmap_scaled.png")), width = 12000, height = 5000, res = 72 * 15)
draw(celltype_heatmap_scaled)
dev.off()


saveRDS(spe_sub_filtered, file.path(output_dir, paste0(tolower(cluster_label), "_celltype_spe.rds")))


if (cluster_label == "Lymphocyte") {
  spe_merged <- readRDS(corrected_spe_file)
} else {
  spe_merged <- readRDS(merged_spe_file)
}

cell_index <- match(colnames(spe_sub), colnames(spe_merged))

if (any(is.na(cell_index))) {
  stop("Some cells in spe_sub were not found in spe_merged.")
}

spe_merged$celltype <- as.character(spe_merged$celltype)
spe_merged$celltype[cell_index] <- as.character(spe_sub$celltype)

spe_merged <- spe_merged[, spe_merged$celltype != "rm"]
spe_merged$celltype <- factor(spe_merged$celltype)

saveRDS(spe_merged, merged_spe_file)

print("Finish!")