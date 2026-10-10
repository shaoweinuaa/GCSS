library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)
library(ggthemes)
library(SingleCellExperiment)
library(SpatialExperiment)
library(tidyverse)
library(dittoSeq)

set.seed(100)

ZScore <- function(x) {
  value_sd <- sd(x, na.rm = TRUE)
  if (!is.finite(value_sd) || value_sd == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / value_sd
}

spe_file <- "./Output/1_Major/major_cluster_spe.rds"
annotation_file <- "./Input/annotation.csv"
output_dir <- "./Output/2_CellType"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

spe <- readRDS(spe_file)
annotation_table <- read.csv(annotation_file, check.names = FALSE)

required_columns <- c("pg_cluster", "celltype")
duplicated_clusters <- annotation_table$pg_cluster[duplicated(annotation_table$pg_cluster)]

if (length(duplicated_clusters) > 0) {
  stop("Duplicated pg_cluster values were found in annotation.csv: ", paste(unique(duplicated_clusters), collapse = ", "))
}

cluster_index <- match(as.character(spe$pg_cluster), as.character(annotation_table$pg_cluster))

if (any(is.na(cluster_index))) {
  missing_clusters <- unique(as.character(spe$pg_cluster)[is.na(cluster_index)])
  stop("The following clusters are missing from annotation.csv: ", paste(missing_clusters, collapse = ", "))
}

spe$celltype <- annotation_table$celltype[cluster_index]
spe <- spe[, spe$celltype != "rm"]

celltype_levels <- annotation_table %>% filter(celltype != "rm") %>% pull(celltype) %>% unique()
spe$celltype <- factor(spe$celltype, levels = celltype_levels)

tableau_colors <- ggthemes_data$tableau$`color-palettes`$regular$`Tableau 20`$value
miller_stone_colors <- ggthemes_data$tableau$`color-palettes`$regular$`Miller Stone`$value

celltype_colors <- colorRampPalette(c(tableau_colors, rev(miller_stone_colors)))(length(celltype_levels))
names(celltype_colors) <- celltype_levels

celltype_umap_plot <- dittoDimPlot(spe, var = "celltype", reduction.use = "UMAP_harmony", size = 0.2, do.label = TRUE, labels.repel = TRUE, labels.highlight = FALSE) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  scale_color_manual(values = celltype_colors) +
  labs(x = "UMAP1", y = "UMAP2", title = "")

ggsave(file.path(output_dir, "celltype_umap.png"), celltype_umap_plot, height = 4, width = 7.5, dpi = 300)


celltype_median_expression <- as.data.frame(
  t(assay(spe, "exprs"))[, rowData(spe)$use_channel, drop = FALSE],
  check.names = FALSE
)

celltype_median_expression$celltype <- spe$celltype

celltype_median_expression <- celltype_median_expression %>%
  group_by(celltype) %>%
  summarise(across(everything(), median), .groups = "drop") %>%
  column_to_rownames("celltype") %>%
  as.matrix()

celltype_median_expression <- celltype_median_expression[celltype_levels, , drop = FALSE]

heatmap_marker_order <- c(
  "E_cadherin",
  "PAN.Keratin",
  "CD3",
  "CD4",
  "CD8",
  "CD20",
  "CD7",
  "CD57",
  "CD45",
  "CD68",
  "CD15",
  "CD14",
  "CD16",
  "CD31",
  "aSMA",
  "Collagen1",
  "VIMENTIN"
)

missing_heatmap_markers <- setdiff(heatmap_marker_order, colnames(celltype_median_expression))

if (length(missing_heatmap_markers) > 0) {
  stop("The following heatmap markers are missing from exprs: ", paste(missing_heatmap_markers, collapse = ", "))
}

celltype_median_expression <- celltype_median_expression[, heatmap_marker_order, drop = FALSE]

heatmap_colors <- colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(101)

celltype_counts <- table(factor(spe$celltype, levels = celltype_levels))

celltype_row_annotation <- rowAnnotation(
  counts = anno_barplot(as.numeric(celltype_counts), gp = grid::gpar(fill = celltype_colors[celltype_levels]))
)

celltype_heatmap_unscaled <- Heatmap(celltype_median_expression,
                                     col = heatmap_colors,
                                     name = "Expression value",
                                     cluster_columns = FALSE,
                                     cluster_rows = FALSE,
                                     right_annotation = celltype_row_annotation,
                                     heatmap_legend_param = list(legend_height = grid::unit(4, "cm"),
                                                                 title_position = "lefttop-rot"))

png(file.path(output_dir, "celltype_heatmap_unscaled.png"), width = 12000, height = 6000, res = 72 * 15)
draw(celltype_heatmap_unscaled)
dev.off()


scaled_celltype_expression <- apply(celltype_median_expression, 2, ZScore)

celltype_heatmap_scaled <- Heatmap(scaled_celltype_expression,
                                   col = heatmap_colors, name = "Z score", 
                                   cluster_columns = FALSE, 
                                   cluster_rows = FALSE, 
                                   right_annotation = celltype_row_annotation,
                                   heatmap_legend_param = list(legend_height = grid::unit(4, "cm"),
                                                               title_position = "lefttop-rot"))

png(file.path(output_dir, "celltype_heatmap_scaled.png"), width = 12000, height = 5000, res = 72 * 15)
draw(celltype_heatmap_scaled)
dev.off()


saveRDS(spe, file.path(output_dir, "celltype_spe.rds"))