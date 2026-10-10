library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)
library(ggthemes)
library(ggpubr)
library(SingleCellExperiment)
library(SpatialExperiment)
library(tidyverse)
library(dittoSeq)
library(scater)
library(harmony)
library(FastPG)

set.seed(100)

ZScore <- function(x) {
  value_sd <- sd(x, na.rm = TRUE)
  if (!is.finite(value_sd) || value_sd == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / value_sd
}

marker_file <- "./Input/marker.csv"
spe_file <- "./Output/2_CellType/celltype_spe.rds"

cluster_label <- "Myeloid" # "Lymphocyte" "Epithelial"

output_dir <- file.path("./Output", cluster_label)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

panel <- read.csv(marker_file, check.names = FALSE, fileEncoding = "GBK")
spe <- readRDS(spe_file)

cell_index <- which(spe$celltype %in% cluster_label)
spe_sub <- spe[, cell_index]

print(cluster_label)
print(length(cell_index))

if (cluster_label == "Myeloid") rowData(spe_sub)$use_channel <- rownames(spe_sub) %in% (panel %>% filter(myeloid_function == "1" | myeloid_lineage == "1") %>% pull(marker))
if (cluster_label == "Lymphocyte") rowData(spe_sub)$use_channel <- rownames(spe_sub) %in% (panel %>% filter(lymphoid_function == "1" | lymphoid_lineage == "1") %>% pull(marker))
if (cluster_label == "Epithelial") rowData(spe_sub)$use_channel <- rownames(spe_sub) %in% (panel %>% filter(tumor_microenvironment_function == "1" | tumor_microenvironment_lineage == "1") %>% pull(marker))

spe_sub <- runUMAP(spe_sub, subset_row = rowData(spe_sub)$use_channel, exprs_values = "exprs", name = "UMAP")

expression_matrix <- t(assay(spe_sub, "exprs"))[, rowData(spe_sub)$use_channel, drop = FALSE]

svd_result <- svd(expression_matrix)
explained_variance <- svd_result$d^2 / sum(svd_result$d^2)
n_pcs <- which(cumsum(explained_variance) >= 0.95)[1]
n_pcs <- max(n_pcs, 2)

if (length(unique(spe_sub$batch_id)) < 2) {
  pca_result <- prcomp(expression_matrix, center = TRUE, scale. = FALSE)
  harmony_embedding <- pca_result$x[, seq_len(n_pcs), drop = FALSE]
} else {
  harmony_embedding <- harmony::HarmonyMatrix(expression_matrix, as.factor(spe_sub$batch_id), do_pca = TRUE, npcs = n_pcs)
}

reducedDim(spe_sub, "harmony") <- harmony_embedding

spe_sub <- runUMAP(spe_sub, dimred = "harmony", name = "UMAP_harmony")

spe_sub$batch_id <- factor(spe_sub$batch_id, levels = unique(spe_sub$batch_id))
legend_columns <- ceiling(length(levels(spe_sub$batch_id)) / 10)

batch_before_plot <- dittoDimPlot(spe_sub,
                                  var = "batch_id",
                                  reduction.use = "UMAP",
                                  size = 0.2) +
  guides(color = guide_legend(ncol = legend_columns,
                              byrow = FALSE,
                              override.aes = list(size = 2))) +
  ggtitle("Batch ID on UMAP before correction")

batch_after_plot <- dittoDimPlot(spe_sub,
                                 var = "batch_id",
                                 reduction.use = "UMAP_harmony",
                                 size = 0.2) +
  guides(color = guide_legend(ncol = legend_columns,
                              byrow = FALSE,
                              override.aes = list(size = 2))) +
  ggtitle("Batch ID on UMAP after correction") +
  labs(x = "UMAP1", y = "UMAP2")

batch_correction_plot <- ggarrange(batch_before_plot, batch_after_plot, ncol = 2)

ggsave(file.path(output_dir, "batch_correction_umap.png"), batch_correction_plot, height = 3.5, width = legend_columns + 7, dpi = 300)


clustering_embedding <- reducedDim(spe_sub, "harmony")
cluster_result <- FastPG::fastCluster(as.matrix(clustering_embedding), k = 50, num_threads = 100)
spe_sub$pg_cluster <- factor(cluster_result$communities)

cluster_colors <- ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value
cluster_colors <- colorRampPalette(cluster_colors)(length(levels(spe_sub$pg_cluster)))
names(cluster_colors) <- levels(spe_sub$pg_cluster)


cluster_median_expression <- as.data.frame(
  t(assay(spe_sub, "exprs"))[, rowData(spe_sub)$use_channel, drop = FALSE],
  check.names = FALSE
)

cluster_median_expression$pg_cluster <- spe_sub$pg_cluster

cluster_median_expression <- cluster_median_expression %>%
  group_by(pg_cluster) %>%
  summarise(across(everything(), median), .groups = "drop") %>%
  column_to_rownames("pg_cluster") %>%
  as.matrix()

cluster_median_expression <- cluster_median_expression[levels(spe_sub$pg_cluster), , drop = FALSE]

heatmap_colors <- colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(101)

cluster_counts <- table(spe_sub$pg_cluster)
cluster_row_annotation <- rowAnnotation(counts = anno_barplot(as.numeric(cluster_counts), gp = grid::gpar(fill = cluster_colors)))

cluster_heatmap_unscaled <- Heatmap(cluster_median_expression,
                                    col = heatmap_colors,
                                    name = "Expression value",
                                    cluster_columns = FALSE,
                                    cluster_rows = FALSE,
                                    right_annotation = cluster_row_annotation,
                                    heatmap_legend_param = list(legend_height = grid::unit(4, "cm"),
                                                                title_position = "lefttop-rot"))

png(file.path(output_dir, paste0(tolower(cluster_label), "_cluster_heatmap_unscaled.png")), width = 9000, height = 10000, res = 72 * 15)
draw(cluster_heatmap_unscaled)
dev.off()


scaled_cluster_expression <- apply(cluster_median_expression, 2, ZScore)

cluster_heatmap_scaled <- Heatmap(scaled_cluster_expression,
                                  col = heatmap_colors,
                                  name = "Z score",
                                  cluster_columns = FALSE,
                                  cluster_rows = FALSE,
                                  right_annotation = cluster_row_annotation,
                                  heatmap_legend_param = list(legend_height = grid::unit(4, "cm"),
                                                              title_position = "lefttop-rot"))

png(file.path(output_dir, paste0(tolower(cluster_label), "_cluster_heatmap_scaled.png")), width = 9000, height = 10000, res = 72 * 15)
draw(cluster_heatmap_scaled)
dev.off()

saveRDS(spe_sub, file.path(output_dir, paste0(tolower(cluster_label), "_cluster_spe.rds")))

print("Finish!")