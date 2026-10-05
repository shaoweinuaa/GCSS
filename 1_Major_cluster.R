library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)
library(ggthemes)
library(ggpubr)
library(plyr)
library(scater)
library(harmony)
library(SingleCellExperiment)
library(SpatialExperiment)
library(tidyverse)
library(dittoSeq)
library(FastPG)

set.seed(220225)

QuantileNormalize <- function(x, n) {
  x[is.na(x)] <- 0
  quantile_values <- quantile(x, probs = (0:(n + 1)) / (n + 1), na.rm = TRUE)
  
  upper_indices <- which(x > quantile_values[n + 1])
  replacement_pool <- x[x > quantile_values[n] & x < quantile_values[n + 1]]
  
  if (length(upper_indices) > 0 && length(replacement_pool) > 0) {
    x[upper_indices] <- sample(replacement_pool, length(upper_indices), replace = TRUE)
  }
  
  as.numeric(scale(x))
}

MinMaxScale <- function(x) {
  value_range <- range(x, na.rm = TRUE)
  if (!all(is.finite(value_range)) || diff(value_range) == 0) {
    return(rep(0, length(x)))
  }
  (x - value_range[1]) / diff(value_range)
}

ZScore <- function(x) {
  value_sd <- sd(x, na.rm = TRUE)
  if (!is.finite(value_sd) || value_sd == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / value_sd
}

raw_csv_dir <- "./Input/csv"
sample_info_file <- "./Input/info.csv"
marker_file <- "./Input/marker.csv"
output_dir <- "./Output/1_Major"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cofactor <- 1
quantile_interval_count <- 101

panel <- read.csv(marker_file, check.names = FALSE, fileEncoding = "GBK")

raw_files <- list.files(raw_csv_dir, pattern = "[.]csv$", full.names = TRUE)
cell_tables <- lapply(raw_files, function(raw_file) {
  cell_data <- read.csv(raw_file)
  signal_columns <- colnames(cell_data)[grep("X[0-9]+", colnames(cell_data))]
  channel_names <- str_split(signal_columns, pattern = "_", simplify = TRUE)[, 1] %>% str_sub(2, nchar(.))
  colnames(cell_data)[grep("X[0-9]+", colnames(cell_data))] <- panel$marker[match(channel_names, panel$channel)]
  cell_data <- cell_data[, !is.na(colnames(cell_data))]
  cell_data
})

names(cell_tables) <- raw_files

for (raw_file in raw_files) {
  current_names <- colnames(cell_tables[[raw_file]])
  duplicated_names <- any(str_detect(current_names, "[.]1$"))
  
  if (duplicated_names) {
    stop(raw_file, " contains duplicated marker names.")
  }
}

cell_table <- do.call("rbind.fill", cell_tables)

sample_info <- read.csv(sample_info_file, check.names = FALSE, fileEncoding = "GBK")
colnames(sample_info)[1] <- "roi_id"

marker_names <- panel %>% pull(marker)
counts_matrix <- cell_table %>% select(all_of(marker_names))

cell_metadata <- cell_table %>%
  select(roi, CellID) %>%
  mutate(roi_id = roi) %>%
  select(-roi) %>%
  left_join(sample_info, by = "roi_id")

roi_prefix <- str_split(cell_metadata$roi_id, "_", simplify = TRUE)[, 1]
cell_metadata$batch_id <- roi_prefix
cell_metadata$sample_id <- roi_prefix

spatial_coordinates <- cell_table %>% select(contains("position"))
colnames(spatial_coordinates) <- c("Pos_X", "Pos_Y")

spe <- SpatialExperiment(
  assays = list(counts = t(counts_matrix)),
  colData = cell_metadata,
  sample_id = as.character(cell_metadata$sample_id),
  image_id = as.character(cell_metadata$roi_id),
  spatialCoords = as.matrix(spatial_coordinates)
)

colnames(spe) <- paste0(spe$roi_id, ".", spe$CellID)

major_markers <- panel %>% filter(as.character(major_annotation) == "1") %>% pull(marker)
rowData(spe)$use_channel <- rownames(spe) %in% major_markers

transformed_expression <- asinh(counts(spe) / cofactor)
normalized_expression <- transformed_expression

for (marker_index in seq_len(nrow(normalized_expression))) {
  normalized_expression[marker_index, ] <- QuantileNormalize(normalized_expression[marker_index, ], quantile_interval_count)
  normalized_expression[marker_index, ] <- MinMaxScale(normalized_expression[marker_index, ])
}

assay(spe, "exprs") <- normalized_expression

spe <- runUMAP(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs", name = "UMAP")

expression_matrix <- t(assay(spe, "exprs"))[, rowData(spe)$use_channel, drop = FALSE]

svd_result <- svd(expression_matrix)
explained_variance <- svd_result$d^2 / sum(svd_result$d^2)
n_pcs <- which(cumsum(explained_variance) >= 0.95)[1]

if (length(unique(spe$batch_id)) < 2) {
  pca_result <- prcomp(expression_matrix, center = TRUE, scale. = FALSE)
  harmony_embedding <- pca_result$x[, seq_len(n_pcs), drop = FALSE]
} else {
  harmony_embedding <- harmony::HarmonyMatrix(expression_matrix, as.factor(spe$batch_id), do_pca = TRUE, npcs = n_pcs)
}

reducedDim(spe, "harmony") <- harmony_embedding

spe <- runUMAP(spe, dimred = "harmony", name = "UMAP_harmony")

spe$batch_id <- factor(spe$batch_id, levels = unique(spe$batch_id))
legend_columns <- ceiling(length(levels(spe$batch_id)) / 10)

batch_before_plot <- dittoDimPlot(spe,
                                  var = "batch_id",
                                  reduction.use = "UMAP",
                                  size = 0.2) + guides(color = guide_legend(
                                    ncol = legend_columns,
                                    byrow = FALSE,
                                    override.aes = list(size = 2)
                                  )) + ggtitle("Batch ID on UMAP before correction")

batch_after_plot <- dittoDimPlot(spe,
                                 var = "batch_id",
                                 reduction.use = "UMAP_harmony",
                                 size = 0.2) + guides(color = guide_legend(
                                   ncol = legend_columns,
                                   byrow = FALSE,
                                   override.aes = list(size = 2)
                                 )) + ggtitle("Batch ID on UMAP after correction") + labs(x = "UMAP1", y = "UMAP2")

batch_correction_plot <- ggarrange(batch_before_plot, batch_after_plot, ncol = 2)

ggsave(file.path(output_dir, "batch_correction_umap.png"), batch_correction_plot, height = 3.5, width = legend_columns + 7, dpi = 300)


clustering_embedding <- reducedDim(spe, "harmony")
cluster_result <- FastPG::fastCluster(as.matrix(clustering_embedding),
                                      k = 30,
                                      num_threads = 100)
spe$pg_cluster <- factor(cluster_result$communities)

cluster_colors <- ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value
cluster_colors <- colorRampPalette(cluster_colors)(length(levels(spe$pg_cluster)))
names(cluster_colors) <- levels(spe$pg_cluster)

cluster_umap_plot <- dittoDimPlot(spe, var = "pg_cluster", reduction.use = "UMAP_harmony", size = 0.2, do.label = TRUE, labels.repel = TRUE, labels.highlight = FALSE) + 
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 4))) + 
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_blank()) + 
  scale_color_manual(values = cluster_colors) + 
  labs(x = "UMAP1", y = "UMAP2", title = "")

ggsave(file.path(output_dir, "major_cluster_umap.png"), cluster_umap_plot, height = 7, width = 8, dpi = 300)

for (marker_name in rownames(spe)) {
  feature_plot <- dittoDimPlot(spe, var = marker_name, reduction.use = "UMAP_harmony", assay = "exprs", size = 0.2) +
    scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
    theme(axis.text = element_text(size = 15), 
          axis.title = element_text(size = 15),
          legend.title = element_blank()) +
    labs(x = "UMAP1", y = "UMAP2")
  
  ggsave(file.path(output_dir, paste0("feature_", marker_name, ".png")), feature_plot, height = 4, width = 4.5, dpi = 300)
}


cluster_median_expression <- as.data.frame(
  t(assay(spe, "exprs"))[, rowData(spe)$use_channel, drop = FALSE],
  check.names = FALSE
)
cluster_median_expression$pg_cluster <- spe$pg_cluster

cluster_median_expression <- cluster_median_expression %>%
  group_by(pg_cluster) %>%
  summarise(across(everything(), median), .groups = "drop") %>%
  column_to_rownames("pg_cluster") %>%
  as.matrix()

cluster_median_expression <- cluster_median_expression[
  levels(spe$pg_cluster),
  ,
  drop = FALSE
]

heatmap_colors <- colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(101)

cluster_counts <- table(spe$pg_cluster)
cluster_row_annotation <- rowAnnotation(counts = anno_barplot(as.numeric(cluster_counts), 
                                                              gp = grid::gpar(fill = cluster_colors)))

unscaled_heatmap <- Heatmap(cluster_median_expression,
                            col = heatmap_colors,
                            name = "Expression value",
                            cluster_columns = FALSE,
                            cluster_rows = FALSE,
                            right_annotation = cluster_row_annotation,
                            heatmap_legend_param = list(legend_height = grid::unit(4, "cm"),
                                                        title_position = "lefttop-rot")
                            )

png(file.path(output_dir, "major_cluster_heatmap_unscaled.png"), width = 9000, height = 10000, res = 72 * 15)
draw(unscaled_heatmap)
dev.off()


scaled_cluster_expression <- apply(cluster_median_expression, 2, ZScore)

scaled_heatmap <- Heatmap(scaled_cluster_expression,
                          col = heatmap_colors,
                          name = "Z score",
                          cluster_columns = FALSE,
                          cluster_rows = FALSE,
                          right_annotation = cluster_row_annotation,
                          heatmap_legend_param = list(legend_height = grid::unit(4, "cm"),
                                                      title_position = "lefttop-rot")
)

png(file.path(output_dir, "major_cluster_heatmap_scaled.png"), width = 9000, height = 10000, res = 72 * 15)
draw(scaled_heatmap)
dev.off()

saveRDS(spe, file.path(output_dir, "major_cluster_spe.rds"))
