suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(fgsea)
  library(pheatmap)
})

deseq2_dir  <- "results/deseq2"
go_kegg_dir <- "results/go"
out_dir     <- "results/gsea"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. 读取DESeq2结果，构建排序基因列表
res_3h  <- read.csv(file.path(deseq2_dir, "DESeq2_3h_JQ1_vs_DMSO.csv"))
res_24h <- read.csv(file.path(deseq2_dir, "DESeq2_24h_JQ1_vs_DMSO.csv"))
res_3h$gene_id  <- sub("\\..*", "", res_3h$gene_id)
res_24h$gene_id <- sub("\\..*", "", res_24h$gene_id)

prepare_ranked_list <- function(res_df, name) {
  res_df <- res_df %>% filter(!is.na(log2FoldChange), is.finite(log2FoldChange))
  rank_list <- res_df$log2FoldChange
  names(rank_list) <- res_df$gene_id
  rank_list <- sort(rank_list, decreasing = TRUE)
  cat(name, "ranked list:", length(rank_list), "genes\n")
  return(rank_list)
}

rank_3h  <- prepare_ranked_list(res_3h,  "3h")
rank_24h <- prepare_ranked_list(res_24h, "24h")

# 2. 从GO结果中提取基因集
get_gene_sets <- function(go_file) {
  if (!file.exists(go_file)) {
    cat("File not found:", go_file, "\n")
    return(NULL)
  }
  go_df <- read.csv(go_file)
  if (nrow(go_df) == 0 || !"geneID" %in% colnames(go_df)) return(list())
  gene_sets <- list()
  for (i in seq_len(nrow(go_df))) {
    gene_sets[[go_df$Description[i]]] <- strsplit(go_df$geneID[i], "/")[[1]]
  }
  return(gene_sets)
}

gs_files <- c("go_3h", "go_24h", "go_combined")
all_gene_sets <- unlist(lapply(gs_files, function(f) {
  get_gene_sets(file.path(go_kegg_dir, paste0(f, ".csv"))) %||% list()
}), recursive = FALSE)

cat("Total gene sets loaded:", length(all_gene_sets), "\n")

# 3. 运行GSEA
run_gsea <- function(rank_list, gene_sets, name) {
  # 基因ID转换，把ENSEMBL转换为ENTREZID
  names_mapped <- tryCatch(
    suppressMessages(bitr(names(rank_list), fromType = "ENSEMBL", toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db, drop = FALSE)),
    error = function(e) NULL
  )
  if (is.null(names_mapped) || nrow(names_mapped) == 0) {
    cat("  WARNING: Failed to map gene IDs\n")
    return(NULL)
  }

  rank_list_eg <- rank_list[names_mapped$ENSEMBL]
  names(rank_list_eg) <- names_mapped$ENTREZID
  rank_list_eg <- rank_list_eg[!is.na(names(rank_list_eg))]
  rank_list_eg <- rank_list_eg[order(abs(rank_list_eg), decreasing = TRUE)]
  rank_list_eg <- rank_list_eg[!duplicated(names(rank_list_eg))]

  gene_sets_eg <- gene_sets[lengths(gene_sets) >= 5 & lengths(gene_sets) <= 500]

  if (length(gene_sets_eg) == 0) {
    cat("  WARNING: No gene sets passed filters\n")
    return(NULL)
  }

  fgsea_results <- as.data.frame(
    fgsea(pathways = gene_sets_eg, stats = rank_list_eg,
          minSize = 5, maxSize = 500, nproc = 4)
  )
  fgsea_results[order(fgsea_results$NES, decreasing = TRUE), ]
}

gsea_3h  <- run_gsea(rank_3h,  all_gene_sets, "3h")
gsea_24h <- run_gsea(rank_24h, all_gene_sets, "24h")

if (!is.null(gsea_3h))  gsea_3h  <- as.data.frame(gsea_3h)
if (!is.null(gsea_24h)) gsea_24h <- as.data.frame(gsea_24h)

# 4. 保存结果
for (tp in c("3h", "24h")) {
  gsea <- if (tp == "3h") gsea_3h else gsea_24h
  if (!is.null(gsea)) {
    gsea_out <- as.data.frame(gsea)
    gsea_out$leadingEdge <- sapply(gsea_out$leadingEdge, paste, collapse = "/")
    write.csv(gsea_out, file.path(out_dir, paste0("GSEA_results_", tp, "h.csv")), row.names = FALSE)
    cat(tp, "h GSEA results saved:", nrow(gsea), "pathways\n")
  }
}

# 5. 合并两个时间的NES热图
if (is.null(gsea_3h) || is.null(gsea_24h)) {
  cat("Skipping combined heatmap - missing GSEA results\n")
  quit(save = "no", status = 0)
}

sig_3h_path  <- gsea_3h  %>% filter(pval < 0.05) %>% pull(pathway)
sig_24h_path <- gsea_24h %>% filter(pval < 0.05) %>% pull(pathway)
all_sig_path <- unique(c(head(sig_3h_path, 15), head(sig_24h_path, 15)))

if (length(all_sig_path) > 0) {
  nes_matrix <- matrix(0, nrow = 2, ncol = length(all_sig_path))
  rownames(nes_matrix) <- c("3h", "24h")
  colnames(nes_matrix) <- all_sig_path
  for (pw in all_sig_path) {
    if (pw %in% gsea_3h$pathway)  nes_matrix["3h", pw]  <- gsea_3h[gsea_3h$pathway == pw, "NES"][1]
    if (pw %in% gsea_24h$pathway) nes_matrix["24h", pw] <- gsea_24h[gsea_24h$pathway == pw, "NES"][1]
  }

  png(file.path(out_dir, "NES_heatmap_combined.png"), width = 8, height = 6, units = "in", res = 300)
  pheatmap(nes_matrix,
           display_numbers = TRUE, number_format = "%.2f",
           breaks = seq(-3, 3, by = 0.1),
           color = colorRampPalette(c("blue", "white", "red"))(100),
           cluster_cols = TRUE, cluster_rows = FALSE,
           main = "Combined NES Matrix (3h vs 24h)",
           fontsize_col = 8, fontsize_row = 12, fontsize_number = 8)
  dev.off()
}

# 6. NES 趋势折线图
gsea_3h$timepoint  <- "3h"
gsea_24h$timepoint <- "24h"
trend_data <- bind_rows(gsea_3h, gsea_24h) %>%
  mutate(sig = pval < 0.05) %>%
  dplyr::select(pathway, timepoint, NES, pval, sig)

top_pathways <- trend_data %>%
  filter(sig) %>%
  group_by(pathway) %>%
  summarise(max_abs = max(abs(NES)), .groups = "drop") %>%
  slice_max(max_abs, n = 15) %>%
  pull(pathway)

if (length(top_pathways) > 0) {
  trend_plot_data <- trend_data %>% filter(pathway %in% top_pathways)

  p <- ggplot(trend_plot_data, aes(x = timepoint, y = NES, color = pathway, group = pathway)) +
    geom_line(size = 1, alpha = 0.7) +
    geom_point(aes(fill = sig), size = 3, shape = 21, color = "black") +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "white")) +
    theme_bw() +
    theme(legend.text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(title = "GSEA NES Trend (3h vs 24h)", x = "Time Point", y = "NES",
         color = "Pathway", fill = "Significant\n(pval < 0.05)") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5)

  ggsave(file.path(out_dir, "NES_trend.png"), p, width = 8, height = 6)
  cat("NES trend plot saved\n")
} else {
  cat("No data for trend plot\n")
}

# 7. 汇总
# 统计GSEA每个时间点显著上调/下调的通路数量
cat("\nGSEA Analysis Summary\n")
for (tp in c("3h", "24h")) {
  gsea <- if (tp == "3h") gsea_3h else gsea_24h
  if (is.null(gsea)) next

# 按 pvalue<0.05 筛选显著通路，按 NES 正负区分上调/下调 
  sig <- nrow(filter(gsea, pval < 0.05))
  up   <- nrow(filter(gsea, pval < 0.05, NES > 0))
  down <- nrow(filter(gsea, pval < 0.05, NES < 0))
  cat("\n", tp, "h (JQ1 vs DMSO):\n  Total:", nrow(gsea),
      "\n  Significant (pval<0.05):", sig,
      "\n    - Up:",   up,
      "\n    - Down:", down, "\n")
}
cat("\nResults saved in:", out_dir)
