suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(apeglm)
})

count_file <- "results/counts/jq1/gene_counts.tsv"
sample_file <- "config/samples.tsv"
out_dir <- "results/deseq2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. 读取counts矩阵
counts <- read.delim(count_file, skip = 1, check.names = FALSE)

# 设置基因ID为行名，删除基因信息列（前6列）
rownames(counts) <- counts[[1]]
counts <- counts[, -(1:6), drop = FALSE]

# 清理列名：提取样本名称
colnames(counts) <- sub("^.*/", "", colnames(counts))  # 删除路径前缀
colnames(counts) <- sub("\\.sorted\\.bam$", "", colnames(counts))  # 删除后缀

cat("Count matrix columns:", colnames(counts), "\n")
cat("Count matrix dimensions:", nrow(counts), "genes ×", ncol(counts), "samples\n")

# 2. 读取样本表型信息
samples <- read.delim(sample_file, stringsAsFactors = FALSE)
samples[] <- lapply(samples, trimws)
samples$time <- gsub("h", "", samples$time)

# 3. 对齐samples和counts
idx <- match(colnames(counts), samples$run_id)
if (anyNA(idx)) {
  miss <- is.na(idx)
  # 从列名推断treatment
  get_trt <- function(x) ifelse(grepl("dmso", tolower(x)), "DMSO", ifelse(grepl("jq1", tolower(x)), "JQ1", NA))
  idx[miss] <- match(tolower(get_trt(colnames(counts)[miss])), tolower(samples$treatment))
}
samples <- samples[idx, ]
samples$run_id <- colnames(counts)
rownames(samples) <- colnames(counts)

# 4. 创建联合分组列（time_treatment）
samples$time   <- gsub("h", "", samples$time)  # 去掉 time 中的 "h"
samples$group <- factor(
  paste(samples$time, samples$treatment, sep = "_")
)
cat("Groups found:", paste(levels(samples$group), collapse = ", "), "\n")

# 5. 构建DESeq2对象
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData = samples,
  design = ~ group
)

cat("Initial gene number:", nrow(dds), "\n")

# 6. 低表达基因过滤：在两个样本中的表达counts>=10
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

cat("After filtering low expression genes:", nrow(dds), "\n")

# 7. DESeq2标准化
dds <- DESeq(dds)

# 8. PCA分析
vsd <- vst(dds, blind = FALSE)

p_pca <- plotPCA(vsd, intgroup = c("group")) +
  theme_bw() +
  labs(title = "PCA Plot")
ggsave(file.path(out_dir, "PCA.png"), p_pca, width = 8, height = 6)

# 9. 差异比较与火山图
run_comparison <- function(dds, contrast_vec, comparison_name, out_dir) {
  
  res <- lfcShrink(dds, contrast = contrast_vec, type = "normal")
  
  res_tbl <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    arrange(padj)
    
  write.csv(res_tbl, file.path(out_dir, paste0("DESeq2_", comparison_name, ".csv")), row.names = FALSE)
  
  # 火山图
  p <- res_tbl %>%
    mutate(sig = ifelse(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1, "DEG", "NS")) %>%
    ggplot(aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("DEG" = "red", "NS" = "gray50")) +
    theme_bw() +
    labs(title = paste0(comparison_name), x = "log2 fold change", y = "-log10 adjusted p-value") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")
  
  ggsave(file.path(out_dir, paste0("volcano_", comparison_name, ".png")), p, width = 8, height = 6)
  
  # 提取差异基因
  deg <- res_tbl %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)
  write.csv(deg, file.path(out_dir, paste0("DEG_", comparison_name, "_all.csv")), row.names = FALSE)
  
  cat(comparison_name, " DEGs: Total:", nrow(deg), "| Up:", sum(deg$log2FoldChange > 0), "| Down:", sum(deg$log2FoldChange < 0), "\n")
}

run_comparison(dds, c("group", "3_JQ1", "3_DMSO"), "3h_JQ1_vs_DMSO", out_dir)
run_comparison(dds, c("group", "24_JQ1", "24_DMSO"), "24h_JQ1_vs_DMSO", out_dir)

# 10. 保存DESeq2对象
saveRDS(dds, file.path(out_dir, "dds.rds"))
saveRDS(vsd, file.path(out_dir, "vsd.rds"))

cat("\nResults saved in:", out_dir, "\n")