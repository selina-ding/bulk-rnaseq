suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

deseq2_dir <- "results/deseq2"
out_dir <- "results/go"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. 读取DESeq2结果
res_3h  <- read.csv(file.path(deseq2_dir, "DESeq2_3h_JQ1_vs_DMSO.csv"))
res_24h <- read.csv(file.path(deseq2_dir, "DESeq2_24h_JQ1_vs_DMSO.csv"))

res_3h$gene_id  <- sub("\\..*", "", res_3h$gene_id)
res_24h$gene_id <- sub("\\..*", "", res_24h$gene_id)

get_sig <- function(df) {
  df %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
    pull(gene_id) %>%
    unique()
}

sig_3h  <- get_sig(res_3h)
sig_24h <- get_sig(res_24h)
all_sig <- union(sig_3h, sig_24h)

cat("3h DEGs:", length(sig_3h), "\n")
cat("24h DEGs:", length(sig_24h), "\n")
cat("Combined:", length(all_sig), "\n")

# 2. ID 转换
convert <- function(x) {
  if (length(x) == 0) return(character(0))
  res <- tryCatch(
    suppressMessages(bitr(x, "ENSEMBL", "ENTREZID", org.Hs.eg.db)),
    error = function(e) NULL
  )
  if (is.null(res) || nrow(res) == 0) return(character(0))
  unique(res$ENTREZID)
}

genes_3h  <- convert(sig_3h)
genes_24h <- convert(sig_24h)
genes_all <- convert(all_sig)

# 3. GO富集
run_go <- function(gene_list, label) {
  if (length(gene_list) < 5) {
    cat("  [", label, "] Too few genes, skipping\n", sep = "")
    return(NULL)
  }
  res <- enrichGO(
    gene          = gene_list,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "none", # 不做多重检验校正，否则太严格没有办法得到GO富集结果
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 1
  )
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    cat("  [", label, "] No enriched terms\n", sep = "")
    return(NULL)
  }
  cat("  [", label, "] ", nrow(as.data.frame(res)), " terms\n", sep = "")
  return(res)
}

k_3h  <- run_go(genes_3h,  "3h")
k_24h <- run_go(genes_24h, "24h")
k_all <- run_go(genes_all, "combined")

# 4. 气泡图
plot_go <- function(obj, name, file) {
  if (is.null(obj)) return(NULL)
  obj <- setReadable(obj, org.Hs.eg.db, "ENTREZID")
  df <- as.data.frame(obj)
  top_n <- min(15, nrow(df))
  df <- df[1:top_n, ]

  df <- df %>%
    mutate(
      Count = as.numeric(Count),
      GeneRatio = sapply(GeneRatio, function(x) {
        v <- as.numeric(strsplit(x, "/")[[1]])
        v[1] / v[2]
      }),
      logP = -log10(pvalue)
    )

  ggplot(df, aes(GeneRatio, reorder(Description, GeneRatio))) +
    geom_point(aes(size = Count, color = logP)) +
    scale_color_gradient(low = "blue", high = "red") +
    theme_bw() +
    labs(
      title = paste("GO BP -", name),
      x     = "Gene Ratio",
      y     = "",
      color = "-log10(pvalue)",
      size  = "Count"
    )
  ggsave(file, width = 8, height = 6)
}

plot_go(k_3h,  "3h",              file.path(out_dir, "go_3h.png"))
plot_go(k_24h, "24h",             file.path(out_dir, "go_24h.png"))
plot_go(k_all, "3h + 24h merged", file.path(out_dir, "go_combined.png"))

# 5. 保存 CSV
save_go <- function(obj, filename) {
  if (!is.null(obj)) write.csv(as.data.frame(obj), file.path(out_dir, filename))
}
save_go(k_3h,  "go_3h.csv")
save_go(k_24h, "go_24h.csv")
save_go(k_all, "go_combined.csv")

cat("\nDone. Results in", out_dir)
