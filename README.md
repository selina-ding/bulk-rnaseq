# Bulk RNA-seq Pipeline — JQ1 vs DMSO (TNBC)

分析 JQ1（BET 抑制剂）处理 SUM159/SUM159R 三阴性乳腺癌细胞系后的转录组变化，时间点 3h 和 24h。

## 目录结构

```
bulk_rnaseq/
├── config/
│   ├── config.sh              # shell 脚本共用配置（路径、线程、参考基因组）
│   └── samples.tsv            # 样本元数据（run_id, cell_line, time, treatment）
├── scripts/
│   ├── 00_check_env.sh        # 环境检查
│   ├── 01_fastqc.sh           # 原始数据质控
│   ├── 02_fastp.sh            # 去接头 + 低质量修剪
│   ├── 03_multiqc.sh          # 汇总质控报告
│   ├── 04_hisat2_index.sh     # 构建 HISAT2 索引
│   ├── 05_hisat2_align.sh     # 比对到参考基因组
│   ├── 06_featurecounts.sh    # 基因计数（featureCounts）
│   ├── 07_deseq2.R            # DESeq2 差异表达分析
│   ├── 08_go.R                # GO 富集分析
│   ├── 09_gsea.R              # GSEA 分析
│   └── lib.sh                 # 公共 bash 函数库
├── env/
│   └── conda_env.yaml         # conda 环境定义
├── tests/
│   └── data/
│       ├── fastq/             # 测试用 FASTQ
│       └── samples_test.tsv   # 测试用样本表
├── data/                      # 输入数据（需自行准备）
│   ├── raw_fastq/
│   ├── trimmed_fastq/
│   └── reference/             # 参考基因组 + HISAT2 索引
├── results/                   # 输出（由脚本生成，不入 git）
│   ├── qc/
│   ├── bam/
│   ├── counts/
│   ├── deseq2/
│   ├── go/
│   └── gsea/
└── logs/
```

## 分析流程

```
fastq  →  fastp  →  hisat2  →  featureCounts  →  DESeq2  →  GO  →  GSEA
(01)      (02)       (05)         (06)             (07)      (08)    (09)
  │                    │
fastqc               multiqc
(01)                  (03)
```

- **比较组**：JQ1 vs DMSO，分 3h 和 24h 两个时间点
- **参考基因组**：GRCh38, gencode v44 注释

## 快速开始

### 1. 环境

```bash
conda env create -f env/conda_env.yaml
conda activate tnbc-rnaseq
bash scripts/00_check_env.sh
```

### 2. 准备参考基因组和输入数据

下载 GRCh38 参考序列和 gencode v44 GTF 放入 `data/reference/`，构建 HISAT2 索引：

```bash
bash scripts/04_hisat2_index.sh
```

原始 FASTQ 放入 `data/raw_fastq/`，路径需与 `config/samples.tsv` 中的 `fq` 列一致。

### 3. 运行

按编号顺序执行：

```bash
# 质控
bash scripts/01_fastqc.sh config/samples.tsv
bash scripts/02_fastp.sh config/samples.tsv
bash scripts/03_multiqc.sh config/samples.tsv

# 比对 + 计数
bash scripts/05_hisat2_align.sh config/samples.tsv
bash scripts/06_featurecounts.sh config/samples.tsv

# 差异表达 + 富集分析
Rscript scripts/07_deseq2.R
Rscript scripts/08_go.R
Rscript scripts/09_gsea.R
```

> **测试**：可先用 `tests/data/samples_test.tsv` 替换 `config/samples.tsv` 跑通流程。

### 4. 使用自己的数据

只需修改 `config/samples.tsv`（保持 5 列不变）：

```
run_id	cell_line	time	treatment	fq
SRRxxx	SUM159	3h	DMSO	data/raw_fastq/SRRxxx.fastq.gz
...
```

如果更换物种，修改 `config/config.sh` 中的 `GENOME_FA`、`GTF`、`HISAT2_INDEX`，以及 R 脚本中的 `org.Hs.eg.db` → 对应物种的 `org.*.eg.db`。

## 输出

### 核心结果文件

| 路径 | 说明 |
|------|------|
| `results/deseq2/PCA.png` | PCA 图 |
| `results/deseq2/volcano_*.png` | 火山图（3h / 24h） |
| `results/deseq2/DESeq2_*.csv` | 全基因差异表达表 |
| `results/deseq2/DEG_*_all.csv` | 显著差异基因（padj<0.05, \|log2FC\|≥1） |
| `results/go/go_*.png` | GO 富集气泡图 |
| `results/go/go_*.csv` | GO 富集结果表 |
| `results/gsea/NES_heatmap_combined.png` | 3h+24h 合并 NES 热图 |
| `results/gsea/NES_trend.png` | NES 时间趋势折线图 |
| `results/gsea/GSEA_results_*.csv` | GSEA 结果表 |
