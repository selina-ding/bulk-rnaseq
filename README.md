# Bulk RNA-seq GSEA 分析完整指南

## 📁 项目文件结构

```
bulk_rnaseq/
├── scripts/
│   ├── 07_deseq2.R                    # DESeq2 分析
│   ├── 08_go_kegg.R                   # GO/KEGG 富集
│   ├── 09_gsea.R                      # ⭐ GSEA 分析（主脚本）
│   ├── install_gsea_packages.R        # 包安装
│   └── check_gsea_prerequisites.R     # 前置条件检查
├── GSEA_QuickRef.md                   # 快速参考
├── GSEA_Guide.md                      # 详细指南
└── results/
    ├── deseq2/                        # DESeq2 输出
    ├── go_kegg/                       # GO/KEGG 输出
    └── gsea/                          # ⭐ GSEA 输出（本脚本生成）
```

## 🚀 快速开始

### 第一步：检查前置条件

```bash
cd /home/dyx/course/linux/bulk_rnaseq
Rscript scripts/check_gsea_prerequisites.R
```

**这个脚本会检查**：
- ✅ 所有必要的输入文件是否存在
- ✅ 所有必要的 R 包是否已安装
- ✅ GO/KEGG 结果是否包含足够的基因集
- ✅ 系统资源（内存、CPU、磁盘）

### 第二步：安装依赖包（如果需要）

```bash
Rscript scripts/install_gsea_packages.R
```

### 第三步：运行 GSEA 分析

```bash
Rscript scripts/09_gsea.R
```

**预计运行时间**：5-15 分钟（取决于基因集数量和系统性能）

## 📊 输出文件说明

| 文件名 | 类型 | 说明 |
|-------|------|------|
| `GSEA_results_3h.csv` | 数据表 | 3h 完整 GSEA 结果（所有基因集） |
| `GSEA_results_24h.csv` | 数据表 | 24h 完整 GSEA 结果（所有基因集） |
| **`01_NES_heatmap_3h.pdf`** | **图表** | **3小时 NES 矩阵热图**（padj<0.05） |
| **`02_NES_heatmap_24h.pdf`** | **图表** | **24小时 NES 矩阵热图**（padj<0.05） |
| **`03_NES_heatmap_combined.pdf`** | **图表** | **3h vs 24h 合并 NES 热图**（关键用途） |
| **`04_NES_trend.pdf`** | **图表** | **GSEA 趋势折线图**（时间动态变化） |
| `05_GSEA_statistics.pdf` | 图表 | 上调/下调通路统计柱状图 |

**⭐ 主要输出**：
1. **01_NES_heatmap_3h.pdf** - 3h 时间点的关键通路
2. **02_NES_heatmap_24h.pdf** - 24h 时间点的关键通路  
3. **04_NES_trend.pdf** - 通路随时间的变化趋势（最重要）

## 📈 三张关键图表详细解读

### 图1: NES 热图（3h）- 01_NES_heatmap_3h.pdf

```
┌─────────────────────────────────────┐
│  3h  │ GO_1 │ GO_2 │ KEGG_1 │ ...  │
├──────┼──────┼──────┼────────┤      │
│ 3h   │ 1.25 │-0.82 │  0.45  │ ...  │  
└─────────────────────────────────────┘
       深红   深蓝   浅蓝
       
深红 = NES > 1 (显著上调)
浅蓝 = NES < 0 (下调)
```

**解读**：
- 红色基因集 = JQ1 处理中上调
- 蓝色基因集 = DMSO 对照中上调
- 显示 3 小时内激活的主要生物过程

### 图2: NES 热图（24h）- 02_NES_heatmap_24h.pdf

类似于 3h 热图，但显示 24 小时时的结果

**对比 3h 热图**：
- 寻找持续出现的基因集（时间一致性）
- 寻找只出现在 24h 的新基因集（时间特异性）

### 图3: NES 趋势折线图 - 04_NES_trend.pdf

```
NES值
 2.0  ┤
      │  ╱╲___   pathway_A (↑激活)
 1.0  │ ╱    ╲___
      │╱         ╲   pathway_B (↓关闭)
 0.0  ├────────────┼─────
      │         ╱╲
-1.0  │    ╱╲__╱  ╲__  pathway_C (波动)
      └─────────────────
      3h        24h
```

**解读趋势**：
- **向上的线** = 生物过程被激活并维持
- **向下的线** = 生物过程被关闭
- **波动的线** = 动态调节的过程
- **红色实心点** = 显著时间点（padj < 0.05）
- **白色空心点** = 不显著时间点

**生物学意义示例**：
- 持续向上 → JQ1 的持续作用机制
- 24h 新增激活 → 延迟响应过程
- 3h 激活但 24h 关闭 → 早期响应

## 📋 CSV 结果表字段说明

### GSEA_results_3h.csv 示例

| pathway | pval | padj | log2err | ES | NES | size | leadingEdge |
|---------|------|------|---------|----|----|------|-------------|
| GO:0008150 | 0.001 | 0.05 | 0.25 | 0.45 | 1.32 | 120 | gene1,gene2,... |

**关键字段**：
- `pathway`: 基因集名称
- `NES`: ⭐ 标准化富集分数（-1 到 1），用于热图的数值
- `padj`: 调整后 p 值，< 0.05 表示显著
- `size`: 该基因集中的基因数
- `leadingEdge`: 驱动富集的关键基因列表

## 🎯 实际应用场景

### 场景 1: 发现 JQ1 的主要作用靶点

1. 查看 **01_NES_heatmap_3h.pdf**，找出红色最深的基因集
2. 查看该基因集的 `leadingEdge` 基因
3. 这些基因就是最直接被 JQ1 激活的基因

### 场景 2: 理解时间依赖的效应

1. 打开 **04_NES_trend.pdf**
2. 寻找从 3h 到 24h 持续上升的线（红色实心点）
3. 这些代表 JQ1 的主要效应通路

### 场景 3: 区分早期和晚期响应

1. 比对 **01** 和 **02** 两张热图
2. 只在 3h 出现的基因集 → 早期响应
3. 只在 24h 出现的基因集 → 晚期响应

## 🔧 常见用法

### 如何提取特定基因集的所有基因

```bash
# 编辑 CSV 文件查看 leadingEdge 字段
# 或使用 R：

gsea_results <- read.csv("results/gsea/GSEA_results_3h.csv")
# 找出最显著的基因集
top_pathway <- gsea_results %>% 
  filter(padj < 0.001) %>%
  slice(1)
# 提取 leading edge 基因
leading_genes <- unlist(strsplit(top_pathway$leadingEdge, ","))
```

### 如何创建自定义热图

```r
# 加载数据
gsea_3h <- read.csv("results/gsea/GSEA_results_3h.csv")
gsea_24h <- read.csv("results/gsea/GSEA_results_24h.csv")

# 过滤想要的通路
sig_pathways <- gsea_3h %>% 
  filter(padj < 0.01, abs(NES) > 1) %>%
  pull(pathway)

# 创建矩阵并绘制
# （详见 GSEA_Guide.md）
```

## ⚠️ 常见问题排查

### Q1: 出现 "No gene sets found" 警告

**检查清单**：
1. 确认 GO/KEGG 分析已运行：`ls results/go_kegg/*.csv`
2. 确认文件不为空：`wc -l results/go_kegg/*.csv`
3. 如果为空，重新运行 `08_go_kegg.R`

### Q2: 热图中数值都是 0

**原因**：基因集转换失败或没有显著基因集

**解决**：
```r
# 检查转换成功的基因集数量
cat("转换成功的基因集数量:", length(all_gene_sets))
```

### Q3: 趋势图中没有线条

**原因**：没有足够显著的基因集

**解决**：
- 修改脚本中的 `padj < 0.05` 为 `padj < 0.1`
- 或修改 `abs(NES) > 0` 来包括弱富集信号

### Q4: 脚本运行很慢

**原因**：基因集太多或系统资源有限

**解决**：
```r
# 修改脚本中的 nproc 参数
nproc = 4  # 改为 2 或 1
```

## 📚 进一步学习资源

### 相关文档
- **GSEA_Guide.md** - 详细原理和参数说明
- **GSEA_QuickRef.md** - 快速参考卡
- **GO/KEGG 分析** - 前置步骤说明

### 官方资源
- [GSEA 官方网站](http://www.broadinstitute.org/gsea)
- [fgsea R 包文档](https://bioconductor.org/packages/fgsea/)
- [clusterProfiler 中文教程](https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html)

## 💾 完整工作流总结

```
1. 数据准备
   └─ samples.tsv, gene_counts.tsv

2. DESeq2 分析 (07_deseq2.R)
   └─ 输出: DESeq2_3h.csv, DESeq2_24h.csv
   
3. GO/KEGG 富集 (08_go_kegg.R)
   └─ 输出: GO_BP_*.csv, KEGG_*.csv

4. ⭐ GSEA 分析 (09_gsea.R) ← 当前
   ├─ 输入: DESeq2 结果 + GO/KEGG 基因集
   └─ 输出: NES 热图 + 趋势图

5. 结果解释和发表
   └─ 整合所有图表讲述生物学故事
```

## 🎓 如何在论文中呈现结果

### 推荐组合
- **Figure A**: DESeq2 火山图（展示显著基因）
- **Figure B**: GO/KEGG 气泡图（识别功能）
- **Figure C**: GSEA NES 合并热图（验证富集趋势）
- **Figure D**: GSEA 趋势图（展示时间动态）

### 图表说明示例

> **Figure 3. Gene set enrichment analysis reveals JQ1-induced pathway dynamics.**
>
> (A) NES heatmap showing significant pathway enrichment at 3h (left) and 24h (right) after JQ1 treatment. Red indicates positive enrichment (JQ1-activated pathways); blue indicates negative enrichment (control-activated pathways).
>
> (B) Temporal trend of pathway NES scores across treatment durations. Each line represents a significantly enriched pathway. Solid dots indicate statistical significance (padj < 0.05). Sustained upward trends reflect persistent pathway activation by JQ1.

## 📞 需要帮助？

如有问题，请查看：
1. **GSEA_QuickRef.md** - 快速问题排查
2. **GSEA_Guide.md** - 详细解释
3. 脚本注释 - 代码级别的说明

祝分析顺利！ 🎉
