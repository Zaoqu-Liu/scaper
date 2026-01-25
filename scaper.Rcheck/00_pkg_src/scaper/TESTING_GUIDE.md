# scaper Package Testing Guide

## 测试前准备 (Pre-Testing Setup)

### 1. 安装依赖包 (Install Dependencies)

```r
# 安装必需的包
install.packages(c("devtools", "microbenchmark", "pheatmap"))

# 确保已安装所有依赖
install.packages(c("magrittr", "xml2", "stringr", "dplyr", "Seurat", 
                   "SeuratObject", "VAM", "knitr", "rmarkdown"))
```

### 2. 构建和检查包 (Build and Check Package)

```bash
# 在终端中运行（Terminal）
cd /Users/liuzaoqu/Desktop/develop/scaper

# R CMD check - 标准检查
R CMD build .
R CMD check scaper_0.2.1.tar.gz

# 或使用 devtools（推荐）
R -e "devtools::check()"
```

---

## 测试步骤 (Testing Steps)

### 步骤 1: 安装更新后的包 (Install Updated Package)

```r
# 方法 1: 从本地安装
devtools::install("/Users/liuzaoqu/Desktop/develop/scaper")

# 方法 2: 从 GitHub 安装（如果已推送）
devtools::install_github("Zaoqu-Liu/scaper")
```

### 步骤 2: 验证安装 (Verify Installation)

```r
# 加载包
library(scaper)

# 检查版本
packageVersion("scaper")  # 应该显示 0.2.1

# 查看帮助
?scaper
?scape
?scapeForSeurat
```

### 步骤 3: 基本功能测试 (Basic Functionality Tests)

```r
library(Seurat)
library(SeuratObject)

# 测试 1: supportedCytokines
cytokines_cytosig <- supportedCytokines(database = "cytosig")
print(length(cytokines_cytosig))  # 应该是 41

cytokines_reactome <- supportedCytokines(database = "reactome")
print(length(cytokines_reactome))  # 应该是 30

# 测试 2: genesetCytoSig
file_cytosig <- system.file("extdata", "IL6_output.csv", package = "scaper")
geneset_il6 <- genesetCytoSig(cytokine.eval = "IL6", file.name = file_cytosig)
print(head(geneset_il6, 10))

# 测试 3: genesetReactome
file_reactome <- system.file("extdata", "IL6_Interleukin6_signaling.xml", package = "scaper")
geneset_il6_reactome <- genesetReactome(cytokine.eval = "IL6", file.name = file_reactome)
print(head(geneset_il6_reactome, 10))

# 测试 4: scape() with CytoSig
data(pbmc_small)
pbmc_small <- NormalizeData(pbmc_small)
counts.matrix <- as.data.frame(t(as.matrix(GetAssayData(pbmc_small, slot = "data"))))

result_scape <- scape(counts.matrix = counts.matrix, 
                      database = "cytosig", 
                      cytokine = c("IL6", "IL4"))
print(dim(result_scape))
print(summary(result_scape))

# 测试 5: scapeForSeurat() with CytoSig
result_seurat <- scapeForSeurat(seurat.object = pbmc_small, 
                                database = "cytosig", 
                                cytokine = c("IL6", "IL4"), 
                                normalize = FALSE)
print(result_seurat)
print(dim(result_seurat@assays$scape@data))
```

### 步骤 4: 科学准确性验证 (Scientific Accuracy Verification)

```r
# 这个测试验证正负权重基因的处理是否正确

# 使用完整的细胞因子列表
result_full <- scape(counts.matrix = counts.matrix, 
                     database = "cytosig", 
                     cytokine = "all")

# 检查输出
print("Output dimensions:")
print(dim(result_full))

print("\nOutput range (should be [0, 1]):")
print(paste("Min:", min(result_full)))
print(paste("Max:", max(result_full)))

print("\nCheck for NA values:")
print(paste("Any NA:", any(is.na(result_full))))

print("\nValues should be between 0 and 1:")
print(paste("All in range:", all(result_full >= 0 & result_full <= 1)))

# 特定细胞因子的活性分数统计
print("\nIL6 activity statistics:")
il6_scores <- result_full[, "IL6"]
print(summary(il6_scores))
```

### 步骤 5: 运行 Benchmark 测试 (Run Benchmark Tests)

```r
# 运行完整的 benchmark 测试
source("/Users/liuzaoqu/Desktop/develop/scaper/benchmark_test.R")

# 输出将显示：
# - 每个函数的执行时间
# - 正确性检查结果
# - 科学准确性验证结果
```

### 步骤 6: Seurat V4/V5 兼容性测试 (Seurat Compatibility Test)

```r
# 检查 Seurat 版本
print(paste("Seurat version:", packageVersion("Seurat")))

# 测试 Seurat V4 对象
if (packageVersion("Seurat") >= "4.0.0") {
  data(pbmc_small)
  pbmc_small <- NormalizeData(pbmc_small)
  
  # 应该正常工作
  result_v4 <- scapeForSeurat(seurat.object = pbmc_small, 
                              database = "cytosig", 
                              cytokine = c("IL6", "IL4"), 
                              normalize = FALSE)
  print("Seurat V4 test: PASS")
}

# 如果有 Seurat V5，也测试一下
if (packageVersion("Seurat") >= "5.0.0") {
  # V5 特定测试
  print("Testing Seurat V5 compatibility...")
  # 应该也能正常工作
}
```

### 步骤 7: 可视化测试 (Visualization Test)

```r
library(pheatmap)

# 创建热图
result <- scapeForSeurat(seurat.object = pbmc_small, 
                         database = "cytosig", 
                         cytokine = "all", 
                         normalize = TRUE)

cytokine_matrix <- as.data.frame(t(as.matrix(result@assays$scape@data)))

# 未聚类热图
pheatmap(cytokine_matrix, 
         fontsize_row = 4, 
         fontsize_col = 7,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Cytokine Activity Scores (Unclustered)")

# 聚类热图
pheatmap(cytokine_matrix, 
         fontsize_row = 4, 
         fontsize_col = 7,
         main = "Cytokine Activity Scores (Clustered)")
```

---

## 预期结果 (Expected Results)

### ✅ 通过标准 (Pass Criteria)

1. **安装成功**: 无错误或警告
2. **版本正确**: 显示 0.2.1
3. **基本功能**: 所有函数都能正常运行
4. **输出范围**: 所有 CDF 值在 [0, 1] 之间
5. **无 NA 值**: 输出不包含 NA
6. **性能提升**: benchmark 显示比旧版本快 15-30%
7. **Seurat 兼容**: V4 和 V5 都能正常工作
8. **可视化**: 热图能正常生成

### ⚠️ 警告标志 (Warning Signs)

- 任何错误消息
- CDF 值超出 [0, 1] 范围
- NA 值出现在输出中
- 性能比旧版本慢
- Seurat 对象访问失败

---

## 故障排除 (Troubleshooting)

### 问题 1: 安装失败

```r
# 解决方案：先安装所有依赖
install.packages(c("magrittr", "xml2", "stringr", "dplyr", 
                   "Seurat", "SeuratObject", "VAM"))
```

### 问题 2: Namespace 冲突

```r
# 解决方案：显式指定包
scaper::scape(...)
```

### 问题 3: Seurat 版本问题

```r
# 检查 Seurat 版本
packageVersion("Seurat")

# 如果版本太旧，升级
install.packages("Seurat")
```

### 问题 4: VAM 包未找到

```r
# VAM 可能需要从 GitHub 安装
devtools::install_github("rwfriedman/VAM")
```

---

## R CMD check 清单 (R CMD Check Checklist)

运行 `R CMD check` 应该通过以下检查：

- [ ] **00install.out**: 安装成功
- [ ] **00check.log**: 无 ERROR
- [ ] **00check.log**: 无 WARNING（或只有可接受的警告）
- [ ] **Examples**: 所有示例都能运行
- [ ] **Tests**: 所有测试通过（如果有）
- [ ] **Vignettes**: 能够成功构建
- [ ] **Namespace**: 无冲突
- [ ] **Documentation**: 所有函数都有文档

---

## 提交前最后检查 (Pre-Submission Final Check)

```r
# 1. 完整检查
devtools::check()

# 2. 构建文档
devtools::document()

# 3. 构建 vignettes
devtools::build_vignettes()

# 4. 拼写检查
devtools::spell_check()

# 5. 检查示例
devtools::run_examples()

# 6. 构建包
devtools::build()
```

---

## 性能基准参考 (Performance Benchmarks Reference)

在 MacBook Pro (M1, 16GB RAM) 上的典型结果：

```
Test 1: scape() CytoSig (5 cytokines, 80 cells, 230 genes)
  Mean: ~150-200 ms

Test 2: scape() Reactome (5 cytokines)
  Mean: ~100-150 ms

Test 3: scapeForSeurat() CytoSig (5 cytokines)
  Mean: ~200-250 ms

Test 4: scapeForSeurat() Reactome (5 cytokines)
  Mean: ~150-200 ms
```

实际时间会根据硬件和数据大小而变化。

---

## 报告问题 (Reporting Issues)

如果测试失败，请收集以下信息：

```r
# 系统信息
sessionInfo()

# 包版本
packageVersion("scaper")
packageVersion("Seurat")
packageVersion("VAM")

# 错误消息
# 复制完整的错误堆栈
```

然后在 GitHub 上创建 issue：
https://github.com/Zaoqu-Liu/scaper/issues

---

## 测试完成 (Testing Complete)

完成所有测试后，您应该确认：

✅ 包安装成功  
✅ 所有示例都能运行  
✅ 科学准确性验证通过  
✅ 性能有所提升  
✅ Seurat 兼容性良好  
✅ 无 linter 错误  
✅ R CMD check 通过  

**恭喜！您的包已准备好发布。** 🎉

---

**测试指南版本**: 1.0  
**最后更新**: 2026年1月25日
