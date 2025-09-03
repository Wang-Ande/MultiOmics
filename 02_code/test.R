# test
library(pak)
pkg_install("mixOmicsTeam/mixOmics")
devtools::install_github("mixOmicsTeam/mixOmics")
library(mixOmics)
data(breast.TCGA)
   
# Extract training data and name each data frame
# Store as list
X <- list(mRNA = breast.TCGA$data.train$mrna, 
          miRNA = breast.TCGA$data.train$mirna, 
          protein = breast.TCGA$data.train$protein)

# Outcome
Y <- breast.TCGA$data.train$subtype
summary(Y)

design <- matrix(0.1, ncol = length(X), nrow = length(X), 
                 dimnames = list(names(X), names(X)))
diag(design) <- 0
design 

# pls
res1.pls.tcga <- pls(X$mRNA, X$protein, ncomp = 1)
cor(res1.pls.tcga$variates$X, res1.pls.tcga$variates$Y)

res2.pls.tcga <- pls(X$mRNA, X$miRNA, ncomp = 1)
cor(res2.pls.tcga$variates$X, res2.pls.tcga$variates$Y)

res3.pls.tcga <- pls(X$protein, X$miRNA, ncomp = 1)
cor(res3.pls.tcga$variates$X, res3.pls.tcga$variates$Y)

# numbers of components
diablo.tcga <- block.plsda(X, Y, ncomp = 5, design = design)

set.seed(123) # For reproducibility, remove for your analyses
perf.diablo.tcga = perf(diablo.tcga, validation = 'Mfold', folds = 10, nrepeat = 10)

#perf.diablo.tcga$error.rate  # Lists the different types of error rates

# Plot of the error rates based on weighted vote
plot(perf.diablo.tcga)
perf.diablo.tcga$choice.ncomp$WeightedVote
ncomp <- perf.diablo.tcga$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]

# chunk takes about 2 min to run
set.seed(123) # for reproducibility
test.keepX <- list(mRNA = c(5:9, seq(10, 25, 5)),
                   miRNA = c(5:9, seq(10, 20, 2)),
                   proteomics = c(seq(5, 25, 5)))

tune.diablo.tcga <- tune.block.splsda(X, Y, ncomp = 2, 
                                      test.keepX = test.keepX, design = design,
                                      validation = 'Mfold', folds = 10, nrepeat = 1, 
                                      BPPARAM = BiocParallel::SnowParam(workers = 2),
                                      dist = "centroids.dist")

list.keepX <- tune.diablo.tcga$choice.keepX

# final model
diablo.tcga <- block.splsda(X, Y, ncomp = ncomp, 
                            keepX = list.keepX, design = design)
diablo.tcga$design
# mRNA variables selected on component 1
selectVar(diablo.tcga, block = 'mRNA', comp = 1)
# plotDiablo
plotDiablo(diablo.tcga, ncomp = 1)
# plotIndiv
plotIndiv(diablo.tcga, ind.names = FALSE, ellipse = TRUE, legend = TRUE, 
          title = 'TCGA, DIABLO comp 1 - 2')
# plotLoadings
plotLoadings(diablo.tcga, comp = 1, contrib = 'max', method = 'median')



# real data ----
# 加载必要的包
library(mixOmics)
library(tidyverse)
library(openxlsx)

# set out path
folder_path <- "./03_result/All/"
if (!dir.exists(folder_path)) {
  dir.create(folder_path)
  print(paste("Folder", folder_path, "created."))
} else {
  print(paste("Folder", folder_path, "already exists."))
}

# Data process ----
## transcriptome ----
tpm_df <- read.csv("./01_data/gene_tpm_matrix.csv", row.names = 1)
colnames(tpm_df)
# 去除非编码RNA
source("../../Code_Box/Transcriptome/run_filter_coding_gene.R")
clean_matrix <- filter_coding_genes(tpm_df)
colnames(clean_matrix)
clean_matrix <- clean_matrix[,!colnames(clean_matrix)%in%c("ense","symbol")]
transcriptome_data <- clean_matrix
transcriptome_data <- clean_matrix[,grep("OCI", colnames(clean_matrix))]
transcriptome_data <- transcriptome_data[,-grep('4W', colnames(transcriptome_data))] # 删除4w样本
# filter low expr
min_sample <- ceiling(ncol(transcriptome_data)/2)
transcriptome_data <- transcriptome_data[rowSums(transcriptome_data > 1) >= min_sample,]
# 提取 gene symbol（假设 Row.names 列有 ENSEMBL|symbol）
transcriptome_data$gene <- sapply(strsplit(rownames(transcriptome_data), "\\|"), `[`, 2)
rownames(transcriptome_data) <- transcriptome_data$gene
# 若提示有重复，进行去重处理
# 1. 处理 NA 值
transcriptome_data <- transcriptome_data[!is.na(transcriptome_data$gene), ]  # remove NA rows
# 2. 处理重复基因
transcriptome_data <- aggregate(. ~ gene, data = transcriptome_data, FUN = max)
# 3. 设定 rownames
rownames(transcriptome_data) <- transcriptome_data$gene
transcriptome_data <- transcriptome_data[,!colnames(transcriptome_data)%in%c("gene")]
transcriptome_data <- transcriptome_data[,order(colnames(transcriptome_data))]

## proteome ----
intensity_df1 <- read.csv("./01_data/report.pg_matrix_fill_norma.csv", row.names = 1) 
colnames(intensity_df1)
colnames(intensity_df1) <- gsub("_11|_M2", "", colnames(intensity_df1))
proteome_data <- intensity_df1
proteome_data <- intensity_df1[,-grep("4W", colnames(intensity_df1))] # 删除4w样本
anno_df1 <- read.xlsx("./01_data/data_anno.xlsx")
anno_df1 <- anno_df1[anno_df1$Protein.Group%in%rownames(proteome_data),]
# 提取 gene symbol（假设 Row.names 列有 ENSEMBL|symbol）
proteome_data$gene <- sapply(strsplit(anno_df1$Genes, "\\;"), `[`, 1)
rownames(proteome_data) <- proteome_data$gene
# 若提示有重复，进行去重处理
# 1. 处理 NA 值
proteome_data <- proteome_data[!is.na(proteome_data$gene), ]  # remove NA rows
# 2. 处理重复基因
proteome_data <- aggregate(. ~ gene, data = proteome_data, FUN = max)
# 3. 设定 rownames
rownames(proteome_data) <- proteome_data$gene
proteome_data <- proteome_data[,!colnames(proteome_data)%in%c("gene")]
proteome_data <- proteome_data[,order(colnames(proteome_data))]

## metabolome ----
intensity_df2 <- read.csv("./01_data/meta_intensity_combined.csv", row.names = 1)
colnames(intensity_df2)
colnames(intensity_df2) <- gsub("_11|_M2", "", colnames(intensity_df2))
anno_df2 <- read.xlsx("./01_data/meta_anno_combined.xlsx")
identical(rownames(intensity_df2), anno_df2$Compound_ID)
rownames(intensity_df2) <- anno_df2$Name
metabolome_data <- intensity_df2
metabolome_data <- intensity_df2[,grep("OCI", colnames(intensity_df2))]
metabolome_data <- metabolome_data[,-grep("4W|QC", colnames(metabolome_data))]
metabolome_data <- metabolome_data[,order(colnames(metabolome_data))]

## X & Y ----
# 检查三个组学矩阵的列名是否完全相同
identical(colnames(transcriptome_data), colnames(proteome_data)) && identical(colnames(proteome_data), colnames(metabolome_data))
X <- list(
  transcriptome = transcriptome_data,
  proteome = proteome_data, 
  metabolome = metabolome_data)
lapply(X, dim)  # 应该显示：若干特征  ×  样本
saveRDS(X, file = file.path(folder_path, "X.rds"))

# 对每个组学矩阵进行 log2(x+1) 转换，并转置
X <- lapply(X, function(x) {
  t(log2(x + 1))
})
# 检查维度确认
lapply(X, dim)  # 应该显示：样本 × 若干特征
# 对特征标准化！
X <- lapply(X, function(x) {
  # 确保数据是：行=样本，列=特征
  # 然后对列（特征）进行标准化
  scale(x, center = TRUE, scale = TRUE)  # 对特征标准化 ✅
})
# 创建响应变量（根据IC50分类）
group_info <- read.xlsx("./01_data/group_info.xlsx")
sample_info <- group_info
sample_info <- group_info[grep("OCI", group_info$id),!colnames(sample_info)%in%c("cell")]
sample_info <- sample_info[-grep("4W", sample_info$id), ]
rownames(sample_info) <- sample_info$id
Y <- factor(sample_info$group, levels = c("WT", "Low", "High"))  # 必须是因子类型
summary(Y)
save(Y, file = file.path(folder_path, "Y.rds"))

# PCA分析每个数据块单独看分组情况
pca_transcriptome <- pca(X$transcriptome, ncomp = 2)
plotIndiv(pca_transcriptome, group = Y, legend = TRUE)
pca_proteome <- pca(X$proteome, ncomp = 2)
plotIndiv(pca_proteome, group = Y, legend = TRUE)
pca_metabolome <- pca(X$metabolome, ncomp = 2)
plotIndiv(pca_metabolome, group = Y, legend = TRUE)

# 计算两两组学pls主成分1的相关性
res1.pls <- pls(X$transcriptome, X$proteome, ncomp = 1)
cor(res1.pls$variates$X, res1.pls$variates$Y)

res2.pls <- pls(X$transcriptome, X$metabolome, ncomp = 1)
cor(res2.pls$variates$X, res2.pls$variates$Y)

res3.pls <- pls(X$proteome, X$metabolome, ncomp = 1)
cor(res3.pls$variates$X, res3.pls$variates$Y)
# 若上述相关性高
# 设置较高的权重（0.8–0.9），就是在告诉模型：
# 这两个 block 的潜在成分应该高度相关，也就是在样本空间里的投影方向尽量一致。
# 若上述相关性较弱，则设计矩阵里设置 0.2
# 算法只会弱弱地考虑它们的共性，更多还是让它们各自解释自己的变化。

# Design Matrix ----
design <- matrix(0.9, ncol = length(X), nrow = length(X),
                 dimnames = list(names(X), names(X)))
diag(design) <- 0  # 对角线设为0
design

# numbers of components
diablo.ncomp <- block.plsda(X, Y, ncomp = 5, design = design)

set.seed(123) # For reproducibility
# validation = 'Mfold', folds = 10, nrepeat = 1 选择交叉验证，k小于最小样本数
# validation = 'loo', nrepeat = 1 ,选择留一法则不需要重复
perf.diablo.ncomp = perf(diablo.ncomp, validation = 'Mfold', folds = 7, nrepeat = 50) 

#perf.diablo.tcga$error.rate  # Lists the different types of error rates

# Plot of the error rates based on weighted vote
pdf(file.path(folder_path, "perf.ncomp.pdf"), width = 7, height = 6)
plot(perf.diablo.ncomp)
dev.off()
perf.diablo.ncomp$choice.ncomp$WeightedVote
ncomp <- 3

# chunk takes about 2 min to run
set.seed(123) # for reproducibility
test.keepX <- list(
  transcriptome = c(5:9, seq(10, 25, 5)),
  proteome = c(5:9, seq(10, 20, 2)),
  metabolome = c(seq(5, 25, 5)))

library(BiocParallel)
registered()
tune.diablo <- tune.block.splsda(X, Y, ncomp = ncomp, 
                                 test.keepX = test.keepX, design = design,
                                 validation = 'Mfold', folds = 7, nrepeat = 50,
                                 BPPARAM = BiocParallel::SnowParam(workers = 16),
                                 dist = "mahalanobis.dist") # mahalanobis.dist,centroids.dist,max,dist

list.keepX.recommended <- tune.diablo$choice.keepX
list.keepX.recommended
# 也可以直接自己设置
list.keepX.manual <- list(
  transcriptome = c(7, 5),   # 成分1多些，成分2少些
  proteome = c(7, 5),
  metabolome = c(6, 10)
)
# final model
set.seed(123) # For reproducibility
diablo.final <- block.splsda(X, Y, ncomp = ncomp, 
                            keepX = list.keepX.recommended, design = design)
diablo.final$design
saveRDS(diablo.final, file = file.path(folder_path, "diablo.final.rds"))
# 模型评估
set.seed(123)
perf.diablo <- perf(diablo.final, validation = 'loo', nrepeat = 1)

# 查看错误率
plot(perf.diablo)

# mRNA variables selected on component 1
selectVar(diablo.final, block = 'transcriptome', comp = 1)
# 变量相关性图
plotVar(diablo.final, cutoff = 0.7)
# plotDiablo 组学相关性图, 手动保存
pdf(file.path(folder_path, "omics_corr_comp3.pdf"), width = 6, height = 6)
plotDiablo(diablo.final, ncomp = 3)
dev.off()
#  circos图显示组学间关系
circosPlot(diablo.final, cutoff = 0.6)
# 网络图
network(diablo.final, block.var.names = TRUE, cutoff = 0.1)
# Not run 保存可以输入cytoscape的格式
library(igraph)
myNetwork <- network(diablo.final, blocks = c(1,2,3), cutoff = 0.4)
write_graph(myNetwork$gR, file = file.path(folder_path, "myNetwork.gml"), format = "gml")
# plotIndiv 手动保存
pdf(file.path(folder_path, "Indiv_PCA.pdf"), width = 7, height = 6)
plotIndiv(diablo.final, ind.names = FALSE, ellipse = TRUE, legend = TRUE, 
          title = 'VR, DIABLO comp 1 - 2')
dev.off()
# plotLoadings
pdf(file.path(folder_path, "Feature_loading_comp3.pdf"), width = 10, height = 7)
plotLoadings(diablo.final, comp = 3, contrib = 'max', method = 'median')
dev.off()
# 提取重要特征
important_features <- selectVar(diablo.final, comp = 3, block = NULL)
saveRDS(important_features, file = file.path(folder_path, "important_features_3.rds"))
write.xlsx(important_features, file = file.path(folder_path, "important_features_3.xlsx"))
transcriptome_important <- important_features$transcriptome$value
