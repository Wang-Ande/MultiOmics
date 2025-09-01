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

# data process ----
## transcriptome ----
tpm_df <- read.csv("./01_data/gene_tpm_matrix.csv", row.names = 1)
colnames(tpm_df)
transcriptome_data <- tpm_df[,grep("OCI", colnames(tpm_df))]
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
intensity_df1 <- read.csv("./01_data/OCI_report.pg_matrix_fill_norma.csv", row.names = 1) 
colnames(intensity_df1)
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
metabolome_data <- intensity_df2[,grep("OCI", colnames(intensity_df2))]
metabolome_data <- metabolome_data[,-grep("4W", colnames(metabolome_data))]
metabolome_data <- metabolome_data[,order(colnames(metabolome_data))]


## merge ----
# 检查三个组学矩阵的列名是否完全相同
identical(colnames(transcriptome_data), colnames(proteome_data)) && identical(colnames(proteome_data), colnames(metabolome_data))
X <- list(
  transcriptome = t(log2(transcriptome_data+1)),  # 行是样本，列是基因
  proteome = t(log2(proteome_data+1)),            # 行是样本，列是蛋白
  metabolome = t(log2(metabolome_data+1))         # 行是样本，列是代谢物
)
lapply(X, dim)
# 数据标准化（重要！）
X <- lapply(X, function(x) {
  scale(x, center = TRUE, scale = TRUE)  # 标准化每个数据块
})
# 移除常数变量
X <- lapply(X, function(x) {
  x[, apply(x, 2, sd) > 0.001]  # 移除标准差接近0的变量
})

# 创建响应变量（根据IC50分类）
group_info <- read.xlsx("./01_data/group_info.xlsx")
sample_info <- group_info[grep("OCI", group_info$id),]
sample_info <- sample_info[-grep("4W", sample_info$id), !colnames(sample_info)%in%c("cell")]
rownames(sample_info) <- sample_info$id
Y <- factor(sample_info$group, levels = c("WT", "Low", "High"))  # 必须是因子类型

# 创建全连接设计矩阵（假设所有组学间都有相关关系）
design <- matrix(0.1, ncol = length(X), nrow = length(X),
                 dimnames = list(names(X), names(X)))
diag(design) <- 0  # 对角线设为0
design

# PCA分析每个数据块单独看分组情况
pca_transcriptome <- pca(X$transcriptome, ncomp = 2)
plotIndiv(pca_transcriptome, group = Y, legend = TRUE)
