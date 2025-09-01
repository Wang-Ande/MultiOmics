library(pak)
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
pak("mixOmics")
library(mixOmics)
library(readr)
library(readxl)

## 2.1 pos --------------------------------------------------

data_pos <- read_csv("01_data/4_KNN_imputed_data_pos.csv")
data_pos <- as.data.frame(data_pos)
data_pos_anno <- data_pos[,1:13]
rownames(data_pos) <- data_pos$Metabolite
data_pos <- data_pos[,-1:-13]
## 2.2 neg --------------------------------------------------

data_neg <- read_csv("01_data/4_KNN_imputed_data_neg.csv")
data_neg <- as.data.frame(data_neg)
data_neg_anno <- data_neg[,1:13]
rownames(data_neg) <- data_neg$Metabolite
data_neg <- data_neg[,-1:-13]
## 2.3 data_group --------------------------------------------------

group_all <- read_excel("01_data/group.xlsx")
group_all <- as.data.frame(group_all)
## 2.4 merge ----
# load data 
# Extract training data and name each data frame
# Store as list
data_merge <- read_csv("./01_data/5_annotated_metabolites.csv")
data_merge <- as.data.frame(data_merge)
data_anno <- data_merge[,1:13]
rownames(data_merge) <- data_merge$Metabolite
data_merge <- data_merge[,-1:-13]
data_merge <- as.data.frame(t(data_merge))
# 行为样本，列为基因
group <- group_all [!group_all $time%in%"M1",]
X <- list(Blood = scale(data_merge[rownames(data_merge)%in%group[group$group%in%c("Blood"),"id"],]), 
          BM = scale(data_merge[rownames(data_merge)%in%group[group$group%in%c("BM"),"id"],]), 
          Brain = scale(data_merge[rownames(data_merge)%in%group[group$group%in%c("Brain"),"id"],]),
          Heart = scale(data_merge[rownames(data_merge)%in%group[group$group%in%c("Heart"),"id"],]),
          Kidney = scale(data_merge[rownames(data_merge)%in%group[group$group%in%c("Kidney"),"id"],]),
          Liver = scale(data_merge[rownames(data_merge)%in%group[group$group%in%c("Liver"),"id"],]),
          Lung = scale(data_merge[rownames(data_merge)%in%group[group$group%in%c("Lung"),"id"],]),
          Spleen = scale(data_merge[rownames(data_merge)%in%group[group$group%in%c("Spleen"),"id"],]))
rownames(X$Blood) <- group[group$id%in%rownames(X$Blood),"sample"]
rownames(X$BM) <- group[group$id%in%rownames(X$BM),"sample"]
rownames(X$Brain) <- group[group$id%in%rownames(X$Brain),"sample"]
rownames(X$Heart) <- group[group$id%in%rownames(X$Heart),"sample"]
rownames(X$Kidney) <- group[group$id%in%rownames(X$Kidney),"sample"]
rownames(X$Liver) <- group[group$id%in%rownames(X$Liver),"sample"]
rownames(X$Lung) <- group[group$id%in%rownames(X$Lung),"sample"]
rownames(X$Spleen) <- group[group$id%in%rownames(X$Spleen),"sample"]
missing_samples <- setdiff(rownames(X$Blood), rownames(X$Spleen))
# 步骤3：创建要添加的NA行
na_row <- rep(NA, ncol(X$Spleen))  # 生成一个长度为列数的NA向量
na_rows <- do.call(rbind, lapply(missing_samples, function(sample) {
  matrix(na_row, nrow = 1, dimnames = list(sample))
}))

# 步骤4：将NA行添加到原始24样本矩阵中
combined_mat <- rbind(X$Spleen, na_rows)

# 步骤5：按原始顺序重新排列行
X$Spleen <- combined_mat[rownames(X$Blood), ]
# Outcome
Y <- factor(group[1:19,"time"],levels = c("M1","M6","M12","M24"))
summary(Y)
design <- matrix(0.1, nrow = 8, ncol = 8, 
                 dimnames = list(names(X), names(X)))
diag(design) <- 0  # 自身关联设为0（DIABLO要求）
# 设置重复交叉验证
set.seed(123)  # 保证可重复性
## setup cluster - use SnowParam() on Widnows
library(BiocParallel)
registered()
BPPARAM <- BiocParallel::SnowParam(workers = 16)
tune.result <- tune.block.splsda(
  X = X, 
  Y = Y,    # 无表型信息时设为NULL（探索性分析）
  design = design,
  ncomp = 2,   # 主成分数（建议从2开始）
  test.keepX = list(
    Blood = c(seq(5,10,2)),     # 测试不同变量数量
    BM = c(seq(5,10,2)),
    Brain = c(seq(5,10,2)), 
    Heart = c(seq(5,10,2)), 
    Kidney = c(seq(5,10,2)), 
    Liver = c(seq(5,10,2)), 
    Lung = c(seq(5,10,2)), 
    Spleen = c(seq(5,10,2))
  ),
  BPPARAM = BPPARAM,
  validation = "Mfold",
  folds = 10,    # 5折交叉验证
  dist = "max.dist"
)

# 查看最优变量数
tune.result$choice.keepX
list.keepX = list(Blood = rep(20,2), 
                  BM = rep(20,2), 
                  Brain = rep(20,2),
                  Heart = rep(20,2), 
                  Kidney = rep(20,2), 
                  Liver = rep(20,2), 
                  Lung = rep(20,2), 
                  Spleen = rep(20,2))
diablo.model <- block.splsda(
  X = X,
  Y = Y,          # 无表型时设为NULL
  design = design,
  ncomp = 2,         # 与tune阶段一致
  keepX = list.keepX  # 使用最优变量数
)

pheatmap::pheatmap(result_cor,fontsize = 5,annotation_col = subset(group,select = c(time,sample)))
data_comp <- diablo.model[["variates"]][["Blood"]]
data_comp <- as.data.frame(data_comp)
data_comp$time <- Y
ggplot(data_comp,aes(x = comp1, y = comp2,color = time)) + 
  geom_point(size = 4) + 
  theme_bw()
saveRDS(diablo.model,file = "03_result/mixOmics/diablo.model_10.rds")
plotIndiv(
  diablo.model,
  legend = TRUE,
  title = "Multi-omics Integration"
)
plotIndiv(diablo.model, legend = TRUE, pch = c(15, 16), ellipse = TRUE)
plotIndiv(diablo.model, legend = TRUE, pch = c(15, 16), ellipse = TRUE)
plotVar(diablo.model, 
        var.names = T, 
        style = "graphics",
        legend = T)

plotDiablo(diablo.model, ncomp = 1)  # 展示第一主成分
plotDiablo(diablo.model, ncomp = 2)  # 展示第一主成分
# 查看swab组学中重要性排名前10的变量
selectVar(diablo.model, block = "platelet", comp = 1)
selectVar(diablo.model, block = "platelet", comp = 2)

selectVar(diablo.model, block = "plasma", comp = 1)
selectVar(diablo.model, block = "plasma", comp = 2)

selectVar(diablo.model, block = "swab", comp = 1)
selectVar(diablo.model, block = "swab", comp = 2)

cor.test(data_swab$WARS1, data_platelet$WARS1, method="spearman")
cor.test(data_plasma$MASP2, data_platelet$ICAM2, method="spearman")

circosPlot(diablo.model, cutoff = 0.7, line = TRUE, size.legend = 0.8,size.variables = 0.8,
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)
network(diablo.model, 
        cutoff = 0.7,graph.scale = 0.1,
        # To save the plot, uncomment below line
        #save = 'png', name.save = 'diablo-network'
)
library(reshape2)
library(dplyr)

net <- network(final_diablo, 
               blocks = c(1,2,3), 
               cutoff = 0.7, graph.scale = 0.2,
               plot.graph = FALSE)

cor_mat <- net$M_swab_plasma
cor_long <- melt(cor_mat, varnames = c("swab", "plasma"), value.name = "cor")
cor_long_filtered <- cor_long %>%
  filter(abs(cor) >= 0.7) %>%
  mutate(color = ifelse(cor > 0, "red", "blue"),
         width = abs(cor) * 3) 
edges <- data.frame(
  from = cor_long_filtered$swab,
  to = cor_long_filtered$plasma,
  value = abs(cor_long_filtered$cor),
  color = cor_long_filtered$color,
  width = cor_long_filtered$width
)
write.csv(edges, "network_edges_swab_plasma.csv", row.names = FALSE)
cor_mat <- net$M_swab_platelet
cor_long <- melt(cor_mat, varnames = c("swab", "platelet"), value.name = "cor")
cor_long_filtered <- cor_long %>%
  filter(abs(cor) >= 0.7) %>%
  mutate(color = ifelse(cor > 0, "red", "blue"),
         width = abs(cor) * 3) 
edges <- data.frame(
  from = cor_long_filtered$swab,
  to = cor_long_filtered$platelet,
  value = abs(cor_long_filtered$cor),
  color = cor_long_filtered$color,
  width = cor_long_filtered$width
)
write.csv(edges, "network_edges_swab_platelet.csv", row.names = FALSE)
cor_mat <- net$M_plasma_platelet
cor_long <- melt(cor_mat, varnames = c("plasma", "platelet"), value.name = "cor")
cor_long_filtered <- cor_long %>%
  filter(abs(cor) >= 0.7) %>%
  mutate(color = ifelse(cor > 0, "red", "blue"),
         width = abs(cor) * 3) 
edges <- data.frame(
  from = cor_long_filtered$plasma,
  to = cor_long_filtered$platelet,
  value = abs(cor_long_filtered$cor),
  color = cor_long_filtered$color,
  width = cor_long_filtered$width
)
write.csv(edges, "network_edges_plasma_platelet.csv", row.names = FALSE)
# loading plot ----
# 提取保存loadings
final_diablo <- diablo.model
load_swab <- final_diablo$loadings$swab[, 1]
load_plasma <- final_diablo$loadings$plasma[, 1]
load_platelet <- final_diablo$loadings$platelet[, 1]
load_protein <- load_protein[load_protein != 0]
load_metab <- load_metab[load_metab != 0]
df_protein <- data.frame(
  feature = names(load_protein),
  loading = load_protein,
  block = "Proteome"
)
df_metab <- data.frame(
  feature = names(load_metab),
  loading = load_metab,
  block = "Metabolome"
)
df_all <- rbind(df_protein, df_metab)
write.csv(df_all, "diablo_loadings.csv", quote = FALSE)

# 绘制DIOBLO loadings图
loadings <- read.csv("table/diablo_loadings.csv", row.names = 1)
loadings$group <- "HSPC"
loadings$group[loadings$loading < 0] <- "Mature_cell"
# 分开数据
loadings_protein <- loadings %>% filter(block == "Proteome")
loadings_metabolite <- loadings %>% filter(block == "Metabolome")

# 图1: Protein
p1 <- ggplot(loadings_protein, aes(x = reorder(feature, abs(loading)), y = loading)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_segment(aes(xend = feature, y = 0, yend = loading), size = 0.5, color = "grey40") +
  geom_point(aes(color = group), size = 3) +
  coord_flip() + scale_x_discrete(expand = c(0.05, 0.05)) +
  theme_bw() + scale_color_nejm() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        axis.text = element_text(color = "black")) +
  labs(title = "Contribution Weights of Proteins on Component 1", x = NULL, y = "Loading")
p1

# 图2: Metabolite
p2 <- ggplot(loadings_metabolite, aes(x = reorder(feature, abs(loading)), y = loading)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_segment(aes(xend = feature, y = 0, yend = loading), size = 0.5, color = "grey40") +
  geom_point(aes(color = group), size = 3) +
  coord_flip() + scale_x_discrete(expand = c(0.2, 0.2)) +
  theme_bw() + scale_color_nejm() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        axis.text = element_text(color = "black")) +
  labs(title = "Contribution Weights of Metabolites on Component 1", x = NULL, y = "Loading")

p1 / p2 + plot_layout(heights = c(4, 1))
ggsave("figure/diablo_loadings_self.pdf", width = 6, height = 6)

plotLoadings(diablo.model, comp = 1, contrib = 'max', method = 'median',size.name = 1,size.legend = 0.8)
plotLoadings(diablo.model, comp = 2, contrib = 'max', method = 'median',size.name = 1,size.legend = 0.8)
cimDiablo(diablo.model, color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
          comp = 1, margin=c(8,20), legend.position = "right")
