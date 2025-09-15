# ------------------------- 2. 基因标注表达热图 --------------------------------
library(readxl)
library(openxlsx)
library(ggplot2)
library(readr)
library(dplyr)
library(circlize)
library(mixOmics)
library(ComplexHeatmap)

# set out path
folder_path <- "./03_result/WT-High/"
if (!dir.exists(folder_path)) {
  dir.create(folder_path)
  print(paste("Folder", folder_path, "created."))
} else {
  print(paste("Folder", folder_path, "already exists."))
}

# data process
diablo.final <- readRDS("./03_result/WT-High/diablo.final.rds")
important_features <- selectVar(diablo.final, comp = 1, block = NULL)
# 提取多组学特征
select_important <- important_features$proteome$name
select_sample <- diablo.final$names$sample

# expr input
expr <- read.csv("./01_data/report.pg_matrix_fill_norma.csv",row.names = 1)
colnames(expr) <- gsub("_11|_M2", "", colnames(expr))
anno <- read.xlsx("./01_data/data_anno.xlsx",rowNames = TRUE)
expr_anno <- merge(expr, anno, by.x = 0, by.y = 0, all.x = TRUE)

# 转换基因名 
y <- expr_anno$Genes
gene1 <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
expr_anno$gene <- gene1

# 选择目标特征
targeted_prote <- expr_anno[expr_anno$gene%in%select_important,]
rownames(targeted_prote) <- targeted_prote$gene
targeted_prote <- targeted_prote[,colnames(targeted_prote)%in%select_sample]
targeted_prote <- targeted_prote[,order(colnames(targeted_prote))]
colnames(targeted_prote)
Sample_order <- c("MOLM13_WT_1", "MOLM13_WT_2", "MOLM13_WT_3", 
                  "MOLM13_6W_1", "MOLM13_6W_2", "MOLM13_6W_3",
                  "MV4_WT_1", "MV4_WT_2", "MV4_WT_3",
                  "MV4_6W_1", "MV4_6W_2", "MV4_6W_3",
                  "OCI_WT_1", "OCI_WT_2", "OCI_WT_3",
                  "OCI_6W_1", "OCI_6W_2", "OCI_6W_3")
# Sample_order <- gsub("OCI","MV4_11",Sample_order)
targeted_prote <- targeted_prote[,Sample_order]

expr_matrix <- targeted_prote

# 定义 Cell line 分组
cell_lines <- case_when(
  grepl("MOLM13", colnames(expr_matrix)) ~ "MOLM13",
  grepl("MV4", colnames(expr_matrix)) ~ "MV4_11",
  grepl("OCI", colnames(expr_matrix)) ~ "OCI_AML2"
)
cell_lines <- factor(cell_lines, levels = c("MOLM13", "MV4_11", "OCI_AML2"))

# 为 Cell line 设置颜色
cell_line_colors <- c(
  MOLM13 = "#66C2A5",
  MV4_11 = "#8DA0CB",
  OCI_AML2 = "#FC8D62"
)

# 创建分组信息
sample_groups <- case_when(
  grepl("WT", colnames(expr_matrix)) ~ "WT",
  grepl("6W", colnames(expr_matrix)) ~ "High"
)

# 转换为因子，设置顺序
sample_groups <- factor(sample_groups, levels = c("WT", "High"))

# 创建颜色映射
group_colors <- c(WT = "#6388B4", High = "#EF6F6A")

# 创建注释列时可以选择要几个层，要么同时线上Cell line和Group 要么只显示其中之一
# 1.创建包含两个分组注释的列注释对象（Cell line 在上面）
ha_col <- HeatmapAnnotation(
  Cell_line = cell_lines,
  Group = sample_groups,
  col = list(
    Cell_line = cell_line_colors,
    Group = group_colors
  ),
  annotation_name_side = "right",
  annotation_legend_param = list(
    Cell_line = list(title = "Cell line", title_gp = gpar(fontface = "bold", fontsize = 10), labels_gp = gpar(fontsize = 8)),
    Group = list(title = "Group", title_gp = gpar(fontface = "bold", fontsize = 10), labels_gp = gpar(fontsize = 8))
  ),
  annotation_height = unit(c(4, 4), "mm")  # 控制两个注释行的高度
)

# 2.只创建 group 的列注释对象
ha_col <- HeatmapAnnotation(
  Group = sample_groups,
  col = list(Group = group_colors),
  annotation_name_side = "right",
  annotation_legend_param = list(title = "Group", 
                                 title_gp = gpar(fontface = "bold",fontsize = 10),
                                 labels_gp = gpar(fontsize = 8))
)

# scale
log_expr <- log2(expr_matrix + 1)
scaled_expr <- t(scale(t(log_expr)))  # 每行（每个基因）标准化

# plot
ht <- Heatmap(
  scaled_expr,
  name = "Z-score",
  top_annotation = ha_col,           # 加上分组注释
  col = circlize::colorRamp2(c(-2, 0, 2), colors = c("#1E90FF", "white", "#FF4500")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(title = "Z-score", 
                              title_gp = gpar(fontface = "bold",fontsize = 10),
                              legend_height = unit(3, "cm"),   # 控制图例高度
                              grid_width = unit(0.35, "cm"),    # 色块宽度
                              labels_gp = gpar(fontsize = 8)   # 标签字体
  ))
ht <- draw(ht,
           heatmap_legend_side = "right",
           annotation_legend_side = "right",
           merge_legend = TRUE,
           padding = unit(c(2, 12, 2, 2), "mm"),  # 调整边距
           legend_gap = unit(6, "mm")) # 图例边距
# res output
cairo_pdf(file.path(folder_path, "Comp1_Proteome_features_heatmap.pdf") , width = 7, height = 2)
ht
dev.off()
