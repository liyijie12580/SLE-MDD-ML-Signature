setwd("C:/Users/sunshiny/Downloads/SLEMDD/")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(dplyr)
load("coDEGgene.RData")
MDD<-read.csv("MDD_expr_combat.csv",row.names = 1)
SLE<-read.csv("SLE_expr.csv",row.names = 1)
MDD_meta<-read.csv("MDD_meta.csv",row.names = 1)
SLE_meta<-read.csv("SLE_meta.csv",row.names = 1)
BiocManager::install("GSVA")
library(ggplot2)
library(tinyarray)
library(GSVA)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(ggpubr)
###########################MDD
#输入分组信息，设置好样本对应的状态
mygroup <-MDD_meta[,c(3)]  #使mygroup最终为character 图4
mygene<-gene

#输入表达矩阵 不要取log（标化后十几内可以），不可以有负值和缺失值---所有基因全表达矩阵
expr_matrix<-as.matrix(MDD)  #注意将表达谱的data.frame转化为matrix，否则后续的gsva分析会报错


#####免疫细胞丰度计算
#导入免疫细胞基因集
load("C:/Users/sunshiny/Documents/ssGSEA/肿瘤浸润免疫细胞基因集.rdata")
immune_signatures <- tisidb_cell

#使用msigdbr包获取基因集
# install.packages("msigdbr")
# library(msigdbr)
# glist <- msigdbr(species = 'Homo sapiens', category = 'H')#选择hallmark
# gene_set<-split(x = glist$gene_symbol, f = glist$gs_name)


#GSVA
gsvaParam<-gsvaParam(as.matrix(expr_matrix), geneSets = immune_signatures, minSize = 1)#参数为表达矩阵和基因集
ssgsea_scores <-gsva(gsvaParam)#GSVA
#ssGSEA
ssgseaParam<-ssgseaParam(as.matrix(expr_matrix), geneSets =immune_signatures, minSize = 1)
my_gsva<-gsva(ssgseaParam)


# 差异分析可视化（示例：Treg细胞）
# treg_data <- data.frame(
#   Score = my_gsva["T_cells_CD8", ],
#   Group = mygroup
# )
# 
# ggplot(treg_data, aes(x=Group, y=Score, fill=Group)) +
#   geom_violin(trim=FALSE, alpha=0.6) +
#   geom_boxplot(width=0.2, fill="white") +
#   scale_fill_manual(values=c("#66C2A5", "#FC8D62")) +
#   labs(title="Treg Cell Infiltration", y="Enrichment Score") +
#   theme_bw() +
#   stat_compare_means(method="wilcox.test", label="p.format")

#免疫细胞浸润分组绘图
draw_boxplot(my_gsva,mygroup,color = c("#1d4a9b","#e5171a") )

# #单基因相关性分析---目的基因与免疫细胞相关性
# mygene <- gene #定义你的目的基因
# nc = t(rbind(my_gsva,exp[mygene,]))  ;#将你的目的基因匹配到表达矩阵---行名匹配--注意大小写
# m = rcorr(nc)$r[1:nrow(my_gsva),(ncol(nc)-length(mygene)+1):ncol(nc)]
# 
# ##计算p值
# p = rcorr(nc)$P[1:nrow(my_gsva),(ncol(nc)-length(mygene)+1):ncol(nc)]
# head(p)
# tmp <- matrix(case_when(as.vector(p) < 0.01 ~ "**",
#                         as.vector(p) < 0.05 ~ "*",
#                         TRUE ~ ""), nrow = nrow(p))
# ##绘制热图
# p1 <- pheatmap(t(m),
#                display_numbers =t(tmp),
#                angle_col =45,
#                color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
#                border_color = "white",
#                cellwidth = 20, 
#                cellheight = 20,
#                width = 7, 
#                height=9.1,
#                treeheight_col = 0,
#                treeheight_row = 0)







###########################SLE
#输入分组信息，设置好样本对应的状态
mygroup <-gsub("case/control: ","",SLE_meta[,c(3)])  #使mygroup最终为character 图4
mygene<-gene

#输入表达矩阵 不要取log（标化后十几内可以），不可以有负值和缺失值---所有基因全表达矩阵
expr_matrix<-as.matrix(SLE)  #注意将表达谱的data.frame转化为matrix，否则后续的gsva分析会报错


#####免疫细胞丰度计算
#导入免疫细胞基因集
load("C:/Users/sunshiny/Documents/ssGSEA/TISIDB肿瘤浸润淋巴细胞基因集.rdata")
immune_signatures <- tisidb_cell

#使用msigdbr包获取基因集
# install.packages("msigdbr")
# library(msigdbr)
# glist <- msigdbr(species = 'Homo sapiens', category = 'H')#选择hallmark
# gene_set<-split(x = glist$gene_symbol, f = glist$gs_name)


#GSVA
gsvaParam<-gsvaParam(as.matrix(expr_matrix), geneSets = immune_signatures, minSize = 1)#参数为表达矩阵和基因集
ssgsea_scores <-gsva(gsvaParam)#GSVA
#ssGSEA
ssgseaParam<-ssgseaParam(as.matrix(expr_matrix), geneSets =immune_signatures, minSize = 1)
my_gsva<-gsva(ssgseaParam)


# 差异分析可视化（示例：Treg细胞）
# treg_data <- data.frame(
#   Score = my_gsva["T_cells_CD8", ],
#   Group = mygroup
# )
# 
# ggplot(treg_data, aes(x=Group, y=Score, fill=Group)) +
#   geom_violin(trim=FALSE, alpha=0.6) +
#   geom_boxplot(width=0.2, fill="white") +
#   scale_fill_manual(values=c("#66C2A5", "#FC8D62")) +
#   labs(title="Treg Cell Infiltration", y="Enrichment Score") +
#   theme_bw() +
#   stat_compare_means(method="wilcox.test", label="p.format")

#免疫细胞浸润分组绘图
draw_boxplot(my_gsva,mygroup,color = c("#1d4a9b","#e5171a") )

# #单基因相关性分析---目的基因与免疫细胞相关性
# mygene <- gene #定义你的目的基因
# nc = t(rbind(my_gsva,exp[mygene,]))  ;#将你的目的基因匹配到表达矩阵---行名匹配--注意大小写
# m = rcorr(nc)$r[1:nrow(my_gsva),(ncol(nc)-length(mygene)+1):ncol(nc)]
# 
# ##计算p值
# p = rcorr(nc)$P[1:nrow(my_gsva),(ncol(nc)-length(mygene)+1):ncol(nc)]
# head(p)
# tmp <- matrix(case_when(as.vector(p) < 0.01 ~ "**",
#                         as.vector(p) < 0.05 ~ "*",
#                         TRUE ~ ""), nrow = nrow(p))
# ##绘制热图
# p1 <- pheatmap(t(m),
#                display_numbers =t(tmp),
#                angle_col =45,
#                color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
#                border_color = "white",
#                cellwidth = 20, 
#                cellheight = 20,
#                width = 7, 
#                height=9.1,
#                treeheight_col = 0,
#                treeheight_row = 0)