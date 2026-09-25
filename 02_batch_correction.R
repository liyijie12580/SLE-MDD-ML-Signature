
#######Figure1########################################################
############	PCA ################################################################################################
setwd("C:/Users/sunshiny/Downloads/SLEMDD/")
MDD<-read.csv("MDD_expr.csv",row.names = 1)
SLE<-read.csv("SLE_expr.csv",row.names = 1)
MDD_meta<-read.csv("MDD_meta.csv",row.names = 1)
SLE_meta<-read.csv("SLE_meta.csv",row.names = 1)
table(colnames(MDD) %in% rownames(MDD_meta))
table(colnames(SLE) %in% rownames(SLE_meta))
#检查是否有批次效应

library(sva)
library(ggplot2)
library(FactoMineR) # PCA函数
library(factoextra) # fviz_pca_ind函数
library(pheatmap) # pheatmap函数
library(sva) # 用于combat和combat_seq
library(limma) # 用于removeBatchEffect
exprSet = as.matrix(MDD)
head(exprSet)[1:4,1:4]
# GSM2612096 GSM2612097 GSM2612098 GSM2612099
# DDR1 /// MIR4640   5.850755   5.577231   5.663056   5.596154
# RFC2               7.092003   6.618856   6.487570   6.565388
# PAX8               5.814709   5.643282   5.363979   5.340978
# PTPN21             4.789368   4.803316   4.612934   4.702126
batch = MDD_meta$batch
head(batch)
# "batch: 1" "batch: 1" "batch: 1" "batch: 1" "batch: 1" "batch: 1"



#原始数据可视化
# PCA分析--------------------------------------------------
# PCA图时要求是行名时样本名，列名时探针名，因此此时需要转换
exp = t(exprSet)
# 将matrix转换为data.frame 
exp = as.data.frame(exp)
dim(exp)
exp[1:4, 1:4]
#              1007_s_at  1053_at   117_at   121_at
# GSM71019.CEL 10.115170 5.345168 6.348024 8.901739
# GSM71020.CEL  8.628044 5.063598 6.663625 9.439977
# GSM71021.CEL  8.779235 5.113116 6.465892 9.540738
# GSM71022.CEL  9.248569 5.179410 6.116422 9.254368
# 添加分组信息
ac <- data.frame(
  row.names = rownames(exp), 
  Group = MDD_meta$batch,
  Case = MDD_meta$case)
# 设置图片标题
proj = "MDDbatch"
pro = proj
this_title <- paste0(pro, '_PCA')

# 绘图
dat.pca <- PCA(exp, graph = FALSE)
p.pca <- fviz_pca_ind(dat.pca,
                      # 只显示点而不显示文本，默认都显示
                      geom.ind = "point",  #c( "point", "text" ) / "point",
                      # geom.ind = 'text',
                      # 设定分类种类
                      col.ind = ac$Group, # color by groups
                      # 设定颜色
                      palette = "jco", # jco/Dark2
                      # 添加椭圆
                      addEllipses = TRUE, # Concentration ellipses
                      # 添加图例标题
                      legend.title = "Groups") +
  ggtitle(this_title) + 
  theme(plot.title = element_text(size=16, hjust = 0.5))
p.pca

p.pca <- fviz_pca_ind(dat.pca,
                      # 只显示点而不显示文本，默认都显示
                      geom.ind = "point",  #c( "point", "text" ) / "point",
                      # geom.ind = 'text',
                      # 设定分类种类
                      col.ind = ac$Case, # color by groups
                      # 设定颜色
                      palette = "jco", # jco/Dark2
                      # 添加椭圆
                      #addEllipses = TRUE, # Concentration ellipses
                      # 添加图例标题
                      legend.title = "Groups") +
  ggtitle(this_title) + 
  theme(plot.title = element_text(size=16, hjust = 0.5))
p.pca

###combat函数是基于贝叶斯框架下的调整数据批次效应的函数，主要适用于已经被过滤和标准化后的数据(芯片数据/非count数据)
#combat去批次-可视化
# 官方展示了三种函数设置方式
# combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
# 参数方法 (par.prior=TRUE): 这种方法假设数据的批次效应可以通过正态分布的参数（均值和方差）来调整。mod=NULL: 这表示在调整批次效应时没有考虑其他协变量。这是一个纯粹的批次效应调整，不考虑如癌症状态等其他因素。
# combat_edata2 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)
# 非参数方法 (par.prior=FALSE): 这种方法不假设批次效应的参数分布，而是直接估计和调整每个批次的效应，更加灵活地适应各种数据分布。只调整均值 (mean.only=TRUE):这种调整只考虑批次间的均值差异，不调整方差。适用于批次间方差相似，但均值有偏差的情况。
# combat_edata3 = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, ref.batch=3)
# 包含协变量 (mod=mod): 这里的模型矩阵 mod 包括了癌症状态（cancer），这意味着调整批次效应的同时也会考虑癌症状态的影响，以防止这些生物变量的影响被误当作批次效应调整掉。指定参考批次 (ref.batch=3): 在进行批次效应调整时，将第三个批次作为参照（基准），调整其他批次的数据以匹配这个参照批次的分布。这适用于当某个批次被认为是质量最高或最标准的情况。
pheno<-MDD_meta
# 设置批次信息
batch <- gsub("batch: ","",pheno$batch) # 批次
batch
# 设置生物学分类，告诉函数不要把生物学差异整没了 
# pheno$cancer <- factor(pheno$studygroup, 
#                        levels = c("whole_blood_case_anxiety_no",
#                                   "whole_blood_case_anxiety_yes",
#                                   "whole_blood_control_anxiety_no"))
pheno$cancer <- factor(pheno$case, 
                       levels = c("CTRL","MDD"))

mod <- model.matrix(~as.factor(cancer), data=pheno)

expr_combat <- ComBat(dat = exprSet, batch = batch, mod = mod,par.prior=TRUE, ref.batch=1)


# PCA分析--------------------------------------------------
# PCA图时要求是行名时样本名，列名时探针名，因此此时需要转换
exp = t(expr_combat)
# 将matrix转换为data.frame 
exp = as.data.frame(exp)
dim(exp)
exp[1:4, 1:4]
# DDR1 /// MIR4640     RFC2     PAX8   PTPN21
# GSM2612096         5.850755 7.092003 5.814709 4.789368
# GSM2612097         5.577231 6.618856 5.643282 4.803316
# GSM2612098         5.663056 6.487570 5.363979 4.612934
# GSM2612099         5.596154 6.565388 5.340978 4.702126
# 添加分组信息
ac <- data.frame(
  row.names = rownames(exp), 
  Group = MDD_meta$batch,
  Case = MDD_meta$case)
# 设置图片标题
proj = "MDDbatch_combat"
pro = proj
this_title <- paste0(pro, '_PCA')

# 绘图
dat.pca <- PCA(exp, graph = FALSE)
p.pca <- fviz_pca_ind(dat.pca,
                      # 只显示点而不显示文本，默认都显示
                      geom.ind = "point",  #c( "point", "text" ) / "point",
                      # geom.ind = 'text',
                      # 设定分类种类
                      col.ind = ac$Group , # color by groups
                      # 设定颜色
                      palette = "jco", # jco/Dark2
                      # 添加椭圆
                      addEllipses = TRUE, # Concentration ellipses
                      # 添加图例标题
                      legend.title = "Groups") +
  ggtitle(this_title) + 
  theme(plot.title = element_text(size=16, hjust = 0.5))
p.pca

p.pca <- fviz_pca_ind(dat.pca,
                      # 只显示点而不显示文本，默认都显示
                      geom.ind = "point",  #c( "point", "text" ) / "point",
                      # geom.ind = 'text',
                      # 设定分类种类
                      col.ind = ac$Case, # color by groups
                      # 设定颜色
                      palette = "jco", # jco/Dark2
                      # 添加椭圆
                      addEllipses = TRUE, # Concentration ellipses
                      # 添加图例标题
                      legend.title = "Groups") +
  ggtitle(this_title) + 
  theme(plot.title = element_text(size=16, hjust = 0.5))
p.pca

write.csv(expr_combat, file = "MDD_expr_combat.csv", row.names = T,quote = F)