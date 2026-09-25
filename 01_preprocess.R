rm(list=ls())  #清空环境内变量
options(stringsAsFactors = F)  #避免自动将字符串转换为R语言因子
library(stringr)
setwd("C:/Users/sunshiny/Downloads/SLEMDD/")


#################################Data loading#c
#############GSE138458
library(GEOquery)
library(limma)
library(affy)
gset <- getGEO('GSE138458', destdir="./",
               AnnotGPL = T,     ## 注释文件
               getGPL = T)       ## 平台文件## 平台文件
gset[[1]]
exp<-exprs(gset[[1]])
load("C:/Users/sunshiny/Downloads/SLEMDD/SLE/GSE138458_eSet.Rdata")
## 获取临床信息
cli<-pData(gset[[1]])
meta2<-data.frame(cli$geo_accession,cli$title)
colnames(meta2)<-c("name","case")
meta2$case<-substr(meta2$case,1,2)
table(meta2$case)
meta2$group<-ifelse(meta2$case == "NC", 0, 1)
table(meta2$group)

######################导入SLE数据GSE138458
SLE<-read.table("C:/Users/sunshiny/Downloads/SLEMDD/GSE138458_series_matrix.txt.gz"
                         ,sep="\t",header=T,row.names = 1,comment.char='!')
dim(SLE)
SLE[1:5,1:10]
#colnames(SLE)[apply(SLE,2,function(x)any(is.na(SLE)))]
SLE[is.na(SLE)]<-0
exprSet<-as.matrix(SLE)


library(data.table)
anno<-fread("C:/Users/sunshiny/Downloads/SLEMDD/GPL10558-50081.txt",data.table=F,header=T)
ids<-data.frame(anno$ID,anno$Symbol)
colnames(ids)<-c("probeset_id","gene")

meta<-fread("C:/Users/sunshiny/Downloads/SLEMDD/meta1.txt",data.table=F,header=T)
meta<-t(meta)
colnames(meta)<-meta[1,]
meta<-meta[-1,]
rownames(meta)<-meta[,1]
colnames(meta)
meta<-meta[,c(1,7,10,11,12)]
colnames(meta)<-c("name","source","case","studygroup","sledaigroup")
head(meta)
write.csv(meta, file = "SLE_meta.csv", row.names = T,quote = F)
#基因和symbol转换
#探针与基因是多对一的关系
#对symbol去重然后查看长度，可以看到总共有多少个基因
length(unique(ids$gene))
#统计gene symbol中每个基因出现的次数，然后排序取其尾部
#所得结果可以看出对应探针数目最多的基因
tail(sort(table(ids$gene)))
#统计对应不同探针数目的基因
table(sort(table(ids$gene)))
plot(table(sort(table(ids$gene))))

#统计exprSet中的probe ID与ids中probe ID对应情况
table(rownames(exprSet) %in% ids$probeset_id)


# %in%为二元操作符，返回其左操作数长度的逻辑向量，指示其中的元素是否匹配
exprSet<-exprSet[(rownames(exprSet) %in% ids$probeset_id),]
dim(exprSet)
#将ids中与exprSet不对应的探针删去
#match返回其第一个参数在第二参数中的（第一个）匹配位置的向量
ids<-ids[match(rownames(exprSet),ids$probeset_id),]
dim(ids)
head(ids)
exprSet[1:4,1:4]
#首先根据ids中的gene symbol，在exprSet中进行匹配，找到对应的行，计算平均表达量
#挑选出平均表达量最大一行，将其探针名保留
tmp<-by(exprSet,ids$gene,function(x) rownames(x) [which.max(rowMeans(x))])
tmp[1:20]
#将tmp转为字符型
probes<-as.character(tmp)
#保留在所有样本里面表达最大的那个探针
exprSet<-exprSet[rownames(exprSet) %in% probes, ]
dim(exprSet)
dim(ids)
#注意此时exprSet和ids的行数不同，需要把ids中多余的行去掉
rownames(exprSet)<-ids[match(rownames(exprSet),ids$probeset_id),2]
head(exprSet)
dim(exprSet)
#31427   336
write.csv(exprSet, file = "SLE_expr.csv", row.names = T,quote = F)

######################导入MDD数据GSE98793
MDD<-read.table("C:/Users/sunshiny/Downloads/SLEMDD/GSE98793_series_matrix.txt.gz"
                ,sep="\t",header=T,row.names = 1,comment.char='!')
dim(MDD)
MDD[1:5,1:10]
MDD[is.na(MDD)]<-0
exprSet<-as.matrix(MDD)

library(data.table)
anno<-fread("C:/Users/sunshiny/Downloads/SLEMDD/GPL570-55999.txt",data.table=F,header=T)
ids<-data.frame(anno$ID,anno$`Gene Symbol`)
colnames(ids)<-c("probeset_id","gene")
meta<-fread("C:/Users/sunshiny/Downloads/SLEMDD/metaMDD.txt",data.table=F,header=T)
meta<-t(meta)
colnames(meta)<-meta[1,]
meta<-meta[-1,]
rownames(meta)<-meta[,1]
colnames(meta)
meta<-meta[,c(1,7,9,10,11,12,13,14)]
colnames(meta)<-c("name","source","case","studygroup","gender","age","tissue","batch")
head(meta)
meta<-data.frame(meta)
meta$studygroup<-paste0(meta$source,"_",gsub(": ","_",meta$studygroup))
write.csv(meta, file = "MDD_meta.csv", row.names = T,quote = F)

#基因和symbol转换
#探针与基因是多对一的关系
#对symbol去重然后查看长度，可以看到总共有多少个基因
length(unique(ids$gene))
#统计gene symbol中每个基因出现的次数，然后排序取其尾部
#所得结果可以看出对应探针数目最多的基因
tail(sort(table(ids$gene)))
#统计对应不同探针数目的基因
table(sort(table(ids$gene)))
plot(table(sort(table(ids$gene))))

#统计exprSet中的probe ID与ids中probe ID对应情况
table(rownames(exprSet) %in% ids$probeset_id)
TRUE 
54715

# %in%为二元操作符，返回其左操作数长度的逻辑向量，指示其中的元素是否匹配
exprSet<-exprSet[(rownames(exprSet) %in% ids$probeset_id),]
dim(exprSet)
#将ids中与exprSet不对应的探针删去
#match返回其第一个参数在第二参数中的（第一个）匹配位置的向量
ids<-ids[match(rownames(exprSet),ids$probeset_id),]
dim(ids)
head(ids)
exprSet[1:4,1:4]
#首先根据ids中的gene symbol，在exprSet中进行匹配，找到对应的行，计算平均表达量
#挑选出平均表达量最大一行，将其探针名保留
tmp<-by(exprSet,ids$gene,function(x) rownames(x) [which.max(rowMeans(x))])
tmp[1:20]
#将tmp转为字符型
probes<-as.character(tmp)
#保留在所有样本里面表达最大的那个探针
exprSet<-exprSet[rownames(exprSet) %in% probes, ]
dim(exprSet)
dim(ids)
#注意此时exprSet和ids的行数不同，需要把ids中多余的行去掉
rownames(exprSet)<-ids[match(rownames(exprSet),ids$probeset_id),2]
head(exprSet)
dim(exprSet)
write.csv(exprSet, file = "MDD_expr.csv", row.names = T,quote = F)




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
