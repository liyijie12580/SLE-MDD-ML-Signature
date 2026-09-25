# install.packages("tidyverse")
# install.packages('e1071')
# install.packages("glmnet")
# install.packages("VennDiagram")
# install.packages("ggplot2")
# 并行计算设置
pacman::p_load(caret, glmnet, plsRglm, xgboost, randomForest, mboost, klaR, MASS, kernlab, doParallel, pROC, DMwR)
# 并行计算设置
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
setwd("C:/Users/sunshiny/Downloads/SLEMDD/")
#加载需要的R包
library(tidyverse)
library(glmnet)
source('C:/Users/sunshiny/Documents/R模板/msvmRFE.R')  
library(VennDiagram)
library(e1071)
library(caret)
#SVM-REF算法输入数据
set.seed(2023)
library(e1071)
data <- data.frame(t(MDD))
group=MDD_meta$case
group<-gsub("MDD",2,group)
group<-gsub("CTRL",1,group)
group<-as.numeric(group)
data_balanced <-data.frame(group,data)
colnames(data_balanced)[1]<-"group"
data_balanced<-as.matrix(data_balanced)
input=data_balanced

nfold = 10 #10倍交叉验证
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))

results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=100)
top.features = WriteFeatures(results, input, save=F)
featsweep = lapply(1:100, FeatSweep.wrap, results, input)
save.image(file = "svm-ref.Rdata")

no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

pdf("svm_rfe.pdf", height = 8, width = 10)
PlotErrors(errors, no.info=no.info)
dev.off()
plot(top.features)#这个图也可以保存


########
stopCluster(cl)




input <- data_balanced
#采用五折交叉验证 (k-fold crossValidation）
svmRFE(input, k = 5, halve.above = 100) #分割数据，分配随机数
nfold = 5
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100) #特征选择
top.features = WriteFeatures(results, input, save=F) #查看主要变量
head(top.features)
#把SVM-REF找到的特征保存到文件
write.csv(top.features,"4feature_svm.csv")
# 选前300个变量进行SVM模型构建，选取的变量越多，运行速度越慢！
featsweep = lapply(1:300, FeatSweep.wrap, results, input) #300个变量
# 画图
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
#绘制基于SVM-REF算法的错误率曲线图
pdf("5B_svm-error.pdf",width = 5,height = 5)
PlotErrors(errors, no.info=no.info) #查看错误率
dev.off()
#绘制基于SVM-REF算法的正确率曲线图
#dev.new(width=4, height=4, bg='white')
pdf("6B_svm-accuracy.pdf",width = 5,height = 5)
Plotaccuracy(1-errors,no.info=no.info) #查看准确率
dev.off()
# 图中红色圆圈所在的位置，即错误率最低点
which.min(errors)
#比较lasso和SVM-REF方法一找出的特征变量，画Venn图
(myoverlap <- intersect(lasso_fea, top.features[1:which.min(errors), "FeatureName"])) #交集
#统计交叉基因有多少个
summary(lasso_fea%in%top.features[1:which.min(errors), "FeatureName"])
#绘制venn图
pdf("7C_lasso_SVM_venn.pdf", width = 15, height = 8)
grid.newpage()
venn.plot<-venn.diagram(list(LASSO=lasso_fea,                         SVM_RFE=as.character(top.features[1:which.min(errors),"FeatureName"])), NULL,
                        fill = c("#E31A1C","#E7B800"),
                        alpha = c(0.5,0.5), cex = 4, cat.fontface=3,
                        category.names = c("LASSO", "SVM_RFE"),
                        main = "Overlap")
grid.draw(venn.plot)
dev.off()



# 关闭并行
stopCluster(cl)