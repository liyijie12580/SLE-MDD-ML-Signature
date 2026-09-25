

#####################################################################Figure5#####################################################################
library(openxlsx)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)
setwd("C:/Users/sunshiny/Downloads/SLEMDD/R代码/Figure5/model/")
setwd("C:/Users/sunshiny/Downloads/SLEMDD/R代码/Figure5/model/")
load("C:/Users/sunshiny/Downloads/SLEMDD/coDEGgene.RData")
##########################################数据预处理########################################################
#####################train

MDDtrain <- read.csv("MDD_expr_combat.csv", row.names = 1,check.names = F,stringsAsFactors = F)
MDDtrain <-MDDtrain[gene,]
MDDtrain_class <- read.csv("MDD_meta.csv", row.names = 1,check.names = F,stringsAsFactors = F)
MDDtrain_class$group<-ifelse(MDDtrain_class$case == "CTRL", 0, 1)
MDDtrain_class<-data.frame( MDDtrain_class$group,MDDtrain_class$gender,
                            MDDtrain_class$studygroup,MDDtrain_class$age,MDDtrain_class$name,
                            row.names = MDDtrain_class$name)
colnames(MDDtrain_class)<-c("outcome","gender","age","studygroup","name")
write.table(MDDtrain_class , file = "Training_class_more_info.txt",sep="\t", row.names = T,quote = F)
MDDtrain_class<-data.frame(MDDtrain_class$outcome,row.names = MDDtrain_class$name)
colnames(MDDtrain_class)<-c("outcome")
MDDtrain_class$outcome<-as.numeric(MDDtrain_class$outcome)
MDDtrain <-MDDtrain[,rownames(MDDtrain_class)]
write.table(MDDtrain, file = "Training_expr.txt", sep="\t",row.names = T,quote = F)
write.table(MDDtrain_class , file = "Training_class.txt",sep="\t", row.names = T,quote = F)


######################test
setwd("C:/Users/sunshiny/Downloads/SLEMDD/R代码/Figure5/model//")
MDDtest <- read.csv("MDDtest_expr.csv", row.names = 1,check.names = F,stringsAsFactors = F)
MDDtest <-MDDtest[intersect(gene,rownames(MDDtest)),]
MDDtest_class <- read.csv("MDDtest_meta.csv", row.names = 1,check.names = F,stringsAsFactors = F)
MDDtest_class<-data.frame(MDDtest_class$group,MDDtest_class$batch,MDDtest_class$gender, 
                          MDDtest_class$Age,
                          row.names = rownames(MDDtest_class))
colnames(MDDtest_class)<-c("outcome","Cohort","gender","age")
MDDtest_class$outcome<-as.numeric(MDDtest_class$outcome)
write.table(MDDtest, file = "Testing_expr.txt", sep="\t",row.names = T,quote = F)
write.table(MDDtest_class , file = "Testing_class.txt",sep="\t", row.names = T,quote = F)



#############加载程序

#################################################正式开始####################################################

code.path="C:/Users/sunshiny/Downloads/SLEMDD/R代码/Figure5/model/Codes"
data.path="C:/Users/sunshiny/Downloads/SLEMDD/R代码/Figure5/model/InputData/"
res.path="C:/Users/sunshiny/Downloads/SLEMDD/R代码/Figure5/model/Result/"
load("C:/Users/sunshiny/Downloads/SLEMDD/coDEGgene.RData")
# 加载模型训练以及模型评估的脚本
source(file.path(code.path, "ML.R"))
# 选择最后生成的模型类型：panML代表生成由不同算法构建的模型； multiLogistic表示抽取其他模型所用到的变量并建立多变量logistic模型
FinalModel <- c("panML", "multiLogistic")[2]

## Training Cohort ---------------------------------------------------------
Train_expr <- read.table(file.path(data.path, "Training_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_class <- read.table(file.path(data.path, "Training_class.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_class$outcome <- as.numeric(Train_class$outcome)
# 剔除训练集标签NA
#Train_class <- Train_class[!is.na(Train_class$outcome), ]
# 提取训练集的共有样本
comsam <- intersect(rownames(Train_class), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_class <- Train_class[comsam,,drop = F]
all(rownames(Train_expr)== gene)

## Validation Cohort -------------------------------------------------------
Test_expr <- read.table(file.path(data.path, "Testing_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Test_class <- read.table(file.path(data.path, "Testing_class.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

# 1. 先过滤匹配样本
comsam <- intersect(rownames(Test_class), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_class <- Test_class[comsam,,drop = F]
Test_expr[Test_expr==0]<-0.001#########不同测序集捕捉基因不同
# 2. 剔除outcome为NA的样本，避免分层后单类别
Test_class <- Test_class[!is.na(Test_class$outcome), ]
keep_sam <- intersect(rownames(Test_class), colnames(Test_expr))
Test_expr <- Test_expr[, keep_sam]
Test_class <- Test_class[keep_sam, , drop = F]

# 3. 拆分队列与性别
library(tidyr)
# Test_class <- Test_class %>%
#   separate(
#     col = Cohort,
#     into = c("Cohort", "gender2"),
#     sep = "_",
#     remove = TRUE,
#     fill = "right"
#   )
table(Test_class$Cohort)
Test_class$gender2 <- Test_class$Cohort
table(Test_class$gender2)
Test_class$outcome <- as.numeric(Test_class$outcome)
all(rownames(Test_expr)== gene)
all(colnames(Test_expr)== rownames(Test_class))

##########scale
# Train_set = scaleData(data = Train_expr, centerFlags = T, scaleFlags = T) 
# names(x = split(as.data.frame(Test_expr), f = Test_class$Cohort))
# Test_set = scaleData(data = Test_expr, cohort = Test_class$Cohort, centerFlags = T, scaleFlags = T)
# summary(apply(Train_set, 2, var))
# summary(apply(Test_set, 2, var))
# # lapply(split(as.data.frame(Test_set), Test_class$Cohort), function(x) summary(apply(x, 2, var))) # 测试scale结果

#============================================================
  # 输入：行=基因，列=样本
  # 输出：Train_set（行=样本，列=基因），Test_set（行=样本，列=基因）
  # ============================================================

# ---------- 0. 原始数据 ----------
# Train_expr: 行=基因，列=样本
# Test_expr:  行=基因，列=样本

# ---------- 1. 取交集基因 ----------
common_gene <- intersect(rownames(Train_expr), rownames(Test_expr))

tr_mat <- Train_expr[common_gene, , drop = FALSE]  # 基因×样本
te_mat <- Test_expr[common_gene, , drop = FALSE]   # 基因×样本

# ---------- 2. 训练集 fit mean/sd（按基因行计算）----------
gene_mean <- apply(tr_mat, 1, mean, na.rm = TRUE)
gene_sd   <- apply(tr_mat, 1, sd,   na.rm = TRUE)

# 过滤 sd≈0
keep <- gene_sd > 1e-6
tr_final <- tr_mat[keep, , drop = FALSE]
te_final <- te_mat[keep, , drop = FALSE]
mean_use <- gene_mean[keep]
sd_use   <- gene_sd[keep]

# ---------- 3. Z-score 标准化（按基因）----------
# 训练集
Train_z <- (tr_final - mean_use) / sd_use

# 测试集：用训练集参数！
Test_z  <- (te_final - mean_use) / sd_use

# ---------- 4. 转置：行=样本，列=基因 ----------
Train_set <- t(Train_z)
Test_set  <- t(Test_z)

# ---------- 5. 确认 ----------
cat("Train_set:", dim(Train_set), "\n")  # 样本×基因
cat("Test_set: ", dim(Test_set),  "\n")
lapply(
  split(as.data.frame(Test_set), Test_class$Cohort),
  function(x) summary(apply(x, 2, var))
)
# 
# # 提取相同基因
# comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
# all(comgene== gene)
# Train_expr <- t(Train_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因
# Test_expr <- t(Test_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因
# Model training and validation -------------------------------------------

## method list --------------------------------------------------------
# 此处记录需要运行的模型，格式为：算法1名称[算法参数]+算法2名称[算法参数]
# 目前仅有Stepglm和Enet支持输入算法参数
methods <- read.xlsx(file.path(code.path, "methods.xlsx"), startRow = 2)
methods <- methods$Model
methods <- gsub("-| ", "", methods)

############### Train the model --------------------------------------------------------

classVar = "outcome" # 设置所要预测的变量名（仅支持[0,1]二元变量格式）
min.selected.var = 3 # 设置模型最少纳入的变量数

## Pre-training 将各方法所用到的变量筛选过程汇总，以减少计算量
Variable = colnames(Train_set)
preTrain.method =  strsplit(methods, "\\+") # 检视所有方法，分析各方法是否需要进行变量预筛选(pre-training)
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1]) # 删除各方法用于构建分类模型的算法，保留用于变量筛选的算法
preTrain.method = unique(unlist(preTrain.method)) # 汇总所有变量筛选算法，去除重复计算
preTrain.method
# "Lasso"             "glmBoost"          "RF"                "Stepglm[both]"     "Stepglm[backward]"

preTrain.var <- list() # 用于保存各算法筛选的变量
set.seed(seed = 666) # 设置建模种子，使得结果可重复
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method, # 变量筛选所需要的机器学习方法
                                 Train_set = Train_set, # 训练集有潜在预测价值的变量
                                 Train_label = Train_class, # 训练集分类标签
                                 mode = "Variable",       # 运行模式，Variable(筛选变量)和Model(获取模型)
                                 classVar = classVar) # 用于训练的分类变量，必须出现在Train_class中
}
preTrain.var[["simple"]] <- colnames(Train_set)# 记录未经筛选的变量集（以便后续代码撰写），可视为使用simple方法（无筛选功能）的变量筛选结果




## Model training
model <- list() # 用于保存各模型的所有信息
set.seed(seed = 666) # 设置建模种子，使得结果可重复
Train_set_bk = Train_set # RunML有一个函数(plsRglm)无法正常传参，需要对训练集数据进行存档备份
for (method in methods){
  cat(match(method, methods), ":", method, "\n")
  method_name = method # 本轮算法名称
  method <- strsplit(method, "\\+")[[1]] # 各步骤算法名称
  
  if (length(method) == 1) method <- c("simple", method) # 如果本方法没有预筛选变量，则认为本方法使用simple方法进行了变量筛选
  
  Variable = preTrain.var[[method[1]]] # 根据方法名称的第一个值，调用先前变量筛选的结果
  Train_set = Train_set_bk[, Variable]   # 对训练集取子集，因为有一个算法原作者写的有点问题，无法正常传参
  Train_label = Train_class            # 所以此处需要修改变量名称，以免函数错误调用对象
  result <- tryCatch( model[[method_name]] <- RunML(method = method[2],        # 根据方法名称第二个值，调用构建的函数分类模型
                                                    Train_set = Train_set,     # 训练集有潜在预测价值的变量
                                                    Train_label = Train_label, # 训练集分类标签
                                                    mode = "Model",            # 运行模式，Variable(筛选变量)和Model(获取模型)
                                                    classVar = classVar)   ,    # 用于训练的分类变量，必须出现在Train_class中
                      error = function(e) {        # 如果发生错误
                        return(NA)              # 返回NA
                      }
  )
  # 如果最终模型纳入的变量数小于预先设定的下限，则认为该算法输出的结果是无意义的
  if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}
Train_set = Train_set_bk; rm(Train_set_bk) # 将数据还原，并移除备份
saveRDS(model, file.path(res.path, "model.rds")) # 报错各模型的所有中间过程

if (FinalModel == "multiLogistic"){
  logisticmodel <- lapply(model, function(fit){ # 根据各算法最终获得的变量，构建多变量Logistic模型，从而以Logistic回归系数和特征表达计算单样本分类概率
    tmp <- glm(formula = Train_class[[classVar]] ~ .,
               family = "binomial", 
               data = as.data.frame(Train_set[, ExtractVar(fit)]))
    tmp$subFeature <- ExtractVar(fit) # 提取当Logistic模型最终使用的预测变量
    return(tmp)
  })
}
saveRDS(logisticmodel, file.path(res.path, "logisticmodel.rds")) # 保存最终以多变量Logistic模型

## Evaluate the model -----------------------------------------------------

# 读取已保存的模型列表
model <- readRDS(file.path(res.path, "model.rds"))
model <- readRDS(file.path(res.path, "logisticmodel.rds")) # 若希望使用多变量保存最终以多变量Logistic模型计算得分，请运行此行
methodsValid <- names(model)

# 根据给定表达量计算样本风险评分
# 预测概率
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalPredictScore(fit = model[[method]], 
                                       new_data = rbind.data.frame(Train_set,Test_set)) # 2.0更新
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat, file.path(res.path, "RS_mat.txt"),sep = "\t", row.names = T, col.names = NA, quote = F) # 输出风险评分文件

# 根据给定表达量预测分类
Class_list <- list()
for (method in methodsValid){
  Class_list[[method]] <- PredictClass(fit = model[[method]], 
                                       new_data = rbind.data.frame(Train_set,Test_set)) # 2.0更新
}
Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
#Class_mat <- cbind.data.frame(Test_class, Class_mat[rownames(Class_mat),]) # 若要合并测试集本身的样本信息文件可运行此行
write.table(Class_mat, file.path(res.path, "Class_mat.txt"), # 测试集经过算法预测出的二分类结果
            sep = "\t", row.names = T, col.names = NA, quote = F)

# 提取所筛选的变量（列表格式）
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}

# 提取所筛选的变量（数据框格式）
fea_df <- lapply(model, function(fit){
  data.frame(ExtractVar(fit))
})
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, file.path(res.path, "fea_df.txt"), # 两列，包含算法以及算法所筛选出的变量
            sep = "\t", row.names = F, col.names = T, quote = F)

# 对各模型计算C-index
AUC_list <- list()

for (method in methodsValid){
  AUC_list[[method]] <- RunEval(fit = model[[method]],     # 分类预测模型
                                Test_set = Test_set,      # 测试集预测变量，应当包含训练集中所有的变量，否则会报错
                                Test_label = Test_class,   # 训练集分类数据，应当包含训练集中所有的变量，否则会报错
                                Train_set = Train_set,    # 若需要同时评估训练集，则给出训练集表达谱，否则置NULL
                                Train_label = Train_class, # 若需要同时评估训练集，则给出训练集分类数据，否则置NULL
                                #Train_name = "Training",       # 若需要同时评估训练集，可给出训练集的标签，否则按“Training”处理
                                cohortVar = "Cohort",      # 重要：用于指定队列的变量，该列必须存在且指定[默认为“Cohort”]，否则会报错
                                classVar = classVar)       # 用于评估的二元分类变量，必须出现在Test_class中
}
AUC_mat <- do.call(rbind, AUC_list)
write.table(AUC_mat, file.path(res.path, "AUC_mat.txt"),
            sep = "\t", row.names = T, col.names = T, quote = F)

# Plot --------------------------------------------------------------------

AUC_mat <- read.table(file.path(res.path, "AUC_mat.txt"),sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
avg_AUC <- apply(AUC_mat, 1, mean)           # 计算每种算法在所有队列中平均AUC
avg_AUC <- sort(avg_AUC, decreasing = T)     # 对各算法AUC由高到低排序
AUC_mat <- AUC_mat[names(avg_AUC), ]      # 对AUC矩阵排序
fea_sel <- fea_list[[rownames(AUC_mat)[1]]] # 最优模型（测试集AUC均值最大）所筛选的特征
avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3)) # 保留三位小数

if(ncol(AUC_mat) < 3) { # 如果用于绘图的队列小于3个
  CohortCol <- c("red","blue") # 则给出两个颜色即可（可自行替换颜色）
} else { # 否则通过brewer.pal赋予超过3个队列的颜色
  CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired") # 设置队列颜色
}
names(CohortCol) <- colnames(AUC_mat)

cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(AUC_mat, # 主矩阵
                    avg_AUC, # 侧边柱状图
                    CohortCol, "steelblue", # 列标签颜色，右侧柱状图颜色
                    cellwidth = cellwidth, cellheight = cellheight, # 热图每个色块的尺寸
                    cluster_columns = F, cluster_rows = F) # 是否对行列进行聚类

pdf(file.path(res.path, "AUC.pdf"), width = cellwidth * ncol(AUC_mat) + 7, height = cellheight * nrow(AUC_mat) * 0.45)
draw(hm)
invisible(dev.off())


#Best
logisticmodel$`glmBoost+Stepglm[both]`
logisticmodel$`glmBoost+Stepglm[backward]`

# 
# > logisticmodel$`glmBoost+Stepglm[both]`
# 
# Call:  glm(formula = Train_class[[classVar]] ~ ., family = "binomial", 
#            data = as.data.frame(Train_set[, ExtractVar(fit)]))
# 
# Coefficients:
#   (Intercept)        KLRB1        IFI27      CEACAM6           HP      C9orf66        WDR41  
# 1.0170      -0.4025       0.4771       0.5279       0.3571       0.2676       0.3762  
# TSHZ3  
# 0.4244  
# 
# Degrees of Freedom: 191 Total (i.e. Null);  184 Residual
# Null Deviance:	    244.4 
# Residual Deviance: 194.7 	AIC: 210.7
# > logisticmodel$`glmBoost+Stepglm[backward]`
# 
# Call:  glm(formula = Train_class[[classVar]] ~ ., family = "binomial", 
#            data = as.data.frame(Train_set[, ExtractVar(fit)]))
# 
# Coefficients:
#   (Intercept)        KLRB1        IFI27      CEACAM6           HP      C9orf66        WDR41  
# 1.0170      -0.4025       0.4771       0.5279       0.3571       0.2676       0.3762  
# TSHZ3  
# 0.4244  
# 
# Degrees of Freedom: 191 Total (i.e. Null);  184 Residual
# Null Deviance:	    244.4 
# Residual Deviance: 194.7 	AIC: 210.7
#####################矫正曲线绘图
final_genes <- names(coef(logisticmodel$`glmBoost+Stepglm[both]`))[-1] #去掉截距
final_genes

save.image("AUC_logisticmodel.RData")
