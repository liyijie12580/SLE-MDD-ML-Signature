setwd("C:/Users/sunshiny/Downloads/SLEMDD/")
# 加载必要包
library(caret)        # 机器学习框架
library(glmnet)       # LASSO回归
library(randomForest) # 随机森林
library(pROC)         # ROC分析
library(DMwR)         # 类别不平衡处理
library(pheatmap)     # 热图可视化

# 读取数据（假设CSV格式，行是样本，列前N为基因，最后一列为分组）
load("coDEGgene.RData")
MDD<-read.csv("MDD_expr_combat.csv",row.names = 1)
MDD_meta<-read.csv("MDD_meta.csv",row.names = 1)
data <- data.frame(t(MDD))
data_balanced <-data
table(MDD_meta$case)
data_balanced$Group <- factor(MDD_meta$case, levels = c("CTRL", "MDD"))
levels(data_balanced$Group)
# 检查数据质量
# sum(is.na(data))  # 查看缺失值
# data <- na.omit(data)  # 删除含缺失值的样本


# Z-score标准化基因表达数据
# preProc <- preProcess(data[, -ncol(data)], method = c("center", "scale"))
# data_norm <- predict(preProc, data)
# 处理类别不平衡（SMOTE过采样）
# set.seed(123)
# dim(data)
# dim(data_balanced)
# data_balanced <- SMOTE(Group ~ ., data_balanced, perc.over = 200)
# table(data_balanced$Group)
# dim(data_balanced)
# 数据拆分
set.seed(123)
# trainIndex <- createDataPartition(data_balanced$Group, p=0.8, list=FALSE)
# table(data_balanced$Group[trainIndex])
# table(data_balanced$Group[-trainIndex])
# train <- data_balanced[trainIndex, ]
# dim(train)
# test <- data_balanced[-trainIndex, ]
# test <-as.matrix(test[, -ncol(test)])
# dim(test)
###############特征基因选取
x <- Train_set
y <-  Train_class$outcome
#trainIndex <- createDataPartition(data_balanced$Group, p=0.8, list=FALSE)
# 方法1：LASSO回归筛选
fit <- glmnet(x, y, family = "binomial");  # 注意binomial代表二分类
plot(fit,label = T,lwd=2)
plot(fit,xvar = "lambda",label = T,lwd=2)
plot(fit, xvar = "dev", label = TRUE)
cvfit <- cv.glmnet(x, y, family = "binomial", nfolds = 10, type.measure = "auc",
                  keep = TRUE)
plot(cvfit)
#换一个type.measure
# cvfit1 <- cv.glmnet(x, y, family = "binomial", type.measure = "auc")
# plot(cvfit1)
# 根据分析结果筛选特征基因
coef.min <- coef(cvfit, s = cvfit$lambda.min);#lambda.1se为8个
coef.min <- coef(cvfit, s = cvfit$lambda.1se);#lambda.1se为8个
# predict(cvfit, newx = test, s = "lambda.min",type = c("class"))
# assess.glmnet(cvfit, x[-trainIndex, ], newy = y[-trainIndex],s = "lambda.min")
selected_genes <- rownames(coef.min)[as.numeric(coef.min) != 0]
print(selected_genes)
lasso_genes <- selected_genes[2:length(selected_genes)]
print(lasso_genes)
write.csv(lasso_genes,"17feature_lasso.csv")
#ROC曲线是非常重要的模型衡量工具
# rocs <- roc.glmnet(cvfit$fit.preval, newy = y)
# 
# best <- cvfit$index["min",] # 提取AUC最大的lambda值
# plot(rocs[[best]], type = "l") # 画出AUC最大的ROC曲线
# invisible(sapply(rocs, lines, col="grey")) # 把所有的ROC都画出来
# lines(rocs[[best]], lwd = 2,col = "red") # 把AUC最大的标红
###
model<-model$RF
rf_model<-model
# 方法2：随机森林特征重要性
model<- randomForest(x, y, ntree = 1000, importance = TRUE)
plot(rf_model, main = "Random forest", lwd = 2)
# 找出误差最小的点
optionTrees <- which.min(rf_model$err.rate[, 1]);

#画图看种多少树合适
oob.error.data <- data.frame(
  Trees=rep(1:nrow(model$err.rate), times=3),
  Type=rep(c("OOB", "GOOD", "BAD"), each=nrow(model$err.rate)),
  Error=c(model$err.rate[,"OOB"],
          model$err.rate[,"GOOD"],
          model$err.rate[,"BAD"]))
ggplot(data=oob.error.data, aes(x=Trees, y=Error)) +
  geom_line(aes(color=Type))

#确定mtry值
oob.values <- vector(length=20)
for(i in 1:20) {
  temp.model <- randomForest(groups ~ ., data = infr, mtry=i, ntree=1000)
  oob.values[i] <- temp.model$err.rate[nrow(temp.model$err.rate),1]
}
oob.values

#最终确定的结果
model <- randomForest(as.factor(groups)~.,data = infr, mtry = 3, ntree = 600,
                      proximity = TRUE, importance = TRUE)
rfi1 <- importance(model, type = 1)
rfi2 <- importance(model, type = 2)

rfi1 <- rfi1[order(rfi1, decreasing = T),]
rfi2 <- rfi2[order(rfi2, decreasing = T),]

rfif <- intersect(names(rfi1[1:30]), names(rfi2[1:30]))
pdf ("randf.pdf", height = 8, width = 10)
varImpPlot(model)
dev.off()

fina_gene <- intersect(top.features[1:68,]$FeatureName, rownames(factors))
fina_gene <- intersect(fina_gene, rfif)
fina_gene

rf_model<- randomForest(x, y, ntree = optionTrees);
importance <- importance(x = rf_model);
varImpPlot(rf_model, sort=TRUE, main="Variable Importance Plot",n.var=30)
rfGenes <- importance[order(importance[, "MeanDecreaseGini"], decreasing = TRUE), ];
# 重要性评分大于2的基因
# rfGenes <- names(rfGenes[rfGenes>0.4])  ;
# rfGenes
# 也可以是重要性评分最高的5个基因
rf_genes <- names(rfGenes[1:20])     



#SVM-REF筛选特征变量







# 取交集作为最终特征
selected_genes1 <- intersect(lasso_genes,gene)
selected_genes1
selected_genes2 <- intersect(lasso_genes,rf_genes)
selected_genes2
selected_genes3 <- intersect(gene,rf_genes)
selected_genes <- intersect(intersect(gene,rf_genes),lasso_genes)
selected_genes<-c(selected_genes1,selected_genes2,selected_genes3)
selected_genes<-selected_genes[!duplicated(selected_genes)]
cat("关键基因:", selected_genes, "\n")

library(VennDiagram)
#：基础VennDiagram包 --------------------------------------------------
venn.plot <- venn.diagram(
  x = list(
    lasso = lasso_genes,
    rf = rf_genes,
    coDEG = gene
  ),
  filename = NULL,  # 不直接保存到文件
  fill = c("#E69F00", "#56B4E9", "#009E73"),  # 自定义颜色
  alpha = 0.5,      # 透明度
  label.col = "black",
  cex = 1.5,        # 标签大小
  fontfamily = "sans",
  cat.cex = 1.2,    # 类别名称大小
  cat.fontfamily = "sans",
  margin = 0.1
)
# 显示图形
grid.newpage()
grid.draw(venn.plot)

###############################配置模型
# 安装必要包
if (!require("pacman")) install.packages("pacman")
pacman::p_load(caret, glmnet, plsRglm, xgboost, randomForest, mboost, klaR, MASS, kernlab, doParallel, pROC, DMwR)

final_features<-gene

# 并行计算设置
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)



# 数据拆分
set.seed(123)
trainIndex <- createDataPartition(data$Group, p=0.8, list=FALSE)
train <- data[trainIndex, ]
test <- data[-trainIndex, ]


# 模型配置 ----------------------------------------------------------------
# 定义12种算法参数网格
 
models <- list(
  Lasso = list(
    method = "glmnet",
    tuneGrid = expand.grid(alpha = 1, lambda = 10^seq(-3, 0, length=50))
  ),
  
  Ridge = list(
    method = "glmnet",
    tuneGrid = expand.grid(alpha = 0, lambda = 10^seq(-3, 1, length=50))
  ),
  
  # Stepglm = list(
  #   method = "glmStepAIC",
  #   tuneGrid = expand.grid(
  #     nvmax = 3:6,
  #     direction = "both"), trControl = trainControl(
  #       method = "cv",
  #       number = 10,
  #       classProbs = TRUE
  #     ), metric = "ROC"
  # ),
  
  XGBoost = list(
    method = "xgbTree",
    tuneGrid = expand.grid(
      nrounds = 100,
      max_depth = 3:6,
      eta = c(0.01, 0.1),
      gamma = 0,
      colsample_bytree = 0.8,
      min_child_weight = 1,
      subsample = 0.8)
  ),
  
  RF = list(
    method = "rf",
    tuneGrid = data.frame(mtry = c(3,5,7))
  ),
  
  Enet = list(
    method = "glmnet",
    tuneGrid = expand.grid(
      alpha = seq(0.1, 0.9, 0.2),
      lambda = 10^seq(-3, 1, length=20))
  ),
  
  # plsRglm = list(
  #   method = "plsRglm",
  #   tuneGrid = data.frame(nt = 1:5)
  # ),
  
  # GBM = list(
  #   method = "gbm",
  #   tuneGrid = expand.grid(
  #     interaction.depth = 3:5,
  #     n.trees = c(100, 200),
  #     shrinkage = 0.1,
  #     n.minobsinnode = 10)
  # ),
  # 
  # NaiveBayes = list(
  #   method = "nb",
  #   tuneGrid = data.frame(usekernel = c(TRUE, FALSE))
  # ),
  # 
  # LDA = list(
  #   method = "lda",
  #   tuneGrid = data.frame()
  # ),
  
  # glmBoost = list(
  #   method = "glmboost",
  #   tuneGrid = expand.grid(mstop = seq(50, 200, 50))
  # ),
  
  SVM = list(
    method = "svmRadial",
    tuneGrid = expand.grid(
      C = 2^seq(-5, 5, 2),
      sigma = 10^seq(-5, 5, length=5))
  )
)
# 模型训练 ----------------------------------------------------------------
ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 3,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  selectionFunction = "best",
  allowParallel = TRUE
)

# 系统化模型训练
results <- list()
for(model_name in names(models)){
  cat("Training", model_name, "...\n")
  
  # 特殊处理逐步回归
  if(model_name == "Stepglm"){
    model <- train(
      x = data_balanced[, final_features],
      y = data_balanced$Group,
      method = models[[model_name]]$method,
      trControl = ctrl,
      tuneGrid = models[[model_name]]$tuneGrid,
      metric = "ROC",
      trace = FALSE
    )
  } else {
    model <- train(
      Group ~ .,
      data = data_balanced[, c(final_features, "Group")],
      method = models[[model_name]]$method,
      trControl = ctrl,
      tuneGrid = models[[model_name]]$tuneGrid,
      metric = "ROC"
    )
  }
  
  results[[model_name]] <- model
}

# 性能评估 ----------------------------------------------------------------
# 收集交叉验证结果
cv_results <- resamples(results)

# 输出统计摘要
summary(cv_results)

# 可视化比较
dotplot(cv_results, metric = "ROC")

# 保存最佳模型
best_model <- results[which.max(sapply(results, function(x) max(x$results$ROC)))]
# saveRDS(best_model, "best_model.rds")
# 特征重要性分析
ggplot(varImp(best_model), top = 10) + 
  labs(title="Top 10关键特征重要性")

# SHAP值解释
library(fastshap)
shap_values <- explain(best_model, X = data_balanced[, final_features])
autoplot(shap_values)


# 关闭并行
stopCluster(cl)
