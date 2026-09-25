# 加载必要包
library(caret)        # 机器学习框架
library(glmnet)       # LASSO回归
library(randomForest) # 随机森林
library(pROC)         # ROC分析
library(DMwR)         # 类别不平衡处理
library(pheatmap)     # 热图可视化
data <- data.frame(t(MDD))
data_balanced <-data.frame(MDD_meta$case,data)
table(data_balanced$case)
data_balanced$Group <- factor(MDD_meta$case, levels = c("CTRL", "MDD"))
levels(data_balanced$Group)
# 读取数据（假设CSV格式，行是样本，列前N为基因，最后一列为分组）
load("coDEGgene.RData")
MDD<-read.csv("MDD_expr_combat.csv",row.names = 1)
MDD_meta<-read.csv("MDD_meta.csv",row.names = 1)
data <- data.frame(t(MDD))
table(MDD_meta$case)
data$Group <- factor(MDD_meta$case, levels = c("CTRL", "MDD"))
levels(data$Group)
# 检查数据质量
sum(is.na(data))  # 查看缺失值
data <- na.omit(data)  # 删除含缺失值的样本

# Z-score标准化基因表达数据
preProc <- preProcess(data[, -ncol(data)], method = c("center", "scale"))
data_norm <- predict(preProc, data)

# 处理类别不平衡（SMOTE过采样）
set.seed(123)
data_balanced <- SMOTE(Group ~ ., data_norm, perc.over = 200)
table(data_balanced$Group)

###############特征基因选取
x <- as.matrix(data_balanced[, -ncol(data_balanced)])
y <- data_balanced$Group
# 方法1：LASSO回归筛选
fit <- glmnet(x, y, family = "binomial");  # 注意binomial代表二分类
cvfit <- cv.glmnet(x, y, family = "binomial", nfolds = 20);# 检查最优λ下的系数
plot(cvfit)
# 根据分析结果筛选特征基因
coef <- coef(fit, s = cvfit$lambda.min);
index <- which(coef != 0);
lasso_genes <- row.names(coef)[index];
lasso_genes <- lassoGene[-1];
###128

# 方法2：随机森林特征重要性
rf_model <- randomForest(x, y, ntree = 1000, importance = TRUE)
plot(rf_model, main = "Random forest", lwd = 2)
# 找出误差最小的点
optionTrees <- which.min(rf_model$err.rate[, 1]);
rf_model<- randomForest(x, y, ntree = optionTrees);
importance <- importance(x = rf_model);
varImpPlot(rf_model, sort=TRUE, main="Variable Importance Plot",n.var=100)
rfGenes <- importance[order(importance[, "MeanDecreaseGini"], decreasing = TRUE), ];
# 重要性评分大于2的基因
# rfGenes <- names(rfGenes[rfGenes>0.4])  ;
# rfGenes
# 也可以是重要性评分最高的5个基因
rf_genes <- names(rfGenes[1:20])          

# 取交集作为最终特征
selected_genes <- intersect(lasso_genes,gene)
selected_genes <- intersect(lasso_genes,rf_genes)
selected_genes <- intersect(gene,rf_genes)
selected_genes <- intersect(intersect(gene,rf_genes),lasso_genes)
cat("关键基因:", selected_genes, "\n")
selected_genes
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


###########################
# LASSO特征初筛
cv_lasso <- cv.glmnet(as.matrix(data_balanced[, -ncol(data_balanced)]), 
                      data_balanced$Group, 
                      family = "binomial", alpha = 1)
lasso_features <- names(which(coef(cv_lasso, s = "lambda.min")[-1] != 0))

# 定义参数网格
param_grids <- list(
  Lasso = expand.grid(alpha = 1, lambda = 10^seq(-5, 0, length=50)),
  #Ridge = expand.grid(alpha = 0, lambda = 10^seq(-3, 2, length=50)),
  XGBoost = expand.grid(
    nrounds = 100,
    max_depth = 3:6,
    eta = c(0.01, 0.1),
    gamma = 0,
    colsample_bytree = 0.8,
    min_child_weight = 1,
    subsample = 0.8
  ),
  RF = data.frame(mtry = c(3,5,7)),
  Enet = expand.grid(alpha = seq(0.1, 0.9, 0.2), 
                     lambda = 10^seq(-3, 1, length=20))
  # 其他算法参数配置...
)

# 定义训练控制
ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 3,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  selectionFunction = "best",
  allowParallel = TRUE
)

# 多算法训练函数
train_models <- function(features) {
  model_list <- list()
  for(algo in names(param_grids)){
    cat("Training", algo, "...\n")
    model <- train(
      Group ~ .,
      data = data_balanced[, c(features, "Group")],
      method = switch(algo,
                      "Lasso" = "glmnet",
                      "XGBoost" = "xgbTree",
                      "RF" = "rf",
                      "Enet" = "glmnet"),
      trControl = ctrl,
      tuneGrid = param_grids[[algo]],
      metric = "ROC"
    )
    model_list[[algo]] <- model
  }
  return(model_list)
}

# 执行特征筛选后的训练
final_models <- train_models(selected_genes)
