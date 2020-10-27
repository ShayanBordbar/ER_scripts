# ER logistic model
# use the created negative and positive set of sequences
# create feature vectors for each sequence
# perform logistic regression with regularization in cross-validation setting

########################################################################################################################
# working directory and libraries
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")

library(aod)
library(ROCR)
library(caret)
library(DMwR)
library(ROSE)
library(VGAM)
library(LiblineaR)
library(MLmetrics)
library(PRROC)
library(purrr)
library(dplyr)
########################################################################################################################
# Functions
dataset_to_numericFeature_factor_label <- function(my_dataset){
  dataset <- my_dataset
  aaclass <- sapply(dataset, class)
  stopifnot(all(aaclass %in% "character"))
  aacolsnum <- colnames(dataset)[1:(ncol(dataset) - 1)]
  dataset[aacolsnum] <- sapply(dataset[aacolsnum], as.numeric)
  colnames(dataset)[ncol(dataset)] <- "label_set"
  dataset <- transform(dataset,  label_set = as.factor(label_set))
  return(dataset)
}
########################################################################################################################
calc_auprc <- function(model, data){
  
  index_Pos <- data$label_set == "Pos"
  index_Neg<- data$label_set == "Neg"
  
  predictions <- predict(model, data, type = "prob")
  
  pr.curve(predictions$Pos[index_Pos], predictions$Pos[index_Neg], curve = TRUE)
  
}
########################################################################################################################

########################################################################################################################
# deciding on features
# processing features
# dividing into test and training
# evaluating
########################################################################################################################
# this block is taken from examples on the internet

# check if any of the features are read as factors using:
# is.factor(data$somefeature)
# For a better understanding of how R is going to deal with the categorical variables, we can
#  use the contrasts() function: contrasts(data$somefeature)

model <- glm(Survived ~.,family=binomial(link='logit'),data=train)
summary(model)
anova(model, test="Chisq")
library(pscl)
pR2(model)
fitted.results <- predict(model,newdata=subset(test,select=c(2,3,4,5,6,7,8)),type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)

misClasificError <- mean(fitted.results != test$Survived)
print(paste('Accuracy',1-misClasificError))

library(ROCR)
p <- predict(model, newdata=subset(test,select=c(2,3,4,5,6,7,8)), type="response")
pr <- prediction(p, test$Survived)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc


xtabs(~admit + rank, data = mydata)
## CIs using profiled log-likelihood
confint(model)
## CIs using standard errors
confint.default(model)

library(aod)
# We can test for an overall effect of rank using the wald.test function of the aod library. 
# The order in which the coefficients are given in the table of coefficients is the same as the
# order of the terms in the model. This is important because the wald.test function refers to the
# coefficients by their order in the model. We use the wald.test function. b supplies the coefficients,
# while Sigma supplies the variance covariance matrix of the error terms, finally Terms tells R which 
# terms in the model are to be tested, in this case, terms 4, 5, and 6, are the three terms for the 
# levels of rank.
wald.test(b = coef(mylogit), Sigma = vcov(mylogit), Terms = 4:6)

## odds ratios and 95% CI
exp(cbind(OR = coef(mylogit), confint(mylogit)))
########################################################################################################################
########################################################################################################################
########################################################################################################################
# deciding on features: TRYING only affinities as sum of LRs
Negative_set_seq_list_chopped_scoreMat_list[[1]]
Positive_set_seq_list_chopped_scoreMat_list[[1]]

aa_alldata <- rbind(cbind(Negative_set_seq_list_chopped_scoreMat_list[[1]], 
                          rep(0, nrow(Negative_set_seq_list_chopped_scoreMat_list[[1]]))),
                    cbind(Positive_set_seq_list_chopped_scoreMat_list[[1]],
                          rep(1, nrow(Positive_set_seq_list_chopped_scoreMat_list[[1]]))))

my_logistic_dataset_list <- list()
my_logistic_dataset_list[[1]] <- list()
# divide data for each positive_negative comb into 10 train and test sets --> in cross-validation scheme
nu_partitionings <- 10
aa_neg_shuffle <- sample(x = c(1:nrow(Negative_set_seq_list_chopped_scoreMat_list[[1]])), 
                         size = nrow(Negative_set_seq_list_chopped_scoreMat_list[[1]]),
                         replace = F)
aa_pos_shuffle <- sample(x = c(1:nrow(Positive_set_seq_list_chopped_scoreMat_list[[1]])),
                         size = nrow(Positive_set_seq_list_chopped_scoreMat_list[[1]]),
                         replace = F)
aa_neg_borders <- ceiling(seq(1, length(aa_neg_shuffle), length.out = nu_partitionings + 1))
aa_pos_borders <- ceiling(seq(1, length(aa_pos_shuffle), length.out = nu_partitionings + 1))
aa_neg_norm <- Negative_set_seq_list_chopped_scoreMat_list[[1]][,1:(ncol(Negative_set_seq_list_chopped_scoreMat_list[[1]]) - 1)]/Negative_set_seq_list_chopped_scoreMat_list[[1]][,ncol(Negative_set_seq_list_chopped_scoreMat_list[[1]])] 
aa_pos_norm <- Positive_set_seq_list_chopped_scoreMat_list[[1]][,1:(ncol(Positive_set_seq_list_chopped_scoreMat_list[[1]]) - 1)]/Positive_set_seq_list_chopped_scoreMat_list[[1]][,ncol(Positive_set_seq_list_chopped_scoreMat_list[[1]])]
for(i in 1:nu_partitionings){
  aa_cur_neg_test_ind <- aa_neg_shuffle[aa_neg_borders[i]:(aa_neg_borders[i+1] - 1*as.integer(i != nu_partitionings))]
  aa_cur_neg_train_ind <- c(1:nrow(aa_neg_norm))[! c(1:nrow(aa_neg_norm)) %in% aa_cur_neg_test_ind]
  
  aa_cur_pos_test_ind <- aa_pos_shuffle[aa_pos_borders[i]:(aa_pos_borders[i+1] - 1*as.integer(i != nu_partitionings))]
  aa_cur_pos_train_ind <- c(1:nrow(aa_pos_norm))[! c(1:nrow(aa_pos_norm)) %in% aa_cur_pos_test_ind]
  
  aa_cur_neg_test <- aa_neg_norm[aa_cur_neg_test_ind, ]
  aa_cur_neg_train <- aa_neg_norm[aa_cur_neg_train_ind,]
  aa_cur_pos_test <- aa_pos_norm[aa_cur_pos_test_ind,]
  aa_cur_pos_train <- aa_pos_norm[aa_cur_pos_train_ind,]
  aa_cur_test <- rbind(cbind(aa_cur_neg_test,
                              rep(0,nrow(aa_cur_neg_test))),
                        cbind(aa_cur_pos_test,
                              rep(1,nrow(aa_cur_pos_test)) ))

  colnames(aa_cur_test)[ncol(aa_cur_test)] <- "label"
  aa_cur_test <- aa_cur_test[sample(x = c(1:nrow(aa_cur_test)),
                                    size = nrow(aa_cur_test), replace = F),]
  aa_cur_train <- rbind(cbind(aa_cur_neg_train,
                              rep(0,nrow(aa_cur_neg_train))),
                        cbind(aa_cur_pos_train,
                              rep(1,nrow(aa_cur_pos_train)) ))
  colnames(aa_cur_train)[ncol(aa_cur_train)] <- "label"
  aa_cur_train <- aa_cur_train[sample(x = c(1:nrow(aa_cur_train)), 
                                      size = nrow(aa_cur_train), replace = F),]
  print(paste0("set ", i ," contains ", nrow(aa_cur_test), 
               " test examples including ", nrow(aa_cur_neg_test), 
               " negative and ", nrow(aa_cur_pos_test), " positives. Set ",
               i, " also contains ",nrow(aa_cur_train) ," training point including ",
               nrow(aa_cur_neg_train), " negatives and ", nrow(aa_cur_pos_train), " positives."
  ))
  my_logistic_dataset_list[[1]][[i]] <- list(train = aa_cur_train,
                                             test = aa_cur_test)
  
}
names(my_logistic_dataset_list[[1]]) <- c(1:length(my_logistic_dataset_list[[1]]))

# running logisitic regression
my_logistic_model_list <- list()
my_logistic_model_list[[1]] <- list()
my_logistic_anova_list <- list()
my_logistic_anova_list[[1]] <- list()
my_logistic_prediction_list <- list()
my_logistic_prediction_list[[1]] <- list()
my_logistic_perf_list <- list()
my_logistic_perf_list[[1]] <- list()

par(mfrow = c(4, 3), mar = c(4,4,4,1))
aaauc <- numeric(length(my_logistic_dataset_list[[1]]))
for(i in 1:length(my_logistic_dataset_list[[1]])){
  aa_cur_train <- as.data.frame(my_logistic_dataset_list[[1]][[i]]$train)
  aa_cur_test <- as.data.frame(my_logistic_dataset_list[[1]][[i]]$test)
  my_logistic_model_list[[1]][[i]] <- glm(label ~.,family=binomial(link='logit'),data=aa_cur_train)
  my_logistic_anova_list[[1]][[i]] <- anova(my_logistic_model_list[[1]][[i]], test="Chisq")
  my_logistic_prediction_list[[1]][[i]] <- predict(my_logistic_model_list[[1]][[i]], newdata=aa_cur_test[, c(1:(ncol(aa_cur_test) - 1))], type="response")
  aapr <- prediction(my_logistic_prediction_list[[1]][[i]], aa_cur_test$label)
  my_logistic_perf_list[[1]][[i]] <- performance(aapr, measure = "tpr", x.measure = "fpr")
  aauc <- performance(aapr, measure = "auc")
  aaauc[i] <- aauc@y.values[[1]]
  plot(my_logistic_perf_list[[1]][[i]], main = paste0("auc: ", (format(aaauc[i], digits = 2))), col = 2)
  abline(a = 0, b = 1, col = 1, lwd = 0.5)
}
print(paste0("average 10 fold CV auc: ",
             mean(aaauc)))
########################################################################################################################
# Using caret package to address imbalanced datasets
aa_neg_norm <- Negative_set_seq_list_chopped_scoreMat_list[[1]][,1:(ncol(Negative_set_seq_list_chopped_scoreMat_list[[1]]) - 1)]/Negative_set_seq_list_chopped_scoreMat_list[[1]][,ncol(Negative_set_seq_list_chopped_scoreMat_list[[1]])] 
aa_pos_norm <- Positive_set_seq_list_chopped_scoreMat_list[[1]][,1:(ncol(Positive_set_seq_list_chopped_scoreMat_list[[1]]) - 1)]/Positive_set_seq_list_chopped_scoreMat_list[[1]][,ncol(Positive_set_seq_list_chopped_scoreMat_list[[1]])]

aa_alldata <- rbind(cbind(aa_neg_norm, 
                          rep(0, nrow(aa_neg_norm))),
                    cbind(aa_pos_norm,
                          rep(1, nrow(aa_pos_norm))))


library(caret)
# detect near zero variance variables
colnames(aa_alldata)[ncol(aa_alldata)] <- "label"
aanzv <- nearZeroVar(aa_alldata, saveMetrics= T)
#  Identifying Correlated Predictors
aadescrCor <- cor(aa_alldata)
summary(aadescrCor[upper.tri(aadescrCor)])
aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
colnames(aa_alldata)[aahighlyCorDescr]
aa_alldata_filtered <- aa_alldata[,-aahighlyCorDescr]
aadescrCor <- cor(aa_alldata_filtered)
summary(aadescrCor[upper.tri(aadescrCor)])
# find Linear Dependencies
aacomboInfo <- findLinearCombos(aa_alldata_filtered)
aacomboInfo

# partitioning to cross validation folds
aatrainIndex <- createFolds(y = aa_alldata_filtered[, ncol(aa_alldata_filtered)], k = 10,list = T)
# simple partition
aatrainpartition <- createDataPartition(y = aa_alldata_filtered[, ncol(aa_alldata_filtered)], times = 1, p = 0.75,list = F)
aa_train_data <- aa_alldata_filtered[aatrainpartition, ]
aa_test_data <-aa_alldata_filtered[-aatrainpartition, ]

aa_fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)
aa_train_data <- as.data.frame(aa_train_data)
aa_train_data$label <- as.factor(aa_train_data$label)

aa_test_data <- as.data.frame(aa_test_data)
aa_test_data$label <- as.factor(aa_test_data$label)


set.seed(825)
aagbmFit1 <- train(label ~ ., data = aa_train_data, 
                 method = "gbm", 
                 trControl = aa_fitControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = FALSE)
aagbmFit1
trellis.par.set(caretTheme())
plot(aagbmFit1)
####
AAfitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 3,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

set.seed(825)
aa_train_data$label <- as.character(levels(aa_train_data$label)[aa_train_data$label])
aa_train_data$label[aa_train_data$label %in% "0"] <- "neg"
aa_train_data$label[aa_train_data$label %in% "1"] <- "pos"

aa_test_data$label <- as.character(levels(aa_test_data$label)[aa_test_data$label])
aa_test_data$label[aa_test_data$label %in% "0"] <- "neg"
aa_test_data$label[aa_test_data$label %in% "1"] <- "pos"
aa_test_data$label <- as.factor(aa_test_data$label)
aa_train_data$label <- as.factor(aa_train_data$label)
# trying a tree based method
AAgbmFit3 <- train(label ~ ., data = aa_train_data, 
                 method = "gbm", 
                 trControl = AAfitControl, 
                 verbose = FALSE, 
                 #tuneGrid = gbmGrid,
                 ## Specify which metric to optimize
                 metric = "ROC")
AAgbmFit3

aapred <- predict(AAgbmFit3, newdata = aa_test_data[, -ncol(aa_test_data)])
head(aapred)
sum(aapred[aa_test_data$label == "pos"] == aa_test_data$label[aa_test_data$label == "pos"])
# trying svm based method
AAsvmFit <- train(label ~ ., data = aa_train_data, 
                method = "svmRadial", 
                trControl = AAfitControl, 
                preProc = c("center", "scale"),
                #tuneLength = 8,
                metric = "ROC")
AAsvmFit


aaresamps <- resamples(list(GBM = AAgbmFit3,
                          SVM = AAsvmFit))
aaresamps
summary(aaresamps)
theme1 <- trellis.par.get()
theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
theme1$plot.symbol$pch = 16
theme1$plot.line$col = rgb(1, 0, 0, .7)
theme1$plot.line$lwd <- 2
trellis.par.set(theme1)
bwplot(aaresamps, layout = c(3, 1))
trellis.par.set(caretTheme())
dotplot(aaresamps, metric = "ROC")
splom(aaresamps)

difValues <- diff(aaresamps)
difValues
summary(difValues)
trellis.par.set(theme1)
bwplot(difValues, layout = c(3, 1))


test_roc <- function(model, data) {
  library(pROC)
  roc_obj <- roc(data$label, 
                 predict(model, data, type = "prob")[, "pos"],
                 levels = c("neg", "pos"))
  ci(roc_obj)
}


aactrl <- trainControl(method = "repeatedcv", repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary)
set.seed(5627)

aaorig_fit <- train(label ~ ., data = aa_train_data, 
                  method = "regLogistic",
                  #nbagg = 50,
                  metric = "ROC",
                  trControl = aactrl)


set.seed(5627)
aactrl <- trainControl(method = "repeatedcv", repeats = 5,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary,
                     ## new option here:
                     sampling = "down")

set.seed(5627)
aadown_inside <- train(label ~ ., data = aa_train_data,
                     method = "regLogistic",
                     #nbagg = 50,
                     metric = "ROC",
                     trControl = aactrl)
summary(aadown_inside)

set.seed(5627)
aadown_inside_svm <- train(label ~ ., data = aa_train_data,
                           method = "svmRadial", 
                           trControl = aactrl, 
                           preProc = c("center", "scale"),
                           metric = "ROC")
summary(aadown_inside)
## now just change that option
aactrl$sampling <- "up"

set.seed(5627)
aaup_inside <- train(label ~ ., data = aa_train_data,
                   method = "regLogistic",
                   #nbagg = 50,
                   metric = "ROC",
                   trControl = aactrl)

aactrl$sampling <- "rose"

set.seed(5627)
aarose_inside <- train(label ~ ., data = aa_train_data,
                     method = "regLogistic",
                     #nbagg = 50,
                     metric = "ROC",
                     trControl = aactrl)

aactrl$sampling <- "smote"

set.seed(5627)
aasmote_inside <- train(label ~ ., data = aa_train_data,
                      method = "regLogistic",
                      #nbagg = 50,
                      metric = "ROC",
                      trControl = aactrl)

aainside_models <- list(original = aaorig_fit,
                      down = aadown_inside,
                      down_svm = aadown_inside_svm,
                      up = aaup_inside,
                      SMOTE = aasmote_inside
                      #,ROSE = aarose_inside
                      )


aainside_resampling <- resamples(aainside_models)

aainside_test <- lapply(aainside_models, test_roc, data = aa_test_data)
aainside_test <- lapply(aainside_test, as.vector)
aainside_test <- do.call("rbind", aainside_test)
colnames(aainside_test) <- c("lower", "ROC", "upper")
aainside_test <- as.data.frame(aainside_test)

summary(aainside_resampling, metric = "ROC")
summary(aainside_resampling)
theme1 <- trellis.par.get()
theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
theme1$plot.symbol$pch = 16
theme1$plot.line$col = rgb(1, 0, 0, .7)
theme1$plot.line$lwd <- 2
trellis.par.set(theme1)
bwplot(aainside_resampling, layout = c(3, 1))
trellis.par.set(caretTheme())
dotplot(aainside_resampling, metric = "ROC")
splom(aainside_resampling)

aainside_test



aapr <- predict(aadown_inside_svm, aa_test_data[,-ncol(aa_test_data)])
sum(aapr[aa_test_data$label == "pos"] == aa_test_data$label[aa_test_data$label == "pos"])

sum(aa_test_data$label == "pos")

######################################################################################################
# Compare logistic reg, svm, RF with and without imbalance tricks

aa_neg_norm <- Negative_set_seq_list_chopped_scoreMat_list[[1]][,1:(ncol(Negative_set_seq_list_chopped_scoreMat_list[[1]]) - 1)]/Negative_set_seq_list_chopped_scoreMat_list[[1]][,ncol(Negative_set_seq_list_chopped_scoreMat_list[[1]])] 
aa_pos_norm <- Positive_set_seq_list_chopped_scoreMat_list[[1]][,1:(ncol(Positive_set_seq_list_chopped_scoreMat_list[[1]]) - 1)]/Positive_set_seq_list_chopped_scoreMat_list[[1]][,ncol(Positive_set_seq_list_chopped_scoreMat_list[[1]])]

aa_alldata <- rbind(cbind(aa_neg_norm, 
                          rep(0, nrow(aa_neg_norm))),
                    cbind(aa_pos_norm,
                          rep(1, nrow(aa_pos_norm))))
colnames(aa_alldata)[ncol(aa_alldata)] <- "label"
aa_alldata <- as.data.frame(aa_alldata)
aa_alldata$label[aa_alldata$label %in% 0] <- "neg"
aa_alldata$label[aa_alldata$label %in% 1] <- "pos"
aa_alldata$label <- as.factor(aa_alldata$label)

#  Identifying Correlated Predictors
aadescrCor <- cor(aa_alldata[,-ncol(aa_alldata)])
summary(aadescrCor[upper.tri(aadescrCor)])
aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
colnames(aa_alldata)[aahighlyCorDescr]
aa_alldata_filtered <- aa_alldata[,-aahighlyCorDescr]

aatrainpartition <- createDataPartition(y = aa_alldata_filtered[, ncol(aa_alldata_filtered)], times = 1, p = 0.75,list = F)
aa_train_data <- aa_alldata_filtered[aatrainpartition, ]
aa_test_data <-aa_alldata_filtered[-aatrainpartition, ]

# logistic regression wo imbalance correction
aactrl <- trainControl(method = "repeatedcv", repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary)

aalogis_orig <- train(label ~ ., data = aa_train_data,
                        method = "regLogistic",
                        #nbagg = 50,
                        metric = "ROC",
                        trControl = aactrl)
# confusion matrix
colnames(aa_test_data)[14] <- "NKX3-1"
aaPred <- predict(aalogis_orig, 
                  aa_test_data[,-ncol(aa_test_data)])
aaTruth <- aa_test_data$label
aaxtab <- table(aaPred, aaTruth)

aa_conf_logreg <- confusionMatrix(aaxtab, positive = "pos")
aa_conf_logreg

# fourfoldplot(aa_conf_logreg$table, color = c("#CC6666", "#99CC99"),
#              conf.level = 0, margin = 1, main = "Confusion Matrix")

# logistic regression with imbalance correction: downsampling
aactrl <- trainControl(method = "repeatedcv", repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       ## new option here:
                       sampling = "down")
colnames(aa_train_data)[14] <- "NKX3_1"
colnames(aa_test_data)[14] <- "NKX3_1"
aalogis_down <- train(label ~ ., data = aa_train_data,
                      method = "regLogistic",
                      #nbagg = 50,
                      metric = "ROC",
                      trControl = aactrl)
colnames(aa_test_data)[14] <- "NKX3_1"
aaPred <- predict(aalogis_down, 
                  aa_test_data[,-ncol(aa_test_data)])
aaTruth <- aa_test_data$label
aaxtab <- table(aaPred, aaTruth)
# confusion matrix
aa_conf_logreg <- confusionMatrix(aaxtab, positive = "pos")
aa_conf_logreg

# logistic regression with imbalance correction: upsampling
aactrl <- trainControl(method = "repeatedcv", repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       ## new option here:
                       sampling = "up")


aalogis_up <- train(label ~ ., data = aa_train_data,
                      method = "regLogistic",
                      #nbagg = 50,
                      metric = "ROC",
                      trControl = aactrl)
aaPred <- predict(aalogis_up, 
                  aa_test_data[,-ncol(aa_test_data)])
aaTruth <- aa_test_data$label
aaxtab <- table(aaPred, aaTruth)
# confusion matrix
aa_conf_logreg <- confusionMatrix(aaxtab, positive = "pos")
aa_conf_logreg

# logistic regression with imbalance correction: rose
aactrl <- trainControl(method = "repeatedcv", repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       ## new option here:
                       sampling = "rose")


colnames(aa_train_data)[14] <- "NKX3_1"
colnames(aa_test_data)[14] <- "NKX3_1"
aalogis_rose <- train(label ~ ., data = aa_train_data,
                    method = "regLogistic",
                    #nbagg = 50,
                    metric = "ROC",
                    trControl = aactrl)

aaPred <- predict(aalogis_rose, 
                  aa_test_data[,-ncol(aa_test_data)])
aaTruth <- aa_test_data$label
aaxtab <- table(aaPred, aaTruth)
# confusion matrix
aa_conf_logreg <- confusionMatrix(aaxtab, positive = "pos")
aa_conf_logreg
# logistic regression with imbalance correction: smote

aactrl <- trainControl(method = "repeatedcv", repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       ## new option here:
                       sampling = "smote")


aalogis_smote <- train(label ~ ., data = aa_train_data,
                      method = "regLogistic",
                      #nbagg = 50,
                      metric = "ROC",
                      trControl = aactrl)
aaPred <- predict(aalogis_smote, 
                  aa_test_data[,-ncol(aa_test_data)])
aaTruth <- aa_test_data$label
aaxtab <- table(aaPred, aaTruth)
# confusion matrix
aa_conf_logreg <- confusionMatrix(aaxtab, positive = "pos")
aa_conf_logreg
######################################################################################################
# using SVM
AAfitControl <- trainControl(method = "repeatedcv",
                             number = 10,
                             repeats = 3,
                             ## Estimate class probabilities
                             classProbs = TRUE,
                             ## Evaluate performance using 
                             ## the following function
                             summaryFunction = twoClassSummary)
AAsvm_orig <- train(label ~ ., data = aa_train_data, 
                  method = "svmRadial", 
                  trControl = AAfitControl, 
                  preProc = c("center", "scale"),
                  #tuneLength = 8,
                  metric = "ROC")

aaPred <- predict(AAsvm_orig, 
                  aa_test_data[,-ncol(aa_test_data)])
aaTruth <- aa_test_data$label
aaxtab <- table(aaPred, aaTruth)
# confusion matrix
aa_conf_logreg <- confusionMatrix(aaxtab, positive = "pos")
aa_conf_logreg

# SVM with downsampling
AAfitControl <- trainControl(method = "repeatedcv",
                             number = 10,
                             repeats = 3,
                             ## Estimate class probabilities
                             classProbs = TRUE,
                             ## Evaluate performance using 
                             ## the following function
                             summaryFunction = twoClassSummary,  
                             sampling = "down")
AAsvm_down <- train(label ~ ., data = aa_train_data, 
                    method = "svmRadial", 
                    trControl = AAfitControl, 
                    preProc = c("center", "scale"),
                    #tuneLength = 8,
                    metric = "ROC")

aaPred <- predict(AAsvm_down, 
                  aa_test_data[,-ncol(aa_test_data)])
aaTruth <- aa_test_data$label
aaxtab <- table(aaPred, aaTruth)
# confusion matrix
aa_conf_logreg <- confusionMatrix(aaxtab, positive = "pos")
aa_conf_logreg

# SVM with downsampling
AAfitControl <- trainControl(method = "repeatedcv",
                             number = 10,
                             repeats = 3,
                             ## Estimate class probabilities
                             classProbs = TRUE,
                             ## Evaluate performance using 
                             ## the following function
                             summaryFunction = twoClassSummary,  
                             sampling = "up")
AAsvm_up <- train(label ~ ., data = aa_train_data, 
                    method = "svmRadial", 
                    trControl = AAfitControl, 
                    preProc = c("center", "scale"),
                    #tuneLength = 8,
                    metric = "ROC")

aaPred <- predict(AAsvm_up, 
                  aa_test_data[,-ncol(aa_test_data)])
aaTruth <- aa_test_data$label
aaxtab <- table(aaPred, aaTruth)
# confusion matrix
aa_conf_logreg <- confusionMatrix(aaxtab, positive = "pos")
aa_conf_logreg
#########################################################################################################
# repeating the logistic regression with length as one factor
# logistic regression with imbalance correction: downsampling
Negative_set_seq_list_chopped_scoreMat_list[[1]]
Positive_set_seq_list_chopped_scoreMat_list[[1]]

aa_alldata <- rbind(cbind(Negative_set_seq_list_chopped_scoreMat_list[[1]][,-ncol(Negative_set_seq_list_chopped_scoreMat_list[[1]])]/100, 
                          rep("neg", nrow(Negative_set_seq_list_chopped_scoreMat_list[[1]]))),
                    cbind(Positive_set_seq_list_chopped_scoreMat_list[[1]][,-ncol(Positive_set_seq_list_chopped_scoreMat_list[[1]])]/100,
                          rep("pos", nrow(Positive_set_seq_list_chopped_scoreMat_list[[1]]))))
colnames(aa_alldata)[ncol(aa_alldata)] <- "label"
aa_alldata <- as.data.frame(aa_alldata)
aa_alldata$label <- as.factor(aa_alldata$label)

boxplot.matrix(rbind(Negative_set_seq_list_chopped_scoreMat_list[[1]], 
      Positive_set_seq_list_chopped_scoreMat_list[[1]])[,-34])
#  Identifying Correlated Predictors
aadescrCor <- cor(as.matrix(aa_alldata[,-ncol(aa_alldata)]))
summary(aadescrCor[upper.tri(aadescrCor)])
aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
colnames(aa_alldata)[aahighlyCorDescr]
aa_alldata_filtered <- aa_alldata[,-aahighlyCorDescr]

aatrainpartition <- createDataPartition(y = aa_alldata_filtered[, ncol(aa_alldata_filtered)], times = 1, p = 0.75,list = F)
aa_train_data2 <- aa_alldata_filtered[aatrainpartition, ]
aa_test_data2 <-aa_alldata_filtered[-aatrainpartition, ]

aa_train_data

aanzv <- nearZeroVar(aa_train_data, saveMetrics= T)

aactrl <- trainControl(method = "repeatedcv", repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       #savePredictions = "final",
                       ## new option here:
                       sampling = "down")


aalogis_down2 <- train(label ~ ., data = aa_train_data2,
                      method = "regLogistic",
                      #preProc = c("center", "scale"),
                      #nbagg = 50,
                      metric = "ROC",
                      trControl = aactrl)

s#colnames(aa_test_data)[14] <- "NKX3-1"
aaPred <- predict(aalogis_down2, 
                  aa_test_data[,-ncol(aa_test_data)])
aaTruth <- aa_test_data$label
aaxtab <- table(aaPred, aaTruth)
# confusion matrix
aa_conf_logreg <- confusionMatrix(aaxtab, positive = "pos")
aa_conf_logreg

# logistic regression with imbalance correction: upsampling
aactrl <- trainControl(method = "repeatedcv",
                       repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       ## new option here:
                       sampling = "up")


aalogis_up <- train(label ~ ., 
                    data = aa_train_data,
                    method = "regLogistic",
                    #nbagg = 50,
                    metric = "ROC",
                    trControl = aactrl)
aaPred <- predict(aalogis_up, 
                  aa_test_data[,-ncol(aa_test_data)])
aaTruth <- aa_test_data$label
aaxtab <- table(aaPred, aaTruth)
# confusion matrix
aa_conf_logreg <- confusionMatrix(aaxtab, positive = "pos")
aa_conf_logreg




aaresamps <- resamples(list(GBM = AAgbmFit3,
                            SVM = AAsvmFit))
aaresamps
summary(aaresamps)
theme1 <- trellis.par.get()
theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
theme1$plot.symbol$pch = 16
theme1$plot.line$col = rgb(1, 0, 0, .7)
theme1$plot.line$lwd <- 2
trellis.par.set(theme1)
bwplot(aaresamps, layout = c(3, 1))
trellis.par.set(caretTheme())
dotplot(aaresamps, metric = "ROC")
splom(aaresamps)

difValues <- diff(aaresamps)
difValues
summary(difValues)
trellis.par.set(theme1)
bwplot(difValues, layout = c(3, 1))

aapr <- prediction(my_logistic_prediction_list[[1]][[i]], aa_cur_test$label)
my_logistic_perf_list[[1]][[i]] <- performance(aapr, measure = "tpr", x.measure = "fpr")
aauc <- performance(aapr, measure = "auc")
aaauc[i] <- aauc@y.values[[1]]
plot(my_logistic_perf_list[[1]][[i]], main = paste0("auc: ", (format(aaauc[i], digits = 2))), col = 2)
abline(a = 0, b = 1, col = 1, lwd = 0.5)
#########################################################################################################
# repeating the logistic regression with no lenngth normalization



Negative_set_seq_list_chopped_notlnnorm_scoreMat_list[[1]]
Positive_set_seq_list_chopped_notlnnorm_scoreMat_list[[1]]
aancol <- ncol(Negative_set_seq_list_chopped_notlnnorm_scoreMat_list[[1]])
aa_alldata <- rbind(cbind(Negative_set_seq_list_chopped_notlnnorm_scoreMat_list[[1]][,-aancol], 
                          rep("neg", nrow(Negative_set_seq_list_chopped_notlnnorm_scoreMat_list[[1]]))),
                    cbind(Positive_set_seq_list_chopped_notlnnorm_scoreMat_list[[1]][,-aancol],
                          rep("pos", nrow(Positive_set_seq_list_chopped_notlnnorm_scoreMat_list[[1]]))))
colnames(aa_alldata)[ncol(aa_alldata)] <- "label"

aa_alldata <- as.data.frame(aa_alldata, stringsAsFactors = F)
for(i in 1:(ncol(aa_alldata) - 1)){
  aa_alldata[, i] <- as.numeric(aa_alldata[,i])
}

aa_alldata$label <- as.factor(aa_alldata$label)
par(mfrow= c(1,1))
boxplot.matrix(rbind(Negative_set_seq_list_chopped_notlnnorm_scoreMat_list[[1]], 
                     Positive_set_seq_list_chopped_notlnnorm_scoreMat_list[[1]])[,-34])
#  Identifying Correlated Predictors
# aadescrCor <- cor(as.matrix(aa_alldata[,-ncol(aa_alldata)]))
# summary(aadescrCor[upper.tri(aadescrCor)])
# aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
# colnames(aa_alldata)[aahighlyCorDescr]
# aa_alldata_filtered <- aa_alldata[,-aahighlyCorDescr]
aanzv <- nearZeroVar(aa_alldata)
aa_alldata <- aa_alldata[, -aanzv]

aa_newpart <- createFolds( y = aa_alldata[, ncol(aa_alldata)],k = 6)

my_logistic_model_list <- list()
my_logistic_model_list[[1]] <- list()
my_logistic_anova_list <- list()
my_logistic_anova_list[[1]] <- list()
my_logistic_prediction_list <- list()
my_logistic_prediction_list[[1]] <- list()
my_logistic_perf_list <- list()
my_logistic_perf_list[[1]] <- list()
par(mfrow = c(4, 3), mar = c(4,4,4,1))
aaauc <- numeric(length(aa_newpart))

for(i in 1:length(aa_newpart)){
  print(i)
  aa_c_train <- aa_alldata[-aa_newpart[[i]], ]
  aa_c_testt <- aa_alldata[aa_newpart[[i]], ]
  aa_weights <- numeric(nrow(aa_c_train))
  for(j in 1:nrow(aa_c_train)){
    aa_weights[j] <- 1/table(aa_c_train$label)[aa_c_train$label[j]]
  }
  my_logistic_model_list[[1]][[i]] <- glm(label ~.,family=binomial(link='logit'),data=aa_c_train, weights = aa_weights)
  my_logistic_anova_list[[1]][[i]] <- anova(my_logistic_model_list[[1]][[i]], test="Chisq")
  my_logistic_prediction_list[[1]][[i]] <- predict(my_logistic_model_list[[1]][[i]], newdata=aa_c_testt[, -ncol(aa_c_testt)], type="response")
  aapr <- prediction(my_logistic_prediction_list[[1]][[i]], aa_c_testt$label)
  my_logistic_perf_list[[1]][[i]] <- performance(aapr, measure = "tpr", x.measure = "fpr")
  aauc <- performance(aapr, measure = "auc")
  aaauc[i] <- aauc@y.values[[1]]
  plot(my_logistic_perf_list[[1]][[i]], main = paste0("auc: ", (format(aaauc[i], digits = 2))), col = 2)
  abline(a = 0, b = 1, col = 1, lwd = 0.5)
  # aaxtab <- table(aaPred, aaTruth)
  # # confusion matrix
  # aa_conf_logreg <- confusionMatrix(aaxtab, positive = "pos")
  
}
print(paste0("average 10 fold CV auc: ",
             mean(aaauc)))

aatrainpartition <- createDataPartition(y = aa_alldata[, ncol(aa_alldata)], times = 1, p = 0.75,list = F)
aa_train_data2 <- aa_alldata[aatrainpartition, ]
aa_test_data2 <-aa_alldata[-aatrainpartition, ]

aactrl <- trainControl(method = "repeatedcv", repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       #savePredictions = "final",
                       ## new option here:
                       sampling = "down")


aalogis_down2 <- train(label ~ ., data = aa_train_data2,
                       method = "regLogistic",
                       #preProc = c("center", "scale"),
                       #nbagg = 50,
                       metric = "ROC",
                       trControl = aactrl)

#s#colnames(aa_test_data2)[14] <- "NKX3-1"
aaPred <- predict(aalogis_down2, 
                  aa_test_data2[,-ncol(aa_test_data2)])
aaTruth <- aa_test_data2$label
aaxtab <- table(aaPred, aaTruth)
# confusion matrix
aa_conf_logreg <- confusionMatrix(aaxtab, positive = "pos")
aa_conf_logreg

# logistic regression with imbalance correction: upsampling
aactrl <- trainControl(method = "repeatedcv",
                       repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       ## new option here:
                       sampling = "up")


aalogis_up2 <- train(label ~ ., 
                    data = aa_train_data2,
                    method = "regLogistic",
                    #nbagg = 50,
                    metric = "ROC",
                    trControl = aactrl)
aaPred <- predict(aalogis_up2, 
                  aa_test_data2[,-ncol(aa_test_data2)])
aaTruth <- aa_test_data2$label
aaxtab <- table(aaPred, aaTruth)
# confusion matrix
aa_conf_logreg <- confusionMatrix(aaxtab, positive = "pos")
aa_conf_logreg

#########################################################################################################
# feature sets:
#sum LR
feature_set_sum_LLR_lengthNormalized <- rbind(Negative_set_seq_list_chopped_scoreMat_list[[2]], 
                                              Positive_set_seq_list_chopped_scoreMat_list[[2]])
colnames(feature_set_sum_LLR_lengthNormalized)[14] <- "NKX3_1"

feature_set_sum_LLR_NotlengthNormalized <- rbind(Negative_set_seq_list_chopped_scoreMat_list[[1]], 
                                                 Positive_set_seq_list_chopped_scoreMat_list[[1]])
colnames(feature_set_sum_LLR_NotlengthNormalized)[14] <- "NKX3_1"
# adjacency
feature_set_Adjmat <- rbind(Negative_set_seq_list_chopped_adjanScore_list[[1]], 
                            Positive_set_seq_list_chopped_adjanScore_list[[1]])
label_set <- c(rep("Neg", nrow(Negative_set_seq_list_chopped_scoreMat_list[[1]])), 
               rep("Pos", nrow(Positive_set_seq_list_chopped_scoreMat_list[[1]])))

# changing the name of NKX3-1 to NKX3_1 so it works in the classification models
aa <- strsplit(colnames(feature_set_Adjmat), split = "_vs_")
for(i in 1:length(aa)){
  aa[[i]][aa[[i]] == "NKX3-1"] <- "NKX3_1"
}
for(i in 1:length(aa)){
  aa[[i]] <- paste(aa[[i]][1], aa[[i]][2], sep = "_vs_")
}
colnames(feature_set_Adjmat) <- unlist(aa)
# overlap
feature_set_Overlapmat <- rbind(Negative_set_seq_list_chopped_OverlapScore_list[[1]], 
                                Positive_set_seq_list_chopped_OverlapScore_list[[1]])
colnames(feature_set_Overlapmat) <- unlist(aa)

my_learning_datasets <- list()
# length normalized sum LR only
my_learning_datasets[[1]] <- as.data.frame(cbind(feature_set_sum_LLR_lengthNormalized,
                                                 label_set), 
                                           stringsAsFactors =F)
my_learning_datasets[[1]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets[[1]])
# Notlength normalized sum LR only
my_learning_datasets[[2]] <- as.data.frame(cbind(feature_set_sum_LLR_NotlengthNormalized,
                                                 label_set), 
                                           stringsAsFactors =F)
my_learning_datasets[[2]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets[[2]])
# adjacancy only
my_learning_datasets[[3]] <- as.data.frame(cbind(feature_set_Adjmat,
                                                 label_set), 
                                           stringsAsFactors =F)
my_learning_datasets[[3]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets[[3]])

# Overlap only
my_learning_datasets[[4]] <- as.data.frame(cbind(feature_set_Overlapmat,
                                                 label_set), 
                                           stringsAsFactors =F)
my_learning_datasets[[4]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets[[4]])
# not length normalized and adj
my_learning_datasets[[5]] <- as.data.frame(cbind(cbind(feature_set_sum_LLR_NotlengthNormalized,
                                                 feature_set_Adjmat),
                                                 label_set), 
                                           stringsAsFactors =F)
my_learning_datasets[[5]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets[[5]])
# not length normalized and overlap
my_learning_datasets[[6]] <- as.data.frame(cbind(cbind(feature_set_sum_LLR_NotlengthNormalized,
                                                       feature_set_Overlapmat),
                                                 label_set), 
                                           stringsAsFactors =F)
my_learning_datasets[[6]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets[[6]])
# not length normalized and adj and overlap
my_learning_datasets[[7]] <- as.data.frame(cbind(cbind(feature_set_sum_LLR_NotlengthNormalized,
                                                       cbind(feature_set_Adjmat, 
                                                             feature_set_Overlapmat)),
                                                 label_set), 
                                           stringsAsFactors =F)
colnames(my_learning_datasets[[7]])[35:595] <- paste(colnames(my_learning_datasets[[7]])[35:595], "adj", sep = "_at_")
colnames(my_learning_datasets[[7]])[596:1156] <- paste(colnames(my_learning_datasets[[7]])[596:1156], "ovl", sep = "_at_")
my_learning_datasets[[7]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets[[7]])

#  length only
my_learning_datasets[[8]] <- as.data.frame(cbind(feature_set_sum_LLR_NotlengthNormalized[, 34], 
                                                 label_set), 
                                           stringsAsFactors =F)
my_learning_datasets[[8]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets[[8]])

names(my_learning_datasets) <- c("Sum_LR_LengthNorm", "Sum_LR", "Adjacency", "Overlap", 
                                 "Sum_LR_plus_Adjacency", "Sum_LR_plus_Overlap", 
                                 "Sum_LR_plus_Adjacency_plus_Overlap", "length_only")
#########################################################################################################
my_data_partitions <- list()
my_learning_datasets_results_logistic <- list()

my_data_partitions[[1]] <- createDataPartition(y = my_learning_datasets[[1]][, ncol(my_learning_datasets[[1]])], times = 1, p = 0.75,list = F)

# running logistic and SVM for all seven feature sets
# examining training and test performance

aanzv <- nearZeroVar(my_learning_datasets[[1]], saveMetrics= F)
aa_all_data <- my_learning_datasets[[1]][,-aanzv]

#  Identifying Correlated Predictors
aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
summary(aadescrCor[upper.tri(aadescrCor)])
aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
colnames(aa_all_data)[aahighlyCorDescr]
aa_all_data <- aa_all_data[,-aahighlyCorDescr]
ncol(aa_all_data)

aa_train_data <- aa_all_data[my_data_partitions[[1]], ]
aa_test_data <-  aa_all_data[-my_data_partitions[[1]], ]

aactrl <- trainControl(method = "repeatedcv", repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       #savePredictions = "final",
                       ## new option here:
                       sampling = "down")


aalogis_down <- train(label_set ~ ., data = aa_train_data,
                       method = "regLogistic",
                       #preProc = c("center", "scale"),
                       #nbagg = 50,
                       metric = "ROC",
                       trControl = aactrl)
my_learning_datasets_results_logistic[[1]] <- aalogis_down
###### test performance
aaPred <- predict(aalogis_down, 
                  aa_test_data[,-ncol(aa_test_data)])
aaTruth <- aa_test_data$label_set
aaxtab <- table(aaPred, aaTruth)
# confusion matrix
aa_conf_logreg <- confusionMatrix(aaxtab, positive = "Pos")
aa_conf_logreg
###### training performance
aaPred2 <- predict(aalogis_down, 
                  aa_train_data[,-ncol(aa_train_data)])
aaTruth2 <- aa_train_data$label_set
aaxtab2 <- table(aaPred2, aaTruth2)
# confusion matrix
aa_conf_logreg2 <- confusionMatrix(aaxtab2, positive = "Pos")
aa_conf_logreg2
#################################################
aanzv <- nearZeroVar(my_learning_datasets[[2]], saveMetrics= F)
aa_all_data <- my_learning_datasets[[2]][,-aanzv]

#  Identifying Correlated Predictors
aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
summary(aadescrCor[upper.tri(aadescrCor)])
aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
colnames(aa_all_data)[aahighlyCorDescr]
aa_all_data <- aa_all_data[,-aahighlyCorDescr]
ncol(aa_all_data)

aa_train_data <- aa_all_data[my_data_partitions[[1]], ]
aa_test_data <-  aa_all_data[-my_data_partitions[[1]], ]

aactrl <- trainControl(method = "repeatedcv", repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       #savePredictions = "final",
                       ## new option here:
                       sampling = "down")


aalogis_down <- train(label_set ~ ., data = aa_train_data,
                      method = "regLogistic",
                      #preProc = c("center", "scale"),
                      #nbagg = 50,
                      metric = "ROC",
                      trControl = aactrl)
my_learning_datasets_results_logistic[[2]] <- aalogis_down
###### test performance
aaPred <- predict(aalogis_down, 
                  aa_test_data[,-ncol(aa_test_data)])
aaTruth <- aa_test_data$label_set
aaxtab <- table(aaPred, aaTruth)
# confusion matrix
aa_conf_logreg <- confusionMatrix(aaxtab, positive = "Pos")
aa_conf_logreg
###### training performance
aaPred2 <- predict(aalogis_down, 
                   aa_train_data[,-ncol(aa_train_data)])
aaTruth2 <- aa_train_data$label_set
aaxtab2 <- table(aaPred2, aaTruth2)
# confusion matrix
aa_conf_logreg2 <- confusionMatrix(aaxtab2, positive = "Pos")
aa_conf_logreg2

#################################################
aanzv <- nearZeroVar(my_learning_datasets[[3]], saveMetrics= F)
aa_all_data <- my_learning_datasets[[3]][,-aanzv]

#  Identifying Correlated Predictors
aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
summary(aadescrCor[upper.tri(aadescrCor)])
aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
colnames(aa_all_data)[aahighlyCorDescr]
aa_all_data <- aa_all_data[,-aahighlyCorDescr]
ncol(aa_all_data)



aa_train_data <- aa_all_data[my_data_partitions[[1]], ]
aa_test_data <-  aa_all_data[-my_data_partitions[[1]], ]

aactrl <- trainControl(method = "repeatedcv", repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       #savePredictions = "final",
                       ## new option here:
                       sampling = "down")


aalogis_down <- train(label_set ~ ., data = aa_train_data,
                      method = "regLogistic",
                      #preProc = c("center", "scale"),
                      #nbagg = 50,
                      metric = "ROC",
                      trControl = aactrl)
my_learning_datasets_results_logistic[[3]] <- aalogis_down
###### test performance
aaPred <- predict(my_learning_datasets_results_logistic[[3]], 
                  aa_test_data[,-ncol(aa_test_data)])
aaTruth <- aa_test_data$label_set
aaxtab <- table(aaPred, aaTruth)
# confusion matrix
aa_conf_logreg <- confusionMatrix(aaxtab, positive = "Pos")
aa_conf_logreg
###### training performance
aaPred2 <- predict(my_learning_datasets_results_logistic[[3]], 
                   aa_train_data[,-ncol(aa_train_data)])
aaTruth2 <- aa_train_data$label_set
aaxtab2 <- table(aaPred2, aaTruth2)
# confusion matrix
aa_conf_logreg2 <- confusionMatrix(aaxtab2, positive = "Pos")
aa_conf_logreg2

################################################# # 4
aanzv <- nearZeroVar(my_learning_datasets[[4]], saveMetrics= F)
aa_all_data <- my_learning_datasets[[4]][,-aanzv]

#  Identifying Correlated Predictors
aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
summary(aadescrCor[upper.tri(aadescrCor)])
aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
colnames(aa_all_data)[aahighlyCorDescr]
aa_all_data <- aa_all_data[,-aahighlyCorDescr]




aa_train_data <- aa_all_data[my_data_partitions[[1]], ]
aa_test_data <-  aa_all_data[-my_data_partitions[[1]], ]

aactrl <- trainControl(method = "repeatedcv", repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       #savePredictions = "final",
                       ## new option here:
                       sampling = "down")


aalogis_down <- train(label_set ~ ., data = aa_train_data,
                      method = "regLogistic",
                      #preProc = c("center", "scale"),
                      #nbagg = 50,
                      metric = "ROC",
                      trControl = aactrl)
my_learning_datasets_results_logistic[[4]] <- aalogis_down
###### test performance
aaPred <- predict(my_learning_datasets_results_logistic[[4]], 
                  aa_test_data[,-ncol(aa_test_data)])
aaTruth <- aa_test_data$label_set
aaxtab <- table(aaPred, aaTruth)
# confusion matrix
aa_conf_logreg <- confusionMatrix(aaxtab, positive = "Pos")
aa_conf_logreg
###### training performance
aaPred2 <- predict(my_learning_datasets_results_logistic[[4]], 
                   aa_train_data[,-ncol(aa_train_data)])
aaTruth2 <- aa_train_data$label_set
aaxtab2 <- table(aaPred2, aaTruth2)
# confusion matrix
aa_conf_logreg2 <- confusionMatrix(aaxtab2, positive = "Pos")
aa_conf_logreg2
ncol(aa_all_data)
################################################# # 5
aanzv <- nearZeroVar(my_learning_datasets[[5]], saveMetrics= F)
aa_all_data <- my_learning_datasets[[5]][,-aanzv]

#  Identifying Correlated Predictors
aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
summary(aadescrCor[upper.tri(aadescrCor)])
aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
colnames(aa_all_data)[aahighlyCorDescr]
aa_all_data <- aa_all_data[,-aahighlyCorDescr]

ncol(aa_all_data)

aa_train_data <- aa_all_data[my_data_partitions[[1]], ]
aa_test_data <-  aa_all_data[-my_data_partitions[[1]], ]

aactrl <- trainControl(method = "repeatedcv", repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = defaultSummary,
                       #savePredictions = "final",
                       ## new option here:
                       sampling = "down")
aalogis_down <- train(label_set ~ ., data = aa_train_data,
                      method = "svmRadial",
                      preProc = c("center", "scale"),
                      #nbagg = 50,
                      metric = "Kappa",
                      trControl = aactrl)
my_learning_datasets_results_logistic[[5]] <- aalogis_down
###### test performance
aaPred <- predict(aalogis_down, 
                  aa_test_data[,-ncol(aa_test_data)])
aaTruth <- aa_test_data$label_set
aaxtab <- table(aaPred, aaTruth)
# confusion matrix
aa_conf_logreg <- confusionMatrix(aaxtab, positive = "Pos")
aa_conf_logreg
###### training performance
aaPred2 <- predict(my_learning_datasets_results_logistic[[5]], 
                   aa_train_data[,-ncol(aa_train_data)])
aaTruth2 <- aa_train_data$label_set
aaxtab2 <- table(aaPred2, aaTruth2)
# confusion matrix
aa_conf_logreg2 <- confusionMatrix(aaxtab2, positive = "Pos")
aa_conf_logreg2

################################################# # 6
aanzv <- nearZeroVar(my_learning_datasets[[6]], saveMetrics= F)
aa_all_data <- my_learning_datasets[[6]][,-aanzv]

#  Identifying Correlated Predictors
aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
summary(aadescrCor[upper.tri(aadescrCor)])
aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
colnames(aa_all_data)[aahighlyCorDescr]
aa_all_data <- aa_all_data[,-aahighlyCorDescr]
ncol(aa_all_data)
aa_train_data <- aa_all_data[my_data_partitions[[1]], ]
aa_test_data <-  aa_all_data[-my_data_partitions[[1]], ]

aactrl <- trainControl(method = "repeatedcv", repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       #savePredictions = "final",
                       ## new option here:
                       sampling = "down")
aalogis_down <- train(label_set ~ ., data = aa_train_data,
                      method = "regLogistic",
                      #preProc = c("center", "scale"),
                      #nbagg = 50,
                      metric = "ROC",
                      trControl = aactrl)
my_learning_datasets_results_logistic[[6]] <- aalogis_down
###### test performance
aaPred <- predict(my_learning_datasets_results_logistic[[6]], 
                  aa_test_data[,-ncol(aa_test_data)])
aaTruth <- aa_test_data$label_set
aaxtab <- table(aaPred, aaTruth)
# confusion matrix
aa_conf_logreg <- confusionMatrix(aaxtab, positive = "Pos")
aa_conf_logreg
###### training performance
aaPred2 <- predict(my_learning_datasets_results_logistic[[6]], 
                   aa_train_data[,-ncol(aa_train_data)])
aaTruth2 <- aa_train_data$label_set
aaxtab2 <- table(aaPred2, aaTruth2)
# confusion matrix
aa_conf_logreg2 <- confusionMatrix(aaxtab2, positive = "Pos")
aa_conf_logreg2

################################################# # 7
aanzv <- nearZeroVar(my_learning_datasets[[7]], saveMetrics= F)
aa_all_data <- my_learning_datasets[[7]][,-aanzv]

#  Identifying Correlated Predictors
aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
summary(aadescrCor[upper.tri(aadescrCor)])
aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
colnames(aa_all_data)[aahighlyCorDescr]
aa_all_data <- aa_all_data[,-aahighlyCorDescr]

ncol(aa_all_data)
ncol(my_learning_datasets[[7]])
aa_train_data <- aa_all_data[my_data_partitions[[1]], ]
aa_test_data <-  aa_all_data[-my_data_partitions[[1]], ]

aactrl <- trainControl(method = "repeatedcv", repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       #savePredictions = "final",
                       ## new option here:
                       sampling = "down")
aalogis_down <- train(label_set ~ ., data = aa_train_data,
                      method = "regLogistic",
                      #preProc = c("center", "scale"),
                      #nbagg = 50,
                      metric = "ROC",
                      trControl = aactrl)
my_learning_datasets_results_logistic[[7]] <- aalogis_down
###### test performance
aaPred <- predict(my_learning_datasets_results_logistic[[7]], 
                  aa_test_data[,-ncol(aa_test_data)])
aaTruth <- aa_test_data$label_set
aaxtab <- table(aaPred, aaTruth)
# confusion matrix
aa_conf_logreg <- confusionMatrix(aaxtab, positive = "Pos")
aa_conf_logreg
###### training performance
aaPred2 <- predict(my_learning_datasets_results_logistic[[7]], 
                   aa_train_data[,-ncol(aa_train_data)])
aaTruth2 <- aa_train_data$label_set
aaxtab2 <- table(aaPred2, aaTruth2)
# confusion matrix
aa_conf_logreg2 <- confusionMatrix(aaxtab2, positive = "Pos")
aa_conf_logreg2
#################################################
#################################################
#################################################
# train a logistic reg classifier with kappa as the metric for hyperparameter choice
# plot the ROC, PR curves for test and training
# do this for all 8 feature settings + have only length as the feature.
# next ~50 lines are commented out to prevent accidental run
# my_learning_datasets_results_logistic_kappa_down <- list()
# my_learning_datasets_results_svmRadial_kappa_down <- list()
# my_learning_datasets_results_logistic_kappa_up <- list()
# my_learning_datasets_results_svmRadial_kappa_up <- list()
# 
# for(i in 1:length(my_learning_datasets)){
#   aanzv <- nearZeroVar(my_learning_datasets[[i]], saveMetrics= F)
#   if(length(aanzv) > 0){
#     aa_all_data <- my_learning_datasets[[i]][,-aanzv]
#   }else{
#     aa_all_data <- my_learning_datasets[[i]]
#   }
#   
#   
#   #  Identifying Correlated Predictors
#   if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
#     aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
#     aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
#     if(length(aahighlyCorDescr) > 0){
#       aa_all_data <- aa_all_data[,-aahighlyCorDescr]
#     }
#   }
#   
#   #summary(aadescrCor[upper.tri(aadescrCor)])
# 
#   print("#########################################################################################")
#   print(paste0("number of features in original feature set number ", i," named ", names(my_learning_datasets)[i] ,":"))
#   print(ncol(my_learning_datasets[[i]]))
#   print("number of features after removing near zero variance and highly correlated features ( > 0.9)")
#   print(ncol(aa_all_data))
#   print("#########################################################################################")
#   aa_train_data <- aa_all_data[my_data_partitions[[1]], ]
#   aa_test_data <-  aa_all_data[-my_data_partitions[[1]], ]
#   
#   
#   aactrl <- trainControl(method = "repeatedcv",
#                          repeats = 5,
#                          classProbs = TRUE,
#                          summaryFunction = defaultSummary,
#                          savePredictions = TRUE,
#                          ## new option here:
#                          sampling = "down")
#   aalogis_down <- train(label_set ~ ., data = aa_train_data,
#                         method = "regLogistic",
#                         #preProc = c("center", "scale"),
#                         #nbagg = 50,
#                         metric = "Kappa",
#                         trControl = aactrl)
#   aaSVM_down <- train(label_set ~ ., data = aa_train_data,
#                         method = "svmRadial",
#                         preProc = c("center", "scale"),
#                         #nbagg = 50,
#                         metric = "Kappa",
#                         trControl = aactrl)
#   aactrl$sampling <- "up"
#   aalogis_up <- train(label_set ~ ., data = aa_train_data,
#                         method = "regLogistic",
#                         #preProc = c("center", "scale"),
#                         #nbagg = 50,
#                         metric = "Kappa",
#                         trControl = aactrl)
#   aaSVM_up <- train(label_set ~ ., data = aa_train_data,
#                       method = "svmRadial",
#                       preProc = c("center", "scale"),
#                       #nbagg = 50,
#                       metric = "Kappa",
#                       trControl = aactrl)
#   my_learning_datasets_results_logistic_kappa_down[[i]] <- aalogis_down
#   my_learning_datasets_results_svmRadial_kappa_down[[i]] <- aaSVM_down
#   my_learning_datasets_results_logistic_kappa_up[[i]] <- aalogis_up
#   my_learning_datasets_results_svmRadial_kappa_up[[i]] <- aaSVM_up
#   
#   
#   
# }

# create an ROC and a PR plot for each feature set, containing four models: logreg (up/down) and SVm (up/down)
aa_gglist <- list()
for(i in 1:length(my_learning_datasets)){
  aanzv <- nearZeroVar(my_learning_datasets[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets[[i]]
  }
  
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  print("#########################################################################################")
  print(paste0("number of features in original feature set number ", i," named ", names(my_learning_datasets)[i] ,":"))
  print(ncol(my_learning_datasets[[i]]))
  print("number of features after removing near zero variance and highly correlated features ( > 0.9)")
  print(ncol(aa_all_data))
  print("#########################################################################################")
  aa_train_data <- aa_all_data[my_data_partitions[[1]], ]
  aa_test_data <-  aa_all_data[-my_data_partitions[[1]], ]
  aaindex_pos <- aa_test_data$label_set == "Pos"
  aaindex_neg <- aa_test_data$label_set == "Neg"
  aapredlogdown <- predict(my_learning_datasets_results_logistic_kappa_down[[i]], 
                           aa_test_data, type = "prob")
  aapredlogup <- predict(my_learning_datasets_results_logistic_kappa_up[[i]], 
                           aa_test_data, type = "prob")
  aapredsvmdown <- predict(my_learning_datasets_results_svmRadial_kappa_down[[i]], 
                           aa_test_data, type = "prob")
  aapredsvmup <- predict(my_learning_datasets_results_svmRadial_kappa_up[[i]], 
                           aa_test_data, type = "prob")
  aapreds_list <- list(aapredlogdown$Pos, aapredlogup$Pos, aapredsvmdown$Pos, aapredsvmup$Pos)
  # List of actual values (same for all)
  m <- length(aapreds_list)
  aactuals_list <- rep(list(aa_test_data$label_set), m)
  # Plot the ROC curves
  aapred <- prediction(aapreds_list, aactuals_list)
  aauc <- performance(aapred, measure = "auc")
  aarocs <- performance(aapred, "tpr", "fpr")
  print(aauc@y.values)
  plot(aarocs, col = as.list(1:m), main = paste0("Test Set ROC Curves for feature set: ", names(my_learning_datasets)[i]))
  abline(a = 0, b = 1, col = 6, lty = 3, lwd = 0.7)
  legend(x = "bottomright", cex = 0.6,
         legend = c(paste0("Log_reg_down ", format(aauc@y.values[[1]], digits = 2)),
                    paste0("Log_reg_up ", format(aauc@y.values[[2]], digits = 2)),
                    paste0("SVM_down ", format(aauc@y.values[[3]], digits = 2)),
                    paste0("SVM_up ", format(aauc@y.values[[4]], digits = 2))),
         fill = 1:m)
  aa_model_list <- list(Log_reg_down = my_learning_datasets_results_logistic_kappa_down[[i]],
                        Log_reg_up = my_learning_datasets_results_logistic_kappa_up[[i]],
                        SVM_down = my_learning_datasets_results_svmRadial_kappa_down[[i]],
                        SVM_up = my_learning_datasets_results_svmRadial_kappa_up[[i]])
  aamodel_list_pr <- aa_model_list %>%
    map(.f = calc_auprc, data = aa_test_data)
  aamodel_list_pr_auc <- aamodel_list_pr %>%
    map(function(the_mod) the_mod$auc.integral)
  aamodel_list_pr_auc <- unlist(aamodel_list_pr_auc)
  
  aaresults_list_pr <- list(NA)
  num_mod <- 1
  
  for(the_pr in aamodel_list_pr){
    
    aaresults_list_pr[[num_mod]] <- data_frame(recall = the_pr$curve[, 1],
                                             precision = the_pr$curve[, 2],
                                             model = names(aamodel_list_pr)[num_mod])
    
    num_mod <- num_mod + 1
    
  }
  
  aaresults_df_pr <- bind_rows(aaresults_list_pr)
  
  aacustom_col <- c("#000000", "#009E73", "#0072B2", "#D55E00")
  
  aa_gglist[[i]] <- ggplot(aes(x = recall, y = precision, group = model), data = aaresults_df_pr) +
    geom_line(aes(color = model), size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = aacustom_col, labels = paste0(names(aa_model_list) , format(aamodel_list_pr_auc, digits = 2))) +
    geom_abline(intercept = sum(aa_test_data$label_set == "Pos")/nrow(aa_test_data),
                slope = 0, color = "gray", size = 1) +
    ggtitle(paste0("PR curve ", names(my_learning_datasets)[i]))+
    theme_bw()
  
}
names(aa_gglist) <- names(my_learning_datasets)

aapred <- predict(my_learning_datasets_results_logistic[[7]], 
                  aa_test_data[,-ncol(aa_test_data)], type = "prob")
aaprcurve <- pr.curve(aapred$Pos[aaindex_pos], aapred$Pos[aaindex_neg], curve = TRUE)
plot(aaprcurve)
line()
plot.roc(aa_test_data$label_set, aapred$Pos, print.auc=T)

####################################################################################################
for(i in 1:length(my_learning_datasets)){
  print(paste("Confusion matrix for test set ", names(my_learning_datasets)[i]))
  aanzv <- nearZeroVar(my_learning_datasets[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets[[i]]
  }
  
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  print("#########################################################################################")
  print(paste0("number of features in original feature set number ", i," named ", names(my_learning_datasets)[i] ,":"))
  print(ncol(my_learning_datasets[[i]]))
  print("number of features after removing near zero variance and highly correlated features ( > 0.9)")
  print(ncol(aa_all_data))
  print("#########################################################################################")
  aa_train_data <- aa_all_data[my_data_partitions[[1]], ]
  aa_test_data <-  aa_all_data[-my_data_partitions[[1]], ]
  
  aa_model_list <- list(Log_reg_down = my_learning_datasets_results_logistic_kappa_down[[i]],
                        Log_reg_up = my_learning_datasets_results_logistic_kappa_up[[i]],
                        SVM_down = my_learning_datasets_results_svmRadial_kappa_down[[i]],
                        SVM_up = my_learning_datasets_results_svmRadial_kappa_up[[i]])
  for(j in 1:length(aa_model_list)){
    aaPred <- predict(aa_model_list[[j]], 
                      aa_test_data)
    aaTruth <- aa_test_data$label_set
    aaxtab <- table(aaPred, aaTruth)
    # confusion matrix
    aa_conf_logreg <- confusionMatrix(aaxtab, positive = "Pos")
    print(names(aa_model_list)[j])
    print(aa_conf_logreg)
  }
}
###################################################################################################
###################################################################################################
###################################################################################################
# repeating with datasets gathered from 1000bp data. --> first positive_negative sets: 
# universe: intersection of (H3k27+) and (union of two ER ChIP replicates)
# positive: intersection of (universe) and (intersection of two eRNA replicates)
# negative: rest of the universe

# feature sets:
#sum LR
feature_set_sum_LLR_1000bp <- rbind(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[1]], 
                                    Positive_set_seq_list_chopped_scoreMat_list_1000bp[[1]])
colnames(feature_set_sum_LLR_1000bp)[14] <- "NKX3_1"
# adjacency
feature_set_Adjmat_1000bp <- rbind(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[1]], 
                                   Positive_set_seq_list_chopped_adjanScore_list_1000bp[[1]])
label_set <- c(rep("Neg", nrow(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[1]])), 
               rep("Pos", nrow(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[1]])))

# changing the name of NKX3-1 to NKX3_1 so it works in the classification models
aa <- strsplit(colnames(feature_set_Adjmat_1000bp), split = "_vs_")
for(i in 1:length(aa)){
  aa[[i]][aa[[i]] == "NKX3-1"] <- "NKX3_1"
}
for(i in 1:length(aa)){
  aa[[i]] <- paste(aa[[i]][1], aa[[i]][2], sep = "_vs_")
}
colnames(feature_set_Adjmat_1000bp) <- unlist(aa)
# overlap
feature_set_Overlapmat_1000bp <- rbind(Negative_set_seq_list_chopped_OverlapScore_list_1000bp[[1]], 
                                       Positive_set_seq_list_chopped_OverlapScore_list_1000bp[[1]])
colnames(feature_set_Overlapmat_1000bp) <- unlist(aa)


my_learning_datasets_1000bp <- list()
# sum LR only
my_learning_datasets_1000bp[[1]] <- as.data.frame(cbind(feature_set_sum_LLR_1000bp,
                                                 label_set), 
                                           stringsAsFactors =F)
my_learning_datasets_1000bp[[1]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp[[1]])
# adjacancy only
my_learning_datasets_1000bp[[2]] <- as.data.frame(cbind(feature_set_Adjmat_1000bp,
                                                 label_set), 
                                           stringsAsFactors =F)
my_learning_datasets_1000bp[[2]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp[[2]])

# Overlap only
my_learning_datasets_1000bp[[3]] <- as.data.frame(cbind(feature_set_Overlapmat_1000bp,
                                                 label_set), 
                                           stringsAsFactors =F)
my_learning_datasets_1000bp[[3]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp[[3]])
# LLR and adj
my_learning_datasets_1000bp[[4]] <- as.data.frame(cbind(cbind(feature_set_sum_LLR_1000bp,
                                                       feature_set_Adjmat_1000bp),
                                                 label_set), 
                                           stringsAsFactors =F)
my_learning_datasets_1000bp[[4]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp[[4]])
# LLR and overlap
my_learning_datasets_1000bp[[5]] <- as.data.frame(cbind(cbind(feature_set_sum_LLR_1000bp,
                                                       feature_set_Overlapmat_1000bp),
                                                 label_set), 
                                           stringsAsFactors =F)
my_learning_datasets_1000bp[[5]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp[[5]])
# LLR and adj and overlap
my_learning_datasets_1000bp[[6]] <- as.data.frame(cbind(cbind(feature_set_sum_LLR_1000bp,
                                                       cbind(feature_set_Adjmat_1000bp, 
                                                             feature_set_Overlapmat_1000bp)),
                                                 label_set), 
                                           stringsAsFactors =F)
colnames(my_learning_datasets_1000bp[[6]])[34:594] <- paste(colnames(my_learning_datasets_1000bp[[6]])[34:594], "adj", sep = "_at_")
colnames(my_learning_datasets_1000bp[[6]])[595:1155] <- paste(colnames(my_learning_datasets_1000bp[[6]])[595:1155], "ovl", sep = "_at_")
my_learning_datasets_1000bp[[6]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp[[6]])

names(my_learning_datasets_1000bp) <- c("Sum_LR","Adjacency", "Overlap", 
                                 "Sum_LR_plus_Adjacency", "Sum_LR_plus_Overlap", 
                                 "Sum_LR_plus_Adjacency_plus_Overlap")
#########################################################################################################
# my_data_partitions <- list()
# my_learning_datasets_results_logistic <- list()

# next ~100 lines are commented out to prevent accidental run
# my_data_partitions[[2]] <- createDataPartition(y = my_learning_datasets_1000bp[[1]][, ncol(my_learning_datasets_1000bp[[1]])], times = 1, p = 0.75,list = F)


# my_learning_datasets_results_logistic_kappa_down_1000bp <- list()
# my_learning_datasets_results_logistic_kappa_up_1000bp <- list()
# my_learning_datasets_results_svmLinear_kappa_down_1000bp <- list()
# my_learning_datasets_results_svmLinear_kappa_up_1000bp <- list()
# my_learning_datasets_results_RF_kappa_down_1000bp <- list()
# my_learning_datasets_results_RF_kappa_up_1000bp <- list()
# 
# # 
# for(i in 3:length(my_learning_datasets_1000bp)){
#   aanzv <- nearZeroVar(my_learning_datasets_1000bp[[i]], saveMetrics= F)
#   if(length(aanzv) > 0){
#     aa_all_data <- my_learning_datasets_1000bp[[i]][,-aanzv]
#   }else{
#     aa_all_data <- my_learning_datasets_1000bp[[i]]
#   }
# 
# 
#   #  Identifying Correlated Predictors
#   if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
#     aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
#     aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
#     if(length(aahighlyCorDescr) > 0){
#       aa_all_data <- aa_all_data[,-aahighlyCorDescr]
#     }
#   }
#   print("#########################################################################################")
#   print(paste0("number of features in original feature set number ", i," named ", names(my_learning_datasets_1000bp)[i] ,":"))
#   print(ncol(my_learning_datasets_1000bp[[i]]))
#   print("number of features after removing near zero variance and highly correlated features ( > 0.9)")
#   print(ncol(aa_all_data))
#   print("#########################################################################################")
#   aa_train_data <- aa_all_data[my_data_partitions[[2]], ]
#   aa_test_data <-  aa_all_data[-my_data_partitions[[2]], ]
# 
# 
#   aactrl <- trainControl(method = "repeatedcv",
#                          repeats = 3,
#                          classProbs = TRUE,
#                          summaryFunction = defaultSummary,
#                          savePredictions = TRUE,
#                          ## new option here:
#                          sampling = "down")
#   print("logis down")
#   aalogis_down <- train(label_set ~ ., data = aa_train_data,
#                         method = "regLogistic",
#                         #preProc = c("center", "scale"),
#                         #nbagg = 50,
#                         metric = "Kappa",
#                         trControl = aactrl)
#   print("svm down")
#   aaSVM_down <- train(label_set ~ ., data = aa_train_data,
#                         method = "svmLinear3",
#                         preProc = c("center", "scale"),
#                         #nbagg = 50,
#                         metric = "Kappa",
#                         trControl = aactrl)
#   print("RF down")
#   aaRF_down <- train(label_set ~ ., data = aa_train_data,
#                       method = "rf",
#                       #preProc = c("center", "scale"),
#                       #nbagg = 50,
#                       metric = "Kappa",
#                       trControl = aactrl)
#   aactrl$sampling <- "up"
#   print("logis up")
#   aalogis_up <- train(label_set ~ ., data = aa_train_data,
#                         method = "regLogistic",
#                         #preProc = c("center", "scale"),
#                         #nbagg = 50,
#                         metric = "Kappa",
#                         trControl = aactrl)
#   print("svm up")
#   aaSVM_up <- train(label_set ~ ., data = aa_train_data,
#                       method = "svmLinear3",
#                       preProc = c("center", "scale"),
#                       #nbagg = 50,
#                       metric = "Kappa",
#                       trControl = aactrl)
#   # aaRF_up <- train(label_set ~ ., data = aa_train_data,
#   #                  method = "rf",
#   #                  #preProc = c("center", "scale"),
#   #                  #nbagg = 50,
#   #                  metric = "Kappa",
#   #                  trControl = aactrl)
#   my_learning_datasets_results_logistic_kappa_down_1000bp[[i]] <- aalogis_down
#   my_learning_datasets_results_svmLinear_kappa_down_1000bp[[i]] <- aaSVM_down
#   my_learning_datasets_results_logistic_kappa_up_1000bp[[i]] <- aalogis_up
#   my_learning_datasets_results_svmLinear_kappa_up_1000bp[[i]] <- aaSVM_up
#   my_learning_datasets_results_RF_kappa_down_1000bp[[i]] <- aaRF_down
#   #my_learning_datasets_results_RF_kappa_up_1000bp[[i]] <- aaRF_up
#   
# }
# save(list = c("my_learning_datasets_results_logistic_kappa_down_1000bp",
#               "my_learning_datasets_results_svmLinear_kappa_down_1000bp", 
#               "my_learning_datasets_results_logistic_kappa_up_1000bp", 
#               "my_learning_datasets_results_svmLinear_kappa_up_1000bp",
#               "my_learning_datasets_results_RF_kappa_down_1000bp", 
#               "my_learning_datasets_results_RF_kappa_up_1000bp"),
#      file = "classification_results_set_1.RData")

load("classification_results_set_1.RData")
# plotting
aa_gglist <- list()

for(i in 2:length(my_learning_datasets_1000bp)){
  aanzv <- nearZeroVar(my_learning_datasets_1000bp[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets_1000bp[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets_1000bp[[i]]
  }
  
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  print("#########################################################################################")
  print(paste0("number of features in original feature set number ", i," named ", names(my_learning_datasets_1000bp)[i] ,":"))
  print(ncol(my_learning_datasets_1000bp[[i]]))
  print("number of features after removing near zero variance and highly correlated features ( > 0.9)")
  print(ncol(aa_all_data))
  print("#########################################################################################")
  aa_train_data <- aa_all_data[my_data_partitions[[2]], ]
  aa_test_data <-  aa_all_data[-my_data_partitions[[2]], ]
  aaindex_pos <- aa_test_data$label_set == "Pos"
  aaindex_neg <- aa_test_data$label_set == "Neg"
  aapredlogdown <- predict(my_learning_datasets_results_logistic_kappa_down_1000bp[[i]], 
                           aa_test_data, type = "prob")
  aapredlogup <- predict(my_learning_datasets_results_logistic_kappa_up_1000bp[[i]], 
                         aa_test_data, type = "prob")
  # aapredsvmdown <- predict(my_learning_datasets_results_svmLinear_kappa_down_1000bp[[i]], 
  #                          aa_test_data, type = "raw")
  # aapredsvmup <- predict(my_learning_datasets_results_svmLinear_kappa_up_1000bp[[i]], 
  #                        aa_test_data, type = "raw")
  aarfdown <- predict(my_learning_datasets_results_RF_kappa_down_1000bp[[i]], 
                      aa_test_data, type = "prob")
  # aarfup <- predict(my_learning_datasets_results_RF_kappa_up_1000bp[[i]], 
  #                   aa_test_data, type = "prob")
  aapreds_list <- list(aapredlogdown$Pos
                       ,aapredlogup$Pos
                       # ,aapredsvmdown$Pos
                       # ,aapredsvmup$Pos
                       ,aarfdown$Pos
                       #,aarfup$Pos
                       )
  # List of actual values (same for all)
  m <- length(aapreds_list)
  aactuals_list <- rep(list(aa_test_data$label_set), m)
  # Plot the ROC curves
  aapred <- prediction(aapreds_list, aactuals_list)
  aauc <- performance(aapred, measure = "auc")
  aarocs <- performance(aapred, "tpr", "fpr")
  print(aauc@y.values)
  plot(aarocs, col = as.list(1:m), main = paste0("Test Set ROC Curves for feature set: ", names(my_learning_datasets_1000bp)[i]))
  abline(a = 0, b = 1, col = 6, lty = 3, lwd = 0.7)
  legend(x = "bottomright", cex = 0.5,
         legend = c(paste0("Log_reg_down ", format(aauc@y.values[[1]], digits = 2)),
                    paste0("Log_reg_up ", format(aauc@y.values[[2]], digits = 2)),
                    #paste0("SVM_down ", format(aauc@y.values[[3]], digits = 2)),
                    #paste0("SVM_up ", format(aauc@y.values[[4]], digits = 2)),
                    paste0("RF_down ", format(aauc@y.values[[3]], digits = 2))
                    #,paste0("RF_up ", format(aauc@y.values[[4]], digits = 2))
                    ),
         
         fill = 1:m)
  aa_model_list <- list(Log_reg_down = my_learning_datasets_results_logistic_kappa_down_1000bp[[i]],
                        Log_reg_up = my_learning_datasets_results_logistic_kappa_up_1000bp[[i]],
                        #SVM_down = my_learning_datasets_results_svmLinear_kappa_down_1000bp[[i]],
                        #SVM_up = my_learning_datasets_results_svmLinear_kappa_up_1000bp[[i]], 
                        RF_down = my_learning_datasets_results_RF_kappa_down_1000bp[[i]]
                        #,RF_up = my_learning_datasets_results_RF_kappa_up_1000bp[[i]]
                        )
  aamodel_list_pr <- aa_model_list %>%
    map(.f = calc_auprc, data = aa_test_data)
  aamodel_list_pr_auc <- aamodel_list_pr %>%
    map(function(the_mod) the_mod$auc.integral)
  aamodel_list_pr_auc <- unlist(aamodel_list_pr_auc)
  
  aaresults_list_pr <- list(NA)
  num_mod <- 1
  
  for(the_pr in aamodel_list_pr){
    
    aaresults_list_pr[[num_mod]] <- data_frame(recall = the_pr$curve[, 1],
                                               precision = the_pr$curve[, 2],
                                               model = names(aamodel_list_pr)[num_mod])
    
    num_mod <- num_mod + 1
    
  }
  
  aaresults_df_pr <- bind_rows(aaresults_list_pr)
  
  aacustom_col <- colorspace::rainbow_hcl(length(aa_model_list))
  
  aa_gglist[[i]] <- ggplot(aes(x = recall, y = precision, group = model), data = aaresults_df_pr) +
    geom_line(aes(color = model), size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = aacustom_col, labels = paste0(names(aa_model_list) , format(aamodel_list_pr_auc, digits = 2))) +
    geom_abline(intercept = sum(aa_test_data$label_set == "Pos")/nrow(aa_test_data),
                slope = 0, color = "gray", size = 1) +
    ggtitle(paste0("PR curve ", names(my_learning_datasets_1000bp)[i]))+
    theme_bw()
  
}
names(aa_gglist) <- names(my_learning_datasets_1000bp)
Prec_recall_curve_posneg_1 <- aa_gglist
aa_gglist[[6]]
# aapred <- predict(my_learning_datasets_results_logistic[[7]], 
#                   aa_test_data[,-ncol(aa_test_data)], type = "prob")
# aaprcurve <- pr.curve(aapred$Pos[aaindex_pos], aapred$Pos[aaindex_neg], curve = TRUE)
# plot(aaprcurve)
# line()
# plot.roc(aa_test_data$label_set, aapred$Pos, print.auc=T)

####################################################################################################
for(i in 1:length(my_learning_datasets_1000bp)){
  print(paste("Confusion matrix for test set ", names(my_learning_datasets_1000bp)[i]))
  aanzv <- nearZeroVar(my_learning_datasets_1000bp[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets_1000bp[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets_1000bp[[i]]
  }
  
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  print("#########################################################################################")
  print(paste0("number of features in original feature set number ", i," named ", names(my_learning_datasets_1000bp)[i] ,":"))
  print(ncol(my_learning_datasets_1000bp[[i]]))
  print("number of features after removing near zero variance and highly correlated features ( > 0.9)")
  print(ncol(aa_all_data))
  print("#########################################################################################")
  aa_train_data <- aa_all_data[my_data_partitions[[2]], ]
  aa_test_data <-  aa_all_data[-my_data_partitions[[2]], ]
  
  aa_model_list <- list(Log_reg_down = my_learning_datasets_results_logistic_kappa_down_1000bp[[i]],
                        Log_reg_up = my_learning_datasets_results_logistic_kappa_up_1000bp[[i]],
                        SVM_down = my_learning_datasets_results_svmLinear_kappa_down_1000bp[[i]],
                        SVM_up = my_learning_datasets_results_svmLinear_kappa_up_1000bp[[i]], 
                        RF_down = my_learning_datasets_results_RF_kappa_down_1000bp[[i]]
                        ,RF_up = my_learning_datasets_results_RF_kappa_up_1000bp[[i]]
                        )
  for(j in 1:length(aa_model_list)){
    aaPred <- predict(aa_model_list[[j]], 
                      aa_test_data)
    aaTruth <- aa_test_data$label_set
    aaxtab <- table(aaPred, aaTruth)
    # confusion matrix
    aa_conf_logreg <- confusionMatrix(aaxtab, positive = "Pos")
    print(names(aa_model_list)[j])
    print(aa_conf_logreg)
  }
}
###################################################################################################
###################################################################################################
###################################################################################################

# repeating with datasets gathered from 1000bp data. --> second positive_negative sets: 
# universe: intersection of (H3k27+) and (union of two ER ChIP replicates)
# positive: intersection of (universe) and (intersection of two eRNA replicates)
# negative: universe minus the union of two eRNA replicates


# feature sets:
#sum LR
feature_set_sum_LLR_1000bp_2 <- rbind(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[2]], 
                                    Positive_set_seq_list_chopped_scoreMat_list_1000bp[[2]])
colnames(feature_set_sum_LLR_1000bp_2)[14] <- "NKX3_1"
# adjacency
feature_set_Adjmat_1000bp_2 <- rbind(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[2]], 
                                   Positive_set_seq_list_chopped_adjanScore_list_1000bp[[2]])
label_set <- c(rep("Neg", nrow(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[2]])), 
               rep("Pos", nrow(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[2]])))

# changing the name of NKX3-1 to NKX3_1 so it works in the classification models
aa <- strsplit(colnames(feature_set_Adjmat_1000bp_2), split = "_vs_")
for(i in 1:length(aa)){
  aa[[i]][aa[[i]] == "NKX3-1"] <- "NKX3_1"
}
for(i in 1:length(aa)){
  aa[[i]] <- paste(aa[[i]][1], aa[[i]][2], sep = "_vs_")
}
colnames(feature_set_Adjmat_1000bp_2) <- unlist(aa)
# overlap
feature_set_Overlapmat_1000bp_2 <- rbind(Negative_set_seq_list_chopped_OverlapScore_list_1000bp[[2]], 
                                       Positive_set_seq_list_chopped_OverlapScore_list_1000bp[[2]])
colnames(feature_set_Overlapmat_1000bp_2) <- unlist(aa)


my_learning_datasets_1000bp_2 <- list()
# sum LR only
my_learning_datasets_1000bp_2[[1]] <- as.data.frame(cbind(feature_set_sum_LLR_1000bp_2,
                                                        label_set), 
                                                  stringsAsFactors =F)
my_learning_datasets_1000bp_2[[1]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_2[[1]])
# adjacancy only
my_learning_datasets_1000bp_2[[2]] <- as.data.frame(cbind(feature_set_Adjmat_1000bp_2,
                                                        label_set), 
                                                  stringsAsFactors =F)
my_learning_datasets_1000bp_2[[2]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_2[[2]])

# Overlap only
my_learning_datasets_1000bp_2[[3]] <- as.data.frame(cbind(feature_set_Overlapmat_1000bp_2,
                                                        label_set), 
                                                  stringsAsFactors =F)
my_learning_datasets_1000bp_2[[3]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_2[[3]])
# LLR and adj
my_learning_datasets_1000bp_2[[4]] <- as.data.frame(cbind(cbind(feature_set_sum_LLR_1000bp_2,
                                                              feature_set_Adjmat_1000bp_2),
                                                        label_set), 
                                                  stringsAsFactors =F)
my_learning_datasets_1000bp_2[[4]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_2[[4]])
# LLR and overlap
my_learning_datasets_1000bp_2[[5]] <- as.data.frame(cbind(cbind(feature_set_sum_LLR_1000bp_2,
                                                              feature_set_Overlapmat_1000bp_2),
                                                        label_set), 
                                                  stringsAsFactors =F)
my_learning_datasets_1000bp_2[[5]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_2[[5]])
# LLR and adj and overlap
my_learning_datasets_1000bp_2[[6]] <- as.data.frame(cbind(cbind(feature_set_sum_LLR_1000bp_2,
                                                              cbind(feature_set_Adjmat_1000bp_2, 
                                                                    feature_set_Overlapmat_1000bp_2)),
                                                        label_set), 
                                                  stringsAsFactors =F)
colnames(my_learning_datasets_1000bp_2[[6]])[34:594] <- paste(colnames(my_learning_datasets_1000bp_2[[6]])[34:594], "adj", sep = "_at_")
colnames(my_learning_datasets_1000bp_2[[6]])[595:1155] <- paste(colnames(my_learning_datasets_1000bp_2[[6]])[595:1155], "ovl", sep = "_at_")
my_learning_datasets_1000bp_2[[6]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_2[[6]])

names(my_learning_datasets_1000bp_2) <- c("Sum_LR","Adjacency", "Overlap", 
                                        "Sum_LR_plus_Adjacency", "Sum_LR_plus_Overlap", 
                                        "Sum_LR_plus_Adjacency_plus_Overlap")


load("classification_results_set_1.RData")

# next ~100 lines are commented out to prevent accidental run
# my_data_partitions[[3]] <- createDataPartition(y = my_learning_datasets_1000bp_2[[1]][, ncol(my_learning_datasets_1000bp_2[[1]])], times = 1, p = 0.75,list = F)

# my_learning_datasets_results_logistic_kappa_down_1000bp_2 <- list()
# my_learning_datasets_results_logistic_kappa_up_1000bp_2 <- list()
# my_learning_datasets_results_svmLinear_kappa_down_1000bp_2 <- list()
# my_learning_datasets_results_svmLinear_kappa_up_1000bp_2 <- list()
# my_learning_datasets_results_RF_kappa_down_1000bp_2 <- list()
# my_learning_datasets_results_RF_kappa_up_1000bp_2 <- list()
# 
# # 
# for(i in 3:length(my_learning_datasets_1000bp_2)){
#   aanzv <- nearZeroVar(my_learning_datasets_1000bp_2[[i]], saveMetrics= F)
#   if(length(aanzv) > 0){
#     aa_all_data <- my_learning_datasets_1000bp_2[[i]][,-aanzv]
#   }else{
#     aa_all_data <- my_learning_datasets_1000bp_2[[i]]
#   }
#   
#   
#   #  Identifying Correlated Predictors
#   if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
#     aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
#     aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
#     if(length(aahighlyCorDescr) > 0){
#       aa_all_data <- aa_all_data[,-aahighlyCorDescr]
#     }
#   }
#   print("#########################################################################################")
#   print(paste0("number of features in original feature set number ", i," named ", names(my_learning_datasets_1000bp_2)[i] ,":"))
#   print(ncol(my_learning_datasets_1000bp_2[[i]]))
#   print("number of features after removing near zero variance and highly correlated features ( > 0.9)")
#   print(ncol(aa_all_data))
#   print("#########################################################################################")
#   aa_train_data <- aa_all_data[my_data_partitions[[3]], ]
#   aa_test_data <-  aa_all_data[-my_data_partitions[[3]], ]
#   
#   
#   aactrl <- trainControl(method = "repeatedcv",
#                          repeats = 3,
#                          classProbs = TRUE,
#                          summaryFunction = defaultSummary,
#                          savePredictions = TRUE,
#                          ## new option here:
#                          sampling = "down")
#   print("aalogis_down")
#   aalogis_down <- train(label_set ~ ., data = aa_train_data,
#                         method = "regLogistic",
#                         #preProc = c("center", "scale"),
#                         #nbagg = 50,
#                         metric = "Kappa",
#                         trControl = aactrl)
#   print("aaSVM_down")
#   aaSVM_down <- train(label_set ~ ., data = aa_train_data,
#                       method = "svmLinear3",
#                       preProc = c("center", "scale"),
#                       #nbagg = 50,
#                       metric = "Kappa",
#                       trControl = aactrl)
#   print("aaRF_down")
#   aaRF_down <- train(label_set ~ ., data = aa_train_data,
#                      method = "rf",
#                      #preProc = c("center", "scale"),
#                      #nbagg = 50,
#                      metric = "Kappa",
#                      trControl = aactrl)
#   aactrl$sampling <- "up"
#   print("aalogis_up")
#   aalogis_up <- train(label_set ~ ., data = aa_train_data,
#                       method = "regLogistic",
#                       #preProc = c("center", "scale"),
#                       #nbagg = 50,
#                       metric = "Kappa",
#                       trControl = aactrl)
#   print("aaSVM_up")
#   aaSVM_up <- train(label_set ~ ., data = aa_train_data,
#                     method = "svmLinear3",
#                     preProc = c("center", "scale"),
#                     #nbagg = 50,
#                     metric = "Kappa",
#                     trControl = aactrl)
#   # print("aaRF_up")
#   # aaRF_up <- train(label_set ~ ., data = aa_train_data,
#   #                  method = "rf",
#   #                  #preProc = c("center", "scale"),
#   #                  #nbagg = 50,
#   #                  metric = "Kappa",
#   #                  trControl = aactrl)
#   my_learning_datasets_results_logistic_kappa_down_1000bp_2[[i]] <- aalogis_down
#   my_learning_datasets_results_svmLinear_kappa_down_1000bp_2[[i]] <- aaSVM_down
#   my_learning_datasets_results_logistic_kappa_up_1000bp_2[[i]] <- aalogis_up
#   my_learning_datasets_results_svmLinear_kappa_up_1000bp_2[[i]] <- aaSVM_up
#   my_learning_datasets_results_RF_kappa_down_1000bp_2[[i]] <- aaRF_down
#   # my_learning_datasets_results_RF_kappa_up_1000bp_2[[i]] <- aaRF_up
#   
# }
# 
# save(list = c("my_learning_datasets_results_logistic_kappa_down_1000bp_2",
#               "my_learning_datasets_results_svmLinear_kappa_down_1000bp_2", 
#               "my_learning_datasets_results_logistic_kappa_up_1000bp_2", 
#               "my_learning_datasets_results_svmLinear_kappa_up_1000bp_2",
#               "my_learning_datasets_results_RF_kappa_down_1000bp_2", 
#               "my_learning_datasets_results_RF_kappa_up_1000bp_2"),
#      file = "classification_results_set_2.RData")

load("classification_results_set_2.RData")
# plotting
aa_gglist_2 <- list()

for(i in 2:length(my_learning_datasets_1000bp_2)){
  aanzv <- nearZeroVar(my_learning_datasets_1000bp_2[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets_1000bp_2[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets_1000bp_2[[i]]
  }
  
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  print("#########################################################################################")
  print(paste0("number of features in original feature set number ", i," named ", names(my_learning_datasets_1000bp)[i] ,":"))
  print(ncol(my_learning_datasets_1000bp_2[[i]]))
  print("number of features after removing near zero variance and highly correlated features ( > 0.9)")
  print(ncol(aa_all_data))
  print("#########################################################################################")
  aa_train_data <- aa_all_data[my_data_partitions[[3]], ]
  aa_test_data <-  aa_all_data[-my_data_partitions[[3]], ]
  aaindex_pos <- aa_test_data$label_set == "Pos"
  aaindex_neg <- aa_test_data$label_set == "Neg"
  aapredlogdown <- predict(my_learning_datasets_results_logistic_kappa_down_1000bp_2[[i]], 
                           aa_test_data, type = "prob")
  aapredlogup <- predict(my_learning_datasets_results_logistic_kappa_up_1000bp_2[[i]], 
                         aa_test_data, type = "prob")
  # aapredsvmdown <- predict(my_learning_datasets_results_svmLinear_kappa_down_1000bp_2[[i]], 
  #                          aa_test_data, type = "prob")
  # aapredsvmup <- predict(my_learning_datasets_results_svmLinear_kappa_up_1000bp_2[[i]], 
  #                        aa_test_data, type = "prob")
  aarfdown <- predict(my_learning_datasets_results_RF_kappa_down_1000bp_2[[i]], 
                      aa_test_data, type = "prob")
  # aarfup <- predict(my_learning_datasets_results_RF_kappa_up_1000bp_2[[i]], 
  #                   aa_test_data, type = "prob")
  aapreds_list <- list(aapredlogdown$Pos
                       ,aapredlogup$Pos
                       #, aapredsvmdown$Pos
                       #,aapredsvmup$Pos
                       , aarfdown$Pos
                       #, aarfup$Pos
                       )
  # List of actual values (same for all)
  m <- length(aapreds_list)
  aactuals_list <- rep(list(aa_test_data$label_set), m)
  # Plot the ROC curves
  aapred <- prediction(aapreds_list, aactuals_list)
  aauc <- performance(aapred, measure = "auc")
  aarocs <- performance(aapred, "tpr", "fpr")
  print(aauc@y.values)
  plot(aarocs, col = as.list(1:m), main = paste0("Test Set ROC Curves for feature set: ", names(my_learning_datasets_1000bp_2)[i]))
  abline(a = 0, b = 1, col = 6, lty = 3, lwd = 0.7)
  legend(x = "bottomright", cex = 0.6,
         legend = c(paste0("Log_reg_down ", format(aauc@y.values[[1]], digits = 2)),
                    paste0("Log_reg_up ", format(aauc@y.values[[2]], digits = 2)),
                    #paste0("SVM_down ", format(aauc@y.values[[3]], digits = 2)),
                    #paste0("SVM_up ", format(aauc@y.values[[4]], digits = 2)),
                    paste0("RF_down ", format(aauc@y.values[[3]], digits = 2))
                    #,paste0("RF_up ", format(aauc@y.values[[4]], digits = 2))
                    ),
         fill = 1:m)
  aa_model_list <- list(Log_reg_down = my_learning_datasets_results_logistic_kappa_down_1000bp_2[[i]],
                        Log_reg_up = my_learning_datasets_results_logistic_kappa_up_1000bp_2[[i]],
                        # SVM_down = my_learning_datasets_results_svmLinear_kappa_down_1000bp_2[[i]],
                        # SVM_up = my_learning_datasets_results_svmLinear_kappa_up_1000bp_2[[i]], 
                        RF_down = my_learning_datasets_results_RF_kappa_down_1000bp_2[[i]]
                        #,RF_up = my_learning_datasets_results_RF_kappa_up_1000bp_2[[i]]
                        )
  aamodel_list_pr <- aa_model_list %>%
    map(.f = calc_auprc, data = aa_test_data)
  aamodel_list_pr_auc <- aamodel_list_pr %>%
    map(function(the_mod) the_mod$auc.integral)
  aamodel_list_pr_auc <- unlist(aamodel_list_pr_auc)
  
  aaresults_list_pr <- list(NA)
  num_mod <- 1
  
  for(the_pr in aamodel_list_pr){
    
    aaresults_list_pr[[num_mod]] <- data_frame(recall = the_pr$curve[, 1],
                                               precision = the_pr$curve[, 2],
                                               model = names(aamodel_list_pr)[num_mod])
    
    num_mod <- num_mod + 1
    
  }
  
  aaresults_df_pr <- bind_rows(aaresults_list_pr)
  
  aacustom_col <- colorspace::rainbow_hcl(length(aa_model_list))
  
  aa_gglist_2[[i]] <- ggplot(aes(x = recall, y = precision, group = model), data = aaresults_df_pr) +
    geom_line(aes(color = model), size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = aacustom_col, labels = paste0(names(aa_model_list) , format(aamodel_list_pr_auc, digits = 2))) +
    geom_abline(intercept = sum(aa_test_data$label_set == "Pos")/nrow(aa_test_data),
                slope = 0, color = "gray", size = 1) +
    ggtitle(paste0("PR curve ", names(my_learning_datasets_1000bp_2)[i]))+
    theme_bw()
  
}
names(aa_gglist_2) <- names(my_learning_datasets_1000bp_2)
Prec_recall_curve_posneg_2 <- aa_gglist_2
aa_gglist_2[[6]]
# aapred <- predict(my_learning_datasets_results_logistic[[7]], 
#                   aa_test_data[,-ncol(aa_test_data)], type = "prob")
# aaprcurve <- pr.curve(aapred$Pos[aaindex_pos], aapred$Pos[aaindex_neg], curve = TRUE)
# plot(aaprcurve)
# line()
# plot.roc(aa_test_data$label_set, aapred$Pos, print.auc=T)

####################################################################################################
for(i in 1:length(my_learning_datasets_1000bp_2)){
  print(paste("Confusion matrix for test set ", names(my_learning_datasets_1000bp_2)[i]))
  aanzv <- nearZeroVar(my_learning_datasets_1000bp_2[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets_1000bp_2[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets_1000bp_2[[i]]
  }
  
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  print("#########################################################################################")
  print(paste0("number of features in original feature set number ", i," named ", names(my_learning_datasets_1000bp_2)[i] ,":"))
  print(ncol(my_learning_datasets_1000bp_2[[i]]))
  print("number of features after removing near zero variance and highly correlated features ( > 0.9)")
  print(ncol(aa_all_data))
  print("#########################################################################################")
  aa_train_data <- aa_all_data[my_data_partitions[[3]], ]
  aa_test_data <-  aa_all_data[-my_data_partitions[[3]], ]
  
  aa_model_list <- list(Log_reg_down = my_learning_datasets_results_logistic_kappa_down_1000bp_2[[i]],
                        Log_reg_up = my_learning_datasets_results_logistic_kappa_up_1000bp_2[[i]],
                        SVM_down = my_learning_datasets_results_svmLinear_kappa_down_1000bp_2[[i]],
                        SVM_up = my_learning_datasets_results_svmLinear_kappa_up_1000bp_2[[i]], 
                        RF_down = my_learning_datasets_results_RF_kappa_down_1000bp_2[[i]],
                        RF_up = my_learning_datasets_results_RF_kappa_up_1000bp_2[[i]])
  for(j in 1:length(aa_model_list)){
    aaPred <- predict(aa_model_list[[j]], 
                      aa_test_data)
    aaTruth <- aa_test_data$label_set
    aaxtab <- table(aaPred, aaTruth)
    # confusion matrix
    aa_conf_logreg <- confusionMatrix(aaxtab, positive = "Pos")
    print(names(aa_model_list)[j])
    print(aa_conf_logreg)
  }
}


####################################################################################################
####################################################################################################
####################################################################################################
# create a common heldout set for pos_neg sets 1 and 2.
# run downsampled random forest with different feature sets trained on sets 1 or 2 on this new common set.
# compare performance

# for second set
# my_data_partitions[[5]] <- createDataPartition(y = my_learning_datasets_1000bp_2[[1]][, ncol(my_learning_datasets_1000bp_2[[1]])], times = 1, p = 0.75,list = F)

aa_all_seq_2 <- c(Negative_set_seq_list_1000bp[[2]],
                  Positive_set_seq_list_1000bp[[2]])
aa_all_seq_1 <- c(Negative_set_seq_list_1000bp[[1]],
                  Positive_set_seq_list_1000bp[[1]])

aamatch <- match(aa_all_seq_2[-my_data_partitions[[5]][, 1]], aa_all_seq_1)

# for first set based on the second set
# my_data_partitions[[4]] <- matrix(setdiff(c(1:nrow(my_learning_datasets_1000bp[[1]])), aamatch), ncol = 1)
# 
# # check if the test set is the same.
# all(aa_all_seq_2[-my_data_partitions[[5]][, 1]] %in% aa_all_seq_1[-my_data_partitions[[4]][, 1]])
# 
# 
# my_learning_datasets_results_RF_kappa_down_1000bp_commonset_1 <- list()
# my_learning_datasets_results_RF_kappa_down_1000bp_commonset_2 <- list()
# 
# for(i in 1:length(my_learning_datasets_1000bp_2)){
# 
#   aanzv1 <- nearZeroVar(my_learning_datasets_1000bp[[i]], saveMetrics= F)
#   if(length(aanzv1) > 0){
#     aa_all_data1 <- my_learning_datasets_1000bp[[i]][,-aanzv1]
#   }else{
#     aa_all_data1 <- my_learning_datasets_1000bp[[i]]
#   }
# 
#   #  Identifying Correlated Predictors
#   if(ncol(aa_all_data1[, 1:(ncol(aa_all_data1) - 1)]) > 1){
#     aadescrCor1 <- cor(aa_all_data1[, 1:(ncol(aa_all_data1) - 1)])
#     aahighlyCorDescr1 <- findCorrelation(aadescrCor1, cutoff = .9)
#     if(length(aahighlyCorDescr1) > 0){
#       aa_all_data1 <- aa_all_data1[,-aahighlyCorDescr1]
#     }
#   }
# 
# 
#   aanzv2 <- nearZeroVar(my_learning_datasets_1000bp_2[[i]], saveMetrics= F)
#   if(length(aanzv2) > 0){
#     aa_all_data2 <- my_learning_datasets_1000bp_2[[i]][,-aanzv2]
#   }else{
#     aa_all_data2 <- my_learning_datasets_1000bp_2[[i]]
#   }
# 
# 
#   #  Identifying Correlated Predictors
#   if(ncol(aa_all_data2[, 1:(ncol(aa_all_data2) - 1)]) > 1){
#     aadescrCor2 <- cor(aa_all_data2[, 1:(ncol(aa_all_data2) - 1)])
#     aahighlyCorDescr2 <- findCorrelation(aadescrCor2, cutoff = .9)
#     if(length(aahighlyCorDescr2) > 0){
#       aa_all_data2 <- aa_all_data2[,-aahighlyCorDescr2]
#     }
#   }
#   aa_train_data1 <- aa_all_data1[my_data_partitions[[4]], ]
#   aa_test_data1 <-  aa_all_data1[-my_data_partitions[[4]], ]
# 
#   aa_train_data2 <- aa_all_data2[my_data_partitions[[5]], ]
#   aa_test_data2 <-  aa_all_data2[-my_data_partitions[[5]], ]
# 
# 
# 
#   aactrl <- trainControl(method = "repeatedcv",
#                          repeats = 3,
#                          classProbs = TRUE,
#                          summaryFunction = defaultSummary,
#                          savePredictions = TRUE,
#                          ## new option here:
#                          sampling = "down")
#   print(i)
#   print("aaRF_down")
#   aaRF_down1 <- train(label_set ~ ., data = aa_train_data1,
#                      method = "rf",
#                      #preProc = c("center", "scale"),
#                      #nbagg = 50,
#                      metric = "Kappa",
#                      trControl = aactrl)
#   aaRF_down2 <- train(label_set ~ ., data = aa_train_data2,
#                       method = "rf",
#                       #preProc = c("center", "scale"),
#                       #nbagg = 50,
#                       metric = "Kappa",
#                       trControl = aactrl)
#   my_learning_datasets_results_RF_kappa_down_1000bp_commonset_1[[i]] <- aaRF_down1
#   my_learning_datasets_results_RF_kappa_down_1000bp_commonset_2[[i]] <- aaRF_down2
# }
# 
# save(list = c("my_learning_datasets_results_RF_kappa_down_1000bp_commonset_1",
#               "my_learning_datasets_results_RF_kappa_down_1000bp_commonset_2"),
#      file = "classification_results_common_test.RData")
# 
# load("classification_results_common_test.RData")
# plotting
aa_gglist_2 <- list()

for(i in 1:length(my_learning_datasets_1000bp_2)){
    aanzv1 <- nearZeroVar(my_learning_datasets_1000bp[[i]], saveMetrics= F)
    if(length(aanzv1) > 0){
      aa_all_data1 <- my_learning_datasets_1000bp[[i]][,-aanzv1]
    }else{
      aa_all_data1 <- my_learning_datasets_1000bp[[i]]
    }

    #  Identifying Correlated Predictors
    if(ncol(aa_all_data1[, 1:(ncol(aa_all_data1) - 1)]) > 1){
      aadescrCor1 <- cor(aa_all_data1[, 1:(ncol(aa_all_data1) - 1)])
      aahighlyCorDescr1 <- findCorrelation(aadescrCor1, cutoff = .9)
      if(length(aahighlyCorDescr1) > 0){
        aa_all_data1 <- aa_all_data1[,-aahighlyCorDescr1]
      }
    }


    aanzv2 <- nearZeroVar(my_learning_datasets_1000bp_2[[i]], saveMetrics= F)
    if(length(aanzv2) > 0){
      aa_all_data2 <- my_learning_datasets_1000bp_2[[i]][,-aanzv2]
    }else{
      aa_all_data2 <- my_learning_datasets_1000bp_2[[i]]
    }


    #  Identifying Correlated Predictors
    if(ncol(aa_all_data2[, 1:(ncol(aa_all_data2) - 1)]) > 1){
      aadescrCor2 <- cor(aa_all_data2[, 1:(ncol(aa_all_data2) - 1)])
      aahighlyCorDescr2 <- findCorrelation(aadescrCor2, cutoff = .9)
      if(length(aahighlyCorDescr2) > 0){
        aa_all_data2 <- aa_all_data2[,-aahighlyCorDescr2]
      }
    }
    aa_train_data1 <- aa_all_data1[my_data_partitions[[4]], ]
    aa_test_data1 <-  aa_all_data1[-my_data_partitions[[4]], ]

    aa_train_data2 <- aa_all_data2[my_data_partitions[[5]], ]
    aa_test_data2 <-  aa_all_data2[-my_data_partitions[[5]], ]
  aaindex_pos1 <- aa_test_data1$label_set == "Pos"
  aaindex_neg1 <- aa_test_data1$label_set == "Neg"
  aaindex_pos2 <- aa_test_data2$label_set == "Pos"
  aaindex_neg2 <- aa_test_data2$label_set == "Neg"

  aarfdown1 <- predict(my_learning_datasets_results_RF_kappa_down_1000bp_commonset_1[[i]], 
                      aa_test_data1, type = "prob")
  aarfdown2 <- predict(my_learning_datasets_results_RF_kappa_down_1000bp_commonset_2[[i]], 
                       aa_test_data2, type = "prob")

  aapreds_list <- list(aarfdown1$Pos
                       ,aarfdown2$Pos)
  # List of actual values (same for all)
  m <- length(aapreds_list)
  stopifnot(all(aa_test_data1$label_set == aa_test_data2$label_set))
  aactuals_list <- rep(list(aa_test_data1$label_set), m)
  # Plot the ROC curves
  aapred <- prediction(aapreds_list, aactuals_list)
  aauc <- performance(aapred, measure = "auc")
  aarocs <- performance(aapred, "tpr", "fpr")
  print(aauc@y.values)
  plot(aarocs, col = as.list(1:m), main = paste0("Test Set ROC Curves for feature set: ", names(my_learning_datasets_1000bp_2)[i]))
  abline(a = 0, b = 1, col = 6, lty = 3, lwd = 0.7)
  legend(x = "bottomright", cex = 0.6,
         legend = c(paste0("RF_down -inters", format(aauc@y.values[[1]], digits = 2)),
                    paste0("RF_down -union", format(aauc@y.values[[2]], digits = 2))),
         fill = 1:m)
  aa_model_list <- list(RF_down1 = my_learning_datasets_results_RF_kappa_down_1000bp_commonset_1[[i]],
                        RF_down2 = my_learning_datasets_results_RF_kappa_down_1000bp_commonset_2[[i]])
  
  
  aadif <- setdiff(colnames(aa_test_data2),
                   colnames(aa_test_data1))
  if(length(aadif) > 0){
    aa_test_data1 <- cbind(aa_test_data1, aa_test_data2[aadif])
  }
  aamodel_list_pr <- aa_model_list %>%
    map(.f = calc_auprc, data = (aa_test_data1))
  aamodel_list_pr_auc <- aamodel_list_pr %>%
    map(function(the_mod) the_mod$auc.integral)
  aamodel_list_pr_auc <- unlist(aamodel_list_pr_auc)
  
  aaresults_list_pr <- list(NA)
  num_mod <- 1
  
  for(the_pr in aamodel_list_pr){
    
    aaresults_list_pr[[num_mod]] <- data_frame(recall = the_pr$curve[, 1],
                                               precision = the_pr$curve[, 2],
                                               model = names(aamodel_list_pr)[num_mod])
    
    num_mod <- num_mod + 1
    
  }
  
  aaresults_df_pr <- bind_rows(aaresults_list_pr)
  
  aacustom_col <- colorspace::rainbow_hcl(length(aa_model_list))
  
  aa_gglist_2[[i]] <- ggplot(aes(x = recall, y = precision, group = model), data = aaresults_df_pr) +
    geom_line(aes(color = model), size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = aacustom_col, labels = paste0(names(aa_model_list) , format(aamodel_list_pr_auc, digits = 2))) +
    geom_abline(intercept = sum(aa_test_data1$label_set == "Pos")/nrow(aa_test_data1),
                slope = 0, color = "gray", size = 1) +
    ggtitle(paste0("PR curve ", names(my_learning_datasets_1000bp_2)[i]))+
    theme_bw()
  
}
names(aa_gglist_2) <- names(my_learning_datasets_1000bp_2)
Prec_recall_curve_posneg_3 <- aa_gglist_2
Prec_recall_curve_posneg_3[[6]]
varImp(my_learning_datasets_results_RF_kappa_down_1000bp_commonset_1[[6]])
varImp(my_learning_datasets_results_RF_kappa_down_1000bp_commonset_2[[6]])


####################################################################################################
# Run logistic, svmRadial and RF with downsampling for each of six feature sets of the third positive_negative set

# universe: (union of two ER ChIP replicates)
# positive: intersection of (universe) and (intersection of two eRNA replicates) and (H3k27ac)
# negative: universe minus the union of two eRNA replicates and minus H3k27 (picked 2500 at random)

feature_set_sum_LLR_1000bp_3
# feature sets:
#sum LR
feature_set_sum_LLR_1000bp_3 <- rbind(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[3]], 
                                      Positive_set_seq_list_chopped_scoreMat_list_1000bp[[3]])
colnames(feature_set_sum_LLR_1000bp_3)[13] <- "NKX3_1"
# adjacency
feature_set_Adjmat_1000bp_3 <- rbind(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[3]], 
                                     Positive_set_seq_list_chopped_adjanScore_list_1000bp[[3]])
label_set <- c(rep("Neg", nrow(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[3]])), 
               rep("Pos", nrow(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[3]])))

# changing the name of NKX3-1 to NKX3_1 so it works in the classification models
aa <- strsplit(colnames(feature_set_Adjmat_1000bp_3), split = "_vs_")
for(i in 1:length(aa)){
  aa[[i]][aa[[i]] == "NKX3-1"] <- "NKX3_1"
}
for(i in 1:length(aa)){
  aa[[i]] <- paste(aa[[i]][1], aa[[i]][2], sep = "_vs_")
}
colnames(feature_set_Adjmat_1000bp_3) <- unlist(aa)
# overlap
feature_set_Overlapmat_1000bp_3 <- rbind(Negative_set_seq_list_chopped_OverlapScore_list_1000bp[[3]], 
                                         Positive_set_seq_list_chopped_OverlapScore_list_1000bp[[3]])
colnames(feature_set_Overlapmat_1000bp_3) <- unlist(aa)


my_learning_datasets_1000bp_3 <- list()
# sum LR only
my_learning_datasets_1000bp_3[[1]] <- as.data.frame(cbind(feature_set_sum_LLR_1000bp_3,
                                                          label_set), 
                                                    stringsAsFactors =F)
my_learning_datasets_1000bp_3[[1]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_3[[1]])
# adjacancy only
my_learning_datasets_1000bp_3[[2]] <- as.data.frame(cbind(feature_set_Adjmat_1000bp_3,
                                                          label_set), 
                                                    stringsAsFactors =F)
my_learning_datasets_1000bp_3[[2]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_3[[2]])

# Overlap only
my_learning_datasets_1000bp_3[[3]] <- as.data.frame(cbind(feature_set_Overlapmat_1000bp_3,
                                                          label_set), 
                                                    stringsAsFactors =F)
my_learning_datasets_1000bp_3[[3]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_3[[3]])
# LLR and adj
my_learning_datasets_1000bp_3[[4]] <- as.data.frame(cbind(cbind(feature_set_sum_LLR_1000bp_3,
                                                                feature_set_Adjmat_1000bp_3),
                                                          label_set), 
                                                    stringsAsFactors =F)
my_learning_datasets_1000bp_3[[4]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_3[[4]])
# LLR and overlap
my_learning_datasets_1000bp_3[[5]] <- as.data.frame(cbind(cbind(feature_set_sum_LLR_1000bp_3,
                                                                feature_set_Overlapmat_1000bp_3),
                                                          label_set), 
                                                    stringsAsFactors =F)
my_learning_datasets_1000bp_3[[5]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_3[[5]])
# LLR and adj and overlap
my_learning_datasets_1000bp_3[[6]] <- as.data.frame(cbind(cbind(feature_set_sum_LLR_1000bp_3,
                                                                cbind(feature_set_Adjmat_1000bp_3, 
                                                                      feature_set_Overlapmat_1000bp_3)),
                                                          label_set), 
                                                    stringsAsFactors =F)
colnames(my_learning_datasets_1000bp_3[[6]])[28:405] <- paste(colnames(my_learning_datasets_1000bp_3[[6]])[28:405], "adj", sep = "_at_")
colnames(my_learning_datasets_1000bp_3[[6]])[406:783] <- paste(colnames(my_learning_datasets_1000bp_3[[6]])[406:783], "ovl", sep = "_at_")
my_learning_datasets_1000bp_3[[6]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_3[[6]])

names(my_learning_datasets_1000bp_3) <- c("Sum_LR","Adjacency", "Overlap", 
                                          "Sum_LR_plus_Adjacency", "Sum_LR_plus_Overlap", 
                                          "Sum_LR_plus_Adjacency_plus_Overlap")

my_data_partitions[[6]] <- createDataPartition(y = my_learning_datasets_1000bp_3[[1]][, ncol(my_learning_datasets_1000bp_3[[1]])], times = 1, p = 0.75,list = F)

# my_learning_datasets_results_logistic_kappa_down_1000bp_2 <- list()
# my_learning_datasets_results_logistic_kappa_up_1000bp_2 <- list()
# my_learning_datasets_results_svmLinear_kappa_down_1000bp_2 <- list()
# my_learning_datasets_results_svmLinear_kappa_up_1000bp_2 <- list()
my_learning_datasets_results_RF_kappa_down_1000bp_3 <- list()
# my_learning_datasets_results_RF_kappa_up_1000bp_2 <- list()

#
for(i in 1:length(my_learning_datasets_1000bp_3)){
  aanzv <- nearZeroVar(my_learning_datasets_1000bp_3[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets_1000bp_3[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets_1000bp_3[[i]]
  }


  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  print("#########################################################################################")
  print(paste0("number of features in original feature set number ", i," named ", names(my_learning_datasets_1000bp_2)[i] ,":"))
  print(ncol(my_learning_datasets_1000bp_3[[i]]))
  print("number of features after removing near zero variance and highly correlated features ( > 0.9)")
  print(ncol(aa_all_data))
  print("#########################################################################################")
  aa_train_data <- aa_all_data[my_data_partitions[[6]], ]
  aa_test_data <-  aa_all_data[-my_data_partitions[[6]], ]


  aactrl <- trainControl(method = "repeatedcv",
                         repeats = 3,
                         classProbs = TRUE,
                         summaryFunction = defaultSummary,
                         savePredictions = TRUE,
                         ## new option here:
                         sampling = "down")
  # print("aalogis_down")
  # aalogis_down <- train(label_set ~ ., data = aa_train_data,
  #                       method = "regLogistic",
  #                       #preProc = c("center", "scale"),
  #                       #nbagg = 50,
  #                       metric = "Kappa",
  #                       trControl = aactrl)
  # print("aaSVM_down")
  # aaSVM_down <- train(label_set ~ ., data = aa_train_data,
  #                     method = "svmLinear3",
  #                     preProc = c("center", "scale"),
  #                     #nbagg = 50,
  #                     metric = "Kappa",
  #                     trControl = aactrl)
  print("aaRF_down")
  aaRF_down <- train(label_set ~ ., data = aa_train_data,
                     method = "rf",
                     #preProc = c("center", "scale"),
                     #nbagg = 50,
                     metric = "Kappa",
                     trControl = aactrl)
  # aactrl$sampling <- "up"
  # print("aalogis_up")
  # aalogis_up <- train(label_set ~ ., data = aa_train_data,
  #                     method = "regLogistic",
  #                     #preProc = c("center", "scale"),
  #                     #nbagg = 50,
  #                     metric = "Kappa",
  #                     trControl = aactrl)
  # print("aaSVM_up")
  # aaSVM_up <- train(label_set ~ ., data = aa_train_data,
  #                   method = "svmLinear3",
  #                   preProc = c("center", "scale"),
  #                   #nbagg = 50,
  #                   metric = "Kappa",
  #                   trControl = aactrl)
  # print("aaRF_up")
  # aaRF_up <- train(label_set ~ ., data = aa_train_data,
  #                  method = "rf",
  #                  #preProc = c("center", "scale"),
  #                  #nbagg = 50,
  #                  metric = "Kappa",
  #                  trControl = aactrl)
  # my_learning_datasets_results_logistic_kappa_down_1000bp_2[[i]] <- aalogis_down
  # my_learning_datasets_results_svmLinear_kappa_down_1000bp_2[[i]] <- aaSVM_down
  # my_learning_datasets_results_logistic_kappa_up_1000bp_2[[i]] <- aalogis_up
  # my_learning_datasets_results_svmLinear_kappa_up_1000bp_2[[i]] <- aaSVM_up
  my_learning_datasets_results_RF_kappa_down_1000bp_3[[i]] <- aaRF_down
  # my_learning_datasets_results_RF_kappa_up_1000bp_2[[i]] <- aaRF_up

}

save(list = c("my_learning_datasets_results_RF_kappa_down_1000bp_3"),
     file = "classification_results_set_3.RData")
load("classification_results_set_3.RData")

aa_gglist_2 <- list()

for(i in 1:length(my_learning_datasets_1000bp_3)){
  aanzv1 <- nearZeroVar(my_learning_datasets_1000bp_3[[i]], saveMetrics= F)
  if(length(aanzv1) > 0){
    aa_all_data1 <- my_learning_datasets_1000bp_3[[i]][,-aanzv1]
  }else{
    aa_all_data1 <- my_learning_datasets_1000bp_3[[i]]
  }
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data1[, 1:(ncol(aa_all_data1) - 1)]) > 1){
    aadescrCor1 <- cor(aa_all_data1[, 1:(ncol(aa_all_data1) - 1)])
    aahighlyCorDescr1 <- findCorrelation(aadescrCor1, cutoff = .9)
    if(length(aahighlyCorDescr1) > 0){
      aa_all_data1 <- aa_all_data1[,-aahighlyCorDescr1]
    }
  }
  
  aa_train_data1 <- aa_all_data1[my_data_partitions[[6]], ]
  aa_test_data1 <-  aa_all_data1[-my_data_partitions[[6]], ]
  
  aaindex_pos1 <- aa_test_data1$label_set == "Pos"
  aaindex_neg1 <- aa_test_data1$label_set == "Neg"

  aarfdown1 <- predict(my_learning_datasets_results_RF_kappa_down_1000bp_3[[i]], 
                       aa_test_data1, type = "prob")
  
  aapreds_list <- list(aarfdown1$Pos)
  # List of actual values (same for all)
  m <- length(aapreds_list)
  aactuals_list <- rep(list(aa_test_data1$label_set), m)
  # Plot the ROC curves
  aapred <- prediction(aapreds_list, aactuals_list)
  aauc <- performance(aapred, measure = "auc")
  aarocs <- performance(aapred, "tpr", "fpr")
  print(aauc@y.values)
  plot(aarocs, col = as.list(1:m), main = paste0("Test Set ROC Curves for feature set: ", names(my_learning_datasets_1000bp_2)[i]))
  abline(a = 0, b = 1, col = 6, lty = 3, lwd = 0.7)
  legend(x = "bottomright", cex = 0.6,
         legend = c(paste0("RF_down -h3k27sets ", "AUC: " , format(aauc@y.values[[1]], digits = 2))),
         fill = 1:m)
  aa_model_list <- list(RF_down1 = my_learning_datasets_results_RF_kappa_down_1000bp_3[[i]])
  
  
  aamodel_list_pr <- aa_model_list %>%
    map(.f = calc_auprc, data = (aa_test_data1))
  aamodel_list_pr_auc <- aamodel_list_pr %>%
    map(function(the_mod) the_mod$auc.integral)
  aamodel_list_pr_auc <- unlist(aamodel_list_pr_auc)
  
  aaresults_list_pr <- list(NA)
  num_mod <- 1
  
  for(the_pr in aamodel_list_pr){
    
    aaresults_list_pr[[num_mod]] <- data_frame(recall = the_pr$curve[, 1],
                                               precision = the_pr$curve[, 2],
                                               model = names(aamodel_list_pr)[num_mod])
    
    num_mod <- num_mod + 1
    
  }
  
  aaresults_df_pr <- bind_rows(aaresults_list_pr)
  
  aacustom_col <- colorspace::rainbow_hcl(length(aa_model_list))
  
  aa_gglist_2[[i]] <- ggplot(aes(x = recall, y = precision, group = model), data = aaresults_df_pr) +
    geom_line(aes(color = model), size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = aacustom_col, labels = paste(names(aa_model_list) ,"AUC: " ,format(aamodel_list_pr_auc, digits = 2))) +
    geom_abline(intercept = sum(aa_test_data1$label_set == "Pos")/nrow(aa_test_data1),
                slope = 0, color = "gray", size = 1) +
    ggtitle(paste0("PR curve ", names(my_learning_datasets_1000bp_3)[i]))+
    theme_bw()
  
}
names(aa_gglist_2) <-  names(my_learning_datasets_1000bp_3)
aa_gglist_2[[6]]

varImp(my_learning_datasets_results_RF_kappa_down_1000bp_3[[6]])

seqLogo::seqLogo(TF.motifs.Expanded_new_pseudo_t[[25]])

####################################################################################################
# Run logistic, svmRadial and RF with downsampling for each of six feature sets of the fourth positive_negative set

# universe: (intesection of H3k27ac and union of two ER ChIP replicates)
# positive: intersection of (universe) and ((intersection of two eRNA replicates) minus (intersection of eRNA after EtOH)) 
# negative: Universe minus union of eRNA peaks after E2

feature_set_sum_LLR_1000bp_4
# feature sets:
#sum LR
feature_set_sum_LLR_1000bp_4 <- rbind(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[4]], 
                                      Positive_set_seq_list_chopped_scoreMat_list_1000bp[[4]])
colnames(feature_set_sum_LLR_1000bp_4)[13] <- "NKX3_1"
# adjacency
feature_set_Adjmat_1000bp_4 <- rbind(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[4]], 
                                     Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]])
label_set <- c(rep("Neg", nrow(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[4]])), 
               rep("Pos", nrow(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[4]])))

# changing the name of NKX3-1 to NKX3_1 so it works in the classification models
aa <- strsplit(colnames(feature_set_Adjmat_1000bp_4), split = "_vs_")
for(i in 1:length(aa)){
  aa[[i]][aa[[i]] == "NKX3-1"] <- "NKX3_1"
}
for(i in 1:length(aa)){
  aa[[i]] <- paste(aa[[i]][1], aa[[i]][2], sep = "_vs_")
}
colnames(feature_set_Adjmat_1000bp_4) <- unlist(aa)
# overlap
feature_set_Overlapmat_1000bp_4 <- rbind(Negative_set_seq_list_chopped_OverlapScore_list_1000bp[[4]], 
                                         Positive_set_seq_list_chopped_OverlapScore_list_1000bp[[4]])
colnames(feature_set_Overlapmat_1000bp_4) <- unlist(aa)


my_learning_datasets_1000bp_4 <- list()
# sum LR only
my_learning_datasets_1000bp_4[[1]] <- as.data.frame(cbind(feature_set_sum_LLR_1000bp_4,
                                                          label_set), 
                                                    stringsAsFactors =F)
my_learning_datasets_1000bp_4[[1]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_4[[1]])
# adjacancy only
my_learning_datasets_1000bp_4[[2]] <- as.data.frame(cbind(feature_set_Adjmat_1000bp_4,
                                                          label_set), 
                                                    stringsAsFactors =F)
my_learning_datasets_1000bp_4[[2]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_4[[2]])

# Overlap only
my_learning_datasets_1000bp_4[[3]] <- as.data.frame(cbind(feature_set_Overlapmat_1000bp_4,
                                                          label_set), 
                                                    stringsAsFactors =F)
my_learning_datasets_1000bp_4[[3]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_4[[3]])
# LLR and adj
my_learning_datasets_1000bp_4[[4]] <- as.data.frame(cbind(cbind(feature_set_sum_LLR_1000bp_4,
                                                                feature_set_Adjmat_1000bp_4),
                                                          label_set), 
                                                    stringsAsFactors =F)
my_learning_datasets_1000bp_4[[4]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_4[[4]])
# LLR and overlap
my_learning_datasets_1000bp_4[[5]] <- as.data.frame(cbind(cbind(feature_set_sum_LLR_1000bp_4,
                                                                feature_set_Overlapmat_1000bp_4),
                                                          label_set), 
                                                    stringsAsFactors =F)
my_learning_datasets_1000bp_4[[5]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_4[[5]])
# LLR and adj and overlap
my_learning_datasets_1000bp_4[[6]] <- as.data.frame(cbind(cbind(feature_set_sum_LLR_1000bp_4,
                                                                cbind(feature_set_Adjmat_1000bp_4, 
                                                                      feature_set_Overlapmat_1000bp_4)),
                                                          label_set), 
                                                    stringsAsFactors =F)
colnames(my_learning_datasets_1000bp_4[[6]])[28:405] <- paste(colnames(my_learning_datasets_1000bp_4[[6]])[28:405], "adj", sep = "_at_")
colnames(my_learning_datasets_1000bp_4[[6]])[406:783] <- paste(colnames(my_learning_datasets_1000bp_4[[6]])[406:783], "ovl", sep = "_at_")
my_learning_datasets_1000bp_4[[6]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_4[[6]])

names(my_learning_datasets_1000bp_4) <- c("Sum_LR","Adjacency", "Overlap", 
                                          "Sum_LR_plus_Adjacency", "Sum_LR_plus_Overlap", 
                                          "Sum_LR_plus_Adjacency_plus_Overlap")

my_data_partitions[[7]] <- createDataPartition(y = my_learning_datasets_1000bp_4[[1]][, ncol(my_learning_datasets_1000bp_4[[1]])], times = 1, p = 0.75,list = F)
############# add 5mer features
# a kmer random forest on 
library(tcR)
#All_5mers <- generate.kmers(.k = 5)
#Sequence_posnegSet_5mers <- list()
aaboth4 <- c(Negative_set_seq_list_char_1000bp[[4]], Positive_set_seq_list_char_1000bp[[4]])

Sequence_posnegSet_5mers[[4]] <- list()
for(i in 1:length(aaboth4)){
  Sequence_posnegSet_5mers[[4]][[i]] <- get.kmers(aaboth4[i], .k = 5)
}
names(Sequence_posnegSet_5mers[[4]]) <- names(aaboth4)

for(i in 4:length(Sequence_posnegSet_5mers)){
  Sequence_posnegSet_5mers_full[[i]] <- list()
  for(j in 1:length(Sequence_posnegSet_5mers[[i]])){
    Sequence_posnegSet_5mers_full[[i]][[j]] <- cbind(All_5mers, integer(length(All_5mers)))
    Sequence_posnegSet_5mers_full[[i]][[j]][, 2] <- as.integer(Sequence_posnegSet_5mers[[i]][[j]]$Count[match(Sequence_posnegSet_5mers_full[[i]][[j]][, 1], 
                                                                                                              Sequence_posnegSet_5mers[[i]][[j]]$Kmers)])
    Sequence_posnegSet_5mers_full[[i]][[j]][is.na(Sequence_posnegSet_5mers_full[[i]][[j]])] <- 0
  }
  names(Sequence_posnegSet_5mers_full[[i]]) <- names(Sequence_posnegSet_5mers[[i]])
}
# create a list of two entries for the two datasets. each: a matrix where each row corresponds to one enhancer, each column corresponds to one kmer. 
# each entry indicates the count of the kmer in that sequence
#Sequence_posnegSet_5mers_full_mat_list <- list()
for(i in 4:length(Sequence_posnegSet_5mers_full)){
  Sequence_posnegSet_5mers_full_mat_list[[i]] <- t(do.call(cbind, Sequence_posnegSet_5mers_full[[i]]))
  Sequence_posnegSet_5mers_full_mat_list[[i]] <- Sequence_posnegSet_5mers_full_mat_list[[i]][seq(2, nrow(Sequence_posnegSet_5mers_full_mat_list[[i]]), 2), ]
  Sequence_posnegSet_5mers_full_mat_list[[i]] <- apply(Sequence_posnegSet_5mers_full_mat_list[[i]], 2, as.integer)
  rownames(Sequence_posnegSet_5mers_full_mat_list[[i]]) <- names(Sequence_posnegSet_5mers_full[[i]])
  colnames(Sequence_posnegSet_5mers_full_mat_list[[i]]) <- Sequence_posnegSet_5mers_full[[i]][[1]][, 1]
}
names(Sequence_posnegSet_5mers_full_mat_list)[4] <- c("Set4")

hist(colSums(Sequence_posnegSet_5mers_full_mat_list[[4]]), breaks = 100)
sum(colSums(Sequence_posnegSet_5mers_full_mat_list[[4]]) == 0)
#

my_learning_datasets_1000bp_4$kmer5 <- as.data.frame(cbind(Sequence_posnegSet_5mers_full_mat_list[[4]],
                                                         c(rep("Neg", length(Negative_set_seq_list_char_1000bp[[4]])),
                                                           rep("Pos", length(Positive_set_seq_list_char_1000bp[[4]])))), 
                                                   stringsAsFactors =F)

my_learning_datasets_1000bp_4$kmer5 <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_4$kmer5)


##############
save(list = c("my_learning_datasets_1000bp_4",
              "my_data_partitions"), file = "Classification_input_set4.RData")
# my_learning_datasets_results_logistic_kappa_down_1000bp_2 <- list()
# my_learning_datasets_results_logistic_kappa_up_1000bp_2 <- list()
# my_learning_datasets_results_svmLinear_kappa_down_1000bp_2 <- list()
# my_learning_datasets_results_svmLinear_kappa_up_1000bp_2 <- list()
my_learning_datasets_results_RF_kappa_down_1000bp_4 <- list()
# my_learning_datasets_results_RF_kappa_up_1000bp_2 <- list()

#
length(my_learning_datasets_1000bp_4)
for(i in 1:1){
  aanzv <- nearZeroVar(my_learning_datasets_1000bp_4[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets_1000bp_4[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets_1000bp_4[[i]]
  }
  
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  print("#########################################################################################")
  print(paste0("number of features in original feature set number ", i," named ", 
               names(my_learning_datasets_1000bp_4)[i] ,":"))
  print(ncol(my_learning_datasets_1000bp_4[[i]]))
  print("number of features after removing near zero variance and highly correlated features ( > 0.9)")
  print(ncol(aa_all_data))
  print("#########################################################################################")
  aa_train_data <- aa_all_data[my_data_partitions[[7]], ]
  aa_test_data <-  aa_all_data[-my_data_partitions[[7]], ]
  
  
  aactrl <- trainControl(method = "repeatedcv",
                         number = 5,
                         repeats = 3,
                         classProbs = TRUE,
                         summaryFunction = mnLogLoss,
                         savePredictions = TRUE,
                         ## new option here:
                         sampling = "down", 
                         search = 'random')
  # aatunegrid <- expand.grid(.mtry = (1:15))
  # print("aalogis_down")
  # aalogis_down <- train(label_set ~ ., data = aa_train_data,
  #                       method = "regLogistic",
  #                       #preProc = c("center", "scale"),
  #                       #nbagg = 50,
  #                       metric = "Kappa",
  #                       trControl = aactrl)
  # print("aaSVM_down")
  # aaSVM_down <- train(label_set ~ ., data = aa_train_data,
  #                     method = "svmLinear3",
  #                     preProc = c("center", "scale"),
  #                     #nbagg = 50,
  #                     metric = "Kappa",
  #                     trControl = aactrl)
  print("aaRF_down")
  aaRF_down <- train(label_set ~ ., data = aa_train_data,
                     method = "svmRadial",
                     preProc = c("center", "scale"),
                     #nbagg = 50,
                     metric = "logLoss",
                     trControl = aactrl,
                     tuneLength  = 10)
  # aactrl$sampling <- "up"
  # print("aalogis_up")
  # aalogis_up <- train(label_set ~ ., data = aa_train_data,
  #                     method = "regLogistic",
  #                     #preProc = c("center", "scale"),
  #                     #nbagg = 50,
  #                     metric = "Kappa",
  #                     trControl = aactrl)
  # print("aaSVM_up")
  # aaSVM_up <- train(label_set ~ ., data = aa_train_data,
  #                   method = "svmLinear3",
  #                   preProc = c("center", "scale"),
  #                   #nbagg = 50,
  #                   metric = "Kappa",
  #                   trControl = aactrl)
  # print("aaRF_up")
  # aaRF_up <- train(label_set ~ ., data = aa_train_data,
  #                  method = "rf",
  #                  #preProc = c("center", "scale"),
  #                  #nbagg = 50,
  #                  metric = "Kappa",
  #                  trControl = aactrl)
  # my_learning_datasets_results_logistic_kappa_down_1000bp_2[[i]] <- aalogis_down
  # my_learning_datasets_results_svmLinear_kappa_down_1000bp_2[[i]] <- aaSVM_down
  # my_learning_datasets_results_logistic_kappa_up_1000bp_2[[i]] <- aalogis_up
  # my_learning_datasets_results_svmLinear_kappa_up_1000bp_2[[i]] <- aaSVM_up
  my_learning_datasets_results_RF_kappa_down_1000bp_4[[i]] <- aaRF_down
  # my_learning_datasets_results_RF_kappa_up_1000bp_2[[i]] <- aaRF_up
  
}

save(list = c("my_learning_datasets_results_RF_kappa_down_1000bp_4"),
     file = "classification_results_set_4.RData")
load("classification_results_set_5.RData")

aa_gglist_2 <- list()

for(i in 1:length(my_learning_datasets_1000bp_4)){
  aanzv1 <- nearZeroVar(my_learning_datasets_1000bp_4[[i]], saveMetrics= F)
  if(length(aanzv1) > 0){
    aa_all_data1 <- my_learning_datasets_1000bp_4[[i]][,-aanzv1]
  }else{
    aa_all_data1 <- my_learning_datasets_1000bp_4[[i]]
  }
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data1[, 1:(ncol(aa_all_data1) - 1)]) > 1){
    aadescrCor1 <- cor(aa_all_data1[, 1:(ncol(aa_all_data1) - 1)])
    aahighlyCorDescr1 <- findCorrelation(aadescrCor1, cutoff = .9)
    if(length(aahighlyCorDescr1) > 0){
      aa_all_data1 <- aa_all_data1[,-aahighlyCorDescr1]
    }
  }
  
  aa_train_data1 <- aa_all_data1[my_data_partitions[[7]], ]
  aa_test_data1 <-  aa_all_data1[-my_data_partitions[[7]], ]
  
  aaindex_pos1 <- aa_test_data1$label_set == "Pos"
  aaindex_neg1 <- aa_test_data1$label_set == "Neg"
  
  aarfdown1 <- predict(my_learning_datasets_results_RF_kappa_down_1000bp_4[[i]], 
                       aa_test_data1, type = "prob")
  
  aapreds_list <- list(aarfdown1$Pos)
  # List of actual values (same for all)
  m <- length(aapreds_list)
  aactuals_list <- rep(list(aa_test_data1$label_set), m)
  # Plot the ROC curves
  aapred <- prediction(aapreds_list, aactuals_list)
  aauc <- performance(aapred, measure = "auc")
  aarocs <- performance(aapred, "tpr", "fpr")
  print(aauc@y.values)
  plot(aarocs, col = as.list(1:m), main = paste0("Test Set ROC Curves for feature set: ",
                                                 names(my_learning_datasets_1000bp_4)[i]))
  abline(a = 0, b = 1, col = 6, lty = 3, lwd = 0.7)
  legend(x = "bottomright", cex = 0.6,
         legend = c(paste0("RF_down -Set4 ", "AUC: " , format(aauc@y.values[[1]], digits = 2))),
         fill = 1:m)
  aa_model_list <- list(RF_down1 = my_learning_datasets_results_RF_kappa_down_1000bp_4[[i]])
  
  
  aamodel_list_pr <- aa_model_list %>%
    map(.f = calc_auprc, data = (aa_test_data1))
  aamodel_list_pr_auc <- aamodel_list_pr %>%
    map(function(the_mod) the_mod$auc.integral)
  aamodel_list_pr_auc <- unlist(aamodel_list_pr_auc)
  
  aaresults_list_pr <- list(NA)
  num_mod <- 1
  
  for(the_pr in aamodel_list_pr){
    
    aaresults_list_pr[[num_mod]] <- data_frame(recall = the_pr$curve[, 1],
                                               precision = the_pr$curve[, 2],
                                               model = names(aamodel_list_pr)[num_mod])
    
    num_mod <- num_mod + 1
    
  }
  
  aaresults_df_pr <- bind_rows(aaresults_list_pr)
  
  aacustom_col <- colorspace::rainbow_hcl(length(aa_model_list))
  
  aa_gglist_2[[i]] <- ggplot(aes(x = recall, y = precision, group = model), data = aaresults_df_pr) +
    geom_line(aes(color = model), size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = aacustom_col, labels = paste(names(aa_model_list) ,"AUC: " ,
                                                             format(aamodel_list_pr_auc, digits = 2))) +
    geom_abline(intercept = sum(aa_test_data1$label_set == "Pos")/nrow(aa_test_data1),
                slope = 0, color = "gray", size = 1) +
    ggtitle(paste0("PR curve ", names(my_learning_datasets_1000bp_4)[i]))+
    theme_bw()
  
}
names(aa_gglist_2) <-  names(my_learning_datasets_1000bp_4)
aa_gglist_2[[7]]

varImp(my_learning_datasets_results_RF_kappa_down_1000bp_5[[4]])
################################################
####### get performance on training
################################################
aa_gglist_3 <- list()

for(i in 1:length(my_learning_datasets_1000bp_4)){
  aanzv1 <- nearZeroVar(my_learning_datasets_1000bp_4[[i]], saveMetrics= F)
  if(length(aanzv1) > 0){
    aa_all_data1 <- my_learning_datasets_1000bp_4[[i]][,-aanzv1]
  }else{
    aa_all_data1 <- my_learning_datasets_1000bp_4[[i]]
  }
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data1[, 1:(ncol(aa_all_data1) - 1)]) > 1){
    aadescrCor1 <- cor(aa_all_data1[, 1:(ncol(aa_all_data1) - 1)])
    aahighlyCorDescr1 <- findCorrelation(aadescrCor1, cutoff = .9)
    if(length(aahighlyCorDescr1) > 0){
      aa_all_data1 <- aa_all_data1[,-aahighlyCorDescr1]
    }
  }
  
  aa_train_data1 <- aa_all_data1[my_data_partitions[[7]], ]
  aa_test_data1 <-  aa_all_data1[-my_data_partitions[[7]], ]
  
  aaindex_pos1 <- aa_test_data1$label_set == "Pos"
  aaindex_neg1 <- aa_test_data1$label_set == "Neg"
  
  aarfdown1 <- predict(my_learning_datasets_results_RF_kappa_down_1000bp_4[[i]], 
                       aa_train_data1, type = "prob")
  
  aapreds_list <- list(aarfdown1$Pos)
  # List of actual values (same for all)
  m <- length(aapreds_list)
  aactuals_list <- rep(list(aa_train_data1$label_set), m)
  # Plot the ROC curves
  aapred <- prediction(aapreds_list, aactuals_list)
  aauc <- performance(aapred, measure = "auc")
  aarocs <- performance(aapred, "tpr", "fpr")
  print(aauc@y.values)
  plot(aarocs, col = as.list(1:m), main = paste0("Train Set ROC Curves for feature set: ", names(my_learning_datasets_1000bp_4)[i]))
  abline(a = 0, b = 1, col = 6, lty = 3, lwd = 0.7)
  legend(x = "bottomright", cex = 0.6,
         legend = c(paste0("RF_down-Set4 ", "AUC: " , format(aauc@y.values[[1]], digits = 2))),
         fill = 1:m)
  aa_model_list <- list(RF_down1 = my_learning_datasets_results_RF_kappa_down_1000bp_4[[i]])
  
  
  aamodel_list_pr <- aa_model_list %>%
    map(.f = calc_auprc, data = (aa_train_data1))
  aamodel_list_pr_auc <- aamodel_list_pr %>%
    map(function(the_mod) the_mod$auc.integral)
  aamodel_list_pr_auc <- unlist(aamodel_list_pr_auc)
  
  aaresults_list_pr <- list(NA)
  num_mod <- 1
  
  for(the_pr in aamodel_list_pr){
    
    aaresults_list_pr[[num_mod]] <- data_frame(recall = the_pr$curve[, 1],
                                               precision = the_pr$curve[, 2],
                                               model = names(aamodel_list_pr)[num_mod])
    
    num_mod <- num_mod + 1
    
  }
  
  aaresults_df_pr <- bind_rows(aaresults_list_pr)
  
  aacustom_col <- colorspace::rainbow_hcl(length(aa_model_list))
  
  aa_gglist_3[[i]] <- ggplot(aes(x = recall, y = precision, group = model), data = aaresults_df_pr) +
    geom_line(aes(color = model), size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = aacustom_col, labels = paste(names(aa_model_list) ,"AUC: " 
                                                             ,format(aamodel_list_pr_auc, digits = 2))) +
    geom_abline(intercept = sum(aa_train_data1$label_set == "Pos")/nrow(aa_train_data1),
                slope = 0, color = "gray", size = 1) +
    ggtitle(paste0("Training PR curve ", names(my_learning_datasets_1000bp_4)[i]))+
    theme_bw()
  
}
names(aa_gglist_3) <-  names(my_learning_datasets_1000bp_4)

aa_gglist_2[[7]]
aa_gglist_3[[7]]
varImp(my_learning_datasets_results_RF_kappa_down_1000bp_5[[7]])

print(my_learning_datasets_results_RF_kappa_down_1000bp_5[[1]])


#########################################################################################
# a kmer random forest on 
library(tcR)
All_5mers <- generate.kmers(.k = 5)
Sequence_posnegSet_5mers <- list()
aaboth1 <- c(Negative_set_seq_list_char_1000bp[[1]], Positive_set_seq_list_char_1000bp[[1]])
aaboth2 <- c(Negative_set_seq_list_char_1000bp[[2]], Positive_set_seq_list_char_1000bp[[2]])

Sequence_posnegSet_5mers[[1]] <- list()
for(i in 1:length(aaboth1)){
  Sequence_posnegSet_5mers[[1]][[i]] <- get.kmers(aaboth1[i], .k = 5)
}
names(Sequence_posnegSet_5mers[[1]]) <- names(aaboth1)

Sequence_posnegSet_5mers[[2]] <- list()
for(i in 1:length(aaboth2)){
  Sequence_posnegSet_5mers[[2]][[i]] <- get.kmers(aaboth2[i], .k = 5)
}
names(Sequence_posnegSet_5mers[[2]]) <- names(aaboth2)

Sequence_posnegSet_5mers_full <- list()
for(i in 1:length(Sequence_posnegSet_5mers)){
  Sequence_posnegSet_5mers_full[[i]] <- list()
  for(j in 1:length(Sequence_posnegSet_5mers[[i]])){
    Sequence_posnegSet_5mers_full[[i]][[j]] <- cbind(All_5mers, integer(length(All_5mers)))
    Sequence_posnegSet_5mers_full[[i]][[j]][, 2] <- as.integer(Sequence_posnegSet_5mers[[i]][[j]]$Count[match(Sequence_posnegSet_5mers_full[[i]][[j]][, 1], 
                                                                                                              Sequence_posnegSet_5mers[[i]][[j]]$Kmers)])
    Sequence_posnegSet_5mers_full[[i]][[j]][is.na(Sequence_posnegSet_5mers_full[[i]][[j]])] <- 0
  }
  names(Sequence_posnegSet_5mers_full[[i]]) <- names(Sequence_posnegSet_5mers[[i]])
}

# create a list of two entries for the two datasets. each: a matrix where each row corresponds to one enhancer, each column corresponds to one kmer. 
# each entry indicates the count of the kmer in that sequence
Sequence_posnegSet_5mers_full_mat_list <- list()
for(i in 1:length(Sequence_posnegSet_5mers_full)){
  Sequence_posnegSet_5mers_full_mat_list[[i]] <- t(do.call(cbind, Sequence_posnegSet_5mers_full[[i]]))
  Sequence_posnegSet_5mers_full_mat_list[[i]] <- Sequence_posnegSet_5mers_full_mat_list[[i]][seq(2, nrow(Sequence_posnegSet_5mers_full_mat_list[[i]]), 2), ]
  Sequence_posnegSet_5mers_full_mat_list[[i]] <- apply(Sequence_posnegSet_5mers_full_mat_list[[i]], 2, as.integer)
  rownames(Sequence_posnegSet_5mers_full_mat_list[[i]]) <- names(Sequence_posnegSet_5mers_full[[i]])
  colnames(Sequence_posnegSet_5mers_full_mat_list[[i]]) <- Sequence_posnegSet_5mers_full[[i]][[1]][, 1]
}
names(Sequence_posnegSet_5mers_full_mat_list) <- c("Set1", "Set2")

hist(colSums(Sequence_posnegSet_5mers_full_mat_list[[1]]), breaks = 100)

sum(colSums(Sequence_posnegSet_5mers_full_mat_list[[2]]) == 0)
#

my_learning_datasets_1000bp$kmer5 <- as.data.frame(cbind(Sequence_posnegSet_5mers_full_mat_list[[1]],
                                       c(rep("Neg", length(Negative_set_seq_list_char_1000bp[[1]])),
                                         rep("Pos", length(Positive_set_seq_list_char_1000bp[[1]])))), 
                                       stringsAsFactors =F)
my_learning_datasets_1000bp_2$kmer5 <- as.data.frame(cbind(Sequence_posnegSet_5mers_full_mat_list[[2]],
                                           c(rep("Neg", length(Negative_set_seq_list_char_1000bp[[2]])),
                                             rep("Pos", length(Positive_set_seq_list_char_1000bp[[2]])))),
                                           stringsAsFactors =F)

my_learning_datasets_1000bp$kmer5 <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp$kmer5)
my_learning_datasets_1000bp_2$kmer5 <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_2$kmer5)

aa_all_seq_2 <- c(Negative_set_seq_list_1000bp[[2]],
                  Positive_set_seq_list_1000bp[[2]])
aa_all_seq_1 <- c(Negative_set_seq_list_1000bp[[1]],
                  Positive_set_seq_list_1000bp[[1]])

aamatch <- match(aa_all_seq_2[-my_data_partitions[[5]][, 1]], aa_all_seq_1)

# for first set based on the second set
#my_data_partitions[[4]] <- matrix(setdiff(c(1:nrow(my_learning_datasets_1000bp[[1]])), aamatch), ncol = 1)

# check if the test set is the same.
all(aa_all_seq_2[-my_data_partitions[[5]][, 1]] %in% aa_all_seq_1[-my_data_partitions[[4]][, 1]])


# my_learning_datasets_results_RF_kappa_down_1000bp_commonset_1 <- list()
# my_learning_datasets_results_RF_kappa_down_1000bp_commonset_2 <- list()

for(i in 7:length(my_learning_datasets_1000bp_2)){
  
  aanzv1 <- nearZeroVar(my_learning_datasets_1000bp[[i]], saveMetrics= F)
  if(length(aanzv1) > 0){
    aa_all_data1 <- my_learning_datasets_1000bp[[i]][,-aanzv1]
  }else{
    aa_all_data1 <- my_learning_datasets_1000bp[[i]]
  }
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data1[, 1:(ncol(aa_all_data1) - 1)]) > 1){
    aadescrCor1 <- cor(aa_all_data1[, 1:(ncol(aa_all_data1) - 1)])
    aahighlyCorDescr1 <- findCorrelation(aadescrCor1, cutoff = .9)
    if(length(aahighlyCorDescr1) > 0){
      aa_all_data1 <- aa_all_data1[,-aahighlyCorDescr1]
    }
  }
  
  
  aanzv2 <- nearZeroVar(my_learning_datasets_1000bp_2[[i]], saveMetrics= F)
  if(length(aanzv2) > 0){
    aa_all_data2 <- my_learning_datasets_1000bp_2[[i]][,-aanzv2]
  }else{
    aa_all_data2 <- my_learning_datasets_1000bp_2[[i]]
  }
  
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data2[, 1:(ncol(aa_all_data2) - 1)]) > 1){
    aadescrCor2 <- cor(aa_all_data2[, 1:(ncol(aa_all_data2) - 1)])
    aahighlyCorDescr2 <- findCorrelation(aadescrCor2, cutoff = .9)
    if(length(aahighlyCorDescr2) > 0){
      aa_all_data2 <- aa_all_data2[,-aahighlyCorDescr2]
    }
  }
  aa_train_data1 <- aa_all_data1[my_data_partitions[[4]], ]
  aa_test_data1 <-  aa_all_data1[-my_data_partitions[[4]], ]
  
  aa_train_data2 <- aa_all_data2[my_data_partitions[[5]], ]
  aa_test_data2 <-  aa_all_data2[-my_data_partitions[[5]], ]
  
  
  
  aactrl <- trainControl(method = "repeatedcv",
                         repeats = 3,
                         classProbs = TRUE,
                         summaryFunction = defaultSummary,
                         savePredictions = TRUE,
                         ## new option here:
                         sampling = "down")
  print(i)
  print("aaRF_down")
  aaRF_down1 <- train(label_set ~ ., data = aa_train_data1,
                     method = "rf",
                     #preProc = c("center", "scale"),
                     #nbagg = 50,
                     metric = "Kappa",
                     trControl = aactrl)
  aaRF_down2 <- train(label_set ~ ., data = aa_train_data2,
                      method = "rf",
                      #preProc = c("center", "scale"),
                      #nbagg = 50,
                      metric = "Kappa",
                      trControl = aactrl)
  
  my_learning_datasets_results_RF_kappa_down_1000bp_commonset_1[[i]] <- aaRF_down1
  my_learning_datasets_results_RF_kappa_down_1000bp_commonset_2[[i]] <- aaRF_down2
}
save(list = c("my_learning_datasets_results_RF_kappa_down_1000bp_commonset_1",
              "my_learning_datasets_results_RF_kappa_down_1000bp_commonset_2"),
     file = "classification_results_common_test.RData")
################################################################################################
################################################################################################
################################################################################################
# create a partitioning that saves a part for test set in all the future modelings.


data_partition_mut_GEMSTAT <- list()
aa_train_valid <- createDataPartition(y = my_learning_datasets_1000bp_4[[1]][, ncol(my_learning_datasets_1000bp_4[[1]])],
                          times = 1, p = 0.80,list = F)
aa_test <- rownames(my_learning_datasets_1000bp_4[[1]])[setdiff(c(1:nrow(my_learning_datasets_1000bp_4[[1]])), aa_train_valid)]
table(unlist(lapply(strsplit(aa_test, split="_"), "[[", 1)))
aa_train_validSet <- my_learning_datasets_1000bp_4[[1]][aa_train_valid,]
aa_trainind <- createDataPartition(y = aa_train_validSet[, ncol(aa_train_validSet)],
                                      times = 1, p = 0.75,list = F)
aa_train <- rownames(aa_train_validSet)[aa_trainind]
table(unlist(lapply(strsplit(aa_train, split="_"), "[[", 1)))

aa_validind <- setdiff(c(1:nrow(aa_train_validSet)), aa_trainind)
aa_valid <- rownames(aa_train_validSet)[aa_validind]
table(unlist(lapply(strsplit(aa_valid, split="_"), "[[", 1)))
data_partition_mut_GEMSTAT[[1]] <- aa_train
data_partition_mut_GEMSTAT[[2]] <- aa_valid
data_partition_mut_GEMSTAT[[3]] <- aa_test
names(data_partition_mut_GEMSTAT) <- c("train", "validation", "test")

table(unlist(lapply(strsplit(data_partition_mut_GEMSTAT[[1]], split="_"), "[[", 1)))
table(unlist(lapply(strsplit(data_partition_mut_GEMSTAT[[2]], split="_"), "[[", 1)))
table(unlist(lapply(strsplit(data_partition_mut_GEMSTAT[[3]], split="_"), "[[", 1)))

# train a RF model on training and validation sets together, test the results on the test set
# create rf train and test data
my_learning_datasets_1000bp_4_ForGEMSTAT <- my_learning_datasets_1000bp_4[c(1, 4, 7)]

save(list = c("my_learning_datasets_1000bp_4_ForGEMSTAT", "data_partition_mut_GEMSTAT"), 
     file = "Classification_set4_GEMSTAT.RData")

my_learning_datasets_1000bp_4_ForGEMSTAT_results <- list()
for(i in 1:length(my_learning_datasets_1000bp_4_ForGEMSTAT)){
  aanzv <- nearZeroVar(my_learning_datasets_1000bp_4_ForGEMSTAT[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]]
  }
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  print("#########################################################################################")
  print(paste0("number of features in original feature set number ", i," named ", 
               names(my_learning_datasets_1000bp_4_ForGEMSTAT)[i] ,":"))
  print(ncol(my_learning_datasets_1000bp_4_ForGEMSTAT[[i]]))
  print("number of features after removing near zero variance and highly correlated features ( > 0.9)")
  print(ncol(aa_all_data))
  print("#########################################################################################")
  aatrainind <- match(c(data_partition_mut_GEMSTAT[[1]], data_partition_mut_GEMSTAT[[2]]),
                      rownames(my_learning_datasets_1000bp_4_ForGEMSTAT[[1]]))
  aagridLen <- c(5, 30, 50)
  aa_train_data <- aa_all_data[aatrainind, ]
  aa_test_data <-  aa_all_data[-aatrainind, ]
  aactrl <- trainControl(method = "repeatedcv",
                         number = 6,
                         repeats = 3,
                         classProbs = TRUE,
                         summaryFunction = mnLogLoss,
                         savePredictions = TRUE,
                         ## new option here:
                         sampling = "down", 
                         search = 'random')
  aaRF_down <- train(label_set ~ ., data = aa_train_data,
                     method = "rf",
                     #preProc = c("center", "scale"),
                     ntree = 1000,
                     metric = "logLoss",
                     trControl = aactrl,
                     tuneLength  = aagridLen[i])
  my_learning_datasets_1000bp_4_ForGEMSTAT_results[[i]] <- aaRF_down
}

save(list = c("my_learning_datasets_1000bp_4_ForGEMSTAT_results"),
     file = "classification_results_set_4_GEMSTAT_logloss.RData")



my_learning_datasets_1000bp_4_ForGEMSTAT_results_2 <- list()
for(i in 1:length(my_learning_datasets_1000bp_4_ForGEMSTAT)){
  aanzv <- nearZeroVar(my_learning_datasets_1000bp_4_ForGEMSTAT[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]]
  }
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  print("#########################################################################################")
  print(paste0("number of features in original feature set number ", i," named ", 
               names(my_learning_datasets_1000bp_4_ForGEMSTAT)[i] ,":"))
  print(ncol(my_learning_datasets_1000bp_4_ForGEMSTAT[[i]]))
  print("number of features after removing near zero variance and highly correlated features ( > 0.9)")
  print(ncol(aa_all_data))
  print("#########################################################################################")
  aatrainind <- match(c(data_partition_mut_GEMSTAT[[1]], data_partition_mut_GEMSTAT[[2]]),
                      rownames(my_learning_datasets_1000bp_4_ForGEMSTAT[[1]]))
  aagridLen <- c(5, 30, 50)
  aa_train_data <- aa_all_data[aatrainind, ]
  aa_test_data <-  aa_all_data[-aatrainind, ]
  aactrl <- trainControl(method = "repeatedcv",
                         number = 6,
                         repeats = 3,
                         classProbs = TRUE,
                         summaryFunction = defaultSummary,
                         savePredictions = TRUE,
                         ## new option here:
                         sampling = "down", 
                         search = 'random')
  aaRF_down <- train(label_set ~ ., data = aa_train_data,
                     method = "rf",
                     #preProc = c("center", "scale"),
                     ntree = 1000,
                     metric = "Kappa",
                     trControl = aactrl,
                     tuneLength  = aagridLen[i])
  my_learning_datasets_1000bp_4_ForGEMSTAT_results_2[[i]] <- aaRF_down
}
save(list = c("my_learning_datasets_1000bp_4_ForGEMSTAT_results_2"),
     file = "classification_results_set_4_GEMSTAT_Kappa.RData")


load("classification_results_set_4_GEMSTAT_logloss.RData")
names(my_learning_datasets_1000bp_4_ForGEMSTAT_results) <- c("Sum_LLR", "Sum_LR_plus_Adjacency", "kmer5")

aa_gglist_2 <- list()

for(i in 1:length(my_learning_datasets_1000bp_4_ForGEMSTAT)){
  aanzv <- nearZeroVar(my_learning_datasets_1000bp_4_ForGEMSTAT[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]]
  }
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  
  aatrainind <- match(c(data_partition_mut_GEMSTAT[[1]], data_partition_mut_GEMSTAT[[2]]),
                      rownames(my_learning_datasets_1000bp_4_ForGEMSTAT[[1]]))
  aagridLen <- c(5, 30, 50)
  aa_train_data <- aa_all_data[aatrainind, ]
  aa_test_data <-  aa_all_data[-aatrainind, ]
  
  aaindex_pos <- aa_test_data$label_set == "Pos"
  aaindex_neg <- aa_test_data$label_set == "Neg"
  
  aarfdown <- predict(my_learning_datasets_1000bp_4_ForGEMSTAT_results[[i]], 
                       aa_test_data, type = "prob")
  
  aapreds_list <- list(aarfdown$Pos)
  # List of actual values (same for all)
  m <- length(aapreds_list)
  aactuals_list <- rep(list(aa_test_data$label_set), m)
  # Plot the ROC curves
  aapred <- prediction(aapreds_list, aactuals_list)
  aauc <- performance(aapred, measure = "auc")
  aarocs <- performance(aapred, "tpr", "fpr")
  print(aauc@y.values)
  plot(aarocs, col = as.list(1:m), main = paste0("Test Set ROC Curves for feature set: ",
                                                 names(my_learning_datasets_1000bp_4_ForGEMSTAT_results)[i]))
  abline(a = 0, b = 1, col = 6, lty = 3, lwd = 0.7)
  legend(x = "bottomright", cex = 0.6,
         legend = c(paste0("RF_down_LOGLOSS ", "AUC: " , format(aauc@y.values[[1]], digits = 2))),
         fill = 1:m)
  aa_model_list <- list(RF_down = my_learning_datasets_1000bp_4_ForGEMSTAT_results[[i]])
  
  
  aamodel_list_pr <- aa_model_list %>%
    map(.f = calc_auprc, data = (aa_test_data))
  aamodel_list_pr_auc <- aamodel_list_pr %>%
    map(function(the_mod) the_mod$auc.integral)
  aamodel_list_pr_auc <- unlist(aamodel_list_pr_auc)
  
  aaresults_list_pr <- list(NA)
  num_mod <- 1
  
  for(the_pr in aamodel_list_pr){
    
    aaresults_list_pr[[num_mod]] <- data_frame(recall = the_pr$curve[, 1],
                                               precision = the_pr$curve[, 2],
                                               model = names(aamodel_list_pr)[num_mod])
    
    num_mod <- num_mod + 1
    
  }
  
  aaresults_df_pr <- bind_rows(aaresults_list_pr)
  
  aacustom_col <- colorspace::rainbow_hcl(length(aa_model_list))
  
  aa_gglist_2[[i]] <- ggplot(aes(x = recall, y = precision, group = model), data = aaresults_df_pr) +
    geom_line(aes(color = model), size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = aacustom_col, labels = paste(names(aa_model_list) ,"AUC: " ,
                                                             format(aamodel_list_pr_auc, digits = 2))) +
    geom_abline(intercept = sum(aa_test_data$label_set == "Pos")/nrow(aa_test_data),
                slope = 0, color = "gray", size = 1) +
    ggtitle(paste0("PR curve ", names(my_learning_datasets_1000bp_4_ForGEMSTAT_results)[i]))+
    theme_bw()
  
}
names(aa_gglist_2) <-  names(my_learning_datasets_1000bp_4_ForGEMSTAT_results)
aa_gglist_2[[3]]

varImp(my_learning_datasets_1000bp_4_ForGEMSTAT_results[[2]])
aa <- varImp(my_learning_datasets_1000bp_4_ForGEMSTAT_results[[2]],scale = F)
aax <- aa$importance
aaxn <- rownames(aax)
aax <- aax$Overall
names(aax) <- aaxn
aax[sort(aax, decreasing = T, index.return=T)$ix[1:20]]


length(sort(c("ESR1_2", "NR5A2", "RARA", "NKX3_1", "JUN_1", "YBX1", 
       "NR3C1", "PBX1", "SP1", "PPARD", "PPARG","NR2F2", "POU5F1", "AR", "FOXA1"))
)
rbind(c())

# TF sets based on my_learning_datasets_1000bp_4_ForGEMSTAT_results 
# set1
c("ESR1_2","FOXA1","JUN_1","NFIB","NR2F2","NR3C1","NR5A2",
  "PBX1","POU5F1","PPARD","RARA","RUNX1","SP1","YBX1")
# set2
c("AR","CEBPB","ESR1_1","FOXA1","JUN_1","NFIB","NKX3_1",
  "NR2C1","NR5A2","PBX1","PPARD","PPARG","RARA","RELA")
# coop set1 based on my_learning_datasets_1000bp_4_ForGEMSTAT_results 

# set1
# "ESR1_2", "FOXA1"
# "ESR1_2", "NR2F2"
# "ESR1_2", "SP1"
# "NR2F2", "PPARD"
# "NR2F2", "YBX1"
# "POU5F1", "SP1"
# "PPARD", "YBX1"
# "RUNX1", "YBX1"



"YBX1_vs_PPARD"
"SP1_vs_POU5F1"
"YBX1_vs_RUNX1"
# "YBX1_vs_NR2F2"
# "PPARD_vs_NR2F2"
# "NR2F2_vs_ESR1_2"



# set2

"ESR1_1_vs_PPARG"
"ESR1_1_vs_FOXA1"


##################################################################################################################
#####################################################################################################################################
load("classification_results_set_4_GEMSTAT_Kappa.RData")
names(my_learning_datasets_1000bp_4_ForGEMSTAT_results_2) <- c("Sum_LLR", "Sum_LR_plus_Adjacency", "kmer5")

aa_gglist_4 <- list()

for(i in 1:length(my_learning_datasets_1000bp_4_ForGEMSTAT)){
  aanzv <- nearZeroVar(my_learning_datasets_1000bp_4_ForGEMSTAT[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]]
  }
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  
  aatrainind <- match(c(data_partition_mut_GEMSTAT[[1]], data_partition_mut_GEMSTAT[[2]]),
                      rownames(my_learning_datasets_1000bp_4_ForGEMSTAT[[1]]))
  aagridLen <- c(5, 30, 50)
  aa_train_data <- aa_all_data[aatrainind, ]
  aa_test_data <-  aa_all_data[-aatrainind, ]
  
  aaindex_pos <- aa_test_data$label_set == "Pos"
  aaindex_neg <- aa_test_data$label_set == "Neg"
  
  aarfdown <- predict(my_learning_datasets_1000bp_4_ForGEMSTAT_results_2[[i]], 
                      aa_test_data, type = "prob")
  
  aapreds_list <- list(aarfdown$Pos)
  # List of actual values (same for all)
  m <- length(aapreds_list)
  aactuals_list <- rep(list(aa_test_data$label_set), m)
  # Plot the ROC curves
  aapred <- prediction(aapreds_list, aactuals_list)
  aauc <- performance(aapred, measure = "auc")
  aarocs <- performance(aapred, "tpr", "fpr")
  print(aauc@y.values)
  plot(aarocs, col = as.list(1:m), main = paste0("Test Set ROC Curves for feature set: ",
                                                 names(my_learning_datasets_1000bp_4_ForGEMSTAT_results_2)[i]))
  abline(a = 0, b = 1, col = 6, lty = 3, lwd = 0.7)
  legend(x = "bottomright", cex = 0.6,
         legend = c(paste0("RF_down_Kappa ", "AUC: " , format(aauc@y.values[[1]], digits = 2))),
         fill = 1:m)
  aa_model_list <- list(RF_down = my_learning_datasets_1000bp_4_ForGEMSTAT_results_2[[i]])
  
  
  aamodel_list_pr <- aa_model_list %>%
    map(.f = calc_auprc, data = (aa_test_data))
  aamodel_list_pr_auc <- aamodel_list_pr %>%
    map(function(the_mod) the_mod$auc.integral)
  aamodel_list_pr_auc <- unlist(aamodel_list_pr_auc)
  
  aaresults_list_pr <- list(NA)
  num_mod <- 1
  
  for(the_pr in aamodel_list_pr){
    
    aaresults_list_pr[[num_mod]] <- data_frame(recall = the_pr$curve[, 1],
                                               precision = the_pr$curve[, 2],
                                               model = names(aamodel_list_pr)[num_mod])
    
    num_mod <- num_mod + 1
    
  }
  
  aaresults_df_pr <- bind_rows(aaresults_list_pr)
  
  aacustom_col <- colorspace::rainbow_hcl(length(aa_model_list))
  
  aa_gglist_4[[i]] <- ggplot(aes(x = recall, y = precision, group = model), data = aaresults_df_pr) +
    geom_line(aes(color = model), size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = aacustom_col, labels = paste(names(aa_model_list) ,"AUC: " ,
                                                             format(aamodel_list_pr_auc, digits = 2))) +
    geom_abline(intercept = sum(aa_test_data$label_set == "Pos")/nrow(aa_test_data),
                slope = 0, color = "gray", size = 1) +
    ggtitle(paste0("PR curve ", names(my_learning_datasets_1000bp_4_ForGEMSTAT_results_2)[i]))+
    theme_bw()
  
}
names(aa_gglist_4) <-  names(my_learning_datasets_1000bp_4_ForGEMSTAT_results_2)
aa_gglist_4[[2]]
print(my_learning_datasets_1000bp_4_ForGEMSTAT_results_2[[2]])
aaimp <- my_learning_datasets_1000bp_4_ForGEMSTAT_results_2[[2]]$finalModel$param$ntree


aa <- varImp(my_learning_datasets_1000bp_4_ForGEMSTAT_results_2[[2]],scale = F)
aax <- aa$importance
aaxn <- rownames(aax)
aax <- aax$Overall
names(aax) <- aaxn
aax[sort(aax, decreasing = T, index.return=T)$ix[1:30]]
# set2
aa_names <- c("AR","ESR1_1","FOXA1","JUN_1","LEF1","NFIB","NKX3_1",
  "NR2C1","NR2F2","PPARD","PPARG","RARA","SP1","YBX1")
# set2

"ESR1_1_vs_ESR1_1"
"ESR_1_vs_FOXA1"
"ESR1_1_vs_NR2F2"
"ESR_1_vs_PPARG"
"ESR_1_vs_YBX1"
"NR2F2_vs_NR2F2"
aa_coop <- rbind(c(2,2), c(2,3) ,c(2, 9), c(2, 11), c(2, 14), c(9, 9))
####################################################################################################
# create a training/test partition using clustering of the feature vectors: clustering based cross validation


ncol(my_learning_datasets_1000bp_4_ForGEMSTAT$Sum_LR)
aaCCV <- scale(as.matrix(my_learning_datasets_1000bp_4_ForGEMSTAT$Sum_LR[, 1:33]))
aakmean <- kmeans(aaCCV, centers = 5)
aalk <- cbind(aakmean$cluster, my_learning_datasets_1000bp_4_ForGEMSTAT$Sum_LR[, 34])
for(i in 1:5){
  print(i)
  print(table(aalk[aalk[, 1] == i, 2]))
}
CCV_training_partition_list <- list()
CCV_training_partition_list[[1]] <- names(aakmean$cluster)[aakmean$cluster %in% c(1,2,3,4)]

save(list = c("CCV_training_partition_list"), file = "CCV_training_partition_list.RData")
# aahc <- hclust(dist(aaCCV))
# plot(aahc)
# aact <- cutree(aahc,h=12)
# table(aact)



my_learning_datasets_1000bp_4_ForGEMSTAT_results_ccv_1 <- list()
for(i in 1:2){
  aanzv <- nearZeroVar(my_learning_datasets_1000bp_4_ForGEMSTAT[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]]
  }
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  print("#########################################################################################")
  print(paste0("number of features in original feature set number ", i," named ", 
               names(my_learning_datasets_1000bp_4_ForGEMSTAT)[i] ,":"))
  print(ncol(my_learning_datasets_1000bp_4_ForGEMSTAT[[i]]))
  print("number of features after removing near zero variance and highly correlated features ( > 0.9)")
  print(ncol(aa_all_data))
  print("#########################################################################################")
  aatrainind <- match(CCV_training_partition_list[[1]] ,
                      rownames(my_learning_datasets_1000bp_4_ForGEMSTAT[[1]]))
  aagridLen <- c(5, 30, 50)
  aa_train_data <- aa_all_data[aatrainind, ]
  aa_test_data <-  aa_all_data[-aatrainind, ]
  aactrl <- trainControl(method = "repeatedcv",
                         number = 6,
                         repeats = 3,
                         classProbs = TRUE,
                         summaryFunction = defaultSummary,
                         savePredictions = TRUE,
                         ## new option here:
                         sampling = "down", 
                         search = 'random')
  aaRF_down <- train(label_set ~ ., data = aa_train_data,
                     method = "rf",
                     #preProc = c("center", "scale"),
                     ntree = 1000,
                     metric = "Kappa",
                     trControl = aactrl,
                     tuneLength  = aagridLen[i])
  my_learning_datasets_1000bp_4_ForGEMSTAT_results_ccv_1[[i]] <- aaRF_down
}
save(list = c("my_learning_datasets_1000bp_4_ForGEMSTAT_results_ccv_1"),
     file = "classification_results_set_4_GEMSTAT_CCV_1_Kappa.RData")

load("classification_results_set_4_GEMSTAT_CCV_1_Kappa.RData")
aa_gglist_4 <- list()

for(i in 1:2){
  aanzv <- nearZeroVar(my_learning_datasets_1000bp_4_ForGEMSTAT[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]]
  }
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  
  aatrainind <- match(CCV_training_partition_list[[1]],
                      rownames(my_learning_datasets_1000bp_4_ForGEMSTAT[[1]]))
  aagridLen <- c(5, 30, 50)
  aa_train_data <- aa_all_data[aatrainind, ]
  aa_test_data <-  aa_all_data[-aatrainind, ]
  
  aaindex_pos <- aa_test_data$label_set == "Pos"
  aaindex_neg <- aa_test_data$label_set == "Neg"
  
  aarfdown <- predict(my_learning_datasets_1000bp_4_ForGEMSTAT_results_ccv_1[[i]], 
                      aa_test_data, type = "prob")
  
  aapreds_list <- list(aarfdown$Pos)
  # List of actual values (same for all)
  m <- length(aapreds_list)
  aactuals_list <- rep(list(aa_test_data$label_set), m)
  # Plot the ROC curves
  aapred <- prediction(aapreds_list, aactuals_list)
  aauc <- performance(aapred, measure = "auc")
  aarocs <- performance(aapred, "tpr", "fpr")
  print(aauc@y.values)
  plot(aarocs, col = as.list(1:m), main = paste0("Test Set ROC Curves for feature set: ",
                                                 names(my_learning_datasets_1000bp_4_ForGEMSTAT_results_ccv_1)[i]))
  abline(a = 0, b = 1, col = 6, lty = 3, lwd = 0.7)
  legend(x = "bottomright", cex = 0.6,
         legend = c(paste0("RF_down_Kappa ", "AUC: " , format(aauc@y.values[[1]], digits = 2))),
         fill = 1:m)
  aa_model_list <- list(RF_down = my_learning_datasets_1000bp_4_ForGEMSTAT_results_ccv_1[[i]])
  
  
  aamodel_list_pr <- aa_model_list %>%
    map(.f = calc_auprc, data = (aa_test_data))
  aamodel_list_pr_auc <- aamodel_list_pr %>%
    map(function(the_mod) the_mod$auc.integral)
  aamodel_list_pr_auc <- unlist(aamodel_list_pr_auc)
  
  aaresults_list_pr <- list(NA)
  num_mod <- 1
  
  for(the_pr in aamodel_list_pr){
    
    aaresults_list_pr[[num_mod]] <- data_frame(recall = the_pr$curve[, 1],
                                               precision = the_pr$curve[, 2],
                                               model = names(aamodel_list_pr)[num_mod])
    
    num_mod <- num_mod + 1
    
  }
  
  aaresults_df_pr <- bind_rows(aaresults_list_pr)
  
  aacustom_col <- colorspace::rainbow_hcl(length(aa_model_list))
  
  aa_gglist_4[[i]] <- ggplot(aes(x = recall, y = precision, group = model), data = aaresults_df_pr) +
    geom_line(aes(color = model), size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = aacustom_col, labels = paste(names(aa_model_list) ,"AUC: " ,
                                                             format(aamodel_list_pr_auc, digits = 2))) +
    geom_abline(intercept = sum(aa_test_data$label_set == "Pos")/nrow(aa_test_data),
                slope = 0, color = "gray", size = 1) +
    ggtitle(paste0("PR curve ", names(my_learning_datasets_1000bp_4_ForGEMSTAT_results_ccv_1)[i]))+
    theme_bw()
  
}
names(aa_gglist_4) <-  names(my_learning_datasets_1000bp_4_ForGEMSTAT_results_ccv_1)
aa_gglist_4[[2]]
print(my_learning_datasets_1000bp_4_ForGEMSTAT_results_ccv_1[[2]])
aa <- varImp(my_learning_datasets_1000bp_4_ForGEMSTAT_results_ccv_1[[2]],scale = F)
aax <- aa$importance
aaxn <- rownames(aax)
aax <- aax$Overall
names(aax) <- aaxn
aax[sort(aax, decreasing = T, index.return=T)$ix[1:30]]

# train perf
for(i in 1:2){
  aanzv <- nearZeroVar(my_learning_datasets_1000bp_4_ForGEMSTAT[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]]
  }
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  
  aatrainind <- match(CCV_training_partition_list[[1]],
                      rownames(my_learning_datasets_1000bp_4_ForGEMSTAT[[1]]))
  aagridLen <- c(5, 30, 50)
  aa_train_data <- aa_all_data[aatrainind, ]
  aa_test_data <-  aa_all_data[-aatrainind, ]
  aarfdown <- predict(my_learning_datasets_1000bp_4_ForGEMSTAT_results_ccv_1[[i]], 
                      aa_train_data, type = "prob")
  aapreds_list <- list(aarfdown$Pos)
  # List of actual values (same for all)
  m <- length(aapreds_list)
  aactuals_list <- rep(list(aa_train_data$label_set), m)
  # Plot the ROC curves
  aapred <- prediction(aapreds_list, aactuals_list)
  aauc <- performance(aapred, measure = "auc")
  aarocs <- performance(aapred, "tpr", "fpr")
  print(aauc@y.values)
  plot(aarocs, col = as.list(1:m), main = paste0("Test Set ROC Curves for feature set: ",
                                                 names(my_learning_datasets_1000bp_4_ForGEMSTAT_results_ccv_1)[i]))
  abline(a = 0, b = 1, col = 6, lty = 3, lwd = 0.7)
  legend(x = "bottomright", cex = 0.6,
         legend = c(paste0("RF_down_Kappa ", "AUC: " , format(aauc@y.values[[1]], digits = 2))),
         fill = 1:m)
  aa_model_list <- list(RF_down = my_learning_datasets_1000bp_4_ForGEMSTAT_results_ccv_1[[i]])
  
  
  aamodel_list_pr <- aa_model_list %>%
    map(.f = calc_auprc, data = (aa_train_data))
  aamodel_list_pr_auc <- aamodel_list_pr %>%
    map(function(the_mod) the_mod$auc.integral)
  aamodel_list_pr_auc <- unlist(aamodel_list_pr_auc)
  
  aaresults_list_pr <- list(NA)
  num_mod <- 1
  
  for(the_pr in aamodel_list_pr){
    
    aaresults_list_pr[[num_mod]] <- data_frame(recall = the_pr$curve[, 1],
                                               precision = the_pr$curve[, 2],
                                               model = names(aamodel_list_pr)[num_mod])
    
    num_mod <- num_mod + 1
    
  }
  
  aaresults_df_pr <- bind_rows(aaresults_list_pr)
  
  aacustom_col <- colorspace::rainbow_hcl(length(aa_model_list))
  
  aa_gglist_4[[i]] <- ggplot(aes(x = recall, y = precision, group = model), data = aaresults_df_pr) +
    geom_line(aes(color = model), size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = aacustom_col, labels = paste(names(aa_model_list) ,"AUC: " ,
                                                             format(aamodel_list_pr_auc, digits = 2))) +
    geom_abline(intercept = sum(aa_train_data$label_set == "Pos")/nrow(aa_train_data),
                slope = 0, color = "gray", size = 1) +
    ggtitle(paste0("PR curve ", names(my_learning_datasets_1000bp_4_ForGEMSTAT_results_ccv_1)[i]))+
    theme_bw()
  
}
aa_gglist_4[[2]]
###############################################################################################
# can I do a full clustered cross validation?
# I can create train/test partitions in different distances and compare RF and GEMSTAT performances: 
# hope to be able to say that GEMSTAT modeling is more generalizable




