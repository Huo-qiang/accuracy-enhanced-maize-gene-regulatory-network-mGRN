load("ML.RData")
ID <- rownames(trainset)
rownames(trainset) <- NULL
trainset <- cbind(ID, trainset)
ID <- rownames(testset)
rownames(testset) <- NULL
testset <- cbind(ID, testset)
class1<-trainset[trainset$CLASS==1,]
class0<-trainset[trainset$CLASS==0,]
train<-rbind(class1,class0)[,-1]
train$CLASS<-factor(train$CLASS)
train.gene<-rbind(class1,class0)[,1]
test.gene<-testset[,1]
test<-testset[,-1]

library(caret)
library(pROC)
folds <- createFolds(train$CLASS, k = 10)

cv_results <- lapply(folds, function(x) {
  train_cv <- train[-x,]
  test_cv <- train[x,]
  model <- randomForest(CLASS ~., data = train_cv)
  pred <- predict(model, newdata = test_cv)
  accuracy <- Accuracy(y_true = test_cv$CLASS, y_pred = pred)
  F1 <- F1_Score(y_true = test_cv$CLASS, y_pred = pred, positive = "1")
  recall <- Recall(y_true = test_cv$CLASS, y_pred = pred, positive = "1")
  auc <- AUC(y_true = test_cv$CLASS, y_pred = pred)
  cm <- confusionMatrix(table(test_cv$CLASS, pred))
  return(list("accuracy" = accuracy, "F1" = F1, "recall" = recall, "auc" = auc, "confusion_matrix" = cm))
})

rnames <- c("accuracy", "recall", "F1","AUC")
cnames <- c(1:10, "mean", "sd")
mertics_10f <- matrix(NA, 4, 12, dimnames = list(rnames, cnames))

com_accuracy <- c(cv_results$Fold01$accuracy,
                  cv_results$Fold02$accuracy,
                  cv_results$Fold03$accuracy,
                  cv_results$Fold04$accuracy,
                  cv_results$Fold05$accuracy,
                  cv_results$Fold06$accuracy,
                  cv_results$Fold07$accuracy,
                  cv_results$Fold08$accuracy,
                  cv_results$Fold09$accuracy,
                  cv_results$Fold10$accuracy)
com_recall <- c(cv_results$Fold01$recall,
                cv_results$Fold02$recall,
                cv_results$Fold03$recall,
                cv_results$Fold04$recall,
                cv_results$Fold05$recall,
                cv_results$Fold06$recall,
                cv_results$Fold07$recall,
                cv_results$Fold08$recall,
                cv_results$Fold09$recall,
                cv_results$Fold10$recall)
com_F1 <- c(cv_results$Fold01$F1,
            cv_results$Fold02$F1,
            cv_results$Fold03$F1,
            cv_results$Fold04$F1,
            cv_results$Fold05$F1,
            cv_results$Fold06$F1,
            cv_results$Fold07$F1,
            cv_results$Fold08$F1,
            cv_results$Fold09$F1,
            cv_results$Fold10$F1)
com_auc <- c(cv_results$Fold01$auc,
             cv_results$Fold02$auc,
             cv_results$Fold03$auc,
             cv_results$Fold04$auc,
             cv_results$Fold05$auc,
             cv_results$Fold06$auc,
             cv_results$Fold07$auc,
             cv_results$Fold08$auc,
             cv_results$Fold09$auc,
             cv_results$Fold10$auc)

for (j in 1:10){
  mertics_10f[1,j] <- com_accuracy[j]
  mertics_10f[2,j] <- com_recall[j]
  mertics_10f[3,j] <- com_F1[j]
  mertics_10f[4,j] <- com_auc[j]
}

mertics_10f[1,11] <- mean(mertics_10f[1,1:10])
mertics_10f[2,11] <- mean(mertics_10f[2,1:10])
mertics_10f[3,11] <- mean(mertics_10f[3,1:10])
mertics_10f[4,11] <- mean(mertics_10f[4,1:10])
mertics_10f[4,12] <- sd(mertics_10f[4,1:10])
mertics_10f[3,12] <- sd(mertics_10f[3,1:10])
mertics_10f[2,12] <- sd(mertics_10f[2,1:10])
mertics_10f[1,12] <- sd(mertics_10f[1,1:10])

write.table(format(mertics_10f,digits=3), "mertics_10f_101cv.txt", quote = F, row.names = T, sep = "\t")
##Adjustment parameters
library(caret)
customRF <- list(type = "Classification", library = "randomForest", loop = NULL,
                 parameters = data.frame(parameter = c("mtry", "ntree"), 
                                         class = rep("numeric", 2), 
                                         label = c("mtry", "ntree")),
                  grid = function(x, y, len = NULL, search = "grid") {
                    if(search == "grid") {
                      out <- expand.grid(mtry = caret::var_seq(p = ncol(x),
                                                              classification = is.factor(y),
                                                              len = len),
                                         ntree = c(500,700,900,1000,1500))
                    } else {
                      out <- data.frame(mtry = unique(sample(1:ncol(x), size = len, replace = TRUE)),
                                        ntree = unique(sample(c(500,700,900,1000,1500), 
                                                              size = len, replace = TRUE)))
                    }
                  },
                 fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
                    randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
                  },
                predict = function(modelFit, newdata, preProc = NULL, submodels = NULL)
                  predict(modelFit, newdata),
                prob = function(modelFit, newdata, preProc = NULL, submodels = NULL)
                  predict(modelFit, newdata, type = "prob"),
                sort = function(x) x[order(x[,1]),],
                levels = function(x) x$classes
                )
trControl <- trainControl(method="repeatedcv", number=10, repeats=3)
tunegrid <- expand.grid(mtry=c(3,10,20,50,100,300,700,1000,2000), 
                        ntree=c(500,700, 800, 1000, 1500, 2000))

rf_custom <- train(CLASS ~., data = train,method=customRF, 
                    metric="Accuracy", tuneGrid=tunegrid, 
                    trControl=trControl)
rf_custom
##predictions
mod <-randomForest(CLASS~.,data=train,importance = T)
#save(w, mod, file="all_rf_101.RData");
rf.pred.test<- predict(mod,type = "prob",test)
rf.pred.value <- ifelse(rf.pred.test[,2]>=0.5,1,0)
table(test$CLASS,rf.pred.value)
test.result<-data.frame(ID=test1$ID,CLASS=rf.pred.value,prob=rf.pred.test[,2])
write.csv(test.result,"seed101.alldata.result.csv",row.names = F)
##
