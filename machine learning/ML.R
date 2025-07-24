# R Code Annotation: Random Forest Classification Model Evaluation and Optimization

```r
# Load data and prepare training and test sets
load("ML.RData")  # Load R data file containing trainset and testset
ID <- rownames(trainset)
rownames(trainset) <- NULL
trainset <- cbind(ID, trainset)  # Add row names as ID column to dataframe
ID <- rownames(testset)
rownames(testset) <- NULL
testset <- cbind(ID, testset)  # Perform the same operation on test set

# Separate data by CLASS and recombine
class1 <- trainset[trainset$CLASS==1,]  # Extract positive class samples
class0 <- trainset[trainset$CLASS==0,]  # Extract negative class samples
train <- rbind(class1,class0)[,-1]  # Merge and remove ID column
train$CLASS <- factor(train$CLASS)  # Convert class to factor type
train.gene <- rbind(class1,class0)[,1]  # Save gene IDs
test.gene <- testset[,1]  # Save test set gene IDs
test <- testset[,-1]  # Test set feature matrix

# Load necessary libraries and set up 10-fold cross-validation
library(caret)
library(pROC)
folds <- createFolds(train$CLASS, k = 10)  # Create indices for 10-fold cross-validation

# Perform cross-validation for each fold
cv_results <- lapply(folds, function(x) {
  train_cv <- train[-x,]  # Training data for current fold
  test_cv <- train[x,]  # Validation data for current fold
  model <- randomForest(CLASS ~., data = train_cv)  # Train random forest model
  pred <- predict(model, newdata = test_cv)  # Predict on validation set
  
  # Calculate various performance metrics
  accuracy <- Accuracy(y_true = test_cv$CLASS, y_pred = pred)
  F1 <- F1_Score(y_true = test_cv$CLASS, y_pred = pred, positive = "1")
  recall <- Recall(y_true = test_cv$CLASS, y_pred = pred, positive = "1")
  auc <- AUC(y_true = test_cv$CLASS, y_pred = pred)
  cm <- confusionMatrix(table(test_cv$CLASS, pred))
  
  return(list("accuracy" = accuracy, "F1" = F1, "recall" = recall, "auc" = auc, "confusion_matrix" = cm))
})

# Create matrix to store cross-validation results
rnames <- c("accuracy", "recall", "F1","AUC")
cnames <- c(1:10, "mean", "sd")
mertics_10f <- matrix(NA, 4, 12, dimnames = list(rnames, cnames))

# Extract metrics from cross-validation results
com_accuracy <- c(cv_results$Fold01$accuracy,
                  cv_results$Fold02$accuracy,
                  # ... omitted code ...
                  cv_results$Fold10$accuracy)
# Similarly extract recall, F1, and auc metrics 

# Fill results into matrix
for (j in 1:10){
  mertics_10f[1,j] <- com_accuracy[j]
  mertics_10f[2,j] <- com_recall[j]
  mertics_10f[3,j] <- com_F1[j]
  mertics_10f[4,j] <- com_auc[j]
}

# Calculate mean and standard deviation for each metric
mertics_10f[1,11] <- mean(mertics_10f[1,1:10])  # accuracy mean
mertics_10f[2,11] <- mean(mertics_10f[2,1:10])  # recall mean
mertics_10f[3,11] <- mean(mertics_10f[3,1:10])  # F1 mean
mertics_10f[4,11] <- mean(mertics_10f[4,1:10])  # AUC mean
mertics_10f[4,12] <- sd(mertics_10f[4,1:10])  # AUC standard deviation
mertics_10f[3,12] <- sd(mertics_10f[3,1:10])  # F1 standard deviation
mertics_10f[2,12] <- sd(mertics_10f[2,1:10])  # Recall standard deviation
mertics_10f[1,12] <- sd(mertics_10f[1,1:10])  # Accuracy standard deviation

# Save cross-validation results to file
write.table(format(mertics_10f,digits=3), "mertics_10f_cv.txt", quote = F, row.names = T, sep = "\t")

# Hyperparameter optimization setup
library(caret)
# Define custom random forest model for parameter tuning with caret package
customRF <- list(
  type = "Classification", 
  library = "randomForest", 
  loop = NULL,
  parameters = data.frame(
    parameter = c("mtry", "ntree"), 
    class = rep("numeric", 2), 
    label = c("mtry", "ntree")
  ),
  # Define parameter grid, fitting function, prediction function, etc.
  grid = function(x, y, len = NULL, search = "grid") {
    # Generate combinations of mtry and ntree parameters
    if(search == "grid") {
      out <- expand.grid(
        mtry = caret::var_seq(p = ncol(x), classification = is.factor(y), len = len),
        ntree = c(500,700,900,1000,1500)
      )
    } else {
      # Random search
      out <- data.frame(
        mtry = unique(sample(1:ncol(x), size = len, replace = TRUE)),
        ntree = unique(sample(c(500,700,900,1000,1500), size = len, replace = TRUE))
      )
    }
  },
  fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
  },
  # Prediction and probability functions
  predict = function(modelFit, newdata, preProc = NULL, submodels = NULL)
    predict(modelFit, newdata),
  prob = function(modelFit, newdata, preProc = NULL, submodels = NULL)
    predict(modelFit, newdata, type = "prob"),
  sort = function(x) x[order(x[,1]),],
  levels = function(x) x$classes
)

# Set cross-validation control parameters: 10-fold cross-validation, repeated 3 times
trControl <- trainControl(method="repeatedcv", number=10, repeats=3)

# Define parameter grid: various combinations of mtry and ntree
tunegrid <- expand.grid(
  mtry = c(3,10,20,50,100,300,700,1000,2000), 
  ntree = c(500,700,800,1000,1500,2000)
)

# Perform parameter tuning
rf_custom <- train(
  CLASS ~., 
  data = train,
  method = customRF, 
  metric = "Accuracy", 
  tuneGrid = tunegrid, 
  trControl = trControl
)
rf_custom  # Display best parameters

# Final model training and prediction
mod <- randomForest(CLASS~., data=train, importance=TRUE)  # Train final model using all training data
#save(w, mod, file="all_rf_101.RData");  # Save model 

# Make predictions on test set
rf.pred.test <- predict(mod, type="prob", test)  # Get prediction probabilities
rf.pred.value <- ifelse(rf.pred.test[,2]>=0.5, 1, 0)  # Convert probabilities to binary classification
table(test$CLASS, rf.pred.value)  # Display confusion matrix

# Save prediction results
test.result <- data.frame(ID=test.gene, CLASS=rf.pred.value, prob=rf.pred.test[,2])
write.csv(test.result, "alldata.result.csv", row.names=FALSE)
```
