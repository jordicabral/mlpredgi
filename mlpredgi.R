#' Building a Function with different Machine Learning algorithms
#'
#' This function runs different machine learning algorithms to predict
#' gene intereaction based on the other variables of the dataset.
#'
#' @param pred as the variable to predict
#' @param data as the dataset which contains the variables
#' @param exclude as optional variables to exclude in the model
#' @return A table with the summary of different models results
#' @author Jordi Cabral
#' @details
#' This function compares different regression models like Linear
#' Regression and some machine learning algorithms: K-NN,
#' Artificial Neural Networks, Support Vector Machine, Decision
#' Tree and Random Forest, and shows the RMSE (Root Mean Square Error)
#' of them to compare the performance with different parameters.
#'
#' The function splits also the dataset in training and test with a
#' proportion of 67-33/100 sample rows, to make the training and
#' validation of each machine learning algorithm.
#' The function allows also the possibility to exclude one or some
#' variables in the model.
#' @seealso \code{lm}
#' @seealso \code{train}
#' @seealso \code{neuralnet}
#' @seealso \code{ksvm}
#' @seealso \code{rpart}
#' @seealso \code{randomForest}
#' @export mlpredgi
#' @importFrom stats lm
#' @import caret
#' @import kableExtra
#' @import neuralnet
#' @import kernlab
#' @import rpart
#' @import randomForest

#' @examples
#' mlpredgi("SGAsco", SGA10_full, exclude="NPC")
#' mlpredgi("SGAsco", SGA10_full, exclude=c("PCC", "fcMF"))
#'
# Building formula
mlpredgi <- function(pred, data, exclude) {
     if(missing(exclude)) {
          f <- as.formula(paste(pred, "~."))
     } else {
          exclude_list <- paste("-", exclude, collapse = " ")
          f <- as.formula(paste(pred, "~.", exclude_list))
     }
     # Linear regression model
     lmod <- lm(f, data=data)
     model_labels <- c("Linear Regression", "KNN default", "KNN param",
                       "ANN default", "ANN 2 HN", "ANN 3HN", "SVM vanilladot",
                       "SVM rbfdot", "SVM laplacedot", "Decision Tree default",
                       "Decision Tree param", "Random Forest default",
                       "Random Forest param")
     results <- data.frame(Model = model_labels, RMSE = NA)
     results[1,2] <- sqrt(mean(lmod$residuals^2))

     # Machine Learning models
     # Split dataset to training and test
     bound <- floor(nrow(data)*0.67)

     seed <- c(12345)

     set.seed(seed)
     row_train <- sample(seq_len(nrow(data)), size = bound)

     data_train <- data[row_train, ]
     data_test <- data[-row_train, ]

     # KNN algorithm
     # Train the model
     set.seed(seed)
     knn_model_1 <- train(f, data=data_train, method="knn")

     # Make the prediction
     pred_knn_1 <- knn_model_1 %>% predict(data_test)

     # Check performance of the model
     results[2,2] <- RMSE(pred_knn_1, data_test[[pred]])

     # Improve performance of the model
     set.seed(seed)
     knn_model2 <- train(f, data=data_train, method="knn",
                         trControl = trainControl("cv", number = 10),
                         preProcess = c("center","scale"),
                         tuneLength = 15)
     pred_knn2 <- knn_model2 %>% predict(data_test)
     results[3,2] <- RMSE(pred_knn2, data_test[[pred]])

     # Artificial Neural Network algorithn
     # Train the model
     set.seed(seed)
     ann_model_1 <- neuralnet(f, data=data_train)

     # Make the prediction
     ann_results_1 <- compute(ann_model_1, data_test)
     ann_pred_1 <- ann_results_1$net.result

     # Check performance of the model
     results[4,2] <- RMSE(ann_pred_1, data_test[[pred]])

     # Improve performance of the model
     set.seed(seed)
     ann_model_2 <- neuralnet(f, data=data_train, hidden = 2)
     ann_results_2 <- compute(ann_model_2, data_test)
     ann_pred_2 <- ann_results_2$net.result
     results[5,2] <- RMSE(ann_pred_2, data_test[[pred]])

     set.seed(seed)
     ann_model_3 <- neuralnet(f, data=data_train, hidden = 3)
     ann_results_3 <- compute(ann_model_3, data_test)
     ann_pred_3 <- ann_results_3$net.result
     results[6,2] <- RMSE(ann_pred_3, data_test[[pred]])

     # Support Vector Machine algorithm
     # Train the model
     set.seed(seed)
     svm_model_1 <- ksvm(f, data=data_train, kernel = "vanilladot")

     # Make the prediction
     svm_pred_1 <- predict(svm_model_1, data_test)

     # Check performance of the model
     results[7,2] <- RMSE(svm_pred_1, data_test[[pred]])

     # Improve performance of the model
     set.seed(seed)
     svm_model_2 <- ksvm(f, data=data_train, kernel = "rbfdot")
     svm_pred_2 <- predict(svm_model_2, data_test)
     results[8,2] <- RMSE(svm_pred_2, data_test[[pred]])

     set.seed(seed)
     svm_model_3 <- ksvm(f, data=data_train, kernel = "laplacedot")
     svm_pred_3 <- predict(svm_model_3, data_test)
     results[9,2] <- RMSE(svm_pred_3, data_test[[pred]])

     # Decision Tree algorithm
     # Train the model
     set.seed(seed)
     dt_model_1 <- rpart(f, data=data_train)

     # Make the prediction
     dt_pred_1 <- predict(dt_model_1, data_test)

     # Check performance of the model
     results[10,2] <- RMSE(dt_pred_1, data_test[[pred]])

     # Improve performance of the model
     set.seed(seed)
     dt_model_2 <- rpart(f, data=data_train, minsplit = 2, minbucket = 1)
     dt_pred_2 <- predict(dt_model_2, data_test)
     results[11,2] <- RMSE(dt_pred_2, data_test[[pred]])

     # Random Forest algorithm
     # Train the model
     set.seed(seed)
     rf_model_1<- randomForest(f, data=data_train)

     # Make the prediction
     rf_pred_1 <- predict(rf_model_1, data_test)

     # Check performance of the model
     results[12,2] <- RMSE(rf_pred_1, data_test[[pred]])

     # Improve performance of the model
     set.seed(seed)
     rf_model_1000<- randomForest(f, data=data_train, ntree=1000)
     rf_pred_1000 <- predict(rf_model_1000, data_test)
     results[13,2] <- RMSE(rf_pred_1000, data_test[[pred]])

     return(results)
}
