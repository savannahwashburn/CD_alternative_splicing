#script for k fold cross validation post op ileum dataset 

library(caret)
library(MASS)
library(pROC)

#------------logistic regression model with k-fold cross validation----------#

#load in Kyle's PC data 
data_k <- read.csv(file = "/Users/swashburn30/Desktop/isoform/top10expressionPCsPheno.csv")

#filter to just PC1-10 and patient id to merge dataframes 
data_k <- data_k[c("RNAid", "consortium_id", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]

#change column names 
colnames(data_k) <- c("sample_ID", "consortium_id", "PC1_k", "PC2_k", "PC3_k", "PC4_k", "PC5_k", "PC6_k", "PC7_k", "PC8_k", "PC9_k", "PC10_k")

#load in dataframe that contains PC scores and metadata
df_full <- read.csv(file = "/Users/swashburn30/Desktop/isoform/validate_kates_results/log_reg_df_2_29_24.csv")

# change race to asian for SB113705_S27
df_full["SB113705_S27", "race"] <- "Asian"

#data_k contains 4 more samples than df_full - need to filter those out 
df_filt_k <- data_k[data_k$sample_ID %in% df_full$sample, ]
df_k_reorder <- df_filt_k[order(match(df_filt_k$sample_ID,df_full$sample)),]
shared_values <- df_k_reorder$sample_ID[df_k_reorder$sample_ID %in% df_full$sample]

not_shared_values <- df_full$sample[!(df_full$sample %in% df_k_reorder$sample_ID)]

df_k_reorder$sample_ID[!(df_k_reorder$sample_ID %in% df_full$sample)]

#change column names again
colnames(df_k_reorder) <- c("sample", "consortium_id", "PC1_k", "PC2_k", "PC3_k", "PC4_k", "PC5_k", "PC6_k", "PC7_k", "PC8_k", "PC9_k", "PC10_k")

#merge the two dataframes  
df_merged <- merge(df_k_reorder, df_full, by = "sample", all = TRUE)

df_merged$affected <- as.factor(df_merged$affected)
df_merged$Batch <- as.factor(df_merged$Batch)


#set seed 
set.seed(5678)

# Create indices for stratified k-fold cross-validation
folds <- createFolds(df_merged$affected, k = 5, list = TRUE, returnTrain = FALSE)

#create variables to store metrics
summary_model <- list()

confusion_matrix <- list()

auc_score <- as.numeric()

roc <- list()

num_folds <- 5

# Perform k-fold cross-validation for logistic regression
cv_results <- for(fold in 1:num_folds) {
  # Extract training and testing sets
  print(fold)
  print(folds[[fold]])
  train_data <- df_merged[-folds[[fold]], ]
  test_data <- df_merged[folds[[fold]], ]
  
  #set R1 as positive level
  
  print("starting to model")
  # Fit logistic regression model
  model <- glm(affected ~ pc2 + PC1 + PC4 + PC1_k + PC2_k + PhenoSex + Batch + Smoking + race + AgeAtSample, data = train_data, family = "binomial")
  
  print("fit model")
  
  summary_model[[fold]] <- summary(model)
  
  # Make predictions on the test set
  predictions <- predict(model, newdata = test_data, type = "response")
  
  print("made prediction with test_new df and model")
  
  #define the accuracy of model
  accuracy <- ifelse(predictions >0.5, "R1", "R0")
  
  #create dataframe with predicted and actual results 
  results <- test_data %>%
    dplyr::select("affected") %>%
    bind_cols(accuracy)
  
  #fix the column names 
  colnames(results) <- c("affected", "predicted_affected")
  
  #make sure results are categorical variable 
  results$predicted_affected <- as.factor(results$predicted_affected)
  
  #make "R1" the postive level
  results$affected <- relevel(results$affected, ref = "R1")
  results$predicted_affected <- relevel(results$predicted_affected, ref = "R1")
  
  #create confusion matrix of results 
  cm <- confusionMatrix(results$affected, results$predicted_affected)
  
  confusion_matrix[[fold]] <- cm
  
  print("created cm")
  
  roc_score <- roc(test_data$affected, predictions)
  roc[[fold]] <- roc_score
  auc_score[fold] <- auc(roc_score)
  
}

#plot roc curve 
plot(roc[[1]], main = "ROC Curve", print.auc = T)

#plot all 5 folds on one graph
plot(roc[[1]], col = 1, lty = 2, main = "ROC")
plot(roc[[2]], col = 4, lty = 3, add = TRUE)
plot(roc[[3]], col = 3, lty = 1, add = TRUE)
plot(roc[[4]], col = 7, lty = 4, add = TRUE)
plot(roc[[5]], col = 6, lty = 5, add = TRUE) 
