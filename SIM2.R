# Set up
#library(parallel)
library(mice)
library(nnet)
library(tidyverse)
#ncores <- detectCores()
M <- 5

#### FUNCTIONS ######################################################################################

# 95% Confidence Intervals.
# Each row contains the lower and upper bound for one outcome category.
CI <- function(proportion, SE){
  lower <- proportion - qnorm(.975)*SE
  upper <- proportion + qnorm(.975)*SE
  return(c(lower, upper))
}

# Distance function that calculates distance between subject i with missing Y in the 
# original data set and subject j with observed Y in the bootstrap sample.
distance <- function(pred.i, pred.j){
  omega <- 1 / length(pred.i)
  sqrt(sum(omega * (pred.i - pred.j)^2))
}

# SE using Rubin's rules with as input one ROW of the proportions matrix:
standard_error <- function(proportions, M){
  sqrt(1/M * sum(proportions*(1-proportions)/300) + (1 + 1/M) / (M-1) * 
         sum((proportions - mean(proportions))^2))
}

#####################################################################################################
#####################################################################################################
##### POPULATION 1: 3 CATEGORIES FOR Y AND FIRST ORDER TERMS ########################################

load("~/SSLBS/Thesis/Thesis/pop1.RData")

### SAMPLING, MISSINGNESS, IMPUTATION, AND ANALYSIS
# Here we first create our samples and introduce missing values, then imputation and evaluation take 
# place. Creating a function allows us to use mc.lapply later on. If responserate80 = TRUE, the 
# response rate is 80%, if responserate80 = FALSE, it is 90%. Including the response rate as an 
# option in the function saves code (i.e. is more efficient).
# Note that the argument 'iterations' is not used within the code of the function, it is just there 
# so we can use it in an lapply or mclapply loop later on.

results_per_samp_first <- function(iterations, population, responserate80 = TRUE){
  # Sample from the population.
  samp <- population[sample(1:50000, 300, replace = FALSE),]
  # Generate response indicators depending on responserate80.
  if(responserate80){
    p_res80 <- exp(samp$x1 + samp$x2 + 0.01*samp$x3) / 
                            (exp(samp$x1 + samp$x2 + 0.01*samp$x3) + 1)
    samp$response_indicators <- rbinom(300, 1, p_res80)}else{
      p_res90 <- exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) / 
                              (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) + 1)
      samp$response_indicators <- rbinom(300, 1, p_res90)
    }
  # Introduce missing values accordingly.
  samp$y[samp$response_indicators == 0] <- NA
  
  # Draw a new sample if the sample does not contain observations from all categories.
  while(length(unique(samp$y))!=4){
    # Sample from the population.
    samp <- population[sample(1:50000, 300, replace = FALSE),]
    # Generate response indicators depending on responserate80.
    if(responserate80){
      p_res80 <- exp(samp$x1 + samp$x2 + 0.01*samp$x3) / 
                              (exp(samp$x1 + samp$x2 + 0.01*samp$x3) + 1)
      samp$response_indicators <- rbinom(300, 1, p_res80)}else{
        p_res90 <- exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) / 
                                (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) + 1)
        samp$response_indicators <- rbinom(300, 1, p_res90)
      }
    # Introduce missing values accordingly.
    samp$y[samp$response_indicators == 0] <- NA
  }
  
  ## Impute using MRNN10 (incorrect missingness model)
  MRNN10 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=3){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2: Using the bootstrap sample, fit one or more multinomial regression models on the 
    # dependent variable using the auxiliary data and save the standardized predicted values.
    mod1.1 <- multinom(y ~ x1 + x2 + x3, data = boot)
    pred1.1 <- fitted(mod1.1)[,-1]
    mod1.2 <- multinom(y ~ x2^2 + x1:x3 + x1 + x2 + x3, data = boot)
    pred1.2 <- fitted(mod1.2)[,-1]
    # If two or more regression models were applied, regress the dependent variable on the predicted
    # values and save the new predicted values. 
    temporary <- cbind(boot[boot$response_indicators == 1, "y"], pred1.1, pred1.2)
    temporary <- as.data.frame(temporary)
    names(temporary) <- c("y", "M1C2", "M1C3", "M2C2", "M2C3")
    mod1.3 <- multinom(y ~ M1C2 + M1C3 + M2C2 + M2C3, temporary)
    pred1.3 <- fitted(mod1.3)[,-1]
    # Note that these models are applied only to respondents.
    # Standardize the predicted values.
    mu1 <- colMeans(pred1.3) # we need to save these for later, that's why we do it by hand 
    sd1 <- apply(pred1.3, 2, function(x) sqrt(var(x)))
    for(i in 1:ncol(pred1.3)){
      pred1.3[,i] <- (pred1.3[,i] - mu1[i]) / sd1[i]
    }
    
    # step 3: Using the bootstrap sample, fit one or more logistic regression models on the response 
    # indicators using the auxiliary data and save the standardized predicted values.
    mod2.1 <- glm(response_indicators ~ x1 + x2, boot, family = binomial)
    pred2.1 <- mod2.1$fitted.values
    mod2.2 <- glm(response_indicators ~ x1*x3, boot, family = binomial)
    pred2.2 <- mod2.2$fitted.values
    # If two or more regression models were applied, regress the response indicators on the predicted 
    # values and save the new predicted values.
    temporary <- cbind(boot$response_indicators, pred2.1, pred2.2)
    temporary <- as.data.frame(temporary)
    mod2.3 <- glm(V1 ~ pred2.1 + pred2.2, temporary, family = binomial)
    pred2.3 <- mod2.3$fitted.values
    # Note that these models are applied to all cases in the bootstrap sample.
    # Standardize the predicted values.
    mu2 <- mean(pred2.3); sd2 <- sqrt(var(pred2.3)) 
    pred2.3 <- (pred2.3 - mu2) / sd2
    
    # step 4: Calculate the predicted outcome probabilities and propensity scores for the
    # nonrespondents in the original sample and standardize these using the parameters obtained in
    # steps 2 and 3.
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    # Predicted outcome probabilities:
    missing_pred1.1 <- predict(mod1.1, newdata = to_be_imputed, type = "probs")[,-1]
    missing_pred1.2 <- predict(mod1.2, newdata = to_be_imputed, type = "probs")[,-1]
    missing_temp <- as.data.frame(cbind(missing_pred1.1, missing_pred1.2))
    names(missing_temp) <- c("M1C2", "M1C3", "M2C2", "M2C3")
    missing_pred1.3 <- predict(mod1.3, newdata = missing_temp, type = "probs")[,-1]
    # Predicted response propensities:
    missing_pred2.1 <- predict(mod2.1, newdata = to_be_imputed, type = "response")
    missing_pred2.2 <- predict(mod2.2, newdata = to_be_imputed, type = "response")
    missing_temp <- as.data.frame(cbind(missing_pred2.1, missing_pred2.2))
    names(missing_temp) <- c("pred2.1", "pred2.2")
    missing_pred2.3 <- predict(mod2.3, newdata = missing_temp, type = "response")
    # Standardize the predicted values.
    for(i in 1:ncol(missing_pred1.3)){
      missing_pred1.3[,i] <- (missing_pred1.3[,i] - mu1[i]) / sd1[i]
    }
    missing_pred2.3 <- (missing_pred2.3 - mu2) / sd2
    
    # step 5: calculate a distance function to define the similarity between subject i with missing Y 
    # in  the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(missing_pred1.3, missing_pred2.3) 
    pred.boot_obs <- cbind(pred1.3, pred2.3[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 6: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # Rank the distances.
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set.
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute.
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  MRNN10_sets <- lapply(1:M, MRNN10, samp)
  
  ## Impute using MRNN01 (incorrect outcome model)
  MRNN01 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=3){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2: Using the bootstrap sample, fit one or more multinomial regression models on the 
    # dependent variable using the auxiliary data and save the predicted values.
    mod1.1 <- multinom(y ~ x1 + x2, data = boot)
    pred1.1 <- fitted(mod1.1)[,-1]
    mod1.2 <- multinom(y ~ x2^2 + x1 + x2, data = boot)
    pred1.2 <- fitted(mod1.2)[,-1]
    # If two or more regression models were applied, regress the dependent variable on the predicted
    # values and save the new predicted values. 
    temporary <- cbind(boot[boot$response_indicators == 1, "y"], pred1.1, pred1.2)
    temporary <- as.data.frame(temporary)
    names(temporary) <- c("y", "M1C2", "M1C3", "M2C2", "M2C3")
    mod1.3 <- multinom(y ~ M1C2 + M1C3 + M2C2 + M2C3, temporary)
    pred1.3 <- fitted(mod1.3)[,-1]
    # Note that these models are applied only to respondents.
    # Standardize the predicted values.
    mu1 <- colMeans(pred1.3) # we need to save these for later, that's why we do it by hand 
    sd1 <- apply(pred1.3, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1.3)){
      pred1.3[,i] <- (pred1.3[,i] - mu1[i]) / sd1[i]
    }
    
    # step 3: Using the bootstrap sample, fit one or more logistic regression models on the response 
    # indicators using the auxiliary data and save the predicted values.
    mod2.1 <- glm(response_indicators ~ x1 + x2 + x3, boot, family = binomial)
    pred2.1 <- mod2.1$fitted.values
    mod2.2 <- glm(response_indicators ~ x1*x3, boot, family = binomial)
    pred2.2 <- mod2.2$fitted.values
    # If two or more regression models were applied, regress the response indicators on the predicted 
    # values and save the new predicted values.
    temporary <- cbind(boot$response_indicators, pred2.1, pred2.2)
    temporary <- as.data.frame(temporary)
    mod2.3 <- glm(V1 ~ pred2.1 + pred2.2, temporary, family = binomial)
    pred2.3 <- mod2.3$fitted.values
    # Note that these models are applied to all cases in the bootstrap sample.
    # Standardize the predicted values.
    mu2 <- mean(pred2.3); sd2 <- sqrt(var(pred2.3)) 
    pred2.3 <- (pred2.3 - mu2) / sd2
    
    # step 4: Calculate the predicted outcome probabilities and propensity scores for the
    # nonrespondents in the original sample and standardize these using the parameters obtained in
    # steps 2 and 3.
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    # Predicted outcome probabilities:
    missing_pred1.1 <- predict(mod1.1, newdata = to_be_imputed, type = "probs")[,-1]
    missing_pred1.2 <- predict(mod1.2, newdata = to_be_imputed, type = "probs")[,-1]
    missing_temp <- as.data.frame(cbind(missing_pred1.1, missing_pred1.2))
    names(missing_temp) <- c("M1C2", "M1C3", "M2C2", "M2C3")
    missing_pred1.3 <- predict(mod1.3, newdata = missing_temp, type = "probs")[,-1]
    # Predicted response propensities:
    missing_pred2.1 <- predict(mod2.1, newdata = to_be_imputed, type = "response")
    missing_pred2.2 <- predict(mod2.2, newdata = to_be_imputed, type = "response")
    missing_temp <- as.data.frame(cbind(missing_pred2.1, missing_pred2.2))
    names(missing_temp) <- c("pred2.1", "pred2.2")
    missing_pred2.3 <- predict(mod2.3, newdata = missing_temp, type = "response")
    # Standardize the predicted values.
    for(i in 1:ncol(missing_pred1.3)){
      missing_pred1.3[,i] <- (missing_pred1.3[,i] - mu1[i]) / sd1[i]
    }
    missing_pred2.3 <- (missing_pred2.3 - mu2) / sd2
    
    # step 5: calculate a distance function to define the similarity between subject i with missing Y 
    # in  the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(missing_pred1.3, missing_pred2.3) 
    pred.boot_obs <- cbind(pred1.3, pred2.3[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 6: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # Rank the distances.
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set.
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute.
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  MRNN01_sets <- lapply(1:M, MRNN01, samp)
  
  ## Impute using MRNN00 (both incorrect)
  MRNN00 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=3){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2: Using the bootstrap sample, fit one or more multinomial regression models on the 
    # dependent variable using the auxiliary data and save the predicted values.
    mod1.1 <- multinom(y ~ x1 + x2, data = boot)
    pred1.1 <- fitted(mod1.1)[,-1]
    mod1.2 <- multinom(y ~ x2^2 + x1 + x2, data = boot)
    pred1.2 <- fitted(mod1.2)[,-1]
    # If two or more regression models were applied, regress the dependent variable on the predicted
    # values and save the new predicted values. 
    temporary <- cbind(boot[boot$response_indicators == 1, "y"], pred1.1, pred1.2)
    temporary <- as.data.frame(temporary)
    names(temporary) <- c("y", "M1C2", "M1C3", "M2C2", "M2C3")
    mod1.3 <- multinom(y ~ M1C2 + M1C3 + M2C2 + M2C3, temporary)
    pred1.3 <- fitted(mod1.3)[,-1]
    # Note that these models are applied only to respondents.
    # Standardize the predicted values.
    mu1 <- colMeans(pred1.3) # we need to save these for later, that's why we do it by hand 
    sd1 <- apply(pred1.3, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1.3)){
      pred1.3[,i] <- (pred1.3[,i] - mu1[i]) / sd1[i]
    }
    
    # step 3: Using the bootstrap sample, fit one or more logistic regression models on the response 
    # indicators using the auxiliary data and save the predicted values.
    mod2.1 <- glm(response_indicators ~ x1 + x2, boot, family = binomial)
    pred2.1 <- mod2.1$fitted.values
    mod2.2 <- glm(response_indicators ~ x1*x3, boot, family = binomial)
    pred2.2 <- mod2.2$fitted.values
    # If two or more regression models were applied, regress the response indicators on the predicted 
    # values and save the new predicted values.
    temporary <- cbind(boot$response_indicators, pred2.1, pred2.2)
    temporary <- as.data.frame(temporary)
    mod2.3 <- glm(V1 ~ pred2.1 + pred2.2, temporary, family = binomial)
    pred2.3 <- mod2.3$fitted.values
    # Note that these models are applied to all cases in the bootstrap sample.
    # Standardize the predicted values.
    mu2 <- mean(pred2.3); sd2 <- sqrt(var(pred2.3)) 
    pred2.3 <- (pred2.3 - mu2) / sd2
    
    # step 4: Calculate the predicted outcome probabilities and propensity scores for the
    # nonrespondents in the original sample and standardize these using the parameters obtained in
    # steps 2 and 3.
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    # Predicted outcome probabilities:
    missing_pred1.1 <- predict(mod1.1, newdata = to_be_imputed, type = "probs")[,-1]
    missing_pred1.2 <- predict(mod1.2, newdata = to_be_imputed, type = "probs")[,-1]
    missing_temp <- as.data.frame(cbind(missing_pred1.1, missing_pred1.2))
    names(missing_temp) <- c("M1C2", "M1C3", "M2C2", "M2C3")
    missing_pred1.3 <- predict(mod1.3, newdata = missing_temp, type = "probs")[,-1]
    # Predicted response propensities:
    missing_pred2.1 <- predict(mod2.1, newdata = to_be_imputed, type = "response")
    missing_pred2.2 <- predict(mod2.2, newdata = to_be_imputed, type = "response")
    missing_temp <- as.data.frame(cbind(missing_pred2.1, missing_pred2.2))
    names(missing_temp) <- c("pred2.1", "pred2.2")
    missing_pred2.3 <- predict(mod2.3, newdata = missing_temp, type = "response")
    # Standardize the predicted values.
    for(i in 1:ncol(missing_pred1.3)){
      missing_pred1.3[,i] <- (missing_pred1.3[,i] - mu1[i]) / sd1[i]
    }
    missing_pred2.3 <- (missing_pred2.3 - mu2) / sd2
    
    # step 5: calculate a distance function to define the similarity between subject i with missing Y 
    # in  the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(missing_pred1.3, missing_pred2.3) 
    pred.boot_obs <- cbind(pred1.3, pred2.3[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 6: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # Rank the distances.
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set.
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute.
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  MRNN00_sets <- lapply(1:M, MRNN00, samp)
  
  ## Impute using DRNN10 (incorrect missingness model)
  DRNN10 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=3){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2a: Using the bootstrap sample, fit a multinomial regression model on the dependent variable
    # using the auxiliary data and save the standardized predicted values.
    mod1 <- multinom(y ~ x1 + x2 + x3, boot)
    pred1 <- fitted(mod1)[,-1]
    # We are manually standardizing because we need to save the means and sds for later.
    mu1 <- colMeans(pred1)
    sd1 <- apply(pred1, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1)){
      pred1[,i] <- (pred1[,i] - mu1[i]) / sd1[i]
    }
    # step 2b: Using the bootstrap sample, fit a logistic regression model on the response indicators 
    # using the auxiliary data and save the standardized predicted values.
    mod2 <- glm(response_indicators ~ x1 + x2, boot, family = binomial)
    pred2 <- mod2$fitted.values
    mu2 <- mean(pred2); sd2 <- sqrt(var(pred2)) 
    pred2 <- (pred2 - mu2) / sd2
    # Note that step 2a only uses respondents and step 2b uses the entire bootstrap sample.
    
    # step 3: calculate a distance function to define the similarity between subject i with missing Y
    # in the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    pred1.2 <- predict(mod1, newdata = to_be_imputed, type = "probs")[,-1]
    pred2.2 <- predict(mod2, newdata = to_be_imputed, type = "response")
    # Standardize the predicted values
    for(i in 1:ncol(pred1.2)){
      pred1.2[,i] <- (pred1.2[,i] - mu1[i]) / sd1[i]
    }
    pred2.2 <- (pred2.2 - mu2) / sd2
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(pred1.2, pred2.2) 
    pred.boot_obs <- cbind(pred1, pred2[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 4: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # rank the distances 
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  DRNN10_sets <- lapply(1:M, DRNN10, samp)
  
  ## Impute using DRNN01 (incorrect outcome model)
  DRNN01 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=3){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2a: Using the bootstrap sample, fit a multinomial regression model on the dependent variable
    # using the auxiliary data and save the standardized predicted values.
    mod1 <- multinom(y ~ x1 + x2, boot)
    pred1 <- fitted(mod1)[,-1]
    # We are manually standardizing because we need to save the means and sds for later.
    mu1 <- colMeans(pred1)
    sd1 <- apply(pred1, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1)){
      pred1[,i] <- (pred1[,i] - mu1[i]) / sd1[i]
    }
    # step 2b: Using the bootstrap sample, fit a logistic regression model on the response indicators 
    # using the auxiliary data and save the standardized predicted values.
    mod2 <- glm(response_indicators ~ x1 + x2 + x3, boot, family = binomial)
    pred2 <- mod2$fitted.values
    mu2 <- mean(pred2); sd2 <- sqrt(var(pred2)) 
    pred2 <- (pred2 - mu2) / sd2
    # Note that step 2a only uses respondents and step 2b uses the entire bootstrap sample.
    
    # step 3: calculate a distance function to define the similarity between subject i with missing Y
    # in the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    pred1.2 <- predict(mod1, newdata = to_be_imputed, type = "probs")[,-1]
    pred2.2 <- predict(mod2, newdata = to_be_imputed, type = "response")
    # Standardize the predicted values
    for(i in 1:ncol(pred1.2)){
      pred1.2[,i] <- (pred1.2[,i] - mu1[i]) / sd1[i]
    }
    pred2.2 <- (pred2.2 - mu2) / sd2
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(pred1.2, pred2.2) 
    pred.boot_obs <- cbind(pred1, pred2[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 4: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # rank the distances 
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  DRNN01_sets <- lapply(1:M, DRNN01, samp)
  
  ## Impute using DRNN00 (both incorrect)
  DRNN00 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=3){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2a: Using the bootstrap sample, fit a multinomial regression model on the dependent variable
    # using the auxiliary data and save the standardized predicted values.
    mod1 <- multinom(y ~ x1 + x2, boot)
    pred1 <- fitted(mod1)[,-1]
    # We are manually standardizing because we need to save the means and sds for later.
    mu1 <- colMeans(pred1)
    sd1 <- apply(pred1, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1)){
      pred1[,i] <- (pred1[,i] - mu1[i]) / sd1[i]
    }
    # step 2b: Using the bootstrap sample, fit a logistic regression model on the response indicators 
    # using the auxiliary data and save the standardized predicted values.
    mod2 <- glm(response_indicators ~ x1 + x2, boot, family = binomial)
    pred2 <- mod2$fitted.values
    mu2 <- mean(pred2); sd2 <- sqrt(var(pred2)) 
    pred2 <- (pred2 - mu2) / sd2
    # Note that step 2a only uses respondents and step 2b uses the entire bootstrap sample.
    
    # step 3: calculate a distance function to define the similarity between subject i with missing Y
    # in the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    pred1.2 <- predict(mod1, newdata = to_be_imputed, type = "probs")[,-1]
    pred2.2 <- predict(mod2, newdata = to_be_imputed, type = "response")
    # Standardize the predicted values
    for(i in 1:ncol(pred1.2)){
      pred1.2[,i] <- (pred1.2[,i] - mu1[i]) / sd1[i]
    }
    pred2.2 <- (pred2.2 - mu2) / sd2
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(pred1.2, pred2.2) 
    pred.boot_obs <- cbind(pred1, pred2[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 4: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # rank the distances 
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  DRNN00_sets <- lapply(1:M, DRNN00, samp)
  
  ## Impute using MICE0 (incorrect model)
  MICE0_sets <- mice(samp[,-4], m = M, maxit = 20)
  
  
  ## EVALUATION
  
  # MRNNI
  # Estimated proportions per data set
  # One column represents one data set
  MRNN10_proportions <- sapply(MRNN10_sets, function(x) table(x$y)/300)
  MRNN01_proportions <- sapply(MRNN01_sets, function(x) table(x$y)/300)
  MRNN00_proportions <- sapply(MRNN00_sets, function(x) table(x$y)/300)
  # Overall proportions:
  MRNN10_avg_props <- rowMeans(MRNN10_proportions)
  MRNN01_avg_props <- rowMeans(MRNN01_proportions)
  MRNN00_avg_props <- rowMeans(MRNN00_proportions)
  # Apply SE function to obtain SEs for every outcome category.
  SE_MRNN10 <- apply(MRNN10_proportions, 1, standard_error, M)
  SE_MRNN01 <- apply(MRNN01_proportions, 1, standard_error, M)
  SE_MRNN00 <- apply(MRNN00_proportions, 1, standard_error, M)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_MRNN10 <- matrix(CI(MRNN10_avg_props, SE_MRNN10), nrow=3)
  CIs_MRNN01 <- matrix(CI(MRNN01_avg_props, SE_MRNN01), nrow=3)
  CIs_MRNN00 <- matrix(CI(MRNN00_avg_props, SE_MRNN00), nrow=3)
  # Then check if the true population proportions are within the CIs (TRUE/FALSE).
  original_proportions <- table(population$y)/nrow(population)
  within_MRNN10 <- CIs_MRNN10[,1]<=original_proportions & original_proportions<=CIs_MRNN10[,2]
  within_MRNN01 <- CIs_MRNN01[,1]<=original_proportions & original_proportions<=CIs_MRNN01[,2]
  within_MRNN00 <- CIs_MRNN00[,1]<=original_proportions & original_proportions<=CIs_MRNN00[,2]
  
  # DRNNI
  # Estimated proportions per data set
  # One column represents one data set
  DRNN10_proportions <- sapply(DRNN10_sets, function(x) table(x$y)/300)
  DRNN01_proportions <- sapply(DRNN01_sets, function(x) table(x$y)/300)
  DRNN00_proportions <- sapply(DRNN00_sets, function(x) table(x$y)/300)
  # Overall proportions:
  DRNN10_avg_props <- rowMeans(DRNN10_proportions)
  DRNN01_avg_props <- rowMeans(DRNN01_proportions)
  DRNN00_avg_props <- rowMeans(DRNN00_proportions)
  # SE using Rubin's rules:
  # Apply function to obtain SEs for every outcome category.
  SE_DRNN10 <- apply(DRNN10_proportions, 1, standard_error, M)
  SE_DRNN01 <- apply(DRNN01_proportions, 1, standard_error, M)
  SE_DRNN00 <- apply(DRNN00_proportions, 1, standard_error, M)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_DRNN10 <- matrix(CI(DRNN10_avg_props, SE_DRNN10), nrow=3)
  CIs_DRNN01 <- matrix(CI(DRNN01_avg_props, SE_DRNN01), nrow=3)
  CIs_DRNN00 <- matrix(CI(DRNN00_avg_props, SE_DRNN00), nrow=3)
  # Then check if the original proportions are within the CIs (TRUE/FALSE).
  within_DRNN10 <- CIs_DRNN10[,1]<=original_proportions & original_proportions<=CIs_DRNN10[,2]
  within_DRNN01 <- CIs_DRNN01[,1]<=original_proportions & original_proportions<=CIs_DRNN01[,2]
  within_DRNN00 <- CIs_DRNN00[,1]<=original_proportions & original_proportions<=CIs_DRNN00[,2]
  
  # MICE
  # Estimated proportions per data set
  # One column represents one data set
  MICE0_proportions <- with(MICE0_sets, table(y)/300)$analyses
  MICE0_proportions <- matrix(unlist(MICE0_proportions), nrow = 3, ncol = M)
  # Overall proportions:
  MICE0_avg_props <- rowMeans(MICE0_proportions)
  # SE using Rubin's rules:
  # Apply function to obtain SEs for every outcome category.
  SE_MICE0 <- apply(MICE0_proportions, 1, standard_error, M)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_MICE0 <- matrix(CI(MICE0_avg_props, SE_MICE0), nrow=3)
  # Then check if the original proportions are within the CIs (TRUE/FALSE).
  within_MICE0 <- CIs_MICE0[,1]<=original_proportions & original_proportions<=CIs_MICE0[,2]
  
  # Output:
  MRNNI_out <- list(MRNN10_avg_props, MRNN01_avg_props, MRNN00_avg_props,
                    SE_MRNN10, SE_MRNN01, SE_MRNN00,
                    within_MRNN10, within_MRNN01, within_MRNN00)
  DRNNI_out <- list(DRNN10_avg_props, DRNN01_avg_props, DRNN00_avg_props,
                    SE_DRNN10, SE_DRNN01, SE_DRNN00,
                    within_DRNN10, within_DRNN01, within_DRNN00)
  MICE_out <- list(MICE0_avg_props, SE_MICE0, within_MICE0)
  
  return(list(MRNNI_out, DRNNI_out, MICE_out))
}

# Repeat 500 times to create, impute, and analyze 500 samples.
# These data need to be saved and exported.
set.seed(1254)
population1_response80 <- lapply(1:500, results_per_samp_first, population = population)
save(population1_response80, file = "S2pop1res80.RData")
set.seed(1255)
population1_response90 <- lapply(1:500, results_per_samp_first, population = population, 
                                 responserate80 = FALSE)
save(population1_response90, file = "S2pop1res90.RData")

#####################################################################################################
#####################################################################################################
##### POPULATION 2: 3 CATEGORIES FOR Y AND SECOND ORDER TERMS #######################################

load("~/SSLBS/Thesis/Thesis/pop2.RData")

### SAMPLING, MISSINGNESS, IMPUTATION, AND ANALYSIS
# Here we first create our samples and introduce missing values, then imputation and evaluation take 
# place. Creating a function allows us to use mc.lapply later on. If responserate80 = TRUE, the 
# response rate is 80%, if responserate80 = FALSE, it is 90%. Including the response rate as an 
# option in the function saves code (i.e. is more efficient).
# Note that the argument 'iterations' is not used within the code of the function, it is just there 
# so we can use it in an lapply or mclapply loop later on.

results_per_samp_second <- function(iterations, population, responserate80 = TRUE){
  # Sample from the population.
  samp <- population[sample(1:50000, 300, replace = FALSE),]
  # Generate response indicators depending on responserate80.
  if(responserate80){
    p_res80 <- exp(samp$x1 + samp$x2 + 0.01*samp$x3) / 
                            (exp(samp$x1 + samp$x2 + 0.01*samp$x3) + 1)
    samp$response_indicators <- rbinom(300, 1, p_res80)}else{
      p_res90 <- exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) / 
                              (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) + 1)
      samp$response_indicators <- rbinom(300, 1, p_res90)
    }
  # Introduce missing values accordingly.
  samp$y[samp$response_indicators == 0] <- NA
  
  # Draw a new sample if the sample does not contain observations from all categories.
  while(length(unique(samp$y))!=4){
    # Sample from the population.
    samp <- population[sample(1:50000, 300, replace = FALSE),]
    # Generate response indicators depending on responserate80.
    if(responserate80){
      p_res80 <- exp(samp$x1 + samp$x2 + 0.01*samp$x3) / 
                              (exp(samp$x1 + samp$x2 + 0.01*samp$x3) + 1)
      samp$response_indicators <- rbinom(300, 1, p_res80)}else{
        p_res90 <- exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) / 
                                (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) + 1)
        samp$response_indicators <- rbinom(300, 1, p_res90)
      }
    # Introduce missing values accordingly.
    samp$y[samp$response_indicators == 0] <- NA
  }
  
  ## Impute using MRNN10 (incorrect missingness model)
  MRNN10 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=3){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2: Using the bootstrap sample, fit one or more multinomial regression models on the 
    # dependent variable using the auxiliary data and save the predicted values.
    mod1.1 <- multinom(y ~ x1 + x2 + x3, data = boot)
    pred1.1 <- fitted(mod1.1)[,-1]
    mod1.2 <- multinom(y ~ x2^2 + x1:x3 + x1 + x2 + x3, data = boot)
    pred1.2 <- fitted(mod1.2)[,-1]
    # If two or more regression models were applied, regress the dependent variable on the predicted
    # values and save the new predicted values. 
    temporary <- cbind(boot[boot$response_indicators == 1, "y"], pred1.1, pred1.2)
    temporary <- as.data.frame(temporary)
    names(temporary) <- c("y", "M1C2", "M1C3", "M2C2", "M2C3")
    mod1.3 <- multinom(y ~ M1C2 + M1C3 + M2C2 + M2C3, temporary)
    pred1.3 <- fitted(mod1.3)[,-1]
    # Note that these models are applied only to respondents.
    # Standardize the predicted values.
    mu1 <- colMeans(pred1.3) # we need to save these for later, that's why we do it by hand 
    sd1 <- apply(pred1.3, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1.3)){
      pred1.3[,i] <- (pred1.3[,i] - mu1[i]) / sd1[i]
    }
    
    # step 3: Using the bootstrap sample, fit one or more logistic regression models on the response 
    # indicators using the auxiliary data and save the predicted values.
    mod2.1 <- glm(response_indicators ~ x1 + x2, boot, family = binomial)
    pred2.1 <- mod2.1$fitted.values
    mod2.2 <- glm(response_indicators ~ x1*x3, boot, family = binomial)
    pred2.2 <- mod2.2$fitted.values
    # If two or more regression models were applied, regress the response indicators on the predicted 
    # values and save the new predicted values.
    temporary <- cbind(boot$response_indicators, pred2.1, pred2.2)
    temporary <- as.data.frame(temporary)
    mod2.3 <- glm(V1 ~ pred2.1 + pred2.2, temporary, family = binomial)
    pred2.3 <- mod2.3$fitted.values
    # Note that these models are applied to all cases in the bootstrap sample.
    # Standardize the predicted values.
    mu2 <- mean(pred2.3); sd2 <- sqrt(var(pred2.3)) 
    pred2.3 <- (pred2.3 - mu2) / sd2
    
    # step 4: Calculate the predicted outcome probabilities and propensity scores for the
    # nonrespondents in the original sample and standardize these using the parameters obtained in
    # steps 2 and 3.
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    # Predicted outcome probabilities:
    missing_pred1.1 <- predict(mod1.1, newdata = to_be_imputed, type = "probs")[,-1]
    missing_pred1.2 <- predict(mod1.2, newdata = to_be_imputed, type = "probs")[,-1]
    missing_temp <- as.data.frame(cbind(missing_pred1.1, missing_pred1.2))
    names(missing_temp) <- c("M1C2", "M1C3", "M2C2", "M2C3")
    missing_pred1.3 <- predict(mod1.3, newdata = missing_temp, type = "probs")[,-1]
    # Predicted response propensities:
    missing_pred2.1 <- predict(mod2.1, newdata = to_be_imputed, type = "response")
    missing_pred2.2 <- predict(mod2.2, newdata = to_be_imputed, type = "response")
    missing_temp <- as.data.frame(cbind(missing_pred2.1, missing_pred2.2))
    names(missing_temp) <- c("pred2.1", "pred2.2")
    missing_pred2.3 <- predict(mod2.3, newdata = missing_temp, type = "response")
    # Standardize the predicted values.
    for(i in 1:ncol(missing_pred1.3)){
      missing_pred1.3[,i] <- (missing_pred1.3[,i] - mu1[i]) / sd1[i]
    }
    missing_pred2.3 <- (missing_pred2.3 - mu2) / sd2
    
    # step 5: calculate a distance function to define the similarity between subject i with missing Y 
    # in  the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(missing_pred1.3, missing_pred2.3) 
    pred.boot_obs <- cbind(pred1.3, pred2.3[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 6: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # Rank the distances.
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set.
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute.
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  MRNN10_sets <- lapply(1:M, MRNN10, samp)
  
  ## Impute using MRNN01 (incorrect outcome model)
  MRNN01 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=3){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2: Using the bootstrap sample, fit one or more multinomial regression models on the 
    # dependent variable using the auxiliary data and save the predicted values.
    mod1.1 <- multinom(y ~ x1 + x2, data = boot)
    pred1.1 <- fitted(mod1.1)[,-1]
    mod1.2 <- multinom(y ~ x2^2 + x1 + x2, data = boot)
    pred1.2 <- fitted(mod1.2)[,-1]
    # If two or more regression models were applied, regress the dependent variable on the predicted
    # values and save the new predicted values. 
    temporary <- cbind(boot[boot$response_indicators == 1, "y"], pred1.1, pred1.2)
    temporary <- as.data.frame(temporary)
    names(temporary) <- c("y", "M1C2", "M1C3", "M2C2", "M2C3")
    mod1.3 <- multinom(y ~ M1C2 + M1C3 + M2C2 + M2C3, temporary)
    pred1.3 <- fitted(mod1.3)[,-1]
    # Note that these models are applied only to respondents.
    # Standardize the predicted values.
    mu1 <- colMeans(pred1.3) # we need to save these for later, that's why we do it by hand 
    sd1 <- apply(pred1.3, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1.3)){
      pred1.3[,i] <- (pred1.3[,i] - mu1[i]) / sd1[i]
    }
    
    # step 3: Using the bootstrap sample, fit one or more logistic regression models on the response 
    # indicators using the auxiliary data and save the predicted values.
    mod2.1 <- glm(response_indicators ~ x1 + x2 + x3, boot, family = binomial)
    pred2.1 <- mod2.1$fitted.values
    mod2.2 <- glm(response_indicators ~ x1*x3, boot, family = binomial)
    pred2.2 <- mod2.2$fitted.values
    # If two or more regression models were applied, regress the response indicators on the predicted 
    # values and save the new predicted values.
    temporary <- cbind(boot$response_indicators, pred2.1, pred2.2)
    temporary <- as.data.frame(temporary)
    mod2.3 <- glm(V1 ~ pred2.1 + pred2.2, temporary, family = binomial)
    pred2.3 <- mod2.3$fitted.values
    # Note that these models are applied to all cases in the bootstrap sample.
    # Standardize the predicted values.
    mu2 <- mean(pred2.3); sd2 <- sqrt(var(pred2.3)) 
    pred2.3 <- (pred2.3 - mu2) / sd2
    
    # step 4: Calculate the predicted outcome probabilities and propensity scores for the
    # nonrespondents in the original sample and standardize these using the parameters obtained in
    # steps 2 and 3.
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    # Predicted outcome probabilities:
    missing_pred1.1 <- predict(mod1.1, newdata = to_be_imputed, type = "probs")[,-1]
    missing_pred1.2 <- predict(mod1.2, newdata = to_be_imputed, type = "probs")[,-1]
    missing_temp <- as.data.frame(cbind(missing_pred1.1, missing_pred1.2))
    names(missing_temp) <- c("M1C2", "M1C3", "M2C2", "M2C3")
    missing_pred1.3 <- predict(mod1.3, newdata = missing_temp, type = "probs")[,-1]
    # Predicted response propensities:
    missing_pred2.1 <- predict(mod2.1, newdata = to_be_imputed, type = "response")
    missing_pred2.2 <- predict(mod2.2, newdata = to_be_imputed, type = "response")
    missing_temp <- as.data.frame(cbind(missing_pred2.1, missing_pred2.2))
    names(missing_temp) <- c("pred2.1", "pred2.2")
    missing_pred2.3 <- predict(mod2.3, newdata = missing_temp, type = "response")
    # Standardize the predicted values.
    for(i in 1:ncol(missing_pred1.3)){
      missing_pred1.3[,i] <- (missing_pred1.3[,i] - mu1[i]) / sd1[i]
    }
    missing_pred2.3 <- (missing_pred2.3 - mu2) / sd2
    
    # step 5: calculate a distance function to define the similarity between subject i with missing Y 
    # in  the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(missing_pred1.3, missing_pred2.3) 
    pred.boot_obs <- cbind(pred1.3, pred2.3[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 6: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # Rank the distances.
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set.
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute.
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  MRNN01_sets <- lapply(1:M, MRNN01, samp)
  
  ## Impute using MRNN00 (both incorrect)
  MRNN00 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=3){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2: Using the bootstrap sample, fit one or more multinomial regression models on the 
    # dependent variable using the auxiliary data and save the predicted values.
    mod1.1 <- multinom(y ~ x1 + x2, data = boot)
    pred1.1 <- fitted(mod1.1)[,-1]
    mod1.2 <- multinom(y ~ x2^2 + x1 + x2, data = boot)
    pred1.2 <- fitted(mod1.2)[,-1]
    # If two or more regression models were applied, regress the dependent variable on the predicted
    # values and save the new predicted values. 
    temporary <- cbind(boot[boot$response_indicators == 1, "y"], pred1.1, pred1.2)
    temporary <- as.data.frame(temporary)
    names(temporary) <- c("y", "M1C2", "M1C3", "M2C2", "M2C3")
    mod1.3 <- multinom(y ~ M1C2 + M1C3 + M2C2 + M2C3, temporary)
    pred1.3 <- fitted(mod1.3)[,-1]
    # Note that these models are applied only to respondents.
    # Standardize the predicted values.
    mu1 <- colMeans(pred1.3) # we need to save these for later, that's why we do it by hand 
    sd1 <- apply(pred1.3, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1.3)){
      pred1.3[,i] <- (pred1.3[,i] - mu1[i]) / sd1[i]
    }
    
    # step 3: Using the bootstrap sample, fit one or more logistic regression models on the response 
    # indicators using the auxiliary data and save the predicted values.
    mod2.1 <- glm(response_indicators ~ x1 + x2, boot, family = binomial)
    pred2.1 <- mod2.1$fitted.values
    mod2.2 <- glm(response_indicators ~ x1*x3, boot, family = binomial)
    pred2.2 <- mod2.2$fitted.values
    # If two or more regression models were applied, regress the response indicators on the predicted 
    # values and save the new predicted values.
    temporary <- cbind(boot$response_indicators, pred2.1, pred2.2)
    temporary <- as.data.frame(temporary)
    mod2.3 <- glm(V1 ~ pred2.1 + pred2.2, temporary, family = binomial)
    pred2.3 <- mod2.3$fitted.values
    # Note that these models are applied to all cases in the bootstrap sample.
    # Standardize the predicted values.
    mu2 <- mean(pred2.3); sd2 <- sqrt(var(pred2.3)) 
    pred2.3 <- (pred2.3 - mu2) / sd2
    
    # step 4: Calculate the predicted outcome probabilities and propensity scores for the
    # nonrespondents in the original sample and standardize these using the parameters obtained in
    # steps 2 and 3.
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    # Predicted outcome probabilities:
    missing_pred1.1 <- predict(mod1.1, newdata = to_be_imputed, type = "probs")[,-1]
    missing_pred1.2 <- predict(mod1.2, newdata = to_be_imputed, type = "probs")[,-1]
    missing_temp <- as.data.frame(cbind(missing_pred1.1, missing_pred1.2))
    names(missing_temp) <- c("M1C2", "M1C3", "M2C2", "M2C3")
    missing_pred1.3 <- predict(mod1.3, newdata = missing_temp, type = "probs")[,-1]
    # Predicted response propensities:
    missing_pred2.1 <- predict(mod2.1, newdata = to_be_imputed, type = "response")
    missing_pred2.2 <- predict(mod2.2, newdata = to_be_imputed, type = "response")
    missing_temp <- as.data.frame(cbind(missing_pred2.1, missing_pred2.2))
    names(missing_temp) <- c("pred2.1", "pred2.2")
    missing_pred2.3 <- predict(mod2.3, newdata = missing_temp, type = "response")
    # Standardize the predicted values.
    for(i in 1:ncol(missing_pred1.3)){
      missing_pred1.3[,i] <- (missing_pred1.3[,i] - mu1[i]) / sd1[i]
    }
    missing_pred2.3 <- (missing_pred2.3 - mu2) / sd2
    
    # step 5: calculate a distance function to define the similarity between subject i with missing Y 
    # in  the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(missing_pred1.3, missing_pred2.3) 
    pred.boot_obs <- cbind(pred1.3, pred2.3[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 6: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # Rank the distances.
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set.
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute.
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  MRNN00_sets <- lapply(1:M, MRNN00, samp)
  
  ## Impute using DRNN10 (incorrect missingness model)
  DRNN10 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=3){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2a: Using the bootstrap sample, fit a multinomial regression model on the dependent variable
    # using the auxiliary data and save the standardized predicted values.
    mod1 <- multinom(y ~ x2^2 + x1:x3 + x1 + x2 + x3, boot)
    pred1 <- fitted(mod1)[,-1]
    # We are manually standardizing because we need to save the means and sds for later.
    mu1 <- colMeans(pred1)
    sd1 <- apply(pred1, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1)){
      pred1[,i] <- (pred1[,i] - mu1[i]) / sd1[i]
    }
    # step 2b: Using the bootstrap sample, fit a logistic regression model on the response indicators 
    # using the auxiliary data and save the standardized predicted values.
    mod2 <- glm(response_indicators ~ x1 + x2, boot, family = binomial)
    pred2 <- mod2$fitted.values
    mu2 <- mean(pred2); sd2 <- sqrt(var(pred2)) 
    pred2 <- (pred2 - mu2) / sd2
    # Note that step 2a only uses respondents and step 2b uses the entire bootstrap sample.
    
    # step 3: calculate a distance function to define the similarity between subject i with missing Y
    # in the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    pred1.2 <- predict(mod1, newdata = to_be_imputed, type = "probs")[,-1]
    pred2.2 <- predict(mod2, newdata = to_be_imputed, type = "response")
    # Standardize the predicted values
    for(i in 1:ncol(pred1.2)){
      pred1.2[,i] <- (pred1.2[,i] - mu1[i]) / sd1[i]
    }
    pred2.2 <- (pred2.2 - mu2) / sd2
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(pred1.2, pred2.2) 
    pred.boot_obs <- cbind(pred1, pred2[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 4: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # rank the distances 
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  DRNN10_sets <- lapply(1:M, DRNN10, samp)
  
  ## Impute using DRNN01 (incorrect outcome model)
  DRNN01 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=3){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2a: Using the bootstrap sample, fit a multinomial regression model on the dependent variable
    # using the auxiliary data and save the standardized predicted values.
    mod1 <- multinom(y ~ x2^2 + x1 + x2, boot)
    pred1 <- fitted(mod1)[,-1]
    # We are manually standardizing because we need to save the means and sds for later.
    mu1 <- colMeans(pred1)
    sd1 <- apply(pred1, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1)){
      pred1[,i] <- (pred1[,i] - mu1[i]) / sd1[i]
    }
    # step 2b: Using the bootstrap sample, fit a logistic regression model on the response indicators 
    # using the auxiliary data and save the standardized predicted values.
    mod2 <- glm(response_indicators ~ x1 + x2 + x3, boot, family = binomial)
    pred2 <- mod2$fitted.values
    mu2 <- mean(pred2); sd2 <- sqrt(var(pred2)) 
    pred2 <- (pred2 - mu2) / sd2
    # Note that step 2a only uses respondents and step 2b uses the entire bootstrap sample.
    
    # step 3: calculate a distance function to define the similarity between subject i with missing Y
    # in the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    pred1.2 <- predict(mod1, newdata = to_be_imputed, type = "probs")[,-1]
    pred2.2 <- predict(mod2, newdata = to_be_imputed, type = "response")
    # Standardize the predicted values
    for(i in 1:ncol(pred1.2)){
      pred1.2[,i] <- (pred1.2[,i] - mu1[i]) / sd1[i]
    }
    pred2.2 <- (pred2.2 - mu2) / sd2
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(pred1.2, pred2.2) 
    pred.boot_obs <- cbind(pred1, pred2[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 4: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # rank the distances 
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  DRNN01_sets <- lapply(1:M, DRNN01, samp)
  
  ## Impute using DRNN00 (both incorrect)
  DRNN00 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=3){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2a: Using the bootstrap sample, fit a multinomial regression model on the dependent variable
    # using the auxiliary data and save the standardized predicted values.
    mod1 <- multinom(y ~ x2^2 + x1 + x2, boot)
    pred1 <- fitted(mod1)[,-1]
    # We are manually standardizing because we need to save the means and sds for later.
    mu1 <- colMeans(pred1)
    sd1 <- apply(pred1, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1)){
      pred1[,i] <- (pred1[,i] - mu1[i]) / sd1[i]
    }
    # step 2b: Using the bootstrap sample, fit a logistic regression model on the response indicators 
    # using the auxiliary data and save the standardized predicted values.
    mod2 <- glm(response_indicators ~ x1 + x2, boot, family = binomial)
    pred2 <- mod2$fitted.values
    mu2 <- mean(pred2); sd2 <- sqrt(var(pred2)) 
    pred2 <- (pred2 - mu2) / sd2
    # Note that step 2a only uses respondents and step 2b uses the entire bootstrap sample.
    
    # step 3: calculate a distance function to define the similarity between subject i with missing Y
    # in the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    pred1.2 <- predict(mod1, newdata = to_be_imputed, type = "probs")[,-1]
    pred2.2 <- predict(mod2, newdata = to_be_imputed, type = "response")
    # Standardize the predicted values
    for(i in 1:ncol(pred1.2)){
      pred1.2[,i] <- (pred1.2[,i] - mu1[i]) / sd1[i]
    }
    pred2.2 <- (pred2.2 - mu2) / sd2
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(pred1.2, pred2.2) 
    pred.boot_obs <- cbind(pred1, pred2[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 4: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # rank the distances 
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  DRNN00_sets <- lapply(1:M, DRNN00, samp)
  
  ## Impute using MICE0 (incorrect model)
  MICE0_sets <- mice(samp[,-4], m = M, maxit = 20)
  
  
  ## EVALUATION
  
  # MRNNI
  # Estimated proportions per data set
  # One column represents one data set
  MRNN10_proportions <- sapply(MRNN10_sets, function(x) table(x$y)/300)
  MRNN01_proportions <- sapply(MRNN01_sets, function(x) table(x$y)/300)
  MRNN00_proportions <- sapply(MRNN00_sets, function(x) table(x$y)/300)
  # Overall proportions:
  MRNN10_avg_props <- rowMeans(MRNN10_proportions)
  MRNN01_avg_props <- rowMeans(MRNN01_proportions)
  MRNN00_avg_props <- rowMeans(MRNN00_proportions)
  # Apply SE function to obtain SEs for every outcome category.
  SE_MRNN10 <- apply(MRNN10_proportions, 1, standard_error, M)
  SE_MRNN01 <- apply(MRNN01_proportions, 1, standard_error, M)
  SE_MRNN00 <- apply(MRNN00_proportions, 1, standard_error, M)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_MRNN10 <- matrix(CI(MRNN10_avg_props, SE_MRNN10), nrow=3)
  CIs_MRNN01 <- matrix(CI(MRNN01_avg_props, SE_MRNN01), nrow=3)
  CIs_MRNN00 <- matrix(CI(MRNN00_avg_props, SE_MRNN00), nrow=3)
  # Then check if the true population proportions are within the CIs (TRUE/FALSE).
  original_proportions <- table(population$y)/nrow(population)
  within_MRNN10 <- CIs_MRNN10[,1]<=original_proportions & original_proportions<=CIs_MRNN10[,2]
  within_MRNN01 <- CIs_MRNN01[,1]<=original_proportions & original_proportions<=CIs_MRNN01[,2]
  within_MRNN00 <- CIs_MRNN00[,1]<=original_proportions & original_proportions<=CIs_MRNN00[,2]
  
  # DRNNI
  # Estimated proportions per data set
  # One column represents one data set
  DRNN10_proportions <- sapply(DRNN10_sets, function(x) table(x$y)/300)
  DRNN01_proportions <- sapply(DRNN01_sets, function(x) table(x$y)/300)
  DRNN00_proportions <- sapply(DRNN00_sets, function(x) table(x$y)/300)
  # Overall proportions:
  DRNN10_avg_props <- rowMeans(DRNN10_proportions)
  DRNN01_avg_props <- rowMeans(DRNN01_proportions)
  DRNN00_avg_props <- rowMeans(DRNN00_proportions)
  # SE using Rubin's rules:
  # Apply function to obtain SEs for every outcome category.
  SE_DRNN10 <- apply(DRNN10_proportions, 1, standard_error, M)
  SE_DRNN01 <- apply(DRNN01_proportions, 1, standard_error, M)
  SE_DRNN00 <- apply(DRNN00_proportions, 1, standard_error, M)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_DRNN10 <- matrix(CI(DRNN10_avg_props, SE_DRNN10), nrow=3)
  CIs_DRNN01 <- matrix(CI(DRNN01_avg_props, SE_DRNN01), nrow=3)
  CIs_DRNN00 <- matrix(CI(DRNN00_avg_props, SE_DRNN00), nrow=3)
  # Then check if the original proportions are within the CIs (TRUE/FALSE).
  within_DRNN10 <- CIs_DRNN10[,1]<=original_proportions & original_proportions<=CIs_DRNN10[,2]
  within_DRNN01 <- CIs_DRNN01[,1]<=original_proportions & original_proportions<=CIs_DRNN01[,2]
  within_DRNN00 <- CIs_DRNN00[,1]<=original_proportions & original_proportions<=CIs_DRNN00[,2]
  
  # MICE
  # Estimated proportions per data set
  # One column represents one data set
  MICE0_proportions <- with(MICE0_sets, table(y)/300)$analyses
  MICE0_proportions <- matrix(unlist(MICE0_proportions), nrow = 3, ncol = M)
  # Overall proportions:
  MICE0_avg_props <- rowMeans(MICE0_proportions)
  # SE using Rubin's rules:
  # Apply function to obtain SEs for every outcome category.
  SE_MICE0 <- apply(MICE0_proportions, 1, standard_error, M)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_MICE0 <- matrix(CI(MICE0_avg_props, SE_MICE0), nrow=3)
  # Then check if the original proportions are within the CIs (TRUE/FALSE).
  within_MICE0 <- CIs_MICE0[,1]<=original_proportions & original_proportions<=CIs_MICE0[,2]
  
  # Output:
  MRNNI_out <- list(MRNN10_avg_props, MRNN01_avg_props, MRNN00_avg_props,
                    SE_MRNN10, SE_MRNN01, SE_MRNN00,
                    within_MRNN10, within_MRNN01, within_MRNN00)
  DRNNI_out <- list(DRNN10_avg_props, DRNN01_avg_props, DRNN00_avg_props,
                    SE_DRNN10, SE_DRNN01, SE_DRNN00,
                    within_DRNN10, within_DRNN01, within_DRNN00)
  MICE_out <- list(MICE0_avg_props, SE_MICE0, within_MICE0)
  
  return(list(MRNNI_out, DRNNI_out, MICE_out))
}

# Repeat 500 times to create, impute, and analyze 500 samples.
# These data need to be saved and exported.
set.seed(1259)
population2_response80 <- lapply(1:500, results_per_samp_second, population = population)
save(population2_response80, file = "S2pop2res80.RData")
set.seed(1300)
population2_response90 <- lapply(1:500, results_per_samp_second, population = population, 
                                 responserate80 = FALSE)
save(population2_response90, file = "S2pop2res90.RData")

#####################################################################################################
#####################################################################################################
##### POPULATION 3: 5 CATEGORIES FOR Y AND FIRST ORDER TERMS ########################################

load("~/SSLBS/Thesis/Thesis/pop3.RData")

### SAMPLING, MISSINGNESS, IMPUTATION, AND ANALYSIS
# Here we first create our samples and introduce missing values, then imputation and evaluation take 
# place. Creating a function allows us to use mc.lapply later on. If responserate80 = TRUE, the 
# response rate is 80%, if responserate80 = FALSE, it is 90%. Including the response rate as an 
# option in the function saves code (i.e. is more efficient).
# Note that the argument 'iterations' is not used within the code of the function, it is just there 
# so we can use it in an lapply or mclapply loop later on.

results_per_samp_third <- function(iterations, population, responserate80 = TRUE){
  # Sample from the population.
  samp <- population[sample(1:50000, 300, replace = FALSE),]
  # Generate response indicators depending on responserate80.
  if(responserate80){
    p_res80 <- exp(samp$x1 + samp$x2 + 0.01*samp$x3) / 
                            (exp(samp$x1 + samp$x2 + 0.01*samp$x3) + 1)
    samp$response_indicators <- rbinom(300, 1, p_res80)}else{
      p_res90 <- exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) / 
                              (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) + 1)
      samp$response_indicators <- rbinom(300, 1, p_res90)
    }
  # Introduce missing values accordingly.
  samp$y[samp$response_indicators == 0] <- NA
  
  # Draw a new sample if the sample does not contain observations from all categories.
  while(length(unique(samp$y))!=6){
    # Sample from the population.
    samp <- population[sample(1:50000, 300, replace = FALSE),]
    # Generate response indicators depending on responserate80.
    if(responserate80){
      p_res80 <- exp(samp$x1 + samp$x2 + 0.01*samp$x3) / 
                              (exp(samp$x1 + samp$x2 + 0.01*samp$x3) + 1)
      samp$response_indicators <- rbinom(300, 1, p_res80)}else{
        p_res90 <- exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) / 
                                (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) + 1)
        samp$response_indicators <- rbinom(300, 1, p_res90)
      }
    # Introduce missing values accordingly.
    samp$y[samp$response_indicators == 0] <- NA
  }
  
  ## Impute using MRNN10 (incorrect missingness model)
  MRNN10 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=5){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2: Using the bootstrap sample, fit one or more multinomial regression models on the 
    # dependent variable using the auxiliary data and save the predicted values.
    mod1.1 <- multinom(y ~ x1 + x2 + x3, data = boot)
    pred1.1 <- fitted(mod1.1)[,-1]
    mod1.2 <- multinom(y ~ x2^2 + x1:x3 + x1 + x2 + x3, data = boot)
    pred1.2 <- fitted(mod1.2)[,-1]
    # If two or more regression models were applied, regress the dependent variable on the predicted
    # values and save the new predicted values. 
    temporary <- cbind(boot[boot$response_indicators == 1, "y"], pred1.1, pred1.2)
    temporary <- as.data.frame(temporary)
    names(temporary) <- c("y", "M1C2", "M1C3", "M1C4", "M1C5", "M2C2", "M2C3", "M2C4", "M2C5")
    mod1.3 <- multinom(y ~ M1C2 + M1C3 + M1C4 + M1C5 + M2C2 + M2C3 + M2C4 + M2C5, temporary)
    pred1.3 <- fitted(mod1.3)[,-1]
    # Note that these models are applied only to respondents.
    # Standardize the predicted values.
    mu1 <- colMeans(pred1.3) # we need to save these for later, that's why we do it by hand 
    sd1 <- apply(pred1.3, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1.3)){
      pred1.3[,i] <- (pred1.3[,i] - mu1[i]) / sd1[i]
    }
    
    # step 3: Using the bootstrap sample, fit one or more logistic regression models on the response 
    # indicators using the auxiliary data and save the predicted values.
    mod2.1 <- glm(response_indicators ~ x1 + x2, boot, family = binomial)
    pred2.1 <- mod2.1$fitted.values
    mod2.2 <- glm(response_indicators ~ x1*x3, boot, family = binomial)
    pred2.2 <- mod2.2$fitted.values
    # If two or more regression models were applied, regress the response indicators on the predicted 
    # values and save the new predicted values.
    temporary <- cbind(boot$response_indicators, pred2.1, pred2.2)
    temporary <- as.data.frame(temporary)
    mod2.3 <- glm(V1 ~ pred2.1 + pred2.2, temporary, family = binomial)
    pred2.3 <- mod2.3$fitted.values
    # Note that these models are applied to all cases in the bootstrap sample.
    # Standardize the predicted values.
    mu2 <- mean(pred2.3); sd2 <- sqrt(var(pred2.3)) 
    pred2.3 <- (pred2.3 - mu2) / sd2
    
    # step 4: Calculate the predicted outcome probabilities and propensity scores for the
    # nonrespondents in the original sample and standardize these using the parameters obtained in
    # steps 2 and 3.
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    # Predicted outcome probabilities:
    missing_pred1.1 <- predict(mod1.1, newdata = to_be_imputed, type = "probs")[,-1]
    missing_pred1.2 <- predict(mod1.2, newdata = to_be_imputed, type = "probs")[,-1]
    missing_temp <- as.data.frame(cbind(missing_pred1.1, missing_pred1.2))
    names(missing_temp) <- c("M1C2", "M1C3", "M1C4", "M1C5", "M2C2", "M2C3", "M2C4", "M2C5")
    missing_pred1.3 <- predict(mod1.3, newdata = missing_temp, type = "probs")[,-1]
    # Predicted response propensities:
    missing_pred2.1 <- predict(mod2.1, newdata = to_be_imputed, type = "response")
    missing_pred2.2 <- predict(mod2.2, newdata = to_be_imputed, type = "response")
    missing_temp <- as.data.frame(cbind(missing_pred2.1, missing_pred2.2))
    names(missing_temp) <- c("pred2.1", "pred2.2")
    missing_pred2.3 <- predict(mod2.3, newdata = missing_temp, type = "response")
    # Standardize the predicted values.
    for(i in 1:ncol(missing_pred1.3)){
      missing_pred1.3[,i] <- (missing_pred1.3[,i] - mu1[i]) / sd1[i]
    }
    missing_pred2.3 <- (missing_pred2.3 - mu2) / sd2
    
    # step 5: calculate a distance function to define the similarity between subject i with missing Y 
    # in  the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(missing_pred1.3, missing_pred2.3) 
    pred.boot_obs <- cbind(pred1.3, pred2.3[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 6: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # Rank the distances.
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set.
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute.
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  MRNN10_sets <- lapply(1:M, MRNN10, samp)
  
  ## Impute using MRNN01 (incorrect outcome model)
  MRNN01 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=5){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2: Using the bootstrap sample, fit one or more multinomial regression models on the 
    # dependent variable using the auxiliary data and save the predicted values.
    mod1.1 <- multinom(y ~ x1 + x2, data = boot)
    pred1.1 <- fitted(mod1.1)[,-1]
    mod1.2 <- multinom(y ~ x2^2 + x1 + x2, data = boot)
    pred1.2 <- fitted(mod1.2)[,-1]
    # If two or more regression models were applied, regress the dependent variable on the predicted
    # values and save the new predicted values. 
    temporary <- cbind(boot[boot$response_indicators == 1, "y"], pred1.1, pred1.2)
    temporary <- as.data.frame(temporary)
    names(temporary) <- c("y", "M1C2", "M1C3", "M1C4", "M1C5", "M2C2", "M2C3", "M2C4", "M2C5")
    mod1.3 <- multinom(y ~ M1C2 + M1C3 + M1C4 + M1C5 + M2C2 + M2C3 + M2C4 + M2C5, temporary)
    pred1.3 <- fitted(mod1.3)[,-1]
    # Note that these models are applied only to respondents.
    # Standardize the predicted values.
    mu1 <- colMeans(pred1.3) # we need to save these for later, that's why we do it by hand 
    sd1 <- apply(pred1.3, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1.3)){
      pred1.3[,i] <- (pred1.3[,i] - mu1[i]) / sd1[i]
    }
    
    # step 3: Using the bootstrap sample, fit one or more logistic regression models on the response 
    # indicators using the auxiliary data and save the predicted values.
    mod2.1 <- glm(response_indicators ~ x1 + x2 + x3, boot, family = binomial)
    pred2.1 <- mod2.1$fitted.values
    mod2.2 <- glm(response_indicators ~ x1*x3, boot, family = binomial)
    pred2.2 <- mod2.2$fitted.values
    # If two or more regression models were applied, regress the response indicators on the predicted 
    # values and save the new predicted values.
    temporary <- cbind(boot$response_indicators, pred2.1, pred2.2)
    temporary <- as.data.frame(temporary)
    mod2.3 <- glm(V1 ~ pred2.1 + pred2.2, temporary, family = binomial)
    pred2.3 <- mod2.3$fitted.values
    # Note that these models are applied to all cases in the bootstrap sample.
    # Standardize the predicted values.
    mu2 <- mean(pred2.3); sd2 <- sqrt(var(pred2.3)) 
    pred2.3 <- (pred2.3 - mu2) / sd2
    
    # step 4: Calculate the predicted outcome probabilities and propensity scores for the
    # nonrespondents in the original sample and standardize these using the parameters obtained in
    # steps 2 and 3.
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    # Predicted outcome probabilities:
    missing_pred1.1 <- predict(mod1.1, newdata = to_be_imputed, type = "probs")[,-1]
    missing_pred1.2 <- predict(mod1.2, newdata = to_be_imputed, type = "probs")[,-1]
    missing_temp <- as.data.frame(cbind(missing_pred1.1, missing_pred1.2))
    names(missing_temp) <- c("M1C2", "M1C3", "M1C4", "M1C5", "M2C2", "M2C3", "M2C4", "M2C5")
    missing_pred1.3 <- predict(mod1.3, newdata = missing_temp, type = "probs")[,-1]
    # Predicted response propensities:
    missing_pred2.1 <- predict(mod2.1, newdata = to_be_imputed, type = "response")
    missing_pred2.2 <- predict(mod2.2, newdata = to_be_imputed, type = "response")
    missing_temp <- as.data.frame(cbind(missing_pred2.1, missing_pred2.2))
    names(missing_temp) <- c("pred2.1", "pred2.2")
    missing_pred2.3 <- predict(mod2.3, newdata = missing_temp, type = "response")
    # Standardize the predicted values.
    for(i in 1:ncol(missing_pred1.3)){
      missing_pred1.3[,i] <- (missing_pred1.3[,i] - mu1[i]) / sd1[i]
    }
    missing_pred2.3 <- (missing_pred2.3 - mu2) / sd2
    
    # step 5: calculate a distance function to define the similarity between subject i with missing Y 
    # in  the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(missing_pred1.3, missing_pred2.3) 
    pred.boot_obs <- cbind(pred1.3, pred2.3[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 6: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # Rank the distances.
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set.
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute.
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  MRNN01_sets <- lapply(1:M, MRNN01, samp)
  
  ## Impute using MRNN10 (both incorrect)
  MRNN00 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=5){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2: Using the bootstrap sample, fit one or more multinomial regression models on the 
    # dependent variable using the auxiliary data and save the predicted values.
    mod1.1 <- multinom(y ~ x1 + x2, data = boot)
    pred1.1 <- fitted(mod1.1)[,-1]
    mod1.2 <- multinom(y ~ x2^2 + x1 + x2, data = boot)
    pred1.2 <- fitted(mod1.2)[,-1]
    # If two or more regression models were applied, regress the dependent variable on the predicted
    # values and save the new predicted values. 
    temporary <- cbind(boot[boot$response_indicators == 1, "y"], pred1.1, pred1.2)
    temporary <- as.data.frame(temporary)
    names(temporary) <- c("y", "M1C2", "M1C3", "M1C4", "M1C5", "M2C2", "M2C3", "M2C4", "M2C5")
    mod1.3 <- multinom(y ~ M1C2 + M1C3 + M1C4 + M1C5 + M2C2 + M2C3 + M2C4 + M2C5, temporary)
    pred1.3 <- fitted(mod1.3)[,-1]
    # Note that these models are applied only to respondents.
    # Standardize the predicted values.
    mu1 <- colMeans(pred1.3) # we need to save these for later, that's why we do it by hand 
    sd1 <- apply(pred1.3, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1.3)){
      pred1.3[,i] <- (pred1.3[,i] - mu1[i]) / sd1[i]
    }
    
    # step 3: Using the bootstrap sample, fit one or more logistic regression models on the response 
    # indicators using the auxiliary data and save the predicted values.
    mod2.1 <- glm(response_indicators ~ x1 + x2, boot, family = binomial)
    pred2.1 <- mod2.1$fitted.values
    mod2.2 <- glm(response_indicators ~ x1*x3, boot, family = binomial)
    pred2.2 <- mod2.2$fitted.values
    # If two or more regression models were applied, regress the response indicators on the predicted 
    # values and save the new predicted values.
    temporary <- cbind(boot$response_indicators, pred2.1, pred2.2)
    temporary <- as.data.frame(temporary)
    mod2.3 <- glm(V1 ~ pred2.1 + pred2.2, temporary, family = binomial)
    pred2.3 <- mod2.3$fitted.values
    # Note that these models are applied to all cases in the bootstrap sample.
    # Standardize the predicted values.
    mu2 <- mean(pred2.3); sd2 <- sqrt(var(pred2.3)) 
    pred2.3 <- (pred2.3 - mu2) / sd2
    
    # step 4: Calculate the predicted outcome probabilities and propensity scores for the
    # nonrespondents in the original sample and standardize these using the parameters obtained in
    # steps 2 and 3.
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    # Predicted outcome probabilities:
    missing_pred1.1 <- predict(mod1.1, newdata = to_be_imputed, type = "probs")[,-1]
    missing_pred1.2 <- predict(mod1.2, newdata = to_be_imputed, type = "probs")[,-1]
    missing_temp <- as.data.frame(cbind(missing_pred1.1, missing_pred1.2))
    names(missing_temp) <- c("M1C2", "M1C3", "M1C4", "M1C5", "M2C2", "M2C3", "M2C4", "M2C5")
    missing_pred1.3 <- predict(mod1.3, newdata = missing_temp, type = "probs")[,-1]
    # Predicted response propensities:
    missing_pred2.1 <- predict(mod2.1, newdata = to_be_imputed, type = "response")
    missing_pred2.2 <- predict(mod2.2, newdata = to_be_imputed, type = "response")
    missing_temp <- as.data.frame(cbind(missing_pred2.1, missing_pred2.2))
    names(missing_temp) <- c("pred2.1", "pred2.2")
    missing_pred2.3 <- predict(mod2.3, newdata = missing_temp, type = "response")
    # Standardize the predicted values.
    for(i in 1:ncol(missing_pred1.3)){
      missing_pred1.3[,i] <- (missing_pred1.3[,i] - mu1[i]) / sd1[i]
    }
    missing_pred2.3 <- (missing_pred2.3 - mu2) / sd2
    
    # step 5: calculate a distance function to define the similarity between subject i with missing Y 
    # in  the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(missing_pred1.3, missing_pred2.3) 
    pred.boot_obs <- cbind(pred1.3, pred2.3[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 6: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # Rank the distances.
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set.
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute.
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  MRNN00_sets <- lapply(1:M, MRNN00, samp)
  
  ## Impute using DRNN10 (incorrect missingness model)
  DRNN10 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=5){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2a: Using the bootstrap sample, fit a multinomial regression model on the dependent variable
    # using the auxiliary data and save the standardized predicted values.
    mod1 <- multinom(y ~ x1 + x2 + x3, boot)
    pred1 <- fitted(mod1)[,-1]
    # We are manually standardizing because we need to save the means and sds for later.
    mu1 <- colMeans(pred1)
    sd1 <- apply(pred1, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1)){
      pred1[,i] <- (pred1[,i] - mu1[i]) / sd1[i]
    }
    # step 2b: Using the bootstrap sample, fit a logistic regression model on the response indicators 
    # using the auxiliary data and save the standardized predicted values.
    mod2 <- glm(response_indicators ~ x1 + x2, boot, family = binomial)
    pred2 <- mod2$fitted.values
    mu2 <- mean(pred2); sd2 <- sqrt(var(pred2)) 
    pred2 <- (pred2 - mu2) / sd2
    # Note that step 2a only uses respondents and step 2b uses the entire bootstrap sample.
    
    # step 3: calculate a distance function to define the similarity between subject i with missing Y
    # in the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    pred1.2 <- predict(mod1, newdata = to_be_imputed, type = "probs")[,-1]
    pred2.2 <- predict(mod2, newdata = to_be_imputed, type = "response")
    # Standardize the predicted values
    for(i in 1:ncol(pred1.2)){
      pred1.2[,i] <- (pred1.2[,i] - mu1[i]) / sd1[i]
    }
    pred2.2 <- (pred2.2 - mu2) / sd2
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(pred1.2, pred2.2) 
    pred.boot_obs <- cbind(pred1, pred2[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 4: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # rank the distances 
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  DRNN10_sets <- lapply(1:5, DRNN10, samp)
  
  ## Impute using DRNN10 (incorrect outcome model)
  DRNN01 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=5){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2a: Using the bootstrap sample, fit a multinomial regression model on the dependent variable
    # using the auxiliary data and save the standardized predicted values.
    mod1 <- multinom(y ~ x1 + x2, boot)
    pred1 <- fitted(mod1)[,-1]
    # We are manually standardizing because we need to save the means and sds for later.
    mu1 <- colMeans(pred1)
    sd1 <- apply(pred1, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1)){
      pred1[,i] <- (pred1[,i] - mu1[i]) / sd1[i]
    }
    # step 2b: Using the bootstrap sample, fit a logistic regression model on the response indicators 
    # using the auxiliary data and save the standardized predicted values.
    mod2 <- glm(response_indicators ~ x1 + x2 + x3, boot, family = binomial)
    pred2 <- mod2$fitted.values
    mu2 <- mean(pred2); sd2 <- sqrt(var(pred2)) 
    pred2 <- (pred2 - mu2) / sd2
    # Note that step 2a only uses respondents and step 2b uses the entire bootstrap sample.
    
    # step 3: calculate a distance function to define the similarity between subject i with missing Y
    # in the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    pred1.2 <- predict(mod1, newdata = to_be_imputed, type = "probs")[,-1]
    pred2.2 <- predict(mod2, newdata = to_be_imputed, type = "response")
    # Standardize the predicted values
    for(i in 1:ncol(pred1.2)){
      pred1.2[,i] <- (pred1.2[,i] - mu1[i]) / sd1[i]
    }
    pred2.2 <- (pred2.2 - mu2) / sd2
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(pred1.2, pred2.2) 
    pred.boot_obs <- cbind(pred1, pred2[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 4: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # rank the distances 
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  DRNN01_sets <- lapply(1:5, DRNN01, samp)
  
  ## Impute using DRNN10 (both incorrect)
  DRNN00 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=5){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2a: Using the bootstrap sample, fit a multinomial regression model on the dependent variable
    # using the auxiliary data and save the standardized predicted values.
    mod1 <- multinom(y ~ x1 + x2, boot)
    pred1 <- fitted(mod1)[,-1]
    # We are manually standardizing because we need to save the means and sds for later.
    mu1 <- colMeans(pred1)
    sd1 <- apply(pred1, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1)){
      pred1[,i] <- (pred1[,i] - mu1[i]) / sd1[i]
    }
    # step 2b: Using the bootstrap sample, fit a logistic regression model on the response indicators 
    # using the auxiliary data and save the standardized predicted values.
    mod2 <- glm(response_indicators ~ x1 + x2, boot, family = binomial)
    pred2 <- mod2$fitted.values
    mu2 <- mean(pred2); sd2 <- sqrt(var(pred2)) 
    pred2 <- (pred2 - mu2) / sd2
    # Note that step 2a only uses respondents and step 2b uses the entire bootstrap sample.
    
    # step 3: calculate a distance function to define the similarity between subject i with missing Y
    # in the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    pred1.2 <- predict(mod1, newdata = to_be_imputed, type = "probs")[,-1]
    pred2.2 <- predict(mod2, newdata = to_be_imputed, type = "response")
    # Standardize the predicted values
    for(i in 1:ncol(pred1.2)){
      pred1.2[,i] <- (pred1.2[,i] - mu1[i]) / sd1[i]
    }
    pred2.2 <- (pred2.2 - mu2) / sd2
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(pred1.2, pred2.2) 
    pred.boot_obs <- cbind(pred1, pred2[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 4: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # rank the distances 
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  DRNN00_sets <- lapply(1:5, DRNN00, samp)
  
  ## Impute using MICE0 (incorrect model)
  MICE0_sets <- mice(samp[,-4], m = M, maxit = 20)
  
  
  ## EVALUATION
  
  
  # MRNNI
  # Estimated proportions per data set
  # One column represents one data set
  MRNN10_proportions <- sapply(MRNN10_sets, function(x) table(x$y)/300)
  MRNN01_proportions <- sapply(MRNN01_sets, function(x) table(x$y)/300)
  MRNN00_proportions <- sapply(MRNN00_sets, function(x) table(x$y)/300)
  # Overall proportions:
  MRNN10_avg_props <- rowMeans(MRNN10_proportions)
  MRNN01_avg_props <- rowMeans(MRNN01_proportions)
  MRNN00_avg_props <- rowMeans(MRNN00_proportions)
  # Apply SE function to obtain SEs for every outcome category.
  SE_MRNN10 <- apply(MRNN10_proportions, 1, standard_error, M)
  SE_MRNN01 <- apply(MRNN01_proportions, 1, standard_error, M)
  SE_MRNN00 <- apply(MRNN00_proportions, 1, standard_error, M)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_MRNN10 <- matrix(CI(MRNN10_avg_props, SE_MRNN10), nrow=5)
  CIs_MRNN01 <- matrix(CI(MRNN01_avg_props, SE_MRNN01), nrow=5)
  CIs_MRNN00 <- matrix(CI(MRNN00_avg_props, SE_MRNN00), nrow=5)
  # Then check if the true population proportions are within the CIs (TRUE/FALSE).
  original_proportions <- table(population$y)/nrow(population)
  within_MRNN10 <- CIs_MRNN10[,1]<=original_proportions & original_proportions<=CIs_MRNN10[,2]
  within_MRNN01 <- CIs_MRNN01[,1]<=original_proportions & original_proportions<=CIs_MRNN01[,2]
  within_MRNN00 <- CIs_MRNN00[,1]<=original_proportions & original_proportions<=CIs_MRNN00[,2]
  
  # DRNNI
  # Estimated proportions per data set
  # One column represents one data set
  DRNN10_proportions <- sapply(DRNN10_sets, function(x) table(x$y)/300)
  DRNN01_proportions <- sapply(DRNN01_sets, function(x) table(x$y)/300)
  DRNN00_proportions <- sapply(DRNN00_sets, function(x) table(x$y)/300)
  # Overall proportions:
  DRNN10_avg_props <- rowMeans(DRNN10_proportions)
  DRNN01_avg_props <- rowMeans(DRNN01_proportions)
  DRNN00_avg_props <- rowMeans(DRNN00_proportions)
  # SE using Rubin's rules:
  # Apply function to obtain SEs for every outcome category.
  SE_DRNN10 <- apply(DRNN10_proportions, 1, standard_error, M)
  SE_DRNN01 <- apply(DRNN01_proportions, 1, standard_error, M)
  SE_DRNN00 <- apply(DRNN00_proportions, 1, standard_error, M)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_DRNN10 <- matrix(CI(DRNN10_avg_props, SE_DRNN10), nrow=5)
  CIs_DRNN01 <- matrix(CI(DRNN01_avg_props, SE_DRNN01), nrow=5)
  CIs_DRNN00 <- matrix(CI(DRNN00_avg_props, SE_DRNN00), nrow=5)
  # Then check if the original proportions are within the CIs (TRUE/FALSE).
  within_DRNN10 <- CIs_DRNN10[,1]<=original_proportions & original_proportions<=CIs_DRNN10[,2]
  within_DRNN01 <- CIs_DRNN01[,1]<=original_proportions & original_proportions<=CIs_DRNN01[,2]
  within_DRNN00 <- CIs_DRNN00[,1]<=original_proportions & original_proportions<=CIs_DRNN00[,2]
  
  # MICE
  # Estimated proportions per data set
  # One column represents one data set
  MICE0_proportions <- with(MICE0_sets, table(y)/300)$analyses
  MICE0_proportions <- matrix(unlist(MICE0_proportions), nrow = 5, ncol = M)
  # Overall proportions:
  MICE0_avg_props <- rowMeans(MICE0_proportions)
  # SE using Rubin's rules:
  # Apply function to obtain SEs for every outcome category.
  SE_MICE0 <- apply(MICE0_proportions, 1, standard_error, M)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_MICE0 <- matrix(CI(MICE0_avg_props, SE_MICE0), nrow=5)
  # Then check if the original proportions are within the CIs (TRUE/FALSE).
  within_MICE0 <- CIs_MICE0[,1]<=original_proportions & original_proportions<=CIs_MICE0[,2]
  
  # Output:
  MRNNI_out <- list(MRNN10_avg_props, MRNN01_avg_props, MRNN00_avg_props,
                    SE_MRNN10, SE_MRNN01, SE_MRNN00,
                    within_MRNN10, within_MRNN01, within_MRNN00)
  DRNNI_out <- list(DRNN10_avg_props, DRNN01_avg_props, DRNN00_avg_props,
                    SE_DRNN10, SE_DRNN01, SE_DRNN00,
                    within_DRNN10, within_DRNN01, within_DRNN00)
  MICE_out <- list(MICE0_avg_props, SE_MICE0, within_MICE0)
  
  return(list(MRNNI_out, DRNNI_out, MICE_out))
}

# Repeat 500 times to create, impute, and analyze 500 samples.
# These data need to be saved and exported.
set.seed(1320)
population3_response80 <- lapply(1:500, results_per_samp_third, population = population)
save(population3_response80, file = "S2pop3res80.RData")
set.seed(1321)
population3_response90 <- lapply(1:500, results_per_samp_third, population = population, 
                                 responserate80 = FALSE)
save(population3_response90, file = "S2pop3res90.RData")

#####################################################################################################
#####################################################################################################
##### POPULATION 4: 5 CATEGORIES FOR 5 AND SECOND ORDER TERMS #######################################

load("~/SSLBS/Thesis/Thesis/pop4.RData")

### SAMPLING, MISSINGNESS, IMPUTATION, AND ANALYSIS
# Here we first create our samples and introduce missing values, then imputation and evaluation take 
# place. Creating a function allows us to use mc.lapply later on. If responserate80 = TRUE, the 
# response rate is 80%, if responserate80 = FALSE, it is 90%. Including the response rate as an 
# option in the function saves code (i.e. is more efficient).
# Note that the argument 'iterations' is not used within the code of the function, it is just there 
# so we can use it in an lapply or mclapply loop later on.

results_per_samp_fourth <- function(iterations, population, responserate80 = TRUE){
  # Sample from the population.
  samp <- population[sample(1:50000, 300, replace = FALSE),]
  # Generate response indicators depending on responserate80.
  if(responserate80){
    p_res80 <- exp(samp$x1 + samp$x2 + 0.01*samp$x3) / 
                            (exp(samp$x1 + samp$x2 + 0.01*samp$x3) + 1)
    samp$response_indicators <- rbinom(300, 1, p_res80)}else{
      p_res90 <- exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) / 
                              (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) + 1)
      samp$response_indicators <- rbinom(300, 1, p_res90)
    }
  # Introduce missing values accordingly.
  samp$y[samp$response_indicators == 0] <- NA
  
  # Draw a new sample if the sample does not contain observations from all categories.
  while(length(unique(samp$y))!=6){
    # Sample from the population.
    samp <- population[sample(1:50000, 300, replace = FALSE),]
    # Generate response indicators depending on responserate80.
    if(responserate80){
      p_res80 <-exp(samp$x1 + samp$x2 + 0.01*samp$x3) / 
                              (exp(samp$x1 + samp$x2 + 0.01*samp$x3) + 1)
      samp$response_indicators <- rbinom(300, 1, p_res80)}else{
        p_res90 <- exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) / 
                                (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) + 1)
        samp$response_indicators <- rbinom(300, 1, p_res90)
      }
    # Introduce missing values accordingly.
    samp$y[samp$response_indicators == 0] <- NA
  }
  
  ## Impute using MRNN10 (incorrect missingness model)
  MRNN10 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=5){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2: Using the bootstrap sample, fit one or more multinomial regression models on the 
    # dependent variable using the auxiliary data and save the predicted values.
    mod1.1 <- multinom(y ~ x1 + x2 + x3, data = boot)
    pred1.1 <- fitted(mod1.1)[,-1]
    mod1.2 <- multinom(y ~ x2^2 + x1:x3 + x1 + x2 + x3, data = boot)
    pred1.2 <- fitted(mod1.2)[,-1]
    # If two or more regression models were applied, regress the dependent variable on the predicted
    # values and save the new predicted values. 
    temporary <- cbind(boot[boot$response_indicators == 1, "y"], pred1.1, pred1.2)
    temporary <- as.data.frame(temporary)
    names(temporary) <- c("y", "M1C2", "M1C3", "M1C4", "M1C5", "M2C2", "M2C3", "M2C4", "M2C5")
    mod1.3 <- multinom(y ~ M1C2 + M1C3 + M1C4 + M1C5 + M2C2 + M2C3 + M2C4 + M2C5, temporary)
    pred1.3 <- fitted(mod1.3)[,-1]
    # Note that these models are applied only to respondents.
    # Standardize the predicted values.
    mu1 <- colMeans(pred1.3) # we need to save these for later, that's why we do it by hand 
    sd1 <- apply(pred1.3, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1.3)){
      pred1.3[,i] <- (pred1.3[,i] - mu1[i]) / sd1[i]
    }
    
    # step 3: Using the bootstrap sample, fit one or more logistic regression models on the response 
    # indicators using the auxiliary data and save the predicted values.
    mod2.1 <- glm(response_indicators ~ x1 + x2, boot, family = binomial)
    pred2.1 <- mod2.1$fitted.values
    mod2.2 <- glm(response_indicators ~ x1*x3, boot, family = binomial)
    pred2.2 <- mod2.2$fitted.values
    # If two or more regression models were applied, regress the response indicators on the predicted 
    # values and save the new predicted values.
    temporary <- cbind(boot$response_indicators, pred2.1, pred2.2)
    temporary <- as.data.frame(temporary)
    mod2.3 <- glm(V1 ~ pred2.1 + pred2.2, temporary, family = binomial)
    pred2.3 <- mod2.3$fitted.values
    # Note that these models are applied to all cases in the bootstrap sample.
    # Standardize the predicted values.
    mu2 <- mean(pred2.3); sd2 <- sqrt(var(pred2.3)) 
    pred2.3 <- (pred2.3 - mu2) / sd2
    
    # step 4: Calculate the predicted outcome probabilities and propensity scores for the
    # nonrespondents in the original sample and standardize these using the parameters obtained in
    # steps 2 and 3.
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    # Predicted outcome probabilities:
    missing_pred1.1 <- predict(mod1.1, newdata = to_be_imputed, type = "probs")[,-1]
    missing_pred1.2 <- predict(mod1.2, newdata = to_be_imputed, type = "probs")[,-1]
    missing_temp <- as.data.frame(cbind(missing_pred1.1, missing_pred1.2))
    names(missing_temp) <- c("M1C2", "M1C3", "M1C4", "M1C5", "M2C2", "M2C3", "M2C4", "M2C5")
    missing_pred1.3 <- predict(mod1.3, newdata = missing_temp, type = "probs")[,-1]
    # Predicted response propensities:
    missing_pred2.1 <- predict(mod2.1, newdata = to_be_imputed, type = "response")
    missing_pred2.2 <- predict(mod2.2, newdata = to_be_imputed, type = "response")
    missing_temp <- as.data.frame(cbind(missing_pred2.1, missing_pred2.2))
    names(missing_temp) <- c("pred2.1", "pred2.2")
    missing_pred2.3 <- predict(mod2.3, newdata = missing_temp, type = "response")
    # Standardize the predicted values.
    for(i in 1:ncol(missing_pred1.3)){
      missing_pred1.3[,i] <- (missing_pred1.3[,i] - mu1[i]) / sd1[i]
    }
    missing_pred2.3 <- (missing_pred2.3 - mu2) / sd2
    
    # step 5: calculate a distance function to define the similarity between subject i with missing Y 
    # in  the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(missing_pred1.3, missing_pred2.3) 
    pred.boot_obs <- cbind(pred1.3, pred2.3[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 6: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # Rank the distances.
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set.
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute.
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  MRNN10_sets <- lapply(1:M, MRNN10, samp)
  
  ## Impute using MRNN01 (incorrect outcome model)
  MRNN01 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=5){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2: Using the bootstrap sample, fit one or more multinomial regression models on the 
    # dependent variable using the auxiliary data and save the predicted values.
    mod1.1 <- multinom(y ~ x1 + x2, data = boot)
    pred1.1 <- fitted(mod1.1)[,-1]
    mod1.2 <- multinom(y ~ x2^2 + x1 + x2, data = boot)
    pred1.2 <- fitted(mod1.2)[,-1]
    # If two or more regression models were applied, regress the dependent variable on the predicted
    # values and save the new predicted values. 
    temporary <- cbind(boot[boot$response_indicators == 1, "y"], pred1.1, pred1.2)
    temporary <- as.data.frame(temporary)
    names(temporary) <- c("y", "M1C2", "M1C3", "M1C4", "M1C5", "M2C2", "M2C3", "M2C4", "M2C5")
    mod1.3 <- multinom(y ~ M1C2 + M1C3 + M1C4 + M1C5 + M2C2 + M2C3 + M2C4 + M2C5, temporary)
    pred1.3 <- fitted(mod1.3)[,-1]
    # Note that these models are applied only to respondents.
    # Standardize the predicted values.
    mu1 <- colMeans(pred1.3) # we need to save these for later, that's why we do it by hand 
    sd1 <- apply(pred1.3, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1.3)){
      pred1.3[,i] <- (pred1.3[,i] - mu1[i]) / sd1[i]
    }
    
    # step 3: Using the bootstrap sample, fit one or more logistic regression models on the response 
    # indicators using the auxiliary data and save the predicted values.
    mod2.1 <- glm(response_indicators ~ x1 + x2 + x3, boot, family = binomial)
    pred2.1 <- mod2.1$fitted.values
    mod2.2 <- glm(response_indicators ~ x1*x3, boot, family = binomial)
    pred2.2 <- mod2.2$fitted.values
    # If two or more regression models were applied, regress the response indicators on the predicted 
    # values and save the new predicted values.
    temporary <- cbind(boot$response_indicators, pred2.1, pred2.2)
    temporary <- as.data.frame(temporary)
    mod2.3 <- glm(V1 ~ pred2.1 + pred2.2, temporary, family = binomial)
    pred2.3 <- mod2.3$fitted.values
    # Note that these models are applied to all cases in the bootstrap sample.
    # Standardize the predicted values.
    mu2 <- mean(pred2.3); sd2 <- sqrt(var(pred2.3)) 
    pred2.3 <- (pred2.3 - mu2) / sd2
    
    # step 4: Calculate the predicted outcome probabilities and propensity scores for the
    # nonrespondents in the original sample and standardize these using the parameters obtained in
    # steps 2 and 3.
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    # Predicted outcome probabilities:
    missing_pred1.1 <- predict(mod1.1, newdata = to_be_imputed, type = "probs")[,-1]
    missing_pred1.2 <- predict(mod1.2, newdata = to_be_imputed, type = "probs")[,-1]
    missing_temp <- as.data.frame(cbind(missing_pred1.1, missing_pred1.2))
    names(missing_temp) <- c("M1C2", "M1C3", "M1C4", "M1C5", "M2C2", "M2C3", "M2C4", "M2C5")
    missing_pred1.3 <- predict(mod1.3, newdata = missing_temp, type = "probs")[,-1]
    # Predicted response propensities:
    missing_pred2.1 <- predict(mod2.1, newdata = to_be_imputed, type = "response")
    missing_pred2.2 <- predict(mod2.2, newdata = to_be_imputed, type = "response")
    missing_temp <- as.data.frame(cbind(missing_pred2.1, missing_pred2.2))
    names(missing_temp) <- c("pred2.1", "pred2.2")
    missing_pred2.3 <- predict(mod2.3, newdata = missing_temp, type = "response")
    # Standardize the predicted values.
    for(i in 1:ncol(missing_pred1.3)){
      missing_pred1.3[,i] <- (missing_pred1.3[,i] - mu1[i]) / sd1[i]
    }
    missing_pred2.3 <- (missing_pred2.3 - mu2) / sd2
    
    # step 5: calculate a distance function to define the similarity between subject i with missing Y 
    # in  the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(missing_pred1.3, missing_pred2.3) 
    pred.boot_obs <- cbind(pred1.3, pred2.3[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 6: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # Rank the distances.
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set.
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute.
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  MRNN01_sets <- lapply(1:M, MRNN01, samp)
  
  ## Impute using MRNN10 (both incorrect)
  MRNN00 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=5){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2: Using the bootstrap sample, fit one or more multinomial regression models on the 
    # dependent variable using the auxiliary data and save the predicted values.
    mod1.1 <- multinom(y ~ x1 + x2, data = boot)
    pred1.1 <- fitted(mod1.1)[,-1]
    mod1.2 <- multinom(y ~ x2^2 + x1 + x2, data = boot)
    pred1.2 <- fitted(mod1.2)[,-1]
    # If two or more regression models were applied, regress the dependent variable on the predicted
    # values and save the new predicted values. 
    temporary <- cbind(boot[boot$response_indicators == 1, "y"], pred1.1, pred1.2)
    temporary <- as.data.frame(temporary)
    names(temporary) <- c("y", "M1C2", "M1C3", "M1C4", "M1C5", "M2C2", "M2C3", "M2C4", "M2C5")
    mod1.3 <- multinom(y ~ M1C2 + M1C3 + M1C4 + M1C5 + M2C2 + M2C3 + M2C4 + M2C5, temporary)
    pred1.3 <- fitted(mod1.3)[,-1]
    # Note that these models are applied only to respondents.
    # Standardize the predicted values.
    mu1 <- colMeans(pred1.3) # we need to save these for later, that's why we do it by hand 
    sd1 <- apply(pred1.3, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1.3)){
      pred1.3[,i] <- (pred1.3[,i] - mu1[i]) / sd1[i]
    }
    
    # step 3: Using the bootstrap sample, fit one or more logistic regression models on the response 
    # indicators using the auxiliary data and save the predicted values.
    mod2.1 <- glm(response_indicators ~ x1 + x2, boot, family = binomial)
    pred2.1 <- mod2.1$fitted.values
    mod2.2 <- glm(response_indicators ~ x1*x3, boot, family = binomial)
    pred2.2 <- mod2.2$fitted.values
    # If two or more regression models were applied, regress the response indicators on the predicted 
    # values and save the new predicted values.
    temporary <- cbind(boot$response_indicators, pred2.1, pred2.2)
    temporary <- as.data.frame(temporary)
    mod2.3 <- glm(V1 ~ pred2.1 + pred2.2, temporary, family = binomial)
    pred2.3 <- mod2.3$fitted.values
    # Note that these models are applied to all cases in the bootstrap sample.
    # Standardize the predicted values.
    mu2 <- mean(pred2.3); sd2 <- sqrt(var(pred2.3)) 
    pred2.3 <- (pred2.3 - mu2) / sd2
    
    # step 4: Calculate the predicted outcome probabilities and propensity scores for the
    # nonrespondents in the original sample and standardize these using the parameters obtained in
    # steps 2 and 3.
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    # Predicted outcome probabilities:
    missing_pred1.1 <- predict(mod1.1, newdata = to_be_imputed, type = "probs")[,-1]
    missing_pred1.2 <- predict(mod1.2, newdata = to_be_imputed, type = "probs")[,-1]
    missing_temp <- as.data.frame(cbind(missing_pred1.1, missing_pred1.2))
    names(missing_temp) <- c("M1C2", "M1C3", "M1C4", "M1C5", "M2C2", "M2C3", "M2C4", "M2C5")
    missing_pred1.3 <- predict(mod1.3, newdata = missing_temp, type = "probs")[,-1]
    # Predicted response propensities:
    missing_pred2.1 <- predict(mod2.1, newdata = to_be_imputed, type = "response")
    missing_pred2.2 <- predict(mod2.2, newdata = to_be_imputed, type = "response")
    missing_temp <- as.data.frame(cbind(missing_pred2.1, missing_pred2.2))
    names(missing_temp) <- c("pred2.1", "pred2.2")
    missing_pred2.3 <- predict(mod2.3, newdata = missing_temp, type = "response")
    # Standardize the predicted values.
    for(i in 1:ncol(missing_pred1.3)){
      missing_pred1.3[,i] <- (missing_pred1.3[,i] - mu1[i]) / sd1[i]
    }
    missing_pred2.3 <- (missing_pred2.3 - mu2) / sd2
    
    # step 5: calculate a distance function to define the similarity between subject i with missing Y 
    # in  the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(missing_pred1.3, missing_pred2.3) 
    pred.boot_obs <- cbind(pred1.3, pred2.3[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 6: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # Rank the distances.
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set.
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute.
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  MRNN00_sets <- lapply(1:M, MRNN00, samp)
  
  ## Impute using DRNN10 (incorrect missingness model)
  DRNN10 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=5){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2a: Using the bootstrap sample, fit a multinomial regression model on the dependent variable
    # using the auxiliary data and save the standardized predicted values.
    mod1 <- multinom(y ~ x2^2 + x1:x3 + x1 + x2 + x3, boot)
    pred1 <- fitted(mod1)[,-1]
    # We are manually standardizing because we need to save the means and sds for later.
    mu1 <- colMeans(pred1)
    sd1 <- apply(pred1, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1)){
      pred1[,i] <- (pred1[,i] - mu1[i]) / sd1[i]
    }
    # step 2b: Using the bootstrap sample, fit a logistic regression model on the response indicators 
    # using the auxiliary data and save the standardized predicted values.
    mod2 <- glm(response_indicators ~ x1 + x2, boot, family = binomial)
    pred2 <- mod2$fitted.values
    mu2 <- mean(pred2); sd2 <- sqrt(var(pred2)) 
    pred2 <- (pred2 - mu2) / sd2
    # Note that step 2a only uses respondents and step 2b uses the entire bootstrap sample.
    
    # step 3: calculate a distance function to define the similarity between subject i with missing Y
    # in the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    pred1.2 <- predict(mod1, newdata = to_be_imputed, type = "probs")[,-1]
    pred2.2 <- predict(mod2, newdata = to_be_imputed, type = "response")
    # Standardize the predicted values
    for(i in 1:ncol(pred1.2)){
      pred1.2[,i] <- (pred1.2[,i] - mu1[i]) / sd1[i]
    }
    pred2.2 <- (pred2.2 - mu2) / sd2
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(pred1.2, pred2.2) 
    pred.boot_obs <- cbind(pred1, pred2[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 4: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # rank the distances 
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  DRNN10_sets <- lapply(1:5, DRNN10, samp)
  
  ## Impute using DRNN10 (incorrect outcome model)
  DRNN01 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=5){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2a: Using the bootstrap sample, fit a multinomial regression model on the dependent variable
    # using the auxiliary data and save the standardized predicted values.
    mod1 <- multinom(y ~ x2^2 + x1 + x2, boot)
    pred1 <- fitted(mod1)[,-1]
    # We are manually standardizing because we need to save the means and sds for later.
    mu1 <- colMeans(pred1)
    sd1 <- apply(pred1, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1)){
      pred1[,i] <- (pred1[,i] - mu1[i]) / sd1[i]
    }
    # step 2b: Using the bootstrap sample, fit a logistic regression model on the response indicators 
    # using the auxiliary data and save the standardized predicted values.
    mod2 <- glm(response_indicators ~ x1 + x2 + x3, boot, family = binomial)
    pred2 <- mod2$fitted.values
    mu2 <- mean(pred2); sd2 <- sqrt(var(pred2)) 
    pred2 <- (pred2 - mu2) / sd2
    # Note that step 2a only uses respondents and step 2b uses the entire bootstrap sample.
    
    # step 3: calculate a distance function to define the similarity between subject i with missing Y
    # in the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    pred1.2 <- predict(mod1, newdata = to_be_imputed, type = "probs")[,-1]
    pred2.2 <- predict(mod2, newdata = to_be_imputed, type = "response")
    # Standardize the predicted values
    for(i in 1:ncol(pred1.2)){
      pred1.2[,i] <- (pred1.2[,i] - mu1[i]) / sd1[i]
    }
    pred2.2 <- (pred2.2 - mu2) / sd2
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(pred1.2, pred2.2) 
    pred.boot_obs <- cbind(pred1, pred2[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 4: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # rank the distances 
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  DRNN01_sets <- lapply(1:5, DRNN01, samp)
  
  ## Impute using DRNN10 (both incorrect)
  DRNN00 <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=5){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2a: Using the bootstrap sample, fit a multinomial regression model on the dependent variable
    # using the auxiliary data and save the standardized predicted values.
    mod1 <- multinom(y ~ x2^2 + x1 + x2, boot)
    pred1 <- fitted(mod1)[,-1]
    # We are manually standardizing because we need to save the means and sds for later.
    mu1 <- colMeans(pred1)
    sd1 <- apply(pred1, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1)){
      pred1[,i] <- (pred1[,i] - mu1[i]) / sd1[i]
    }
    # step 2b: Using the bootstrap sample, fit a logistic regression model on the response indicators 
    # using the auxiliary data and save the standardized predicted values.
    mod2 <- glm(response_indicators ~ x1 + x2, boot, family = binomial)
    pred2 <- mod2$fitted.values
    mu2 <- mean(pred2); sd2 <- sqrt(var(pred2)) 
    pred2 <- (pred2 - mu2) / sd2
    # Note that step 2a only uses respondents and step 2b uses the entire bootstrap sample.
    
    # step 3: calculate a distance function to define the similarity between subject i with missing Y
    # in the original data set and subject j with observed Y in the bootstrap sample based on the  
    # predictive scores.
    
    # Select all cases in original data set with missing Y.
    to_be_imputed <- samp[samp$response_indicators == 0, -1]
    # Use function predict() to get predicted values for these cases.
    pred1.2 <- predict(mod1, newdata = to_be_imputed, type = "probs")[,-1]
    pred2.2 <- predict(mod2, newdata = to_be_imputed, type = "response")
    # Standardize the predicted values
    for(i in 1:ncol(pred1.2)){
      pred1.2[,i] <- (pred1.2[,i] - mu1[i]) / sd1[i]
    }
    pred2.2 <- (pred2.2 - mu2) / sd2
    # Create matrices with predicted values for subjects with missing Y in the original data set and 
    # subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(pred1.2, pred2.2) 
    pred.boot_obs <- cbind(pred1, pred2[boot$response_indicators == 1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot$response_indicators == 1)
    # Apply distance function.
    # Here, each column contains the distances from a subject i with missing Y in the original data
    # set to all subjects j with observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 4: For each incomplete subject in the original data set, define the imputation set as the 5 
    # nearest neighbors and impute by randomly drawing from the imputation set.
    # rank the distances 
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the indices of the 
    # subjects with observed Y in the bootstrap sample with the smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with missing Y in the original
    # data set.
    donors <- apply(indices, 2, function(x) boot[x, "y"])
    # Randomly draw one donor per subject with missing Y in the original data set
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute
    imputed_set <- samp
    imputed_set$y[imputed_set$response_indicators == 0] <- new.values
    return(imputed_set)
  }
  
  # Apply the function M times to obtain M completed data sets.
  DRNN00_sets <- lapply(1:5, DRNN00, samp)
  
  ## Impute using MICE0 (incorrect model)
  MICE0_sets <- mice(samp[,-4], m = M, maxit = 20)
  
  
  ## EVALUATION
  
  
  # MRNNI
  # Estimated proportions per data set
  # One column represents one data set
  MRNN10_proportions <- sapply(MRNN10_sets, function(x) table(x$y)/300)
  MRNN01_proportions <- sapply(MRNN01_sets, function(x) table(x$y)/300)
  MRNN00_proportions <- sapply(MRNN00_sets, function(x) table(x$y)/300)
  # Overall proportions:
  MRNN10_avg_props <- rowMeans(MRNN10_proportions)
  MRNN01_avg_props <- rowMeans(MRNN01_proportions)
  MRNN00_avg_props <- rowMeans(MRNN00_proportions)
  # Apply SE function to obtain SEs for every outcome category.
  SE_MRNN10 <- apply(MRNN10_proportions, 1, standard_error, M)
  SE_MRNN01 <- apply(MRNN01_proportions, 1, standard_error, M)
  SE_MRNN00 <- apply(MRNN00_proportions, 1, standard_error, M)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_MRNN10 <- matrix(CI(MRNN10_avg_props, SE_MRNN10), nrow=5)
  CIs_MRNN01 <- matrix(CI(MRNN01_avg_props, SE_MRNN01), nrow=5)
  CIs_MRNN00 <- matrix(CI(MRNN00_avg_props, SE_MRNN00), nrow=5)
  # Then check if the true population proportions are within the CIs (TRUE/FALSE).
  original_proportions <- table(population$y)/nrow(population)
  within_MRNN10 <- CIs_MRNN10[,1]<=original_proportions & original_proportions<=CIs_MRNN10[,2]
  within_MRNN01 <- CIs_MRNN01[,1]<=original_proportions & original_proportions<=CIs_MRNN01[,2]
  within_MRNN00 <- CIs_MRNN00[,1]<=original_proportions & original_proportions<=CIs_MRNN00[,2]
  
  # DRNNI
  # Estimated proportions per data set
  # One column represents one data set
  DRNN10_proportions <- sapply(DRNN10_sets, function(x) table(x$y)/300)
  DRNN01_proportions <- sapply(DRNN01_sets, function(x) table(x$y)/300)
  DRNN00_proportions <- sapply(DRNN00_sets, function(x) table(x$y)/300)
  # Overall proportions:
  DRNN10_avg_props <- rowMeans(DRNN10_proportions)
  DRNN01_avg_props <- rowMeans(DRNN01_proportions)
  DRNN00_avg_props <- rowMeans(DRNN00_proportions)
  # SE using Rubin's rules:
  # Apply function to obtain SEs for every outcome category.
  SE_DRNN10 <- apply(DRNN10_proportions, 1, standard_error, M)
  SE_DRNN01 <- apply(DRNN01_proportions, 1, standard_error, M)
  SE_DRNN00 <- apply(DRNN00_proportions, 1, standard_error, M)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_DRNN10 <- matrix(CI(DRNN10_avg_props, SE_DRNN10), nrow=5)
  CIs_DRNN01 <- matrix(CI(DRNN01_avg_props, SE_DRNN01), nrow=5)
  CIs_DRNN00 <- matrix(CI(DRNN00_avg_props, SE_DRNN00), nrow=5)
  # Then check if the original proportions are within the CIs (TRUE/FALSE).
  within_DRNN10 <- CIs_DRNN10[,1]<=original_proportions & original_proportions<=CIs_DRNN10[,2]
  within_DRNN01 <- CIs_DRNN01[,1]<=original_proportions & original_proportions<=CIs_DRNN01[,2]
  within_DRNN00 <- CIs_DRNN00[,1]<=original_proportions & original_proportions<=CIs_DRNN00[,2]
  
  # MICE
  # Estimated proportions per data set
  # One column represents one data set
  MICE0_proportions <- with(MICE0_sets, table(y)/300)$analyses
  MICE0_proportions <- matrix(unlist(MICE0_proportions), nrow = 5, ncol = M)
  # Overall proportions:
  MICE0_avg_props <- rowMeans(MICE0_proportions)
  # SE using Rubin's rules:
  # Apply function to obtain SEs for every outcome category.
  SE_MICE0 <- apply(MICE0_proportions, 1, standard_error, M)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_MICE0 <- matrix(CI(MICE0_avg_props, SE_MICE0), nrow=5)
  # Then check if the original proportions are within the CIs (TRUE/FALSE).
  within_MICE0 <- CIs_MICE0[,1]<=original_proportions & original_proportions<=CIs_MICE0[,2]
  
  # Output:
  MRNNI_out <- list(MRNN10_avg_props, MRNN01_avg_props, MRNN00_avg_props,
                    SE_MRNN10, SE_MRNN01, SE_MRNN00,
                    within_MRNN10, within_MRNN01, within_MRNN00)
  DRNNI_out <- list(DRNN10_avg_props, DRNN01_avg_props, DRNN00_avg_props,
                    SE_DRNN10, SE_DRNN01, SE_DRNN00,
                    within_DRNN10, within_DRNN01, within_DRNN00)
  MICE_out <- list(MICE0_avg_props, SE_MICE0, within_MICE0)
  
  return(list(MRNNI_out, DRNNI_out, MICE_out))
}

# Repeat 500 times to create, impute, and analyze 500 samples.
# These data need to be saved and exported.
set.seed(1323)
population4_response80 <- lapply(1:500, results_per_samp_fourth, population = population)
save(population4_response80, file = "S2pop4res80.RData")
set.seed(1325)
population4_response90 <- lapply(1:500, results_per_samp_fourth, population = population, 
                                 responserate80 = FALSE)
save(population4_response90, file = "S2pop4res90.RData")
