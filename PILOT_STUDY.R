##### PILOT STUDY ###################################################################################
# To determine the appropriate number of imputations, the performance of the three MI procedures with
# correctly specified models (DRNN11, MICE1, and MRNN11) will be compared for 5, 10, and 20 
# imputations.

# Set up
#library(parallel)
library(mice)
library(nnet)
library(tidyverse)
set.seed(2122)
#ncores <- detectCores()

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

# Generate auxiliary data
N <- 50000 # population size
x1 <- rbinom(N, 1, .5)
x2 <- rnorm(N)
x3 <- rnorm(N, 100, 15)

# Generate dependent variable:
# Generate outcome probabilities for each category
p2 <- exp(x1 - 5*x2 + .05*x3)/ (1 + exp(x1 - 5*x2 + .05*x3) + exp(2*x1 + 2*x2 - .03*x3))
p3 <- exp(2*x1 + 2*x2 - .03*x3)/ (1 + exp(x1 - 5*x2 + .05*x3) + exp(2*x1 + 2*x2 - .03*x3))
p1 <- 1 - p2 - p3
probs <- cbind(p1, p2, p3)
# Generate multinomial samples using the outcome probabilities 
y <- apply(probs, 1, function(x){rmultinom(1, 1, x)})
# Combine/Transform into a usable y variable
new1 <- ifelse(y[1,] == 1, 1, 0)
new2 <- ifelse(y[2,] == 1, 2, 0)
new3 <- ifelse(y[3,] == 1, 3, 0)
y <- as.factor(new1 + new2 + new3)

# Create a dataframe with the dependent and auxiliary variables.
population <- cbind.data.frame(y, x1, x2, x3)
save(population, file = "pop1.RData")

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
  p_res80 <- (exp(samp$x1 + samp$x2 + 0.01*samp$x3) / 
                        (exp(samp$x1 + samp$x2 + 0.01*samp$x3) + 1))
samp$response_indicators <- rbinom(300, 1, p_res80)}else{
  p_res90 <- (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) / 
                          (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) + 1))
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
    p_res80 <- (exp(samp$x1 + samp$x2 + 0.01*samp$x3) / 
                            (exp(samp$x1 + samp$x2 + 0.01*samp$x3) + 1))
    samp$response_indicators <- rbinom(300, 1, p_res80)}else{
      p_res90 <- (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) / 
                              (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) + 1))
      samp$response_indicators <- rbinom(300, 1, p_res90)
    }
  # Introduce missing values accordingly.
  samp$y[samp$response_indicators == 0] <- NA
}

## Impute using MRNNI
MRNNI <- function(iterations, samp){
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

# Apply the function M times to obtain M completed data sets.
MRNNI_sets5 <- lapply(1:5, MRNNI, samp)
MRNNI_sets10 <- lapply(1:10, MRNNI, samp)
MRNNI_sets20 <- lapply(1:20, MRNNI, samp)

## Impute using DRNNI
DRNNI <- function(iterations, samp){
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
  # rank the distances. 
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
DRNNI_sets5 <- lapply(1:5, DRNNI, samp)
DRNNI_sets10 <- lapply(1:10, DRNNI, samp)
DRNNI_sets20 <- lapply(1:20, DRNNI, samp)

## Impute using MICE
MICE_sets5 <- mice(samp, m = 5, maxit = 20)
MICE_sets10 <- mice(samp, m = 10, maxit = 20)
MICE_sets20 <- mice(samp, m = 20, mait = 20)


## EVALUATION

# MRNNI
# Estimated proportions per data set
# One column represents one data set
MRNNI5_proportions <- sapply(MRNNI_sets5, function(x) table(x$y)/300)
MRNNI10_proportions <- sapply(MRNNI_sets10, function(x) table(x$y)/300)
MRNNI20_proportions <- sapply(MRNNI_sets20, function(x) table(x$y)/300)
# Overall proportions:
MRNNI5_avg_props <- rowMeans(MRNNI5_proportions)
MRNNI10_avg_props <- rowMeans(MRNNI10_proportions)
MRNNI20_avg_props <- rowMeans(MRNNI20_proportions)
# Apply SE function to obtain SEs for every outcome category.
SE_MRNNI5 <- apply(MRNNI5_proportions, 1, standard_error, 5)
SE_MRNNI10 <- apply(MRNNI10_proportions, 1, standard_error, 10)
SE_MRNNI20 <- apply(MRNNI20_proportions, 1, standard_error, 20)
# Coverage rates:
# First calculate the 95% CI for the average proportions of each sample using the SEs.
# Each row contains the lower and upper bound for one outcome category.
CIs_MRNNI5 <- matrix(CI(MRNNI5_avg_props, SE_MRNNI5), nrow=3)
CIs_MRNNI10 <- matrix(CI(MRNNI10_avg_props, SE_MRNNI10), nrow=3)
CIs_MRNNI20 <- matrix(CI(MRNNI20_avg_props, SE_MRNNI20), nrow=3)
# Then check if the true population proportions are within the CIs (TRUE/FALSE).
original_proportions <- table(population$y)/nrow(population)
within_MRNNI5 <- CIs_MRNNI5[,1]<=original_proportions & original_proportions<=CIs_MRNNI5[,2]
within_MRNNI10 <- CIs_MRNNI10[,1]<=original_proportions & original_proportions<=CIs_MRNNI10[,2]
within_MRNNI20 <- CIs_MRNNI20[,1]<=original_proportions & original_proportions<=CIs_MRNNI20[,2]

# DRNNI
# Estimated proportions per data set
# One column represents one data set
DRNNI5_proportions <- sapply(DRNNI_sets5, function(x) table(x$y)/300)
DRNNI10_proportions <- sapply(DRNNI_sets10, function(x) table(x$y)/300)
DRNNI20_proportions <- sapply(DRNNI_sets20, function(x) table(x$y)/300)
# Overall proportions:
DRNNI5_avg_props <- rowMeans(DRNNI5_proportions)
DRNNI10_avg_props <- rowMeans(DRNNI10_proportions)
DRNNI20_avg_props <- rowMeans(DRNNI20_proportions)
# SE using Rubin's rules:
# Apply function to obtain SEs for every outcome category.
SE_DRNNI5 <- apply(DRNNI5_proportions, 1, standard_error, 5)
SE_DRNNI10 <- apply(DRNNI10_proportions, 1, standard_error, 10)
SE_DRNNI20 <- apply(DRNNI20_proportions, 1, standard_error, 20)
# Coverage rates:
# First calculate the 95% CI for the average proportions of each sample using the SEs.
# Each row contains the lower and upper bound for one outcome category.
CIs_DRNNI5 <- matrix(CI(DRNNI5_avg_props, SE_DRNNI5), nrow=3)
CIs_DRNNI10 <- matrix(CI(DRNNI10_avg_props, SE_DRNNI10), nrow=3)
CIs_DRNNI20 <- matrix(CI(DRNNI20_avg_props, SE_DRNNI20), nrow=3)
# Then check if the original proportions are within the CIs (TRUE/FALSE).
within_DRNNI5 <- CIs_DRNNI5[,1]<=original_proportions & original_proportions<=CIs_DRNNI5[,2]
within_DRNNI10 <- CIs_DRNNI10[,1]<=original_proportions & original_proportions<=CIs_DRNNI10[,2]
within_DRNNI20 <- CIs_DRNNI20[,1]<=original_proportions & original_proportions<=CIs_DRNNI20[,2]

# MICE
# Estimated proportions per data set
# One column represents one data set
MICE5_proportions <- with(MICE_sets5, table(y)/300)$analyses
MICE5_proportions <- matrix(unlist(MICE5_proportions), nrow = 3, ncol = 5)
MICE10_proportions <- with(MICE_sets10, table(y)/300)$analyses
MICE10_proportions <- matrix(unlist(MICE10_proportions), nrow = 3, ncol = 10)
MICE20_proportions <- with(MICE_sets20, table(y)/300)$analyses
MICE20_proportions <- matrix(unlist(MICE20_proportions), nrow = 3, ncol = 20)
# Overall proportions:
MICE5_avg_props <- rowMeans(MICE5_proportions)
MICE10_avg_props <- rowMeans(MICE10_proportions)
MICE20_avg_props <- rowMeans(MICE20_proportions)
# SE using Rubin's rules:
# Apply function to obtain SEs for every outcome category.
SE_MICE5 <- apply(MICE5_proportions, 1, standard_error, 5)
SE_MICE10 <- apply(MICE10_proportions, 1, standard_error, 10)
SE_MICE20 <- apply(MICE20_proportions, 1, standard_error, 20)
# Coverage rates:
# First calculate the 95% CI for the average proportions of each sample using the SEs.
# Each row contains the lower and upper bound for one outcome category.
CIs_MICE5 <- matrix(CI(MICE5_avg_props, SE_MICE5), nrow=3)
CIs_MICE10 <- matrix(CI(MICE10_avg_props, SE_MICE10), nrow=3)
CIs_MICE20 <- matrix(CI(MICE20_avg_props, SE_MICE20), nrow=3)
# Then check if the original proportions are within the CIs (TRUE/FALSE).
within_MICE5 <- CIs_MICE5[,1]<=original_proportions & original_proportions<=CIs_MICE5[,2]
within_MICE10 <- CIs_MICE10[,1]<=original_proportions & original_proportions<=CIs_MICE10[,2]
within_MICE20 <- CIs_MICE20[,1]<=original_proportions & original_proportions<=CIs_MICE20[,2]

# Output:
MRNNI_out <- list(MRNNI5_avg_props, MRNNI10_avg_props, MRNNI20_avg_props, 
              SE_MRNNI5, SE_MRNNI10, SE_MRNNI20,
              within_MRNNI5, within_MRNNI10, within_MRNNI20)
DRNNI_out <- list(DRNNI5_avg_props, DRNNI10_avg_props, DRNNI20_avg_props, 
              SE_DRNNI5, SE_DRNNI10, SE_DRNNI20,
              within_DRNNI5, within_DRNNI10, within_DRNNI20)
MICE_out <- list(MICE5_avg_props, MICE10_avg_props, MICE20_avg_props,
              SE_MICE5, SE_MICE10, SE_MICE20,
              within_MICE5, within_MICE10, within_MICE20)

return(list(MRNNI_out, DRNNI_out, MICE_out))
}

# Repeat 500 times to create, impute, and analyze 500 samples.
# These data need to be saved and exported.
population1_response80 <- lapply(1:500, results_per_samp_first, population = population)
save(population1_response80, file = "pop1res80.RData")
# Because I am running this over different sessions, I need to set a seed for every session
set.seed(22022022)
population1_response90 <- lapply(1:500, results_per_samp_first, population = population, 
                                   responserate80 = FALSE)
save(population1_response90, file = "pop1res90.RData")

#####################################################################################################
#####################################################################################################
##### POPULATION 2: 3 CATEGORIES FOR Y AND SECOND ORDER TERMS #######################################

set.seed(1234)
# Generate auxiliary data
N <- 50000 # population size
x1 <- rbinom(N, 1, .5)
x2 <- rnorm(N)
x3 <- rnorm(N, 100, 15)

# Generate dependent variable:
# Generate outcome probabilities for each category
p2 <- exp(.4*x2^2 - .02*x1*x3 + x1 + x2 + .01*x3)/ 
  (1 + exp(.4*x2^2 - .02*x1*x3 + x1 + x2 + .01*x3) + exp(x2^2 + .04*x1*x3 + x1 + .2*x2 - .05*x3))
p3 <- exp(x2^2 + .04*x1*x3 + x1 + .2*x2 - .05*x3)/ 
  (1 + exp(.4*x2^2 - .02*x1*x3 + x1 + x2 + .01*x3) + exp(x2^2 + .04*x1*x3 + x1 + .2*x2 - .05*x3))
p1 <- 1 - p2 - p3
probs <- cbind(p1, p2, p3)
# Generate multinomial samples using the outcome probabilities
y <- apply(probs, 1, function(x){rmultinom(1, 1, x)})
# Combine/Transform into a usable y variable
new1 <- ifelse(y[1,] == 1, 1, 0)
new2 <- ifelse(y[2,] == 1, 2, 0)
new3 <- ifelse(y[3,] == 1, 3, 0)
y <- as.factor(new1 + new2 + new3)

# Create a dataframe with the dependent and auxiliary variables.
population <- cbind.data.frame(y, x1, x2, x3)
save(population, file = "pop2.RData")

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
    p_res80 <- (exp(samp$x1 + samp$x2 + 0.01*samp$x3) / 
                            (exp(samp$x1 + samp$x2 + 0.01*samp$x3) + 1))
    samp$response_indicators <- rbinom(300, 1, p_res80)}else{
      p_res90 <- (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) / 
                              (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) + 1))
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
      p_res80 <- (exp(samp$x1 + samp$x2 + 0.01*samp$x3) / 
                              (exp(samp$x1 + samp$x2 + 0.01*samp$x3) + 1))
      samp$response_indicators <- rbinom(300, 1, p_res80)}else{
        p_res90 <- (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) / 
                                (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) + 1))
        samp$response_indicators <- rbinom(300, 1, p_res90)
      }
    # Introduce missing values accordingly.
    samp$y[samp$response_indicators == 0] <- NA
  }
  
  ## Impute using MRNNI
  MRNNI <- function(iterations, samp){
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
    
    # step 5: calculate a distance function to define the similarity between subject i with missing 
    # Y  in  the original data set and subject j with observed Y in the bootstrap sample based on the
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
    
    # step 6: For each incomplete subject in the original data set, define the imputation set as the 
    # 5 nearest neighbors and impute by randomly drawing from the imputation set.
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
  MRNNI_sets5 <- lapply(1:5, MRNNI, samp)
  MRNNI_sets10 <- lapply(1:10, MRNNI, samp)
  MRNNI_sets20 <- lapply(1:20, MRNNI, samp)
  
  ## Impute using DRNNI
  DRNNI <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=3){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2a: Using the bootstrap sample, fit a multinomial regression model on the dependent 
    # variable using the auxiliary data and save the standardized predicted values.
    mod1 <- multinom(y ~ x2^2 + x1:x3 + x1 + x2 + x3, boot)
    pred1 <- fitted(mod1)[,-1]
    # We are manually standardizing because we need to save the means and sds for later.
    mu1 <- colMeans(pred1)
    sd1 <- apply(pred1, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1)){
      pred1[,i] <- (pred1[,i] - mu1[i]) / sd1[i]
    }
    # step 2b: Using the bootstrap sample, fit a logistic regression model on the response 
    # indicators using the auxiliary data and save the standardized predicted values.
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
    
    # step 4: For each incomplete subject in the original data set, define the imputation set as the 
    # 5  nearest neighbors and impute by randomly drawing from the imputation set.
    # rank the distances. 
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
  DRNNI_sets5 <- lapply(1:5, DRNNI, samp)
  DRNNI_sets10 <- lapply(1:10, DRNNI, samp)
  DRNNI_sets20 <- lapply(1:20, DRNNI, samp)
  
  ## Impute using MICE
  MICE_sets5 <- mice(samp, m = 5, maxit = 20)
  MICE_sets10 <- mice(samp, m = 10, maxit = 20)
  MICE_sets20 <- mice(samp, m = 20, maxit = 20)
  
  
  ## EVALUATION
  
  # MRNNI
  # Estimated proportions per data set
  # One column represents one data set
  MRNNI5_proportions <- sapply(MRNNI_sets5, function(x) table(x$y)/300)
  MRNNI10_proportions <- sapply(MRNNI_sets10, function(x) table(x$y)/300)
  MRNNI20_proportions <- sapply(MRNNI_sets20, function(x) table(x$y)/300)
  # Overall proportions:
  MRNNI5_avg_props <- rowMeans(MRNNI5_proportions)
  MRNNI10_avg_props <- rowMeans(MRNNI10_proportions)
  MRNNI20_avg_props <- rowMeans(MRNNI20_proportions)
  # Apply SE function to obtain SEs for every outcome category.
  SE_MRNNI5 <- apply(MRNNI5_proportions, 1, standard_error, 5)
  SE_MRNNI10 <- apply(MRNNI10_proportions, 1, standard_error, 10)
  SE_MRNNI20 <- apply(MRNNI20_proportions, 1, standard_error, 20)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_MRNNI5 <- matrix(CI(MRNNI5_avg_props, SE_MRNNI5), nrow=3)
  CIs_MRNNI10 <- matrix(CI(MRNNI10_avg_props, SE_MRNNI10), nrow=3)
  CIs_MRNNI20 <- matrix(CI(MRNNI20_avg_props, SE_MRNNI20), nrow=3)
  # Then check if the true population proportions are within the CIs (TRUE/FALSE).
  original_proportions <- table(population$y)/nrow(population)
  within_MRNNI5 <- CIs_MRNNI5[,1]<=original_proportions & original_proportions<=CIs_MRNNI5[,2]
  within_MRNNI10 <- CIs_MRNNI10[,1]<=original_proportions & original_proportions<=CIs_MRNNI10[,2]
  within_MRNNI20 <- CIs_MRNNI20[,1]<=original_proportions & original_proportions<=CIs_MRNNI20[,2]
  
  # DRNNI
  # Estimated proportions per data set
  # One column represents one data set
  DRNNI5_proportions <- sapply(DRNNI_sets5, function(x) table(x$y)/300)
  DRNNI10_proportions <- sapply(DRNNI_sets10, function(x) table(x$y)/300)
  DRNNI20_proportions <- sapply(DRNNI_sets20, function(x) table(x$y)/300)
  # Overall proportions:
  DRNNI5_avg_props <- rowMeans(DRNNI5_proportions)
  DRNNI10_avg_props <- rowMeans(DRNNI10_proportions)
  DRNNI20_avg_props <- rowMeans(DRNNI20_proportions)
  # SE using Rubin's rules:
  # Apply function to obtain SEs for every outcome category.
  SE_DRNNI5 <- apply(DRNNI5_proportions, 1, standard_error, 5)
  SE_DRNNI10 <- apply(DRNNI10_proportions, 1, standard_error, 10)
  SE_DRNNI20 <- apply(DRNNI20_proportions, 1, standard_error, 20)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_DRNNI5 <- matrix(CI(DRNNI5_avg_props, SE_DRNNI5), nrow=3)
  CIs_DRNNI10 <- matrix(CI(DRNNI10_avg_props, SE_DRNNI10), nrow=3)
  CIs_DRNNI20 <- matrix(CI(DRNNI20_avg_props, SE_DRNNI20), nrow=3)
  # Then check if the original proportions are within the CIs (TRUE/FALSE).
  within_DRNNI5 <- CIs_DRNNI5[,1]<=original_proportions & original_proportions<=CIs_DRNNI5[,2]
  within_DRNNI10 <- CIs_DRNNI10[,1]<=original_proportions & original_proportions<=CIs_DRNNI10[,2]
  within_DRNNI20 <- CIs_DRNNI20[,1]<=original_proportions & original_proportions<=CIs_DRNNI20[,2]
  
  # MICE
  # Estimated proportions per data set
  # One column represents one data set
  MICE5_proportions <- with(MICE_sets5, table(y)/300)$analyses
  MICE5_proportions <- matrix(unlist(MICE5_proportions), nrow = 3, ncol = 5)
  MICE10_proportions <- with(MICE_sets10, table(y)/300)$analyses
  MICE10_proportions <- matrix(unlist(MICE10_proportions), nrow = 3, ncol = 10)
  MICE20_proportions <- with(MICE_sets20, table(y)/300)$analyses
  MICE20_proportions <- matrix(unlist(MICE20_proportions), nrow = 3, ncol = 20)
  # Overall proportions:
  MICE5_avg_props <- rowMeans(MICE5_proportions)
  MICE10_avg_props <- rowMeans(MICE10_proportions)
  MICE20_avg_props <- rowMeans(MICE20_proportions)
  # SE using Rubin's rules:
  # Apply function to obtain SEs for every outcome category.
  SE_MICE5 <- apply(MICE5_proportions, 1, standard_error, 5)
  SE_MICE10 <- apply(MICE10_proportions, 1, standard_error, 10)
  SE_MICE20 <- apply(MICE20_proportions, 1, standard_error, 20)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_MICE5 <- matrix(CI(MICE5_avg_props, SE_MICE5), nrow=3)
  CIs_MICE10 <- matrix(CI(MICE10_avg_props, SE_MICE10), nrow=3)
  CIs_MICE20 <- matrix(CI(MICE20_avg_props, SE_MICE20), nrow=3)
  # Then check if the original proportions are within the CIs (TRUE/FALSE).
  within_MICE5 <- CIs_MICE5[,1]<=original_proportions & original_proportions<=CIs_MICE5[,2]
  within_MICE10 <- CIs_MICE10[,1]<=original_proportions & original_proportions<=CIs_MICE10[,2]
  within_MICE20 <- CIs_MICE20[,1]<=original_proportions & original_proportions<=CIs_MICE20[,2]
  
  # Output:
  MRNNI_out <- list(MRNNI5_avg_props, MRNNI10_avg_props, MRNNI20_avg_props, 
                    SE_MRNNI5, SE_MRNNI10, SE_MRNNI20,
                    within_MRNNI5, within_MRNNI10, within_MRNNI20)
  DRNNI_out <- list(DRNNI5_avg_props, DRNNI10_avg_props, DRNNI20_avg_props, 
                    SE_DRNNI5, SE_DRNNI10, SE_DRNNI20,
                    within_DRNNI5, within_DRNNI10, within_DRNNI20)
  MICE_out <- list(MICE5_avg_props, MICE10_avg_props, MICE20_avg_props,
                   SE_MICE5, SE_MICE10, SE_MICE20,
                   within_MICE5, within_MICE10, within_MICE20)
  
  return(list(MRNNI_out, DRNNI_out, MICE_out))
}

# Repeat 500 times to create, impute, and analyze 500 samples.
# These data need to be saved and exported.
population2_response80 <- lapply(1:500, results_per_samp_second, population = population)
save(population2_response80, file = "pop2res80.RData")
set.seed(23022022)
population2_response90 <- lapply(1:500, results_per_samp_second, population = population, 
                                   responserate80 = FALSE)
save(population2_response90, file = "pop2res90.RData")

#####################################################################################################
#####################################################################################################
##### POPULATION 3: 5 CATEGORIES FOR Y AND FIRST ORDER TERMS ########################################

set.seed(2402022)
# Generate auxiliary data
N <- 50000 # population size
x1 <- rbinom(N, 1, .5)
x2 <- rnorm(N)
x3 <- rnorm(N, 100, 15)

# Generate dependent variable:
# Generate outcome probabilities for each category
p2 <- exp(.2*x1 + .1*x2 - .01*x3)/
  (1 + exp(.2*x1 + .1*x2 - .01*x3) + exp(-.1*x1 - .01*x2 + .01*x3) + 
     exp(.1*x1 + .02*x2 - .01*x3) + exp(.57*x1 - .05*x2 - .01*x3))
p3 <- exp(-.1*x1 - .01*x2 + .01*x3)/
  (1 + exp(.2*x1 + .1*x2 - .01*x3) + exp(-.1*x1 - .01*x2 + .01*x3) + 
     exp(.1*x1 + .02*x2 - .01*x3) + exp(.57*x1 - .05*x2 - .01*x3))
p4 <- exp(.1*x1 + .02*x2 - .01*x3)/
  (1 + exp(.2*x1 + .1*x2 - .01*x3) + exp(-.1*x1 - .01*x2 + .01*x3) + 
     exp(.1*x1 + .02*x2 - .01*x3) + exp(.57*x1 - .05*x2 - .01*x3))
p5 <- exp(.57*x1 - .05*x2 - .01*x3)/
  (1 + exp(.2*x1 + .1*x2 - .01*x3) + exp(-.1*x1 - .01*x2 + .01*x3) + 
     exp(.1*x1 + .02*x2 - .01*x3) + exp(.57*x1 - .05*x2 - .01*x3))
p1 <- 1 - p2 - p3 - p4 - p5
probs <- cbind(p1, p2, p3, p4, p5)
# Generate multinomial samples using the outcome probabilities
y <- apply(probs, 1, function(x){rmultinom(1, 1, x)})
# Combine/Transform into a usable y variable
new1 <- ifelse(y[1,] == 1, 1, 0)
new2 <- ifelse(y[2,] == 1, 2, 0)
new3 <- ifelse(y[3,] == 1, 3, 0)
new4 <- ifelse(y[4,] == 1, 4, 0)
new5 <- ifelse(y[5,] == 1, 5, 0)
y <- as.factor(new1 + new2 + new3 + new4 + new5)

# Create a dataframe with the dependent and auxiliary variables.
population <- cbind.data.frame(y, x1, x2, x3)
save(population, file = "pop3.RData")

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
    p_res80 <- (exp(samp$x1 + samp$x2 + 0.01*samp$x3) / 
                            (exp(samp$x1 + samp$x2 + 0.01*samp$x3) + 1))
    samp$response_indicators <- rbinom(300, 1, p_res80)}else{
      p_res90 <- (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) / 
                              (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) + 1))
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
      p_res80 <- (exp(samp$x1 + samp$x2 + 0.01*samp$x3) / 
                              (exp(samp$x1 + samp$x2 + 0.01*samp$x3) + 1))
      samp$response_indicators <- rbinom(300, 1, p_res80)}else{
        p_res90 <- (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) / 
                                (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) + 1))
        samp$response_indicators <- rbinom(300, 1, p_res90)
      }
    # Introduce missing values accordingly.
    samp$y[samp$response_indicators == 0] <- NA
  }
  
  ## Impute using MRNNI
  MRNNI <- function(iterations, samp){
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
  MRNNI_sets5 <- lapply(1:5, MRNNI, samp)
  MRNNI_sets10 <- lapply(1:10, MRNNI, samp)
  MRNNI_sets20 <- lapply(1:20, MRNNI, samp)
  
  ## Impute using DRNNI
  DRNNI <- function(iterations, samp){
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
    # rank the distances.
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
  DRNNI_sets5 <- lapply(1:5, DRNNI, samp)
  DRNNI_sets10 <- lapply(1:10, DRNNI, samp)
  DRNNI_sets20 <- lapply(1:20, DRNNI, samp)
  
  ## Impute using MICE
  MICE_sets5 <- mice(samp, m = 5, maxit = 20)
  MICE_sets10 <- mice(samp, m = 10, maxit = 20)
  MICE_sets20 <- mice(samp, m = 20, maxit = 20)
  
  
  ## EVALUATION
  
  # MRNNI
  # Estimated proportions per data set
  # One column represents one data set
  MRNNI5_proportions <- sapply(MRNNI_sets5, function(x) table(x$y)/300)
  MRNNI10_proportions <- sapply(MRNNI_sets10, function(x) table(x$y)/300)
  MRNNI20_proportions <- sapply(MRNNI_sets20, function(x) table(x$y)/300)
  # Overall proportions:
  MRNNI5_avg_props <- rowMeans(MRNNI5_proportions)
  MRNNI10_avg_props <- rowMeans(MRNNI10_proportions)
  MRNNI20_avg_props <- rowMeans(MRNNI20_proportions)
  # Apply SE function to obtain SEs for every outcome category.
  SE_MRNNI5 <- apply(MRNNI5_proportions, 1, standard_error, 5)
  SE_MRNNI10 <- apply(MRNNI10_proportions, 1, standard_error, 10)
  SE_MRNNI20 <- apply(MRNNI20_proportions, 1, standard_error, 20)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_MRNNI5 <- matrix(CI(MRNNI5_avg_props, SE_MRNNI5), nrow=5)
  CIs_MRNNI10 <- matrix(CI(MRNNI10_avg_props, SE_MRNNI10), nrow=5)
  CIs_MRNNI20 <- matrix(CI(MRNNI20_avg_props, SE_MRNNI20), nrow=5)
  # Then check if the true population proportions are within the CIs (TRUE/FALSE).
  original_proportions <- table(population$y)/nrow(population)
  within_MRNNI5 <- CIs_MRNNI5[,1]<=original_proportions & original_proportions<=CIs_MRNNI5[,2]
  within_MRNNI10 <- CIs_MRNNI10[,1]<=original_proportions & original_proportions<=CIs_MRNNI10[,2]
  within_MRNNI20 <- CIs_MRNNI20[,1]<=original_proportions & original_proportions<=CIs_MRNNI20[,2]
  
  # DRNNI
  # Estimated proportions per data set
  # One column represents one data set
  DRNNI5_proportions <- sapply(DRNNI_sets5, function(x) table(x$y)/300)
  DRNNI10_proportions <- sapply(DRNNI_sets10, function(x) table(x$y)/300)
  DRNNI20_proportions <- sapply(DRNNI_sets20, function(x) table(x$y)/300)
  # Overall proportions:
  DRNNI5_avg_props <- rowMeans(DRNNI5_proportions)
  DRNNI10_avg_props <- rowMeans(DRNNI10_proportions)
  DRNNI20_avg_props <- rowMeans(DRNNI20_proportions)
  # SE using Rubin's rules:
  # Apply function to obtain SEs for every outcome category.
  SE_DRNNI5 <- apply(DRNNI5_proportions, 1, standard_error, 5)
  SE_DRNNI10 <- apply(DRNNI10_proportions, 1, standard_error, 10)
  SE_DRNNI20 <- apply(DRNNI20_proportions, 1, standard_error, 20)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_DRNNI5 <- matrix(CI(DRNNI5_avg_props, SE_DRNNI5), nrow=5)
  CIs_DRNNI10 <- matrix(CI(DRNNI10_avg_props, SE_DRNNI10), nrow=5)
  CIs_DRNNI20 <- matrix(CI(DRNNI20_avg_props, SE_DRNNI20), nrow=5)
  # Then check if the original proportions are within the CIs (TRUE/FALSE).
  within_DRNNI5 <- CIs_DRNNI5[,1]<=original_proportions & original_proportions<=CIs_DRNNI5[,2]
  within_DRNNI10 <- CIs_DRNNI10[,1]<=original_proportions & original_proportions<=CIs_DRNNI10[,2]
  within_DRNNI20 <- CIs_DRNNI20[,1]<=original_proportions & original_proportions<=CIs_DRNNI20[,2]
  
  # MICE
  # Estimated proportions per data set
  # One column represents oen data set
  MICE5_proportions <- with(MICE_sets5, table(y)/300)$analyses
  MICE5_proportions <- matrix(unlist(MICE5_proportions), nrow = 5, ncol = 5)
  MICE10_proportions <- with(MICE_sets10, table(y)/300)$analyses
  MICE10_proportions <- matrix(unlist(MICE10_proportions), nrow = 5, ncol = 10)
  MICE20_proportions <- with(MICE_sets20, table(y)/300)$analyses
  MICE20_proportions <- matrix(unlist(MICE20_proportions), nrow = 5, ncol = 20)
  # Overall proportions:
  MICE5_avg_props <- rowMeans(MICE5_proportions)
  MICE10_avg_props <- rowMeans(MICE10_proportions)
  MICE20_avg_props <- rowMeans(MICE20_proportions)
  # SE using Rubin's rules:
  # Apply function to obtain SEs for every outcome category.
  SE_MICE5 <- apply(MICE5_proportions, 1, standard_error, 5)
  SE_MICE10 <- apply(MICE10_proportions, 1, standard_error, 10)
  SE_MICE20 <- apply(MICE20_proportions, 1, standard_error, 20)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_MICE5 <- matrix(CI(MICE5_avg_props, SE_MICE5), nrow=5)
  CIs_MICE10 <- matrix(CI(MICE10_avg_props, SE_MICE10), nrow=5)
  CIs_MICE20 <- matrix(CI(MICE20_avg_props, SE_MICE20), nrow=5)
  # Then check if the original proportions are within the CIs (TRUE/FALSE).
  within_MICE5 <- CIs_MICE5[,1]<=original_proportions & original_proportions<=CIs_MICE5[,2]
  within_MICE10 <- CIs_MICE10[,1]<=original_proportions & original_proportions<=CIs_MICE10[,2]
  within_MICE20 <- CIs_MICE20[,1]<=original_proportions & original_proportions<=CIs_MICE20[,2]
  
  # Output:
  MRNNI_out <- list(MRNNI5_avg_props, MRNNI10_avg_props, MRNNI20_avg_props, 
                    SE_MRNNI5, SE_MRNNI10, SE_MRNNI20,
                    within_MRNNI5, within_MRNNI10, within_MRNNI20)
  DRNNI_out <- list(DRNNI5_avg_props, DRNNI10_avg_props, DRNNI20_avg_props, 
                    SE_DRNNI5, SE_DRNNI10, SE_DRNNI20,
                    within_DRNNI5, within_DRNNI10, within_DRNNI20)
  MICE_out <- list(MICE5_avg_props, MICE10_avg_props, MICE20_avg_props,
                   SE_MICE5, SE_MICE10, SE_MICE20,
                   within_MICE5, within_MICE10, within_MICE20)
  
  return(list(MRNNI_out, DRNNI_out, MICE_out))
}

population3_response80 <- lapply(1:500, results_per_samp_third, population = population)
save(population3_response80, file = "pop3res80.RData")
set.seed(1456)
population3_response90 <- lapply(1:500, results_per_samp_third, population = population, 
                                   responserate80 = FALSE)
save(population3_response90, file = "pop3res90.RData")

#####################################################################################################
#####################################################################################################
##### POPULATION 4: 5 CATEGORIES FOR 5 AND SECOND ORDER TERMS #######################################

set.seed(25022022)
# Generate auxiliary data
N <- 50000 # population size
x1 <- rbinom(N, 1, .5)
x2 <- rnorm(N)
x3 <- rnorm(N, 100, 15)

# Generate dependent variable:
# Generate outcome probabilities for each category
p2 <- exp(-.1*x2^2 - .01*x1*x3 - .2*x1 + .1*x2 + .01*x3) /
  (1 + exp(-.1*x2^2 - .01*x1*x3 - .2*x1 + .1*x2 + .01*x3) + 
     exp(-.2*x2^2 - .01*x1*x3 + x1 + .2*x2 - .01*x3) + 
     exp(.15*x2^2 + .01*x1*x3 + .3*x1 + .08*x2 - .03*x3) +
     exp(-.08*x2^2 - .01*x1*x3 + .7*x1 - .02*x2 - .01*x3))
p3 <- exp(-.2*x2^2 - .01*x1*x3 + x1 + .2*x2 - .01*x3) /
  (1 + exp(-.1*x2^2 - .01*x1*x3 - .2*x1 + .1*x2 + .01*x3) + 
     exp(-.2*x2^2 - .01*x1*x3 + x1 + .2*x2 - .01*x3) + 
     exp(.15*x2^2 + .01*x1*x3 + .3*x1 + .08*x2 - .03*x3) +
     exp(-.08*x2^2 - .01*x1*x3 + .7*x1 - .02*x2 - .01*x3))
p4 <- exp(.15*x2^2 + .01*x1*x3 + .3*x1 + .08*x2 - .03*x3) /
  (1 + exp(-.1*x2^2 - .01*x1*x3 - .2*x1 + .1*x2 + .01*x3) + 
     exp(-.2*x2^2 - .01*x1*x3 + x1 + .2*x2 - .01*x3) + 
     exp(.15*x2^2 + .01*x1*x3 + .3*x1 + .08*x2 - .03*x3) +
     exp(-.08*x2^2 - .01*x1*x3 + .7*x1 - .02*x2 - .01*x3))
p5 <- exp(-.08*x2^2 - .01*x1*x3 + .7*x1 - .02*x2 - .01*x3) /
  (1 + exp(-.1*x2^2 - .01*x1*x3 - .2*x1 + .1*x2 + .01*x3) + 
     exp(-.2*x2^2 - .01*x1*x3 + x1 + .2*x2 - .01*x3) + 
     exp(.15*x2^2 + .01*x1*x3 + .3*x1 + .08*x2 - .03*x3) +
     exp(-.08*x2^2 - .01*x1*x3 + .7*x1 - .02*x2 - .01*x3))
p1 <- 1 - p2 - p3 - p4 - p5
probs <- cbind(p1, p2, p3, p4, p5)
# Generate multinomial samples using the outcome probabilities
y <- apply(probs, 1, function(x){rmultinom(1, 1, x)})
# Combine/Transform into a usable y variable
new1 <- ifelse(y[1,] == 1, 1, 0)
new2 <- ifelse(y[2,] == 1, 2, 0)
new3 <- ifelse(y[3,] == 1, 3, 0)
new4 <- ifelse(y[4,] == 1, 4, 0)
new5 <- ifelse(y[5,] == 1, 5, 0)
y <- as.factor(new1 + new2 + new3 + new4 + new5)

# Create a dataframe with the dependent and auxiliary variables.
population <- cbind.data.frame(y, x1, x2, x3)
save(population, file = "pop4.RData")

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
    p_res80 <- (exp(samp$x1 + samp$x2 + 0.01*samp$x3) / 
                            (exp(samp$x1 + samp$x2 + 0.01*samp$x3) + 1))
    samp$response_indicators <- rbinom(300, 1, p_res80)}else{
      p_res90 <- (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) / 
                              (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) + 1))
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
      p_res80 <- (exp(samp$x1 + samp$x2 + 0.01*samp$x3) / 
                              (exp(samp$x1 + samp$x2 + 0.01*samp$x3) + 1))
      samp$response_indicators <- rbinom(300, 1, p_res80)}else{
        p_res90 <- (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) / 
                                (exp(1 + samp$x1 + samp$x2 + 0.01*samp$x3) + 1))
        samp$response_indicators <- rbinom(300, 1, p_res90)
      }
    # Introduce missing values accordingly.
    samp$y[samp$response_indicators == 0] <- NA
  }
  
  ## Impute using MRNNI
  MRNNI <- function(iterations, samp){
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
    
    # step 5: calculate a distance function to define the similarity between subject i with missing 
    # Y  in  the original data set and subject j with observed Y in the bootstrap sample based on the
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
    
    # step 6: For each incomplete subject in the original data set, define the imputation set as the 
    # 5 nearest neighbors and impute by randomly drawing from the imputation set.
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
  MRNNI_sets5 <- lapply(1:5, MRNNI, samp)
  MRNNI_sets10 <- lapply(1:10, MRNNI, samp)
  MRNNI_sets20 <- lapply(1:20, MRNNI, samp)
  
  ## Impute using DRNNI
  DRNNI <- function(iterations, samp){
    # step 1: Bootstrap the original data set.
    boot <- samp[sample(1:300, 300, replace = TRUE),]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot$response_indicators == 1, "y"]))!=5){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # step 2a: Using the bootstrap sample, fit a multinomial regression model on the dependent 
    # variable using the auxiliary data and save the standardized predicted values.
    mod1 <- multinom(y ~ x2^2 + x1:x3 + x1 + x2 + x3, boot)
    pred1 <- fitted(mod1)[,-1]
    # We are manually standardizing because we need to save the means and sds for later.
    mu1 <- colMeans(pred1)
    sd1 <- apply(pred1, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred1)){
      pred1[,i] <- (pred1[,i] - mu1[i]) / sd1[i]
    }
    # step 2b: Using the bootstrap sample, fit a logistic regression model on the response 
    # indicators using the auxiliary data and save the standardized predicted values.
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
    
    # step 4: For each incomplete subject in the original data set, define the imputation set as the 
    # 5  nearest neighbors and impute by randomly drawing from the imputation set.
    # rank the distances. 
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
  DRNNI_sets5 <- lapply(1:5, DRNNI, samp)
  DRNNI_sets10 <- lapply(1:10, DRNNI, samp)
  DRNNI_sets20 <- lapply(1:20, DRNNI, samp)
  
  ## Impute using MICE
  MICE_sets5 <- mice(samp, m = 5, maxit = 20)
  MICE_sets10 <- mice(samp, m = 10, maxit = 20)
  MICE_sets20 <- mice(samp, m = 20, maxit = 20)
  
  
  ## EVALUATION
  
  # MRNNI
  # Estimated proportions per data set
  # One column represents one data set
  MRNNI5_proportions <- sapply(MRNNI_sets5, function(x) table(x$y)/300)
  MRNNI10_proportions <- sapply(MRNNI_sets10, function(x) table(x$y)/300)
  MRNNI20_proportions <- sapply(MRNNI_sets20, function(x) table(x$y)/300)
  # Overall proportions:
  MRNNI5_avg_props <- rowMeans(MRNNI5_proportions)
  MRNNI10_avg_props <- rowMeans(MRNNI10_proportions)
  MRNNI20_avg_props <- rowMeans(MRNNI20_proportions)
  # Apply SE function to obtain SEs for every outcome category.
  SE_MRNNI5 <- apply(MRNNI5_proportions, 1, standard_error, 5)
  SE_MRNNI10 <- apply(MRNNI10_proportions, 1, standard_error, 10)
  SE_MRNNI20 <- apply(MRNNI20_proportions, 1, standard_error, 20)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_MRNNI5 <- matrix(CI(MRNNI5_avg_props, SE_MRNNI5), nrow=5)
  CIs_MRNNI10 <- matrix(CI(MRNNI10_avg_props, SE_MRNNI10), nrow=5)
  CIs_MRNNI20 <- matrix(CI(MRNNI20_avg_props, SE_MRNNI20), nrow=5)
  # Then check if the true population proportions are within the CIs (TRUE/FALSE).
  original_proportions <- table(population$y)/nrow(population)
  within_MRNNI5 <- CIs_MRNNI5[,1]<=original_proportions & original_proportions<=CIs_MRNNI5[,2]
  within_MRNNI10 <- CIs_MRNNI10[,1]<=original_proportions & original_proportions<=CIs_MRNNI10[,2]
  within_MRNNI20 <- CIs_MRNNI20[,1]<=original_proportions & original_proportions<=CIs_MRNNI20[,2]
  
  # DRNNI
  # Estimated proportions per data set
  # One column represents one data set
  DRNNI5_proportions <- sapply(DRNNI_sets5, function(x) table(x$y)/300)
  DRNNI10_proportions <- sapply(DRNNI_sets10, function(x) table(x$y)/300)
  DRNNI20_proportions <- sapply(DRNNI_sets20, function(x) table(x$y)/300)
  # Overall proportions:
  DRNNI5_avg_props <- rowMeans(DRNNI5_proportions)
  DRNNI10_avg_props <- rowMeans(DRNNI10_proportions)
  DRNNI20_avg_props <- rowMeans(DRNNI20_proportions)
  # SE using Rubin's rules:
  # Apply function to obtain SEs for every outcome category.
  SE_DRNNI5 <- apply(DRNNI5_proportions, 1, standard_error, 5)
  SE_DRNNI10 <- apply(DRNNI10_proportions, 1, standard_error, 10)
  SE_DRNNI20 <- apply(DRNNI20_proportions, 1, standard_error, 20)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_DRNNI5 <- matrix(CI(DRNNI5_avg_props, SE_DRNNI5), nrow=5)
  CIs_DRNNI10 <- matrix(CI(DRNNI10_avg_props, SE_DRNNI10), nrow=5)
  CIs_DRNNI20 <- matrix(CI(DRNNI20_avg_props, SE_DRNNI20), nrow=5)
  # Then check if the original proportions are within the CIs (TRUE/FALSE).
  within_DRNNI5 <- CIs_DRNNI5[,1]<=original_proportions & original_proportions<=CIs_DRNNI5[,2]
  within_DRNNI10 <- CIs_DRNNI10[,1]<=original_proportions & original_proportions<=CIs_DRNNI10[,2]
  within_DRNNI20 <- CIs_DRNNI20[,1]<=original_proportions & original_proportions<=CIs_DRNNI20[,2]
  
  # MICE
  # Estimated proportions per data set
  # One column represents one data set
  MICE5_proportions <- with(MICE_sets5, table(y)/300)$analyses
  MICE5_proportions <- matrix(unlist(MICE5_proportions), nrow = 5, ncol = 5)
  MICE10_proportions <- with(MICE_sets10, table(y)/300)$analyses
  MICE10_proportions <- matrix(unlist(MICE10_proportions), nrow = 5, ncol = 10)
  MICE20_proportions <- with(MICE_sets20, table(y)/300)$analyses
  MICE20_proportions <- matrix(unlist(MICE20_proportions), nrow = 5, ncol = 20)
  # Overall proportions:
  MICE5_avg_props <- rowMeans(MICE5_proportions)
  MICE10_avg_props <- rowMeans(MICE10_proportions)
  MICE20_avg_props <- rowMeans(MICE20_proportions)
  # SE using Rubin's rules:
  # Apply function to obtain SEs for every outcome category.
  SE_MICE5 <- apply(MICE5_proportions, 1, standard_error, 5)
  SE_MICE10 <- apply(MICE10_proportions, 1, standard_error, 10)
  SE_MICE20 <- apply(MICE20_proportions, 1, standard_error, 20)
  # Coverage rates:
  # First calculate the 95% CI for the average proportions of each sample using the SEs.
  # Each row contains the lower and upper bound for one outcome category.
  CIs_MICE5 <- matrix(CI(MICE5_avg_props, SE_MICE5), nrow=5)
  CIs_MICE10 <- matrix(CI(MICE10_avg_props, SE_MICE10), nrow=5)
  CIs_MICE20 <- matrix(CI(MICE20_avg_props, SE_MICE20), nrow=5)
  # Then check if the original proportions are within the CIs (TRUE/FALSE).
  within_MICE5 <- CIs_MICE5[,1]<=original_proportions & original_proportions<=CIs_MICE5[,2]
  within_MICE10 <- CIs_MICE10[,1]<=original_proportions & original_proportions<=CIs_MICE10[,2]
  within_MICE20 <- CIs_MICE20[,1]<=original_proportions & original_proportions<=CIs_MICE20[,2]
  
  # Output:
  MRNNI_out <- list(MRNNI5_avg_props, MRNNI10_avg_props, MRNNI20_avg_props, 
                    SE_MRNNI5, SE_MRNNI10, SE_MRNNI20,
                    within_MRNNI5, within_MRNNI10, within_MRNNI20)
  DRNNI_out <- list(DRNNI5_avg_props, DRNNI10_avg_props, DRNNI20_avg_props, 
                    SE_DRNNI5, SE_DRNNI10, SE_DRNNI20,
                    within_DRNNI5, within_DRNNI10, within_DRNNI20)
  MICE_out <- list(MICE5_avg_props, MICE10_avg_props, MICE20_avg_props,
                   SE_MICE5, SE_MICE10, SE_MICE20,
                   within_MICE5, within_MICE10, within_MICE20)
  
  return(list(MRNNI_out, DRNNI_out, MICE_out))
}

population4_response80 <- lapply(1:500, results_per_samp_fourth, population = population)
save(population4_response80, file = "pop4res80.RData")
set.seed(1607)
population4_response90 <- lapply(1:500, results_per_samp_fourth, population = population, 
                                   responserate80 = FALSE)
save(population4_response90, file = "pop4res90.RData")
