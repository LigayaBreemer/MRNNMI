### MRNNMI new code??
library(nnet)

MRNNMI <- function(outcome_models, missing_models, data, response_indicator, 
                   imputations = 5){
  # Outcome_models: a list containing the multinom model objects for the outcome
  # models.
  # Missing_models: a list containing the glm model objects for the missingness 
  # models.
  # Data: the data set to be imputed (as a dataframe).
  # Response_indicator: vector of response indicators for the response
  # variable.
  # Imputations: integer specifying the number of imputed data sets to be 
  # created.
  # Imputed_data (output): a list containing the imputed data sets.
  
  # Checks
  if(is.null(data))stop("no data set given to be imputed")
  if(is.null(outcome_models))stop("argument 'outcome_models' missing")
  if(is.null(missing_models))stop("argument 'missing_models' missing")
  if(is.null(response_indicator))stop("argument 'response_indicator' missing")
  if(is.null(imputations))stop(" 'imputations' must be a positive integer")
  if(!inherits(imputations, "integer"))
    stop(" 'imputations' must be a positive integer")
  if(length(imputations) != 1)stop(" 'imputations' must be a positive integer")
  if(imputations < 1)stop(" 'imputations' must be a positive integer")
  if(class(outcome_models) != "list") stop("outcome_models must be a list")
  if(class(missing_models) != "list") stop("missing_models must be a list")
  for(mod in outcome_models){
    if(!inherit(mod, "multinom"))
      stop("the outcome models must be of class 'multinom' ")
    }
  for(mod in missing_models){
    if(!inherit(mod, "glm"))
      stop("the missingness models must be of class 'glm' ")
  }
  if(class(data) != "data.frame")stop("the data must be of class 'data.frame' ")
  
  
  N <- nrow(data)
  y.name <- attr(attr(outcome_models[[1]]$terms, "factors"), "dimnames")[[1]][1]
  C <- length(unique(data[, y.name]))
  imputed_data <- list()
  
  for(i in 1:imputations){
    # Step 1: Bootstrap the original data set.
    indices <- sample(1:N, N, replace = TRUE)
    boot <- data[indices,]
    boot_Rind <- response_indicator[indices]
    # Check that the bootstrap does not miss a category. 
    while(length(unique(boot[boot_Rind == 1, y.name]))!=3){
      boot <- samp[sample(1:300, 300, replace = TRUE),]
    }
    
    # Step 2: outcome models
    # If two or more regression models were applied, regress the dependent 
    # variable on the predicted values and save the new predicted values. 
    predlist <- lapply(outcome_models, FUN = function(x){fitted(x)[,-1]})
    pred <- matrix(unlist(predlist), ncol = length(predlist)*(C - 1))
    temporary <- as.data.frame(cbind(boot[boot_Rind == 1, y.name], pred))
    l.pred <- ncol(pred)
    names(temporary) <- c(y.name, paste0("V", 1:l.pred))
    mod.o <- multinom(y.name~., temporary)
    pred.outcome <- fitted(mod)[,-1]
    # Note that these models are applied only to respondents.
    # Standardize the predicted values. We need to save these for later, that's 
    # why we do it by hand.
    mu1 <- colMeans(pred.outcome) 
    sd1 <- apply(pred.outcome, 2, function(x) sqrt(var(x))) 
    for(i in 1:ncol(pred.outcome)){
      pred.outcome[,i] <- (pred.outcome[,i] - mu1[i]) / sd1[i]
    }
    
    # Step 3: missingness models
    # If two or more regression models were applied, regress the dependent 
    # variable on the predicted values and save the new predicted values. 
    predlist <- lapply(missing_models, FUN = function(x){x$fitted.values})
    pred <- matrix(unlist(predlist), ncol = length(predlist))
    temporary2 <- as.data.frame(cbind(boot_Rind, pred))
    l.pred <- ncol(pred)
    names(temporary2) <- c(y.name, paste0("V", 1:l.pred))
    mod.m <- glm(y.name~., temporary2, family = binomial)
    pred.response <- mod$fitted.values
    # Note that these models are applied to all cases in the bootstrap sample.
    # Standardize the predicted values.
    mu2 <- mean(pred.response)
    sd2 <- sqrt(var(pred.response))
    pred.response <- (pred.response - mu2) / sd2
    
    # step 4: Calculate the predicted outcome probabilities and propensity 
    # scores for the nonrespondents in the original sample and standardize these
    # using the parameters obtained in steps 2 and 3.
    # Select all cases in original data set with missing Y.
    to_be_imputed <- data[response_indicators == 0, -y.name]
    # Use function predict() to get predicted values for these cases.
    # Predicted outcome probabilities:
    missing_predlist <- lapply(outcome_models, FUN = function(x){
      matrix(matrix(predict(x, newdata = to_be_imputed, type = "probs"), 
                    ncol = C)[, -1], ncol = C - 1)})
    missing_pred <- matrix(unlist(missing_predlist), 
                           ncol = length(missing_predlist)*(C - 1))
    names(missing_pred) <- names(temporary)[-1]
    missing_pred.outcome <- matrix(matrix(predict(mod.o, newdata = missing_pred, 
                                                  type = "probs"), 
                                          ncol = C)[, -1], ncol = C - 1)
    # Predicted response propensities:
    missing_predlist <- lapply(missing_models, 
                               FUN = function(x){x$fitted.values})
    missing_pred <- matrix(unlist(missing_predlist), 
                           ncol = length(missing_predlist))
    names(missing_pred) <- names(temporary2)[-1]
    missing_pred.response <- predict(mod.m, newdata = missing_pred, 
                                     type = "response")
    # Standardize the predicted values
    for(i in 1:ncol(missing_pred.outcome)){
      missing_pred.outcome[,i] <- (missing_pred.outcome[,i] - mu1[i]) / sd1[i]
    }
    missing_pred.response <- (missing_pred.response - mu2) / sd2
    
    # step 5: calculate a distance function to define the similarity between 
    # subject i with missing Y in  the original data set and subject j with 
    # observed Y in the bootstrap sample based on the predictive scores.
    # Create matrices with predicted values for subjects with missing Y in the 
    # original data set and subjects with observed Y in the bootstrap sample.
    pred.to_be_imputed <- cbind(missing_pred.outcome, missing_pred.response)
    pred.boot_obs <- cbind(pred.outcome, pred.response[boot_Rind==1])
    # save the indices of the bootstrap subjects with observed Y.
    indx <- which(boot_Rind == 1)
    # Apply distance function. Here, each column contains the distances from a 
    # subject i with missing Y in the original data set to all subjects j with 
    # observed Y in the bootstrap sample.
    distances <- apply(pred.to_be_imputed, 1, 
                       function(x) apply(pred.boot_obs, 1, distance, x))
    
    # step 6: For each incomplete subject in the original data set, define the 
    # imputation set as the 5 nearest neighbors and impute by randomly drawing 
    # from the imputation set.
    # Rank the distances.
    ranks <- apply(distances, 2, order)
    # For each subject (column) with missing Y in the original data set, get the 
    # indices of the subjects with observed Y in the bootstrap sample with the 
    # smallest distance.
    indices <- apply(ranks, 2, function(x) indx[x<=5])
    # Using the indices, get the list of donor values for each subject with 
    # missing Y in the original data set.
    donors <- apply(indices, 2, function(x) boot[x, y.name])
    # Randomly draw one donor per subject with missing Y in the original data 
    # set.
    new.values <- apply(donors, 2, function(x) sample(x,1))
    # Impute.
    imputed_set <- data
    imputed_set[response_indicators == 0, y.name] <- new.values
    imputed_data[[i]] <- imputed_set
  }
  
  return(imputed_data)
}


# Distance function that calculates distance between subject i with missing Y in
# the original data set and subject j with observed Y in the bootstrap sample.
distance <- function(pred.i, pred.j){
  omega <- 1 / length(pred.i)
  sqrt(sum(omega * (pred.i - pred.j)^2))
}

