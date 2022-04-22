library(tidyverse)
library(nnet)
library(car)
#library(haven)
#AD000091N01V01P2016ANAV1 <- read_sav("\\\\CBSP.NL\\Productie\\Secundair\\MPOnderzoek_SEC1\\Input\\AD000091N01V01P2016ANAV1.SAV")
#save(AD000091N01V01P2016ANAV1, file = "\\\\CBSP.NL\\Productie\\Secundair\\MPOnderzoek_SEC1\\Werk\\Ligaya Breemer\\EBB2016")
load("\\\\CBSP.NL\\Productie\\Secundair\\MPOnderzoek_SEC1\\Werk\\Ligaya Breemer\\EBB2016")
set.seed(0604)

summary(AD000091N01V01P2016ANAV1)

# variables I want to use are: 
# Age (29)
# Gender (27)
# Marital status (26)
# Education (52 or 54)
# Maybe average working hours (43)
# Ethnicity (42)
# Survey weight (21)

my_data <- AD000091N01V01P2016ANAV1[,c(29, 27, 26, 42, 43, 21, 54)]

# Remove people under the age of 15 and above the age of 90.
my_data <- filter(my_data, LFT1SEPT >= 15 & LFT1SEPT <= 90)
# Remove people with missing values.
my_data <- na.omit(my_data)

# Draw a SRS of size 500.
my_sample <- my_data[sample(1:nrow(my_data), 500, replace = FALSE),]
original_sample <- my_sample
save(original_sample, file = "\\\\CBSP.NL\\Productie\\Secundair\\MPOnderzoek_SEC1\\Werk\\Ligaya Breemer\\original_sample")

# Generate response probabilities.
p <- exp(log(.5) + .19*my_sample$LFT1SEPT - .15*my_sample$EBBPB8URENWERK + 
           0.0004*my_sample$LFT1SEPT*my_sample$EBBPB8URENWERK) / 
  (exp(log(.5) + .19*my_sample$LFT1SEPT - .15*my_sample$EBBPB8URENWERK +
         0.0004*my_sample$LFT1SEPT*my_sample$EBBPB8URENWERK) + 1)
summary(p)
hist(p)
# Generate response indicators.
my_sample$REDU <- rbinom(500, 1, p)
# Introduce missing values.
my_sample$OPLNIVSOI2016AGG2HB[my_sample$REDU == 0] <- NA
sum(my_sample$REDU==1)/500
save(my_sample, file = "\\\\CBSP.NL\\Productie\\Secundair\\MPOnderzoek_SEC1\\Werk\\Ligaya Breemer\\my_sample")
load("\\\\CBSP.NL\\Productie\\Secundair\\MPOnderzoek_SEC1\\Werk\\Ligaya Breemer\\my_sample")

# Model selection: outcome model
set.seed(1104)
mod1 <- multinom(OPLNIVSOI2016AGG2HB ~ LFT1SEPT + EBBHHBGESLACHT + EBBHHBBURGST + 
                   ETNGRP + EBBPB8URENWERK, 
                 weights = EBBGEWJAARGEWICHTA,
                 data = my_sample)
summary(mod1)
Anova(mod1)
Anova(mod1)$`Pr(>Chisq)`*5 # Bonferroni correction

# Interactions

mod2 <- multinom(OPLNIVSOI2016AGG2HB ~ LFT1SEPT + EBBHHBGESLACHT + EBBHHBBURGST  
                   + ETNGRP + EBBPB8URENWERK + LFT1SEPT:EBBHHBGESLACHT + 
                   LFT1SEPT:EBBHHBBURGST + LFT1SEPT:ETNGRP + 
                   LFT1SEPT:EBBPB8URENWERK + EBBHHBGESLACHT:EBBHHBBURGST + 
                   EBBHHBGESLACHT:ETNGRP + EBBHHBGESLACHT:EBBPB8URENWERK + 
                   EBBHHBBURGST:ETNGRP + EBBHHBBURGST:EBBPB8URENWERK + 
                   ETNGRP:EBBPB8URENWERK, 
                 weights = EBBGEWJAARGEWICHTA,
                 data = my_sample)


summary(mod2)
Anova(mod2)
Anova(mod2)$`Pr(>Chisq)`*15 # Bonferroni correction

mod3 <- multinom(OPLNIVSOI2016AGG2HB ~ LFT1SEPT + EBBHHBGESLACHT + EBBHHBBURGST  
                 + ETNGRP + EBBPB8URENWERK + LFT1SEPT:EBBPB8URENWERK + 
                   EBBHHBGESLACHT:EBBHHBBURGST + EBBHHBGESLACHT:ETNGRP + 
                   EBBHHBGESLACHT:EBBPB8URENWERK + EBBHHBBURGST:ETNGRP, 
                 weights = EBBGEWJAARGEWICHTA,
                 data = my_sample)


summary(mod3)
Anova(mod3)
Anova(mod3)$`Pr(>Chisq)`*10 # Bonferroni correction

mod4 <- step(mod2)
step(mod2, k = log(500*.844))
mod5 <- step(mod1, OPLNIVSOI2016AGG2HB ~ LFT1SEPT + EBBHHBGESLACHT + 
               EBBHHBBURGST + ETNGRP + EBBPB8URENWERK + LFT1SEPT:EBBHHBGESLACHT
             + LFT1SEPT:EBBHHBBURGST + LFT1SEPT:ETNGRP + LFT1SEPT:EBBPB8URENWERK
             + EBBHHBGESLACHT:EBBHHBBURGST + EBBHHBGESLACHT:ETNGRP + 
               EBBHHBGESLACHT:EBBPB8URENWERK + EBBHHBBURGST:ETNGRP + 
               EBBHHBBURGST:EBBPB8URENWERK + ETNGRP:EBBPB8URENWERK)
# Models 3, 4, and 5 are selected

# Model selection: missingness model

mod6 <- glm(REDU ~ LFT1SEPT + EBBHHBGESLACHT + EBBHHBBURGST + ETNGRP + 
              EBBPB8URENWERK + LFT1SEPT:EBBHHBGESLACHT + LFT1SEPT:EBBHHBBURGST +
              LFT1SEPT:ETNGRP + LFT1SEPT:EBBPB8URENWERK + 
              EBBHHBGESLACHT:EBBHHBBURGST + EBBHHBGESLACHT:ETNGRP + 
              EBBHHBGESLACHT:EBBPB8URENWERK + EBBHHBBURGST:ETNGRP + 
              EBBHHBBURGST:EBBPB8URENWERK + ETNGRP:EBBPB8URENWERK,
            family = binomial,
            data = my_sample,
            weights = EBBGEWJAARGEWICHTA)
summary(mod6)
Anova(mod6)
Anova(mod6)$`Pr(>Chisq)`*15 # Bonferroni correction

mod7 <- glm(REDU ~ LFT1SEPT + EBBHHBGESLACHT + EBBHHBBURGST + ETNGRP + 
              EBBPB8URENWERK + LFT1SEPT:EBBHHBGESLACHT + EBBHHBBURGST:ETNGRP + 
              ETNGRP:EBBPB8URENWERK,
            family = binomial,
            data = my_sample,
            weights = EBBGEWJAARGEWICHTA)
summary(mod7)
Anova(mod7)
Anova(mod7)$`Pr(>Chisq)`*8 # Bonferroni correction

mod8 <- step(mod6)
Anova(mod8)

mod9 <- glm(REDU ~ LFT1SEPT + EBBHHBGESLACHT + ETNGRP + EBBPB8URENWERK,
            family = binomial,
            data = my_sample,
            weights = EBBGEWJAARGEWICHTA)
summary(mod9)
Anova(mod9)

mod10 <- step(mod9, REDU ~ LFT1SEPT + EBBHHBGESLACHT + EBBHHBBURGST + ETNGRP + 
                EBBPB8URENWERK + LFT1SEPT:EBBHHBGESLACHT + LFT1SEPT:EBBHHBBURGST +
                LFT1SEPT:ETNGRP + LFT1SEPT:EBBPB8URENWERK + 
                EBBHHBGESLACHT:EBBHHBBURGST + EBBHHBGESLACHT:ETNGRP + 
                EBBHHBGESLACHT:EBBPB8URENWERK + EBBHHBBURGST:ETNGRP + 
                EBBHHBBURGST:EBBPB8URENWERK + ETNGRP:EBBPB8URENWERK)
mod10$call
mod10$aic

mod11 <- step(mod8, REDU ~ LFT1SEPT + EBBHHBGESLACHT + EBBHHBBURGST + ETNGRP + 
                EBBPB8URENWERK + LFT1SEPT:EBBHHBGESLACHT + LFT1SEPT:EBBHHBBURGST +
                LFT1SEPT:ETNGRP + LFT1SEPT:EBBPB8URENWERK + 
                EBBHHBGESLACHT:EBBHHBBURGST + EBBHHBGESLACHT:ETNGRP + 
                EBBHHBGESLACHT:EBBPB8URENWERK + EBBHHBBURGST:ETNGRP + 
                EBBHHBBURGST:EBBPB8URENWERK + ETNGRP:EBBPB8URENWERK)
mod11$call
mod11$aic
# Models 7 and 11 are selected for imputation

# Imputation

# function to calculate the 95% confidence intervals.
# each row contains the lower and upper bound for one outcome category.
CI <- function(proportion, SE){
  lower <- proportion - qnorm(.975)*SE
  upper <- proportion + qnorm(.975)*SE
  return(c(lower, upper))
}

# distance function that calculates the distance between subject i with missing
# Y in the original data set and subject j with observed Y in the bootstrap 
# sample.
distance <- function(pred.i, pred.j){
  omega <- 1 / length(pred.i)
  sqrt(sum(omega * (pred.i - pred.j)^2))
}

# SE using Rubin's rules with as input one ROW of the proportions matrix.
standard_error <- function(proportions, M, n){
  sqrt(1/M * sum(proportions*(1 - proportions)/n) + (1 + 1/M) / (M - 1) * 
         sum((proportions - mean(proportions))^2))
}

MRNNMI <- function(iterations, samp){
  # step 1: draw a boostrap sample
  boot <- samp[sample(1:500, 500, replace = TRUE),]
  # check that the target variable still has the same number of categories
  temp <- length(unique(samp[samp$REDU == 1, 7]))
  while(length(unique(boot[boot$REDU == 1, 7])) != temp){
    boot <- samp[sample(1:500, 500, replace = TRUE),]
  }
  
  # step 2: using the bootstrap sample, fit one or more multinomial regression
  # models regressing the target variable on the auxiliary variables and save 
  # the standardized predicted values.
  mod1.1 <- multinom(OPLNIVSOI2016AGG2HB ~ LFT1SEPT + EBBHHBGESLACHT + 
                       EBBHHBBURGST + ETNGRP + EBBPB8URENWERK + 
                       LFT1SEPT:EBBPB8URENWERK + EBBHHBGESLACHT:EBBHHBBURGST + 
                       EBBHHBGESLACHT:ETNGRP + EBBHHBGESLACHT:EBBPB8URENWERK + 
                       EBBHHBBURGST:ETNGRP, 
                     weights = EBBGEWJAARGEWICHTA,
                     data = boot)
  pred1.1 <- fitted(mod1.1)[,-1]
  mod1.2 <- multinom(formula = OPLNIVSOI2016AGG2HB ~ LFT1SEPT + EBBHHBGESLACHT + 
                       EBBHHBBURGST + ETNGRP + EBBPB8URENWERK + 
                       LFT1SEPT:EBBHHBGESLACHT + LFT1SEPT:EBBPB8URENWERK + 
                       EBBHHBGESLACHT:EBBHHBBURGST + EBBHHBGESLACHT:ETNGRP + 
                       EBBHHBGESLACHT:EBBPB8URENWERK + EBBHHBBURGST:ETNGRP + 
                       EBBHHBBURGST:EBBPB8URENWERK + ETNGRP:EBBPB8URENWERK, 
                     data = boot, 
                     weights = EBBGEWJAARGEWICHTA)
  pred1.2 <- fitted(mod1.2)[,-1]
  mod1.3 <- multinom(formula = OPLNIVSOI2016AGG2HB ~ LFT1SEPT + EBBHHBGESLACHT + 
                       EBBHHBBURGST + ETNGRP + EBBPB8URENWERK + LFT1SEPT:ETNGRP 
                     + EBBHHBBURGST:ETNGRP + EBBHHBGESLACHT:EBBPB8URENWERK + 
                       EBBHHBGESLACHT:EBBHHBBURGST + EBBHHBBURGST:EBBPB8URENWERK
                     + EBBHHBGESLACHT:ETNGRP + LFT1SEPT:EBBPB8URENWERK, 
                     data = boot, 
                     weights = EBBGEWJAARGEWICHTA)
  pred1.3 <- fitted(mod1.3)[,-1]
  # regress the dependent variable on the predicted values and save the new 
  # predicted values.
  temporary <- cbind(boot[boot$REDU == 1, 7], pred1.1, 
                     pred1.2, pred1.3)
  temporary <- as.data.frame(temporary)
  names(temporary) <- c("y", "M1C2", "M1C3", "M1C4", "M1C5", "M1C6", "M2C2", 
                        "M2C3", "M2C4", "M2C5", "M2C6", "M3C2", "M3C3", "M3C4",
                        "M3C5", "M3C6")
  mod1.4 <- multinom(y ~ M1C2 + M1C3 + M1C4 + M1C5 + M1C6 + M2C2 + M2C3 + M2C4 +
                       M2C5 + M2C6 + M2C6 + M3C2 + M3C3 + M3C4 + M3C5 + M3C6,
                     data = temporary)
  pred1.4 <- fitted(mod1.3)[,-1]
  # standardize the predicted values.
  mu1 <- colMeans(pred1.4)
  sd1 <- apply(pred1.4, 2, function(x) sqrt(var(x)))
  for(i in 1:ncol(pred1.4)){
    pred1.4[,i] <- (pred1.4[,i] - mu1[i]) / sd1[i]
  }
  # step 3: using the boostrap sample, fiiiiit one or mroe logistic regression 
  # models regressing the respons indicator on the auxiliary variables and save
  # the standardized predicted values.
  mod2.1 <- glm(REDU ~ LFT1SEPT + EBBHHBGESLACHT + EBBHHBBURGST + ETNGRP + 
                  EBBPB8URENWERK + LFT1SEPT:EBBHHBGESLACHT + EBBHHBBURGST:ETNGRP
                + ETNGRP:EBBPB8URENWERK,
                family = binomial,
                data = boot,
                weights = EBBGEWJAARGEWICHTA)
  pred2.1 <- mod2.1$fitted.values
  mod2.2 <- glm(formula = REDU ~ LFT1SEPT + EBBHHBGESLACHT + EBBHHBBURGST + 
                  ETNGRP + EBBPB8URENWERK + EBBHHBBURGST:ETNGRP + 
                  ETNGRP:EBBPB8URENWERK + EBBHHBGESLACHT:ETNGRP + 
                  LFT1SEPT:ETNGRP + EBBHHBBURGST:EBBPB8URENWERK + 
                  LFT1SEPT:EBBPB8URENWERK + LFT1SEPT:EBBHHBGESLACHT + 
                  EBBHHBGESLACHT:EBBHHBBURGST + EBBHHBGESLACHT:EBBPB8URENWERK, 
                family = binomial, 
                data = boot, 
                weights = EBBGEWJAARGEWICHTA)
  pred2.2 <- mod2.2$fitted.values
  # regress the response indicator on the predicted values and save the new 
  # predicted values.
  temporary <- cbind(boot$REDU, pred2.1, 
                     pred2.2)
  temporary <- as.data.frame(temporary)
  names(temporary) <- c("R", "p1", "p2")
  mod2.3 <- glm(R ~ p1 + p2,
                family = binomial,
                data = temporary)
  pred2.3 <- mod2.3$fitted.values
  # standardize the predicted values.
  mu2 <- mean(pred2.3); sd2 <- sqrt(var(pred2.3))
  pred2.3 <- (pred2.3 - mu2) / sd2
  
  # step 4: calculate the predicted outcome probabilities and propensity scores
  # for the nonrespondents in the original sample and standardize using the 
  # parameters obtained in steps 2 and 3.
  # select all nonrespondents in the original data set.
  to_be_imputed <- samp[samp$REDU == 0, -7]
  # use predict() to get the predicted values.
  missing_pred1.1 <- predict(mod1.1, newdata = to_be_imputed, 
                             type = "probs")[, -1]
  missing_pred1.2 <- predict(mod1.2, newdata = to_be_imputed,
                             type = "probs")[, -1]
  missing_pred1.3 <- predict(mod1.3, newdata = to_be_imputed,
                             type = "probs")[, -1]
  missing_temp <- as.data.frame(cbind(missing_pred1.1, missing_pred1.2, 
                                      missing_pred1.3))
  names(missing_temp) <- c("M1C2", "M1C3", "M1C4", "M1C5", "M1C6", "M2C2", 
                           "M2C3", "M2C4", "M2C5", "M2C6", "M3C2", "M3C3",
                           "M3C4", "M3C5", "M3C6")
  missing_pred1.4 <- predict(mod1.4, newdata = missing_temp, 
                             type = "probs")[, -1]
  missing_pred2.1 <- predict(mod2.1, newdata = to_be_imputed)
  missing_pred2.2 <- predict(mod2.2, newdata = to_be_imputed)
  missing_temp <- as.data.frame(cbind(missing_pred2.1, missing_pred2.2))
  names(missing_temp) <- c("p1", "p2")
  missing_pred2.3 <- predict(mod2.3, newdata = missing_temp)
  # standardize the predicted values.
  for(i in 1:ncol(missing_pred1.4)){
    missing_pred1.4[,i] <- (missing_pred1.4[,i] - mu1[i]) / sd1[i]
  }
  missing_pred2.3 <- (missing_pred2.3 - mu2) / sd2
  
  # step 5: calculate the distances between the nonrespondents i in the original
  # sample and the respondents j in the bootstrap sample.
  pred.to_be_imputed <- cbind(missing_pred1.4, missing_pred2.3)
  pred.bootstrap <- cbind(pred1.4, pred2.3[boot$REDU == 1])
  # save the indices of the boostrap respondents.
  indx <- which(boot$REDU == 1)
  # apply distance function. Here, each column contains the distances from a
  # subject i with misisng Y in the original sample to subject j with
  # observed Y in the boostrap sample.
  distances <- apply(pred.to_be_imputed, 1, function(x) apply(pred.bootstrap, 1, 
                                                              distance, x))
  # step 6: for each nonrespondent in the original sample, define the imputation
  # set as the k-nearest neighborhood and impute by randomly drawing a donor.
  # rank the distances
  ranks <- apply(distances, 2, order)
  # for each nonrespondent in the original data set (column), get the indices of
  # the bootstrap respondents with the smallest distance to them.
  indices <- apply(ranks, 2, function(x) indx[x <= 20])
  # get the list of possible donor values for each nonrespondent in the original
  # sample.
  donors <- apply(indices, 2, function(x) boot[x, 7])
  # randomly draw one donor per nonrespondent.
  new.values <- unlist(lapply(donors, function(x)sample(unlist(x), 1)))
  # impute.
  imputed_set <- samp
  imputed_set[imputed_set$REDU == 0, 7] <- new.values
  
  return(imputed_set)
}

set.seed(1229)
imputed_data <- lapply(1:5, MRNNMI, my_sample)
save(imputed_data, file = "\\\\CBSP.NL\\Productie\\Secundair\\MPOnderzoek_SEC1\\Werk\\Ligaya Breemer\\imputed_data")

# Analysis
# Restart the R session here, else tidyverse will change the data frames into
# tibbles and the function unique() will output a tibble instead of a vector.
# Calculate the aggregate-level proportions for the original sample
categories <- unique(original_sample[, 7])
original_aggregate <- sapply(categories, FUN = function(x){
  indicator <- original_sample[,7] == x
  weights <- original_sample[indicator, 6]
  sum(weights)/sum(original_sample[,6])
})

# Calculate the aggregate-level proportions for the imputed data sets

imputed_aggregates <- sapply(imputed_data, function(x){
  categories <- unique(x[, 7])
  sapply(categories, function(y){
    indicator <- x[,7] == y
    weights <- x[indicator, 6]
    sum(weights)/sum(x[,6])
  })
})
imputed_aggregate <- rowMeans(imputed_aggregates)

rbind(original_aggregate, imputed_aggregate)
sum(abs(original_aggregate-imputed_aggregate))/6


# Calculate the domain-level proportions for the original sample.

# Filter domains
original_dom1 <- original_sample[original_sample$EBBHHBGESLACHT == "1",]
original_dom2 <- original_sample[original_sample$EBBHHBGESLACHT == "2",]
# Domain 1
categories <- unique(original_dom1[, 7])
original_men <- sapply(categories, FUN = function(x){
  indicator <- original_dom1[,7] == x
  weights <- original_dom1[indicator, 6]
  sum(weights)/sum(original_dom1[,6])
})
# Domain 2
categories <- unique(original_dom2[, 7])
original_women <- sapply(categories, FUN = function(x){
  indicator <- original_dom2[,7] == x
  weights <- original_dom2[indicator, 6]
  sum(weights)/sum(original_dom2[,6])
})

# Calculate the domain-level proportions for the imputed data.

# Filter domains.
imputed_dom1 <- lapply(imputed_data, function(x) x[x$EBBHHBGESLACHT == "1",])
imputed_dom2 <- lapply(imputed_data, function(x) x[x$EBBHHBGESLACHT == "2",])
# Domain 1 (category 99 is empty)
imputed_mens <- sapply(imputed_dom1, function(x){
  categories <- unique(x[, 7])
  sapply(categories, FUN = function(y){
    indicator <- x[,7] == y
    weights <- x[indicator, 6]
    sum(weights)/sum(x[,6])
  })
})
imputed_men <- rowMeans(imputed_mens)
imputed_men <- c(imputed_men, 0)
# Domain 2 
imputed_womens <- sapply(imputed_dom2, function(x){
  categories <- unique(x[, 7])
  sapply(categories, FUN = function(y){
    indicator <- x[,7] == y
    weights <- x[indicator, 6]
    sum(weights)/sum(x[,6])
  })
})
imputed_women <- rowMeans(imputed_womens)

rbind(original_men, imputed_men, original_women, imputed_women)
sum(abs(original_men - imputed_men))/6
sum(abs(original_women - imputed_women))/6

# SEs
SE_aggregate <- apply(imputed_aggregates, 1, standard_error, 5, 500)
SE_men <- apply(imputed_mens, 1, standard_error, 5, 258)
SE_women <- apply(imputed_womens, 1, standard_error, 5, 242)
