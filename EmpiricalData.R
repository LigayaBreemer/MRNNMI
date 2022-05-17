library(tidyverse)
library(nnet)
library(car)
library(boot)
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

# Frequencies

mean(my_sample$LFT1SEPT)
mean(my_sample$EBBPB8URENWERK)

table(my_sample$EBBHHBBURGST)
attr(my_sample$EBBHHBBURGST, "labels")
table(my_sample$OPLNIVSOI2016AGG2HB)
attr(my_sample$OPLNIVSOI2016AGG2HB, "labels")
table(my_sample$ETNGRP)
attr(my_sample$ETNGRP, "labels")


### Model selection: outcome model #############################################
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

### Model selection: missingness model #########################################

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

##### Imputation ###############################################################

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

#### Analysis ##################################################################
# Restart the R session here, else tidyverse will change the data frames into
# tibbles and the function unique() will output a tibble instead of a vector.

### Proportions ################################################################
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

### Standard Errors ############################################################

# SE using Rubin's rules with as input one ROW of the proportions matrix.
standard_error <- function(proportions, M, n){
  sqrt(1/M * sum(proportions*(1-proportions)/n) + (1 + 1/M) / (M - 1) * 
         sum((proportions - mean(proportions))^2))
}

# SE using Rubin's rules + bootstrap within variance with as input one ROW of
# the proportions matrix and the within variance of each imputed data set.
SE_bootW <- function(proportions, within_var, M){
  sqrt(1/M * sum(within_var) + (1 + 1/M) / (M - 1) * 
         sum((proportions - mean(proportions))^2))
}

# SE using Rubin's rules + bootstrap variances with as input the within variance
# of each imputed data set and the between variance.
SE_boot <- function(within_var, between_var, M){
  sqrt(1/M * sum(within_var) + (1 + 1/M) / (M - 1) * between_var)
}

# SE calculated using Rubin's rules, no bootstrap #
SE_aggregate <- apply(imputed_aggregates, 1, standard_error, 5, 500)
SE_men <- apply(imputed_mens, 1, standard_error, 5, 258)
SE_women <- apply(imputed_womens, 1, standard_error, 5, 242)


# SE calculated using Rubin's rules, bootstrap within-variance #

# Aggregate level
# First, the number of bootstrap repetitions is determined by increasing the
# number of repetitions until convergence.
set.seed(1105)
b200 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 200, stype = "i")
mean(b200$t)

b300 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 300, stype = "i")
mean(b300$t)

b400 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 400, stype = "i")
mean(b400$t)

b500 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 500, stype = "i")
mean(b500$t)

b600 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 600, stype = "i")
mean(b600$t)

b700 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 700, stype = "i")
mean(b700$t)

b800 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 800, stype = "i")
mean(b800$t)

b900 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 900, stype = "i")
mean(b900$t)

b1000 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 1000, stype = "i")
mean(b1000$t)

b1100 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 1100, stype = "i")
mean(b1100$t)

b1200 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 1200, stype = "i")
mean(b1200$t) # no longer accuracte up to 5 decimal points?

b1300 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 1300, stype = "i")
mean(b1300$t)

b1400 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 1400, stype = "i")
mean(b1400$t) # new stability?

b1500 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 1500, stype = "i")
mean(b1500$t)

b1600 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 1600, stype = "i")
mean(b1600$t)

b1700 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 1700, stype = "i")
mean(b1700$t)

b1800 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 1800, stype = "i")
mean(b1800$t)

b1900 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 1900, stype = "i")
mean(b1900$t)

b2000 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  c <- unique(new_dat[,7])[1]
  indicator <- new_dat[,7] == c
  weights <- new_dat[indicator, 6]
  p <- sum(weights)/sum(new_dat[,6])
  p*(1-p)/500
}, R = 2000, stype = "i")
mean(b2000$t)
# After R= 1500, all estimate are appr. .000364-.000366

# Now calculate the bootstrap estimates of the variances for each category and
# each imputed data set
B1 <- boot(imputed_data[[1]], function(x, i){
  new_dat <- x[i,]
  categories <- unique(x[, 7])
  sapply(categories, function(y){
    indicator <- new_dat[,7] == y
    weights <- new_dat[indicator, 6]
    p <- sum(weights)/sum(new_dat[,6])
    p*(1-p)/500
  })
}, R=1500, stype = "i")
within_var1 <- colMeans(B1$t)

B2 <- boot(imputed_data[[2]], function(x, i){
  new_dat <- x[i,]
  categories <- unique(x[, 7])
  sapply(categories, function(y){
    indicator <- new_dat[,7] == y
    weights <- new_dat[indicator, 6]
    p <- sum(weights)/sum(new_dat[,6])
    p*(1-p)/500
  })
}, R=1500, stype = "i")
within_var2 <- colMeans(B2$t)

B3 <- boot(imputed_data[[3]], function(x, i){
  new_dat <- x[i,]
  categories <- unique(x[, 7])
  sapply(categories, function(y){
    indicator <- new_dat[,7] == y
    weights <- new_dat[indicator, 6]
    p <- sum(weights)/sum(new_dat[,6])
    p*(1-p)/500
  })
}, R=1500, stype = "i")
within_var3 <- colMeans(B3$t)

B4 <- boot(imputed_data[[4]], function(x, i){
  new_dat <- x[i,]
  categories <- unique(x[, 7])
  sapply(categories, function(y){
    indicator <- new_dat[,7] == y
    weights <- new_dat[indicator, 6]
    p <- sum(weights)/sum(new_dat[,6])
    p*(1-p)/500
  })
}, R=1500, stype = "i")
within_var4 <- colMeans(B4$t)

B5 <- boot(imputed_data[[5]], function(x, i){
  new_dat <- x[i,]
  categories <- unique(x[, 7])
  sapply(categories, function(y){
    indicator <- new_dat[,7] == y
    weights <- new_dat[indicator, 6]
    p <- sum(weights)/sum(new_dat[,6])
    p*(1-p)/500
  })
}, R=1500, stype = "i")
within_var5 <- colMeans(B5$t)
within_var <- rbind(within_var1, within_var2, within_var3, within_var4, 
                    within_var5)
SEbootW_aggregate <- numeric(6)
for(i in 1:6){
  SEbootW_aggregate[i] <- SE_bootW(imputed_aggregates[i,], within_var[,i], 5)
}

# Domain level
# Men
set.seed(1605)
B1_men <- boot(imputed_dom1[[1]], function(x, i){
  new_dat <- x[i,]
  categories <- unique(x[, 7])
  sapply(categories, function(y){
    indicator <- new_dat[,7] == y
    weights <- new_dat[indicator, 6]
    p <- sum(weights)/sum(new_dat[,6])
    p*(1-p)/500
  })
}, R=1500, stype = "i")
within_var1_men <- colMeans(B1_men$t)

B2_men <- boot(imputed_dom1[[2]], function(x, i){
  new_dat <- x[i,]
  categories <- unique(x[, 7])
  sapply(categories, function(y){
    indicator <- new_dat[,7] == y
    weights <- new_dat[indicator, 6]
    p <- sum(weights)/sum(new_dat[,6])
    p*(1-p)/500
  })
}, R=1500, stype = "i")
within_var2_men <- colMeans(B2_men$t)

B3_men <- boot(imputed_dom1[[3]], function(x, i){
  new_dat <- x[i,]
  categories <- unique(x[, 7])
  sapply(categories, function(y){
    indicator <- new_dat[,7] == y
    weights <- new_dat[indicator, 6]
    p <- sum(weights)/sum(new_dat[,6])
    p*(1-p)/500
  })
}, R=1500, stype = "i")
within_var3_men <- colMeans(B3_men$t)

B4_men <- boot(imputed_dom1[[4]], function(x, i){
  new_dat <- x[i,]
  categories <- unique(x[, 7])
  sapply(categories, function(y){
    indicator <- new_dat[,7] == y
    weights <- new_dat[indicator, 6]
    p <- sum(weights)/sum(new_dat[,6])
    p*(1-p)/500
  })
}, R=1500, stype = "i")
within_var4_men <- colMeans(B4_men$t)

B5_men <- boot(imputed_dom1[[5]], function(x, i){
  new_dat <- x[i,]
  categories <- unique(x[, 7])
  sapply(categories, function(y){
    indicator <- new_dat[,7] == y
    weights <- new_dat[indicator, 6]
    p <- sum(weights)/sum(new_dat[,6])
    p*(1-p)/500
  })
}, R=1500, stype = "i")
within_var5_men <- colMeans(B5_men$t)
within_var_men <- rbind(within_var1_men, within_var2_men, within_var3_men, 
                        within_var4_men, within_var5_men)
SEbootW_men <- numeric(6)
for(i in 1:6){
  SEbootW_men[i] <- SE_bootW(imputed_mens[i,], within_var_men[,i], 5)
}

# Women
B1_women <- boot(imputed_dom2[[1]], function(x, i){
  new_dat <- x[i,]
  categories <- unique(x[, 7])
  sapply(categories, function(y){
    indicator <- new_dat[,7] == y
    weights <- new_dat[indicator, 6]
    p <- sum(weights)/sum(new_dat[,6])
    p*(1-p)/500
  })
}, R=1500, stype = "i")
within_var1_women <- colMeans(B1_women$t)

B2_women <- boot(imputed_dom2[[2]], function(x, i){
  new_dat <- x[i,]
  categories <- unique(x[, 7])
  sapply(categories, function(y){
    indicator <- new_dat[,7] == y
    weights <- new_dat[indicator, 6]
    p <- sum(weights)/sum(new_dat[,6])
    p*(1-p)/500
  })
}, R=1500, stype = "i")
within_var2_women <- colMeans(B2_women$t)

B3_women <- boot(imputed_dom2[[3]], function(x, i){
  new_dat <- x[i,]
  categories <- unique(x[, 7])
  sapply(categories, function(y){
    indicator <- new_dat[,7] == y
    weights <- new_dat[indicator, 6]
    p <- sum(weights)/sum(new_dat[,6])
    p*(1-p)/500
  })
}, R=1500, stype = "i")
within_var3_women <- colMeans(B3_women$t)

B4_women <- boot(imputed_dom2[[4]], function(x, i){
  new_dat <- x[i,]
  categories <- unique(x[, 7])
  sapply(categories, function(y){
    indicator <- new_dat[,7] == y
    weights <- new_dat[indicator, 6]
    p <- sum(weights)/sum(new_dat[,6])
    p*(1-p)/500
  })
}, R=1500, stype = "i")
within_var4_women <- colMeans(B4_women$t)

B5_women <- boot(imputed_dom2[[5]], function(x, i){
  new_dat <- x[i,]
  categories <- unique(x[, 7])
  sapply(categories, function(y){
    indicator <- new_dat[,7] == y
    weights <- new_dat[indicator, 6]
    p <- sum(weights)/sum(new_dat[,6])
    p*(1-p)/500
  })
}, R=1500, stype = "i")
within_var5_women <- colMeans(B5_women$t)
within_var_women <- rbind(within_var1_women, within_var2_women, within_var3_women, 
                        within_var4_women, within_var5_women)
SEbootW_women <- numeric(6)
for(i in 1:6){
  SEbootW_women[i] <- SE_bootW(imputed_womens[i,], within_var_women[,i], 5)
}

# SE using Rubin's rules with bootstrap variances.
# Aggregate level
# First, the number of bootstraps has to be determined again.
set.seed(1257)
b200 <- boot(imputed_aggregates[1,], function(x, i){sum((x[i] - mean(x[i]))^2)}, 
          R=200, stype = "i")
mean(b200$t)

b300 <- boot(imputed_aggregates[1,], function(x, i){sum((x[i] - mean(x[i]))^2)}, 
             R=300, stype = "i")
mean(b300$t)

b400 <- boot(imputed_aggregates[1,], function(x, i){sum((x[i] - mean(x[i]))^2)}, 
             R=400, stype = "i")
mean(b400$t)

b500 <- boot(imputed_aggregates[1,], function(x, i){sum((x[i] - mean(x[i]))^2)}, 
             R=500, stype = "i")
mean(b500$t)

b600 <- boot(imputed_aggregates[1,], function(x, i){sum((x[i] - mean(x[i]))^2)}, 
             R=600, stype = "i")
mean(b600$t)

b700 <- boot(imputed_aggregates[1,], function(x, i){sum((x[i] - mean(x[i]))^2)}, 
             R=700, stype = "i")
mean(b700$t)

b800 <- boot(imputed_aggregates[1,], function(x, i){sum((x[i] - mean(x[i]))^2)}, 
             R=800, stype = "i")
mean(b800$t)
# After 500 bootstraps the estimate seems to converge around .000340.

B <- boot(t(imputed_aggregates), function(x, i){
  new_dat <- x[i,]
  apply(new_dat, 2, function(y){
    sum((y - mean(y))^2)
  })
}, R=500, stype = "i")

between_var <- colMeans(B$t)

SEboot_aggregate <- numeric(6)
for(i in 1:6){
  SEboot_aggregate[i] <- SE_boot(within_var[,i], between_var[i], 5)
}

# Domain level
# Men

B_men <- boot(t(imputed_mens), function(x, i){
  new_dat <- x[i,]
  apply(new_dat, 2, function(y){
    sum((y - mean(y))^2)
  })
}, R=500, stype = "i")

between_var_men <- colMeans(B_men$t)

SEboot_men <- numeric(6)
for(i in 1:6){
  SEboot_men[i] <- SE_boot(within_var_men[,i], between_var_men[i], 5)
}

# Women

B_women <- boot(t(imputed_womens), function(x, i){
  new_dat <- x[i,]
  apply(new_dat, 2, function(y){
    sum((y - mean(y))^2)
  })
}, R=500, stype = "i")

between_var_women <- colMeans(B_women$t)

SEboot_women <- numeric(6)
for(i in 1:6){
  SEboot_women[i] <- SE_boot(within_var_women[,i], between_var_women[i], 5)
}

# Compare 

SE_aggregate
SEboot_aggregate

SE_men
SEboot_men

SE_women
SEboot_women
