# For different populations simply change the 1 in "pop1res**" and "pop1" and 
# "population1_response**" to the number of the population.
# Also note that nrow = 3 should be changed to nrow = 5 for pop 3 and 4 and that the CRs
# should be calculated by dividing by 2500 instead of 1500.

### POPULATION 1 RESPONSE RATE 80%
load("~/SSLBS/Thesis/Thesis/pop1res80.RData")
load("~/SSLBS/Thesis/Thesis/pop1.RData")

original_proportions <- table(population$y)/nrow(population)

# function for calculating confidence intervals of the coverage rates
CI <- function(x) c(x - 1.96*sqrt(x*(1-x)/1500), x + 1.96*sqrt(x*(1-x)/1500))

### MRNNI ###
MRNNI5_avg_props <- unlist(lapply(population1_response80, function(x){unlist(x[[1]][[1]])}))
MRNNI5_avg_props <- t(matrix(MRNNI5_avg_props, nrow = 3)) # every column is one outcome category
MRNNI5_AEP <- colMeans(MRNNI5_avg_props) # average estimated proportions
SE_MRNNI5 <- unlist(lapply(population1_response80, function(x){unlist(x[[1]][[4]])}))
SE_MRNNI5 <- t(matrix(SE_MRNNI5, nrow = 3)) # every column is one outcome category
MRNNI5_ASE <- colMeans(SE_MRNNI5) # average standard errors
within_MRNNI5 <- unlist(lapply(population1_response80, function(x){unlist(x[[1]][[7]])}))
MRNNI5_CR <- sum(within_MRNNI5)/1500 # coverage rate
MRNNI5_CI <- CI(MRNNI5_CR)

MRNNI10_avg_props <- unlist(lapply(population1_response80, function(x){unlist(x[[1]][[2]])}))
MRNNI10_avg_props <- t(matrix(MRNNI10_avg_props, nrow = 3)) # every column is one outcome category
MRNNI10_AEP <- colMeans(MRNNI10_avg_props) # average estimated proportions
SE_MRNNI10 <- unlist(lapply(population1_response80, function(x){unlist(x[[1]][[5]])}))
SE_MRNNI10 <- t(matrix(SE_MRNNI10, nrow = 3)) # every column is one outcome category
MRNNI10_ASE <- colMeans(SE_MRNNI10) # average standard errors
within_MRNNI10 <- unlist(lapply(population1_response80, function(x){unlist(x[[1]][[8]])}))
MRNNI10_CR <- sum(within_MRNNI10)/1500 # coverage rate
MRNNI10_CI <- CI(MRNNI10_CR)

MRNNI20_avg_props <- unlist(lapply(population1_response80, function(x){unlist(x[[1]][[3]])}))
MRNNI20_avg_props <- t(matrix(MRNNI20_avg_props, nrow = 3)) # every column is one outcome category
MRNNI20_AEP <- colMeans(MRNNI20_avg_props) # average estimated proportions
SE_MRNNI20 <- unlist(lapply(population1_response80, function(x){unlist(x[[1]][[6]])}))
SE_MRNNI20 <- t(matrix(SE_MRNNI20, nrow = 3)) # every column is one outcome category
MRNNI20_ASE <- colMeans(SE_MRNNI20) # average standard errors
within_MRNNI20 <- unlist(lapply(population1_response80, function(x){unlist(x[[1]][[9]])}))
MRNNI20_CR <- sum(within_MRNNI20)/1500 # coverage rate
MRNNI20_CI <- CI(MRNNI20_CR)

### DRNNI ###
DRNNI5_avg_props <- unlist(lapply(population1_response80, function(x){unlist(x[[2]][[1]])}))
DRNNI5_avg_props <- t(matrix(DRNNI5_avg_props, nrow = 3)) # every column is one outcome category
DRNNI5_AEP <- colMeans(DRNNI5_avg_props) # average estimated proportions
SE_DRNNI5 <- unlist(lapply(population1_response80, function(x){unlist(x[[2]][[4]])}))
SE_DRNNI5 <- t(matrix(SE_DRNNI5, nrow = 3)) # every column is one outcome category
DRNNI5_ASE <- colMeans(SE_DRNNI5) # average standard errors
within_DRNNI5 <- unlist(lapply(population1_response80, function(x){unlist(x[[2]][[7]])}))
DRNNI5_CR <- sum(within_DRNNI5)/1500 # coverage rate
DRNNI5_CI <- CI(DRNNI5_CR)

DRNNI10_avg_props <- unlist(lapply(population1_response80, function(x){unlist(x[[2]][[2]])}))
DRNNI10_avg_props <- t(matrix(DRNNI10_avg_props, nrow = 3)) # every column is one outcome category
DRNNI10_AEP <- colMeans(DRNNI10_avg_props) # average estimated proportions
SE_DRNNI10 <- unlist(lapply(population1_response80, function(x){unlist(x[[2]][[5]])}))
SE_DRNNI10 <- t(matrix(SE_DRNNI10, nrow = 3)) # every column is one outcome category
DRNNI10_ASE <- colMeans(SE_DRNNI10) # average standard errors
within_DRNNI10 <- unlist(lapply(population1_response80, function(x){unlist(x[[2]][[8]])}))
DRNNI10_CR <- sum(within_DRNNI10)/1500 # coverage rate
DRNNI10_CI <- CI(DRNNI10_CR)

DRNNI20_avg_props <- unlist(lapply(population1_response80, function(x){unlist(x[[2]][[3]])}))
DRNNI20_avg_props <- t(matrix(DRNNI20_avg_props, nrow = 3)) # every column is one outcome category
DRNNI20_AEP <- colMeans(DRNNI20_avg_props) # average estimated proportions
SE_DRNNI20 <- unlist(lapply(population1_response80, function(x){unlist(x[[2]][[6]])}))
SE_DRNNI20 <- t(matrix(SE_DRNNI20, nrow = 3)) # every column is one outcome category
DRNNI20_ASE <- colMeans(SE_DRNNI20) # average standard errors
within_DRNNI20 <- unlist(lapply(population1_response80, function(x){unlist(x[[2]][[9]])}))
DRNNI20_CR <- sum(within_DRNNI20)/1500 # coverage rate
DRNNI20_CI <- CI(DRNNI20_CR)

### MICE ###
MICE5_avg_props <- unlist(lapply(population1_response80, function(x){unlist(x[[3]][[1]])}))
MICE5_avg_props <- t(matrix(MICE5_avg_props, nrow = 3)) # every column is one outcome category
MICE5_AEP <- colMeans(MICE5_avg_props) # average estimated proportions
SE_MICE5 <- unlist(lapply(population1_response80, function(x){unlist(x[[3]][[4]])}))
SE_MICE5 <- t(matrix(SE_MICE5, nrow = 3)) # every column is one outcome category
MICE5_ASE <- colMeans(SE_MICE5) # average standard errors
within_MICE5 <- unlist(lapply(population1_response80, function(x){unlist(x[[3]][[7]])}))
MICE5_CR <- sum(within_MICE5)/1500 # coverage rate
MICE5_CI <- CI(MICE5_CR)

MICE10_avg_props <- unlist(lapply(population1_response80, function(x){unlist(x[[3]][[2]])}))
MICE10_avg_props <- t(matrix(MICE10_avg_props, nrow = 3)) # every column is one outcome category
MICE10_AEP <- colMeans(MICE10_avg_props) # average estimated proportions
SE_MICE10 <- unlist(lapply(population1_response80, function(x){unlist(x[[3]][[5]])}))
SE_MICE10 <- t(matrix(SE_MICE10, nrow = 3)) # every column is one outcome category
MICE10_ASE <- colMeans(SE_MICE10) # average standard errors
within_MICE10 <- unlist(lapply(population1_response80, function(x){unlist(x[[3]][[8]])}))
MICE10_CR <- sum(within_MICE10)/1500 # coverage rate
MICE10_CI <- CI(MICE10_CR)

MICE20_avg_props <- unlist(lapply(population1_response80, function(x){unlist(x[[3]][[3]])}))
MICE20_avg_props <- t(matrix(MICE20_avg_props, nrow = 3)) # every column is one outcome category
MICE20_AEP <- colMeans(MICE20_avg_props) # average estimated proportions
SE_MICE20 <- unlist(lapply(population1_response80, function(x){unlist(x[[3]][[6]])}))
SE_MICE20 <- t(matrix(SE_MICE20, nrow = 3)) # every column is one outcome category
MICE20_ASE <- colMeans(SE_MICE20) # average standard errors
within_MICE20 <- unlist(lapply(population1_response80, function(x){unlist(x[[3]][[9]])}))
MICE20_CR <- sum(within_MICE20)/1500 # coverage rate
MICE20_CI <- CI(MICE20_CR)

results_pop1res80 <- list(original_proportions, 
                          MRNNI5_AEP, MRNNI5_ASE, MRNNI5_CR, MRNNI5_CI,
                          MRNNI10_AEP, MRNNI10_ASE, MRNNI10_CR, MRNNI10_CI,
                          MRNNI20_AEP, MRNNI20_ASE, MRNNI20_CR, MRNNI20_CI,
                          DRNNI5_AEP, DRNNI5_ASE, DRNNI5_CR, DRNNI5_CI,
                          DRNNI10_AEP, DRNNI10_ASE, DRNNI10_CR, DRNNI10_CI,
                          DRNNI20_AEP, DRNNI20_ASE, DRNNI20_CR, DRNNI20_CI,
                          MICE5_AEP, MICE5_ASE, MICE5_CR, MICE5_CI,
                          MICE10_AEP, MICE10_ASE, MICE10_CR, MICE10_CI,
                          MICE20_AEP, MICE20_ASE, MICE20_CR, MICE20_CI)
save(results_pop1res80, file = "results_pop1res80.RData")

### POPULATION 1 RESPONSE RATE 90%
load("~/SSLBS/Thesis/Thesis/pop1res90.RData")

### MRNNI ###
MRNNI5_avg_props <- unlist(lapply(population1_response90, function(x){unlist(x[[1]][[1]])}))
MRNNI5_avg_props <- t(matrix(MRNNI5_avg_props, nrow = 3)) # every column is one outcome category
MRNNI5_AEP <- colMeans(MRNNI5_avg_props) # average estimated proportions
SE_MRNNI5 <- unlist(lapply(population1_response90, function(x){unlist(x[[1]][[4]])}))
SE_MRNNI5 <- t(matrix(SE_MRNNI5, nrow = 3)) # every column is one outcome category
MRNNI5_ASE <- colMeans(SE_MRNNI5) # average standard errors
within_MRNNI5 <- unlist(lapply(population1_response90, function(x){unlist(x[[1]][[7]])}))
MRNNI5_CR <- sum(within_MRNNI5)/1500 # coverage rate
MRNNI5_CI <- CI(MRNNI5_CR)

MRNNI10_avg_props <- unlist(lapply(population1_response90, function(x){unlist(x[[1]][[2]])}))
MRNNI10_avg_props <- t(matrix(MRNNI10_avg_props, nrow = 3)) # every column is one outcome category
MRNNI10_AEP <- colMeans(MRNNI10_avg_props) # average estimated proportions
SE_MRNNI10 <- unlist(lapply(population1_response90, function(x){unlist(x[[1]][[5]])}))
SE_MRNNI10 <- t(matrix(SE_MRNNI10, nrow = 3)) # every column is one outcome category
MRNNI10_ASE <- colMeans(SE_MRNNI10) # average standard errors
within_MRNNI10 <- unlist(lapply(population1_response90, function(x){unlist(x[[1]][[8]])}))
MRNNI10_CR <- sum(within_MRNNI10)/1500 # coverage rate
MRNNI10_CI <- CI(MRNNI10_CR)

MRNNI20_avg_props <- unlist(lapply(population1_response90, function(x){unlist(x[[1]][[3]])}))
MRNNI20_avg_props <- t(matrix(MRNNI20_avg_props, nrow = 3)) # every column is one outcome category
MRNNI20_AEP <- colMeans(MRNNI20_avg_props) # average estimated proportions
SE_MRNNI20 <- unlist(lapply(population1_response90, function(x){unlist(x[[1]][[6]])}))
SE_MRNNI20 <- t(matrix(SE_MRNNI20, nrow = 3)) # every column is one outcome category
MRNNI20_ASE <- colMeans(SE_MRNNI20) # average standard errors
within_MRNNI20 <- unlist(lapply(population1_response90, function(x){unlist(x[[1]][[9]])}))
MRNNI20_CR <- sum(within_MRNNI20)/1500 # coverage rate
MRNNI20_CI <- CI(MRNNI20_CR)

### DRNNI ###
DRNNI5_avg_props <- unlist(lapply(population1_response90, function(x){unlist(x[[2]][[1]])}))
DRNNI5_avg_props <- t(matrix(DRNNI5_avg_props, nrow = 3)) # every column is one outcome category
DRNNI5_AEP <- colMeans(DRNNI5_avg_props) # average estimated proportions
SE_DRNNI5 <- unlist(lapply(population1_response90, function(x){unlist(x[[2]][[4]])}))
SE_DRNNI5 <- t(matrix(SE_DRNNI5, nrow = 3)) # every column is one outcome category
DRNNI5_ASE <- colMeans(SE_DRNNI5) # average standard errors
within_DRNNI5 <- unlist(lapply(population1_response90, function(x){unlist(x[[2]][[7]])}))
DRNNI5_CR <- sum(within_DRNNI5)/1500 # coverage rate
DRNNI5_CI <- CI(DRNNI5_CR)

DRNNI10_avg_props <- unlist(lapply(population1_response90, function(x){unlist(x[[2]][[2]])}))
DRNNI10_avg_props <- t(matrix(DRNNI10_avg_props, nrow = 3)) # every column is one outcome category
DRNNI10_AEP <- colMeans(DRNNI10_avg_props) # average estimated proportions
SE_DRNNI10 <- unlist(lapply(population1_response90, function(x){unlist(x[[2]][[5]])}))
SE_DRNNI10 <- t(matrix(SE_DRNNI10, nrow = 3)) # every column is one outcome category
DRNNI10_ASE <- colMeans(SE_DRNNI10) # average standard errors
within_DRNNI10 <- unlist(lapply(population1_response90, function(x){unlist(x[[2]][[8]])}))
DRNNI10_CR <- sum(within_DRNNI10)/1500 # coverage rate
DRNNI10_CI <- CI(DRNNI10_CR)

DRNNI20_avg_props <- unlist(lapply(population1_response90, function(x){unlist(x[[2]][[3]])}))
DRNNI20_avg_props <- t(matrix(DRNNI20_avg_props, nrow = 3)) # every column is one outcome category
DRNNI20_AEP <- colMeans(DRNNI20_avg_props) # average estimated proportions
SE_DRNNI20 <- unlist(lapply(population1_response90, function(x){unlist(x[[2]][[6]])}))
SE_DRNNI20 <- t(matrix(SE_DRNNI20, nrow = 3)) # every column is one outcome category
DRNNI20_ASE <- colMeans(SE_DRNNI20) # average standard errors
within_DRNNI20 <- unlist(lapply(population1_response90, function(x){unlist(x[[2]][[9]])}))
DRNNI20_CR <- sum(within_DRNNI20)/1500 # coverage rate
DRNNI20_CI <- CI(DRNNI20_CR)

### MICE ###
MICE5_avg_props <- unlist(lapply(population1_response90, function(x){unlist(x[[3]][[1]])}))
MICE5_avg_props <- t(matrix(MICE5_avg_props, nrow = 3)) # every column is one outcome category
MICE5_AEP <- colMeans(MICE5_avg_props) # average estimated proportions
SE_MICE5 <- unlist(lapply(population1_response90, function(x){unlist(x[[3]][[4]])}))
SE_MICE5 <- t(matrix(SE_MICE5, nrow = 3)) # every column is one outcome category
MICE5_ASE <- colMeans(SE_MICE5) # average standard errors
within_MICE5 <- unlist(lapply(population1_response90, function(x){unlist(x[[3]][[7]])}))
MICE5_CR <- sum(within_MICE5)/1500 # coverage rate
MICE5_CI <- CI(MICE5_CR)

MICE10_avg_props <- unlist(lapply(population1_response90, function(x){unlist(x[[3]][[2]])}))
MICE10_avg_props <- t(matrix(MICE10_avg_props, nrow = 3)) # every column is one outcome category
MICE10_AEP <- colMeans(MICE10_avg_props) # average estimated proportions
SE_MICE10 <- unlist(lapply(population1_response90, function(x){unlist(x[[3]][[5]])}))
SE_MICE10 <- t(matrix(SE_MICE10, nrow = 3)) # every column is one outcome category
MICE10_ASE <- colMeans(SE_MICE10) # average standard errors
within_MICE10 <- unlist(lapply(population1_response90, function(x){unlist(x[[3]][[8]])}))
MICE10_CR <- sum(within_MICE10)/1500 # coverage rate
MICE10_CI <- CI(MICE10_CR)

MICE20_avg_props <- unlist(lapply(population1_response90, function(x){unlist(x[[3]][[3]])}))
MICE20_avg_props <- t(matrix(MICE20_avg_props, nrow = 3)) # every column is one outcome category
MICE20_AEP <- colMeans(MICE20_avg_props) # average estimated proportions
SE_MICE20 <- unlist(lapply(population1_response90, function(x){unlist(x[[3]][[6]])}))
SE_MICE20 <- t(matrix(SE_MICE20, nrow = 3)) # every column is one outcome category
MICE20_ASE <- colMeans(SE_MICE20) # average standard errors
within_MICE20 <- unlist(lapply(population1_response90, function(x){unlist(x[[3]][[9]])}))
MICE20_CR <- sum(within_MICE20)/1500 # coverage rate
MICE20_CI <- CI(MICE20_CR)

results_pop1res90 <- list(original_proportions, 
                          MRNNI5_AEP, MRNNI5_ASE, MRNNI5_CR, MRNNI5_CI,
                          MRNNI10_AEP, MRNNI10_ASE, MRNNI10_CR, MRNNI10_CI,
                          MRNNI20_AEP, MRNNI20_ASE, MRNNI20_CR, MRNNI20_CI,
                          DRNNI5_AEP, DRNNI5_ASE, DRNNI5_CR, DRNNI5_CI,
                          DRNNI10_AEP, DRNNI10_ASE, DRNNI10_CR, DRNNI10_CI,
                          DRNNI20_AEP, DRNNI20_ASE, DRNNI20_CR, DRNNI20_CI,
                          MICE5_AEP, MICE5_ASE, MICE5_CR, MICE5_CI,
                          MICE10_AEP, MICE10_ASE, MICE10_CR, MICE10_CI,
                          MICE20_AEP, MICE20_ASE, MICE20_CR, MICE20_CI)
save(results_pop1res90, file = "results_pop1res90.RData")
