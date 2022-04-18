# For different populations simply change the 1 in "pop1res**" and "pop1" and 
# "population1_response**" to the number of the population.
# Also note that nrow = 3 should be changed to nrow = 5 for pop 3 and 4 and that the CRs
# should be calculated by dividing by 2500 instead of 1500.

### POPULATION 4 RESPONSE RATE 80%
load("~/SSLBS/Thesis/Thesis/S2pop4res80.RData")
load("~/SSLBS/Thesis/Thesis/pop4.RData")

original_proportions <- table(population$y)/nrow(population)

### MRNNI ###
MRNN10_avg_props <- unlist(lapply(population4_response80, function(x){unlist(x[[1]][[1]])}))
MRNN10_avg_props <- t(matrix(MRNN10_avg_props, nrow = 5)) # every column is one outcome category
MRNN10_AEP <- colMeans(MRNN10_avg_props) # average estimated proportions
SE_MRNN10 <- unlist(lapply(population4_response80, function(x){unlist(x[[1]][[4]])}))
SE_MRNN10 <- t(matrix(SE_MRNN10, nrow = 5)) # every column is one outcome category
MRNN10_ASE <- colMeans(SE_MRNN10) # average standard errors
within_MRNN10 <- unlist(lapply(population4_response80, function(x){unlist(x[[1]][[7]])}))
MRNN10_CR <- sum(within_MRNN10)/2500 # coverage rate

MRNN01_avg_props <- unlist(lapply(population4_response80, function(x){unlist(x[[1]][[2]])}))
MRNN01_avg_props <- t(matrix(MRNN01_avg_props, nrow = 5)) # every column is one outcome category
MRNN01_AEP <- colMeans(MRNN01_avg_props) # average estimated proportions
SE_MRNN01 <- unlist(lapply(population4_response80, function(x){unlist(x[[1]][[5]])}))
SE_MRNN01 <- t(matrix(SE_MRNN01, nrow = 5)) # every column is one outcome category
MRNN01_ASE <- colMeans(SE_MRNN01) # average standard errors
within_MRNN01 <- unlist(lapply(population4_response80, function(x){unlist(x[[1]][[8]])}))
MRNN01_CR <- sum(within_MRNN01)/2500 # coverage rate

MRNN00_avg_props <- unlist(lapply(population4_response80, function(x){unlist(x[[1]][[3]])}))
MRNN00_avg_props <- t(matrix(MRNN00_avg_props, nrow = 5)) # every column is one outcome category
MRNN00_AEP <- colMeans(MRNN00_avg_props) # average estimated proportions
SE_MRNN00 <- unlist(lapply(population4_response80, function(x){unlist(x[[1]][[6]])}))
SE_MRNN00 <- t(matrix(SE_MRNN00, nrow = 5)) # every column is one outcome category
MRNN00_ASE <- colMeans(SE_MRNN00) # average standard errors
within_MRNN00 <- unlist(lapply(population4_response80, function(x){unlist(x[[1]][[9]])}))
MRNN00_CR <- sum(within_MRNN00)/2500 # coverage rate

### DRNNI ###
DRNN10_avg_props <- unlist(lapply(population4_response80, function(x){unlist(x[[2]][[1]])}))
DRNN10_avg_props <- t(matrix(DRNN10_avg_props, nrow = 5)) # every column is one outcome category
DRNN10_AEP <- colMeans(DRNN10_avg_props) # average estimated proportions
SE_DRNN10 <- unlist(lapply(population4_response80, function(x){unlist(x[[2]][[4]])}))
SE_DRNN10 <- t(matrix(SE_DRNN10, nrow = 5)) # every column is one outcome category
DRNN10_ASE <- colMeans(SE_DRNN10) # average standard errors
within_DRNN10 <- unlist(lapply(population4_response80, function(x){unlist(x[[2]][[7]])}))
DRNN10_CR <- sum(within_DRNN10)/2500 # coverage rate

DRNN01_avg_props <- unlist(lapply(population4_response80, function(x){unlist(x[[2]][[2]])}))
DRNN01_avg_props <- t(matrix(DRNN01_avg_props, nrow = 5)) # every column is one outcome category
DRNN01_AEP <- colMeans(DRNN01_avg_props) # average estimated proportions
SE_DRNN01 <- unlist(lapply(population4_response80, function(x){unlist(x[[2]][[5]])}))
SE_DRNN01 <- t(matrix(SE_DRNN01, nrow = 5)) # every column is one outcome category
DRNN01_ASE <- colMeans(SE_DRNN01) # average standard errors
within_DRNN01 <- unlist(lapply(population4_response80, function(x){unlist(x[[2]][[8]])}))
DRNN01_CR <- sum(within_DRNN01)/2500 # coverage rate

DRNN00_avg_props <- unlist(lapply(population4_response80, function(x){unlist(x[[2]][[3]])}))
DRNN00_avg_props <- t(matrix(DRNN00_avg_props, nrow = 5)) # every column is one outcome category
DRNN00_AEP <- colMeans(DRNN00_avg_props) # average estimated proportions
SE_DRNN00 <- unlist(lapply(population4_response80, function(x){unlist(x[[2]][[6]])}))
SE_DRNN00 <- t(matrix(SE_DRNN00, nrow = 5)) # every column is one outcome category
DRNN00_ASE <- colMeans(SE_DRNN00) # average standard errors
within_DRNN00 <- unlist(lapply(population4_response80, function(x){unlist(x[[2]][[9]])}))
DRNN00_CR <- sum(within_DRNN00)/2500 # coverage rate

### MICE ###

MICE0_avg_props <- unlist(lapply(population4_response80, function(x){unlist(x[[3]][[1]])}))
MICE0_avg_props <- t(matrix(MICE0_avg_props, nrow = 5)) # every column is one outcome category
MICE0_AEP <- colMeans(MICE0_avg_props) # average estimated proportions
SE_MICE0 <- unlist(lapply(population4_response80, function(x){unlist(x[[3]][[2]])}))
SE_MICE0 <- t(matrix(SE_MICE0, nrow = 5)) # every column is one outcome category
MICE0_ASE <- colMeans(SE_MICE0) # average standard errors
within_MICE0 <- unlist(lapply(population4_response80, function(x){unlist(x[[3]][[3]])}))
MICE0_CR <- sum(within_MICE0)/2500 # coverage rate

results_S2pop4res80 <- list(original_proportions, 
                          MRNN10_AEP, MRNN10_ASE, MRNN10_CR,
                          MRNN01_AEP, MRNN01_ASE, MRNN01_CR,
                          MRNN00_AEP, MRNN00_ASE, MRNN00_CR,
                          DRNN10_AEP, DRNN10_ASE, DRNN10_CR,
                          DRNN01_AEP, DRNN01_ASE, DRNN01_CR,
                          DRNN00_AEP, DRNN00_ASE, DRNN00_CR,
                          MICE0_AEP, MICE0_ASE, MICE0_CR)
save(results_S2pop4res80, file = "results_S2pop4res80.RData")

### POPULATION 4 RESPONSE RATE 90%
load("~/SSLBS/Thesis/Thesis/S2pop4res90.RData")

### MRNNI ###
MRNN10_avg_props <- unlist(lapply(population4_response90, function(x){unlist(x[[1]][[1]])}))
MRNN10_avg_props <- t(matrix(MRNN10_avg_props, nrow = 5)) # every column is one outcome category
MRNN10_AEP <- colMeans(MRNN10_avg_props) # average estimated proportions
SE_MRNN10 <- unlist(lapply(population4_response90, function(x){unlist(x[[1]][[4]])}))
SE_MRNN10 <- t(matrix(SE_MRNN10, nrow = 5)) # every column is one outcome category
MRNN10_ASE <- colMeans(SE_MRNN10) # average standard errors
within_MRNN10 <- unlist(lapply(population4_response90, function(x){unlist(x[[1]][[7]])}))
MRNN10_CR <- sum(within_MRNN10)/2500 # coverage rate

MRNN01_avg_props <- unlist(lapply(population4_response90, function(x){unlist(x[[1]][[2]])}))
MRNN01_avg_props <- t(matrix(MRNN01_avg_props, nrow = 5)) # every column is one outcome category
MRNN01_AEP <- colMeans(MRNN01_avg_props) # average estimated proportions
SE_MRNN01 <- unlist(lapply(population4_response90, function(x){unlist(x[[1]][[5]])}))
SE_MRNN01 <- t(matrix(SE_MRNN01, nrow = 5)) # every column is one outcome category
MRNN01_ASE <- colMeans(SE_MRNN01) # average standard errors
within_MRNN01 <- unlist(lapply(population4_response90, function(x){unlist(x[[1]][[8]])}))
MRNN01_CR <- sum(within_MRNN01)/2500 # coverage rate

MRNN00_avg_props <- unlist(lapply(population4_response90, function(x){unlist(x[[1]][[3]])}))
MRNN00_avg_props <- t(matrix(MRNN00_avg_props, nrow = 5)) # every column is one outcome category
MRNN00_AEP <- colMeans(MRNN00_avg_props) # average estimated proportions
SE_MRNN00 <- unlist(lapply(population4_response90, function(x){unlist(x[[1]][[6]])}))
SE_MRNN00 <- t(matrix(SE_MRNN00, nrow = 5)) # every column is one outcome category
MRNN00_ASE <- colMeans(SE_MRNN00) # average standard errors
within_MRNN00 <- unlist(lapply(population4_response90, function(x){unlist(x[[1]][[9]])}))
MRNN00_CR <- sum(within_MRNN00)/2500 # coverage rate

### DRNNI ###
DRNN10_avg_props <- unlist(lapply(population4_response90, function(x){unlist(x[[2]][[1]])}))
DRNN10_avg_props <- t(matrix(DRNN10_avg_props, nrow = 5)) # every column is one outcome category
DRNN10_AEP <- colMeans(DRNN10_avg_props) # average estimated proportions
SE_DRNN10 <- unlist(lapply(population4_response90, function(x){unlist(x[[2]][[4]])}))
SE_DRNN10 <- t(matrix(SE_DRNN10, nrow = 5)) # every column is one outcome category
DRNN10_ASE <- colMeans(SE_DRNN10) # average standard errors
within_DRNN10 <- unlist(lapply(population4_response90, function(x){unlist(x[[2]][[7]])}))
DRNN10_CR <- sum(within_DRNN10)/2500 # coverage rate

DRNN01_avg_props <- unlist(lapply(population4_response90, function(x){unlist(x[[2]][[2]])}))
DRNN01_avg_props <- t(matrix(DRNN01_avg_props, nrow = 5)) # every column is one outcome category
DRNN01_AEP <- colMeans(DRNN01_avg_props) # average estimated proportions
SE_DRNN01 <- unlist(lapply(population4_response90, function(x){unlist(x[[2]][[5]])}))
SE_DRNN01 <- t(matrix(SE_DRNN01, nrow = 5)) # every column is one outcome category
DRNN01_ASE <- colMeans(SE_DRNN01) # average standard errors
within_DRNN01 <- unlist(lapply(population4_response90, function(x){unlist(x[[2]][[8]])}))
DRNN01_CR <- sum(within_DRNN01)/2500 # coverage rate

DRNN00_avg_props <- unlist(lapply(population4_response90, function(x){unlist(x[[2]][[3]])}))
DRNN00_avg_props <- t(matrix(DRNN00_avg_props, nrow = 5)) # every column is one outcome category
DRNN00_AEP <- colMeans(DRNN00_avg_props) # average estimated proportions
SE_DRNN00 <- unlist(lapply(population4_response90, function(x){unlist(x[[2]][[6]])}))
SE_DRNN00 <- t(matrix(SE_DRNN00, nrow = 5)) # every column is one outcome category
DRNN00_ASE <- colMeans(SE_DRNN00) # average standard errors
within_DRNN00 <- unlist(lapply(population4_response90, function(x){unlist(x[[2]][[9]])}))
DRNN00_CR <- sum(within_DRNN00)/2500 # coverage rate

### MICE ###

MICE0_avg_props <- unlist(lapply(population4_response90, function(x){unlist(x[[3]][[1]])}))
MICE0_avg_props <- t(matrix(MICE0_avg_props, nrow = 5)) # every column is one outcome category
MICE0_AEP <- colMeans(MICE0_avg_props) # average estimated proportions
SE_MICE0 <- unlist(lapply(population4_response90, function(x){unlist(x[[3]][[2]])}))
SE_MICE0 <- t(matrix(SE_MICE0, nrow = 5)) # every column is one outcome category
MICE0_ASE <- colMeans(SE_MICE0) # average standard errors
within_MICE0 <- unlist(lapply(population4_response90, function(x){unlist(x[[3]][[3]])}))
MICE0_CR <- sum(within_MICE0)/2500 # coverage rate

results_S2pop4res90 <- list(original_proportions, 
                            MRNN10_AEP, MRNN10_ASE, MRNN10_CR,
                            MRNN01_AEP, MRNN01_ASE, MRNN01_CR,
                            MRNN00_AEP, MRNN00_ASE, MRNN00_CR,
                            DRNN10_AEP, DRNN10_ASE, DRNN10_CR,
                            DRNN01_AEP, DRNN01_ASE, DRNN01_CR,
                            DRNN00_AEP, DRNN00_ASE, DRNN00_CR,
                            MICE0_AEP, MICE0_ASE, MICE0_CR)
save(results_S2pop4res90, file = "results_S2pop4res90.RData")