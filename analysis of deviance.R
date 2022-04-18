library(nnet)
library(car)
set.seed(1052)

load("~/SSLBS/Thesis/Thesis/pop1.RData")
g <- replicate(500, {
  samp <- population[sample(1:50000, 300, replace = FALSE),]
  mod <- multinom(y ~ x1 + x2 + x3, samp)
  c(Anova(mod, 3)$`Pr(>Chisq)`, Anova(mod, 3)$`LR Chisq`)
})
p_g <- rowMeans(g[1:3,]) # average p-values
var_g <- apply(g[1:3,], 1, function(x)var(x)*(500-1)/500) # variances of the p-values
sigrate_g <- apply(g[1:3,], 1, function(x){sum(x<.05)/500}) # significance rate
LRT_g <- rowMeans(g[4:6,]) # average LRT statistics

load("~/SSLBS/Thesis/Thesis/pop2.RData")
f <- replicate(500, {
  samp <- population[sample(1:50000, 300, replace = FALSE),]
  mod <- multinom(y ~ x2^2 + x1:x3 + x1 + x2 + x3, samp)
  c(Anova(mod, 3)$`Pr(>Chisq)`, Anova(mod, 3)$`LR Chisq`)
})
p_f <- rowMeans(f[1:4,])
var_f <- apply(f[1:4,], 1, function(x)var(x)*(500-1)/500)
sigrate_f <- apply(f[1:4,], 1, function(x){sum(x<.05)/500})
LRT_f <- rowMeans(f[5:8,])

load("~/SSLBS/Thesis/Thesis/pop3.RData")
h <- replicate(500, {
  samp <- population[sample(1:50000, 300, replace = FALSE),]
  mod <- multinom(y ~ x1 + x2 + x3, samp)
  c(Anova(mod, 3)$`Pr(>Chisq)`, Anova(mod, 3)$`LR Chisq`)
})
p_h <- rowMeans(h[1:3,])
var_h <- apply(h[1:3,], 1, function(x)var(x)*(500-1)/500)
sigrate_h <- apply(h[1:3,], 1, function(x){sum(x<.05)/500})
LRT_h <- rowMeans(h[4:6,])

load("~/SSLBS/Thesis/Thesis/pop4.RData")
i <- replicate(500, {
  samp <- population[sample(1:50000, 300, replace = FALSE),]
  mod <- multinom(y ~ x2^2 + x1:x3 + x1 + x2 + x3, samp)
  c(Anova(mod, 3)$`Pr(>Chisq)`, Anova(mod, 3)$`LR Chisq`)
})
p_i <- rowMeans(i[1:4,])
var_i <- apply(i[1:4,], 1, function(x)var(x)*(500-1)/500)
sigrate_i <- apply(i[1:4,], 1, function(x){sum(x<.05)/500})
LRT_i <- rowMeans(i[5:8,])


