# Simulate Categorical Data
# IRT framework
# 5 items, each with 4 categories
# Composite reliability ~ 0.75
# Higher proportion in higher categories

library(dplyr)
library(mirt)
library(lavaan)
library(semTools)
library(mnormt)

N = 1e5

set.seed(1234)

# Set higher composite reliability for continuous items
# Composite reliability = 0.84 = 0.84/(0.84 + 0.16)
# sum(lambda)^2 = 0.75
# sum(Error variance) = 0.25

# Simulate latent factor with N(0,1)
eta_scores <- rnorm(N)
# Simulate loadings
# Fixed at normal ogive: 1.7 * c(.6, .7, .8, .5, .3)
lambda <- 1.7 * c(.6, .7, .8, .5, .3)

# Simulate error variances
# Sum(error_var) = 8.101633
error_vars <- sum(1.7 * c(.6, .7, .8, .5, .3))^2 * 0.25 / 0.75
# Check composite reliability
sum(lambda)^2/(sum(lambda)^2 + error_var)

# Simulate latent continuous responses
item_star <- eta_scores %*% t(lambda) + rnorm(N * length(lambda), 
                                              mean = 0, 
                                              sd = sqrt(error_var/5))

# Thresholds
item_thres <- rbind(apply(item_star, 2, function(x) qnorm(0.05, mean = mean(x), sd = sd(x))),
                    apply(item_star, 2, function(x) qnorm(0.15, mean = mean(x), sd = sd(x))),
                    apply(item_star, 2, function(x) qnorm(0.4, mean = mean(x), sd = sd(x))))
rownames(item_thres) <- c("Thres1", "Thres2", "Thres3")

item_scores <- vapply(
  seq_along(lambda),
  FUN = \(i) {
    findInterval(item_star[, i], item_thres[, i])
  },
  FUN.VALUE = numeric(N))

# Check proportions
apply(item_scores, 2, table)/N

colnames(item_scores) <- paste0("Item", 1:5)
colnames(item_star) <- paste0("Item", 1:5)
item_summary <- as.data.frame(rbind(mean = apply(item_scores, 2, mean),
                      variance = apply(item_scores, 2, var)))
colnames(item_summary) <- paste0("Item", 1:5)

# Check lambdas 
test <- cfa(model = "F1 =~ Item1 + Item2 + Item3 + Item4 + Item5", item_star, std.lv = T)
coef(test)
lambda
# Check composite reliability using function
compRelSEM(test)

#---------------------------Check simulated data----------------------------
# Check proportions
apply(item_scores, 2, table)/nrow(item_scores)

# Fit IRT model
fit <- mirt(data = item_scores, 
               model = "F1 = Item1-5",
               itemtype = "graded", 
               verbose = FALSE)
summary(fit)
coef(fit)

# Fit categorical CFA
test_cat <- cfa(model = "F1 =~ Item1 + Item2 + Item3 + Item4 + Item5", 
                item_scores, 
                ordered = c("Item1", "Item2", "Item3", "Item4", "Item5"),
                std.lv = T)
coef(test_cat)
compRelSEM(test_cat)
