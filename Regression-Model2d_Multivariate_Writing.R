## Documentation

# Collects TF activity matrix and performs (A) univariate regression, (B) multiple regression, 
# (C) variable selection with lasso and (D) variable selection with non-local prior
# TODO: For the actual pipeline, move the modeling script to another file


library(filenamer)
library(io)
library(Matrix)
library(glmnet)

library(mvtnorm)
library(ncvreg)
library(nlme)
library(mgcv)
library(mombf)
library(caret)

library(plyr)
library(tidyverse)


# Input expression data matrix (rownames: gene symbols)
expr <- qread("data/expr_tissue-median_gtex.rds");
expr <- expr$data;



# input TF-target list of targets

# library(OmnipathR)
# library(decoupleR)
# tf.targets <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
# write_rds(tf.targets, "data/net_TRI.rds")
tf.targets <- read_rds("data/net_TRI.rds");



# input target query
target <- "CCND1"


# input TFs of the target (e.g. from our manually-curated literature database)
tfs.query <- c("CTNNB1", "JUNB", "SNIP1", "STAT3", "STAT5",
               "OCT1", "E2F1", "E2F4", "RELA", "NFKB1")


# input candidate TFs of the target (e.g. from CollecTRI, motif binding database, ChIP-seq database, etc.)
tfs.candidates <- subset(tf.targets, tf.targets$target == "CCND1")$source




# Infer TF activity based on mean expression of target genes --------------

# 1 Infer TF activity matrix (all samples by all TFs)


# @param target.tfs   a list of TFs whose activity you want to obtain from the input expression data matrix
get_activity <- function(target.tfs){
  tf_activity <- c()
  for(i in 1:length(target.tfs)){
    tf.activity.targets <- subset(tf.targets, tf.targets$source == target.tfs[i])$target
    tf.activity.targets.present <- tf.activity.targets[tf.activity.targets %in% rownames(expr)] # Some targets in stat3_targets are not present in expr
    mean_exp <- colMeans(matrix(unlist(expr[tf.activity.targets.present, ]), ncol = dim(expr)[2]), na.rm = TRUE) # Collect the mean expression of tfs across all tissues (mean of each column). colMeans only takes a matrix with 2 dimensions as input.
    tf_activity <- rbind2(tf_activity, mean_exp) 
  }
  rownames(tf_activity) <- target.tfs[1:length(target.tfs)]
  colnames(tf_activity) <- colnames(expr)
  return(tf_activity)
}


query.tf.activities <- get_activity(tfs.query)
sum(is.na(query.tf.activities[, "Adipose - Subcutaneous"]))
query.tf.activities <- query.tf.activities[!is.na(query.tf.activities[, "Adipose - Subcutaneous"]), ] # Remove TFs with Na activities
query.tf.activities |> view()
tfs.query <- rownames(query.tf.activities) # Update with only TFs with known activities


tfs.candidates <- tfs.candidates[tfs.candidates %in% rownames(expr)] # Remove TFs whose expression is not in expr
candidate.tf.activities <- get_activity(tfs.candidates)
candidate.tf.activities |> view()






# Investigation: NA values  -----------------------------------------------------------

sum(is.na(query.tf.activities[, "Adipose - Subcutaneous"])) # 3 of the 10 input TFs have NA values in tf.activity. Why?
# Investigate why there are NA values and how to interpret them
# Each TF should have at lease 1 known target, which is the query TF ("CCND1")
# Targets whose expression were not found in expr should be filtered out.
# Some possibilities investigated:
rownames(expr)[str_detect(rownames(expr), "CCND1")] # expr contains "CCND1" expression
rownames(expr)[str_detect(rownames(expr), "SNIP1")] # expr contains "SNIP1" expression
tf.targets[str_detect(tf.targets$source, "SNIP1"), ] # CollecTRI does not contain "SNIP1". May be another gene symbol

# TODO: Investigate the gene symbols of "SNIP1" is used in collecTRI, which uses HGNC symbols (http://bioconductor.org/packages/release/data/experiment/manuals/dorothea/man/dorothea.pdf)
# (A) try to use the correct gene symbols.
# (B) if you cannot find the expression of the tf.targets, remove the TF from the equation

# -------------------------------------------------------------------------




# # 2 Obtaining the plotting data matrix: TFs and Target
plotting.data.query <- rbind2(expr[target, ], query.tf.activities)
rownames(plotting.data.query)[1] <- target # replace rowname with target gene symbol
plotting.data.query |> view()

plotting.data.candidates <- rbind2(expr[target, ], candidate.tf.activities)
rownames(plotting.data.candidates)[1] <- target # replace rowname with target gene symbol
plotting.data.candidates |> view()



# 3 Normalize expression level of target genes (divide each TF expression by its 
# mean before colMeans) before averaging to obtain TF activity

# NOT NEEDED to do this currently






# -------------------------------------------------------------------------
# Regression models (uni, multi, lasso, non-local prior) to infer regulatory effects (betas) ------------------




# Univariate regression model ---------------------------------------------

# # @param target.expj    expression vector of target across N samples
# # @param tf.activities  TF activity matrix (N samples by K TFs)
# # @param target         gene name of target
# # @param tf.candidates  gene names of candidate TFs
# ulinear_model <- function(expr, tf.activities, target, tf.candidates){
#   unlist(lapply(tf.candidates, 
#                 function(tf){
#                   fit <- lm(expr[target, ] ~ -1 + tf.activities[tf, ]);
#                   coef(fit)
#                   summary(fit)
#                   plot(fit$model)
#                 }
#   ))
# }
# TODO: Error - this function only produces 1 univariate regression out of 10


# # @param tf.activities  TF activity matrix containing activities of the TFs 
# # @param expr           An expression data matrix
# univariate_regression_inference <- function(tf.activities, expr){
#   lm_uni_tf_list <- list()
#   tf.activities.present <- tf.activities[!is.na(tf.activities[, "Adipose - Subcutaneous"]), ] # Remove TFs with Na activities
#   
#   for (i in tf.activities.present){
#     tf_activity <- tf.activities.present[i,]
#     target_expr <- expr[target, ]
#     
#     lm_uni_tf_list[[i]] <-  lm(target_expr ~ -1 + tf_activity)
#     
#     pdf(file = paste0("./output/",i,"rm.pdf"),width=7,height=7) # Save plots as pdf in output folder
#     print(plot(target_expr ~ -1 + tf_activity))
#     print(plot(lm_uni_tf_list[[i]]))
#     dev.off()
#   }
# }
# 
# univariate_regression_inference(query.tf.activities, expr)







# Multivariate regression model -------------------------------------------


# # @param plotting.data  Plotting data matrix of target expression + TF activity (N samples by 1 Target + K TFs)
# # @param target         gene symbol of target
# # @param tfs            gene symbols of candidate TFs
mlinear_model <- function(plotting.data, target, tfs){
  fit <- lm(plotting.data[target, ] ~ -1 + t(plotting.data[tfs, ]));
  return(fit)
}



# Lasso multivariable variable selection model ----------------------------
# Uses Gaussian linear model (GLM) or least-squares


lasso_model <- function(plotting.data, target, tfs){
  require(glmnet)
  fit <- glmnet(t(plotting.data[tfs, ]), plotting.data[target, ], alpha = 1, intercept = FALSE);
  return(fit)
}



# mombf multivariable variable selection model ----------------------------
# uses Bayesian information criterion (BIC) for model selection. ie. Bayesian model selection (BMS)


mombf_model <- function(plotting.data, target, tfs){
  require(mombf)
  fit <- bestBIC(plotting.data[target, ] ~ -1 + t(plotting.data[tfs, ]))
  return(fit)
}



# Run regression models -------------------------------------------------------

# For query tfs
hist(expr[tfs.query, ]) # Inspect data distribution of TF expression and TF activity respectively
hist(plotting.data.query[tfs.query, ])


fit.expr.b <- mlinear_model(expr, target, tfs.query) # Using TF expression to predict Target expression
fit.activity.b <- mlinear_model(plotting.data.query, target, tfs.query) # Using TF activity to predict Target expression

fit.expr.c <- lasso_model(expr, target, tfs.query)
fit.activity.c <- lasso_model(plotting.data.query, target, tfs.query)

fit.expr.d <- mombf_model(expr, target, tfs.query)
fit.activity.d <- mombf_model(plotting.data.query, target, tfs.query)


# For all candidate tfs
hist(expr[tfs.candidates, ])
hist(plotting.data.candidates[tfs.candidates, ])


fit.expr.b.candidates <- mlinear_model(expr, target, tfs.candidates)
fit.activity.b.candidates <- mlinear_model(plotting.data.candidates, target, tfs.candidates)

fit.expr.c.candidates <- lasso_model(expr, target, tfs.candidates)
fit.activity.c.candidates <- lasso_model(plotting.data.candidates, target, tfs.candidates)

fit.expr.d.candidates <- mombf_model(expr, target, tfs.candidates)
fit.activity.d.candidates <- mombf_model(plotting.data.candidates, target, tfs.candidates)




# Compare model performances ----------------------------------------------

# Metrics
summary(fit.expr.b)
coef(fit.expr.b)

fit.expr.d$topmodel.fit
confint(fit.activity.d)
fit.activity.d$topmodel.fit


# RMSE of predicted target expression (Y hat)

plot_prediction_calc_nrmse <- function(model, data, y_test){
  prediction <- predict(model, as.data.frame(data))
  plot(y_test, prediction, xlab = "Actual Values", ylab = "Predicted Values", main = "Predicted vs Actual Target expression")
  
  nrmse <- sqrt(mean((prediction - y_test)^2)) / (max(y_test) - min(y_test))
  
  output <- list(prediction = prediction, nrmse = nrmse)
  return(output)
}

calc_nrmse <- function(prediction, y_test){
  output <- sqrt(mean((prediction - y_test)^2)) / (max(y_test) - min(y_test))
  return(output)
}



pred.expr.b <- plot_prediction_calc_nrmse(fit.expr.b, expr, expr[target, ])
pred.activity.b <- plot_prediction_calc_nrmse(fit.activity.b, plotting.data.query, expr[target, ])


pred.expr.b.candidates <- plot_prediction_calc_nrmse(fit.expr.b.candidates, expr, expr[target, ])
pred.activity.b.candidates <- plot_prediction_calc_nrmse(fit.activity.b.candidates, plotting.data.candidates, expr[target, ])



fit.expr.c
pred.expr.c <- predict(fit.expr.c, 
                       newx = t(expr[tfs.query, ]), type = "response", s = 0.05) # Gives the fitted values at lambda = 0.05
nrmse.expr.c <- calc_nrmse(pred.expr.c, expr[target, ])
pred.activity.c <- predict(fit.activity.c, 
                           newx = t(plotting.data.query[tfs.query, ]), type = "response", s = 0.05) # Gives the fitted values at lambda = 0.05
nrmse.activity.c <- calc_nrmse(pred.activity.c, expr[target, ])




fit.expr.c.candidates
pred.expr.c.candidates <- predict(fit.expr.c.candidates, 
                       newx = t(expr[tfs.candidates, ]), type = "response", s = 0.05) # Gives the fitted values at lambda = 0.05
nrmse.expr.c.candidates <- calc_nrmse(pred.expr.c.candidates, expr[target, ])




# Cross validation and choosing lambda value for Lasso model
#
x <- t(expr[tfs.query,])
y <- expr[target,]
cvfit <- cv.glmnet(x, y, type.measure = "mse", nfolds = 20)
cvfit
# print(cvfit)
# cvfit$lambda.min # Lambda at which smallest MSE is achieved
# pred.expr.c <- predict(fit.expr.c, newx = x, type = "response", s = cvfit$lambda.min) # Uses the lambda.min to fit lasso
# nrmse.expr.c <- calc_nrmse(pred.expr.c, expr[target, ])






# Plotting figures

pred.expr.b$nrmse
pred.activity.b$nrmse

pred.expr.b.candidates$nrmse
pred.activity.b.candidates$nrmse


nrmse.expr.c
nrmse.activity.c

plot(fit.expr.c, xvar = "lambda", label = TRUE)
plot(fit.expr.c, xvar = "dev", label = TRUE)



nrmse.expr.c.candidates

plot(fit.expr.c.candidates, xvar = "lambda", label = TRUE)
plot(fit.expr.c.candidates, xvar = "dev", label = TRUE)







# Create a data frame with the metric values
metric_df <- data.frame(
  model = c("TF expression", "TF activity", "CollecTRI database"),
  NRMSE = c(0.1417954, 0.1124271, 1.659121e-14),
  color = c("red", "blue", "green")
)

# Create a plot object with ggplot()
plot <- ggplot(data = metric_df, aes(x = model, y = NRMSE, fill = color))

# Customize the plot with additional layers and options
plot + 
  geom_bar(stat = "identity", width = .4) +
  geom_text(aes(label = sprintf("%.2e", NRMSE)), position = position_stack(vjust = 0.9), hjust = 1) +
  labs(x = "Model", 
       y = "Normalized root mean squared error (NRMSE) (log scale)", 
       title = "Model Performance") +
  theme_minimal() + 
  scale_y_log10() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  coord_flip()












# Obtain predicted regulatory coefficients (Beta hat) from regression models

beta.hat.b <- coef(fit.expr.b)
beta.hat.c <- coef(fit.expr.c)
beta.hat.d <- coef(fit.expr.d)

beta.hat.b
beta.hat.c
beta.hat.d

hist(beta.hat.b)
# hist(beta.hat.c) # Error: presence of non-number values
hist(beta.hat.d)





# TODO Construct a simulation model with known Beta values

RMSE <- (beta.hat - beta)^2

