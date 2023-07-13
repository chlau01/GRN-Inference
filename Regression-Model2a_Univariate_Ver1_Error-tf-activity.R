# if (!require("BiocManager")){
#   install.packages("BiocManager")
# }
  
# BiocManager::install("dorothea")

# install.packages("remotes")
# library(remotes)
# remotes::install_github('saezlab/OmnipathR')
# remotes::install_github('saezlab/decoupleR')

# if(!require(installr)){
#   install.packages("installr");
#   require(installr)
# }
# updateR()

library(dorothea)
library(OmnipathR)
library(decoupleR)
library(tidyverse)


# CollecTRI ----------------------------------------------

# CollecTRI (TF-based search) ---------------------------------------------

net_TRI <- decoupleR::get_collectri(split_complexes = FALSE)
head(net_TRI)

TF_query <- "SOX9"
tf_targets <- net_TRI |> 
  filter(grepl(TF_query, source))
head(tf_targets)


# Retrieve targets transcript levels -----------------------------------------------

gtex_df <- read_rds("data/expr_tissue-median_gtex.rds")
gtex_df <- as.data.frame(gtex_df$data)
head(gtex_df)


targets_transcripts <- gtex_df[tf_targets$target,]
head(targets_transcripts)


if (nrow(tf_targets) == nrow(targets_transcripts)) {
  print("All targets found")
} else {
  print("check missing rows/redundant rows")
}


# TF activity = Targets transcript level  ------------------------------------------

########### Literature search: How to approximate TF activity from transcript level?


tf_activity <- 
tf_activity




# univariate regressions ------------------------------------------------------


target_exp <- gtex_df["COL2A1",]
target_exp <- t(target_exp) #turn into column vector
target_exp


lm_target = lm(tf_activity~target_exp)
summary(lm_target)
plot(lm_target$model)


