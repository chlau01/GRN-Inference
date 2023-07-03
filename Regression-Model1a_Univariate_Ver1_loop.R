library(tidyverse)


# import ------------------------------------------------------------------


gtex_df <- read_rds("data/expr_tissue-median_gtex.rds")
glimpse(gtex_df)
view(gtex_df)

gtex_df <- as.data.frame(gtex_df$data)
# gtex_df <- as.matrix(gtex_df$data)


# Convert data.frame to tibble. But rownames lost
# gtex_df <- as_tibble(gtex_df)


head(gtex_df)

# gtex_df <- t(gtex_df)



# univariate regressions (not looped) ------------------------------------------------------


target_exp <- gtex_df["COL2A1",]
target_exp <- t(target_exp) #turn into column vector
target_exp

tf_exp <- gtex_df["SOX9",]
tf_exp <- t(tf_exp)
tf_exp

lm_target = lm(tf_exp~target_exp)
summary(lm_target)
plot(lm_target$model)



tf_exp <- gtex_df["KLF4",]
tf_exp <- t(tf_exp)
tf_exp

lm_target = lm(tf_exp~target_exp)
summary(lm_target)
plot(lm_target$model)



tf_exp <- gtex_df["ARID5B",]
tf_exp <- t(tf_exp)
tf_exp

lm_target = lm(tf_exp~target_exp)
summary(lm_target)
plot(lm_target$model)



tf_exp <- gtex_df["SOX5",]
tf_exp <- t(tf_exp)
tf_exp

lm_target = lm(tf_exp~target_exp)
summary(lm_target)
plot(lm_target$model)



tf_exp <- gtex_df["SOX6",]
tf_exp <- t(tf_exp)
tf_exp

lm_target = lm(tf_exp~target_exp)
summary(lm_target)
plot(lm_target$model)



tf_exp <- gtex_df["NKX3-2",]
tf_exp <- t(tf_exp)
tf_exp

lm_target = lm(tf_exp~target_exp)
summary(lm_target)
plot(lm_target$model)



tf_exp <- gtex_df["PAX9",]
tf_exp <- t(tf_exp)
tf_exp

lm_target = lm(tf_exp~target_exp)
summary(lm_target)
plot(lm_target$model)



tf_exp <- gtex_df["NFATC1",]
tf_exp <- t(tf_exp)
tf_exp

lm_target = lm(tf_exp~target_exp)
summary(lm_target)
plot(lm_target$model)



tf_exp <- gtex_df["SP1",]
tf_exp <- t(tf_exp)
tf_exp

lm_target = lm(tf_exp~target_exp)
summary(lm_target)
plot(lm_target$model)



tf_exp <- gtex_df["SP3",]
tf_exp <- t(tf_exp)
tf_exp

lm_target = lm(tf_exp~target_exp)
summary(lm_target)
plot(lm_target$model)



# univariate regression (loop) --------------------------------------------


target_exp <- gtex_df["COL2A1",]
target_exp <- t(target_exp)
target_exp

tf_list <- c("SOX9", "KLF4", "ARID5B", "SOX5", "SOX6",
               "NKX3-2", "PAX9", "NFATC1", "SP1", "SP3")
lm_uni_tf_list=list()


tf_exp <- gtex_df["SOX9",]
tf_exp


for (i in tf_list){
  tf_exp <- t(gtex_df[i,])
  tf_exp <- t(tf_exp)
  print(tf_exp)
  
  lm_uni_tf_list[[i]] = lm(tf_exp~target_exp)
  summary(lm_uni_tf_list[i])
  
  pdf(file = paste0("./output/",i,"rm.pdf"),width=7,height=7)
  print(plot(lm_uni_tf_list[[i]]))
  dev.off()
}



