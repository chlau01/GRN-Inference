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

# CollecTRI (target-based search) -----------------------------------------

# Searching and recording literature-supported TFs for a target (not looped) -------------------------------------------------------------------------


net_TRI <- decoupleR::get_collectri(split_complexes = FALSE)
head(net_TRI)


target_query <- "COL2A1"
target_tf <- net_TRI |> 
  filter(grepl(target_query, target))
view(target_tf)

literature_found_Jun23 <- c("SOX9", "KLF4", "ARID5B", "SOX5", "SOX6",
                            "NKX3-2", "PAX9", "NFATC1", "SP1", "SP3")

literature_found_Jun30 <- c("")

literature_found_Jul07 <- c("")


# Loopable (start)

(literature_found_Jun23 <- tibble(
  source = literature_found_Jun23, 
  target = target_query,
  literature_Jun23 = "found"))

target_tf=as.data.frame(target_tf)

target_tf_rec <- target_tf |> 
  full_join(literature_found_Jun23, join_by(source == source, target == target))
view(target_tf_rec)
target_tf_rec=as.data.frame(target_tf_rec)

target_tf_rec <- target_tf_rec |> 
  mutate(literature_in_collecTRI = ifelse(is.na(mor), FALSE, ifelse(is.na(target_tf_rec$literature_Jun23), NA, TRUE))) |> 
  arrange(literature_Jun23, desc(literature_in_collecTRI))
  
view(target_tf_rec)

# Loopable (end)



# Searching and recording literature-supported TFs for a target (loop) -------------------------------------------------------------------------


