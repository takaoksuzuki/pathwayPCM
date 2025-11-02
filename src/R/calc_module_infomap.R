library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(igraph)

# === Load file names for z-score results ===
#-----------------------------------------------------------------------------------------
fin_z <- list.files(getwd(), pattern ="^zscore_ch_GTDB207_cellshape_motility_sporulation_", full.names = F)
#-----------------------------------------------------------------------------------------

for (cfin in fin_z) {
  # === Read z-score data ===
  # Note: do NOT use read_tsv() here, otherwise color codes will not be read correctly!
  mydata <- read.delim(cfin)
  # Rename columns
  colnames(mydata) <- c("parameter","zscore")
  
  # === Create edge list ===
  # Split transition states:
  # qAB -> A, B
  mybuf <- mydata %>% mutate(state_1=str_sub(parameter, 2, 2)) 
  mydf <- mybuf %>% mutate(state_2=str_sub(parameter, 3, 3)) %>% select(state_1, state_2, zscore)
  
  # === Convert z-score to “non-z-score” ===
  mydf <- mydf %>% mutate(nonzscore = 100 - zscore) %>% select(-zscore)
  
  # === Set threshold ===
  #-------------------------------
  mythreshold_lower <- 30
  #-------------------------------
  # Filter edges based on the threshold
  myedge <- mydf %>% filter(nonzscore > mythreshold_lower)
  
  if (nrow(mydf) > 0){
    # === Convert to graph ===
    g <- graph_from_data_frame(myedge, directed = TRUE)
    # Assign edge weights
    E(g)$weight <- myedge$nonzscore
    
    #モジュール検出：infomap法
    #------------------------------------------------------
    aaa <- cluster_infomap(g, e.weights = E(g)$weight)
    myfout <- "infomap"
    #------------------------------------------------------
    
    #モジュールの保存
    mymodule <- tibble(label = names(membership(aaa)), module = membership(aaa)) 
    
    #ファイル名
    fout <- str_replace(cfin, "zscore_ch", "modulePathway")
    write_tsv(mymodule, fout)
  }
}
