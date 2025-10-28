library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(igraph)

# === Load file names for z-score results ===
#-----------------------------------------------------------------------------------------
fin_z <- list.files(getwd(), pattern ="^zscore_ch_GTDB207_cellshape_motility_sporulation_", full.names = F)
#-----------------------------------------------------------------------------------------

# === qrate result files will be loaded inside the for-loop ===

for (cfin in fin_z) {
  # Generate the corresponding qrate file name
  fin_q <- str_replace(cfin, "zscore", "qrate")

  # === Read z-score data ===
  # Note: do NOT use read_tsv() here, otherwise color codes will not be read correctly!
  mydata_z <- read.delim(cfin)
  # Rename columns
  colnames(mydata_z) <- c("parameter","zscore")

  # === Read qrate data ===
  # Note: do NOT use read_tsv() here, otherwise color codes will not be read correctly!
  mydata_q <- read.delim(fin_q)
  # Rename columns
  colnames(mydata_q) <- c("parameter","qrates")
  
  # === Merge z-score and qrate data ===
  mydata <- left_join(mydata_z, mydata_q, by="parameter")
  
  # === Create edge list ===
  # Split transition states:
  # qAB -> A, B
  mybuf <- mydata %>% mutate(state_1=str_sub(parameter, 2, 2)) 
  mydf <- mybuf %>% mutate(state_2=str_sub(parameter, 3, 3)) %>% select(state_1, state_2, zscore, qrates)
  
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
    E(g)$weight <- myedge$qrates

    # === Plot settings ===
    #----------------------------------------------------------
    # Define layout for visualization
    mylayout=layout_with_graphopt

    # output file name
    fout <- str_replace(cfin, "zscore_ch", "plot")
    fout <- str_replace(fout, "tsv", "png")
    
    # Save plot as PNG
    png(fout, width=3000, height=2000)
    
    par(mar = c(0, 0, 0, 0)) #  Set margins: bottom, left, top, right
    # ---- Plot graph ----
    plot.igraph(g,
                layout = mylayout,
                vertex.color = "white",
                vertex.size = 10,
                vertex.shape = "circle",
                vertex.label.cex = 5,
                vertex.frame.width=5,
                edge.lty = E(g)$lty,
                edge.arrow.size = 4,
                edge.curved = 0.2)
    dev.off()
  }
}
