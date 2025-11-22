library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(igraph)

# === Set directories ===
data_dir <- "data"
results_dir <- "results"

# === Load z-score result files from results/ ===
fin_z <- list.files(results_dir, 
                    pattern = "^zscore_ch_GTDB207_cellshape_motility_sporulation_", 
                    full.names = TRUE)

# === Loop over each z-score file ===
for (cfin in fin_z) {

  # --- Read z-score data ---
  mydata <- read.delim(cfin)
  colnames(mydata) <- c("parameter","zscore")

  # --- Create edge list (qAB → A, B) ---
  mybuf <- mydata %>% mutate(state_1 = str_sub(parameter, 2, 2))
  mydf  <- mybuf %>% mutate(state_2 = str_sub(parameter, 3, 3)) %>%
           select(state_1, state_2, zscore)

  # --- Convert z-score to non-z-score ---
  mydf <- mydf %>% mutate(nonzscore = 100 - zscore) %>% select(-zscore)

  # --- Threshold filtering ---
  mythreshold_lower <- 30
  myedge <- mydf %>% filter(nonzscore > mythreshold_lower)

  if (nrow(mydf) > 0) {

    # --- Build weighted graph ---
    g <- graph_from_data_frame(myedge, directed = TRUE)
    E(g)$weight <- myedge$nonzscore

    # --- Infomap clustering ---
    aaa <- cluster_infomap(g, e.weights = E(g)$weight)

    # --- Extract modules ---
    mymodule <- tibble(
      label  = names(membership(aaa)),
      module = membership(aaa)
    )

    # --- Output file name ---
    fout_name <- basename(cfin) %>%
                 str_replace("zscore_ch", "modulePathway")
    fout <- file.path(results_dir, fout_name)

    # --- Save output ---
    write_tsv(mymodule, fout)
  }
}
