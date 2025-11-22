library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(igraph)

# === Directories ===
results_dir <- file.path("results")   # BayesTraits output
output_dir  <- file.path("output")    # Pre-computed alternative output

# Use results/ as primary search path
fin_z <- list.files(results_dir,
                    pattern = "^zscore_ch_GTDB207_cellshape_motility_sporulation_",
                    full.names = TRUE)

if (length(fin_z) == 0) {
  # fallback: use pre-computed example output
  fin_z <- list.files(output_dir,
                      pattern = "^zscore_ch_GTDB207_cellshape_motility_sporulation_",
                      full.names = TRUE)
}

# Warn if still nothing
if (length(fin_z) == 0) {
  stop("No zscore_* files found in results/ or output/.")
}

for (cfin in fin_z) {

  # Matching qrate file
  fin_q <- str_replace(cfin, "zscore", "qrate")

  # === Read z-score ===
  mydata_z <- read.delim(cfin)
  colnames(mydata_z) <- c("parameter", "zscore")

  # === Read qrate ===
  mydata_q <- read.delim(fin_q)
  colnames(mydata_q) <- c("parameter", "qrates")

  # === Merge ===
  mydata <- left_join(mydata_z, mydata_q, by = "parameter")

  # === Extract states ===
  mydf <- mydata %>%
    mutate(state_1 = str_sub(parameter, 2, 2),
           state_2 = str_sub(parameter, 3, 3)) %>%
    select(state_1, state_2, zscore, qrates)

  # nonzscore
  mydf <- mydf %>%
    mutate(nonzscore = 100 - zscore) %>%
    select(-zscore)

  # Threshold
  mythreshold_lower <- 30
  myedge <- mydf %>% filter(nonzscore > mythreshold_lower)

  if (nrow(myedge) > 0) {

    g <- graph_from_data_frame(myedge, directed = TRUE)
    E(g)$weight <- myedge$qrates

    # Layout
    mylayout <- layout_with_graphopt

    # Output filename → results/
    fout <- basename(cfin)
    fout <- str_replace(fout, "zscore_ch", "plot")
    fout <- str_replace(fout, ".tsv", ".png")
    fout <- file.path(results_dir, fout)

    png(fout, width = 3000, height = 2000)
    par(mar = c(0, 0, 0, 0))

    plot.igraph(g,
                layout = mylayout,
                vertex.color = "white",
                vertex.size = 10,
                vertex.shape = "circle",
                vertex.label.cex = 5,
                vertex.frame.width = 5,
                edge.arrow.size = 4,
                edge.curved = 0.2)

    dev.off()
  }
}
