library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(igraph)

# === Directories ===
results_dir <- "results"   # BayesTraits results (priority)
output_dir  <- "output"    # Pre-computed example results

# === Try reading zscore_* from results/ first ===
fin_z <- list.files(
  results_dir,
  pattern = "^zscore_ch_GTDB207_cellshape_motility_sporulation_",
  full.names = TRUE
)

# === Fallback to output/ if results/ has no zscore files ===
if (length(fin_z) == 0) {
  fin_z <- list.files(
    output_dir,
    pattern = "^zscore_ch_GTDB207_cellshape_motility_sporulation_",
    full.names = TRUE
  )
  input_dir <- output_dir
} else {
  input_dir <- results_dir
}

# === If still no files → stop ===
if (length(fin_z) == 0) {
  stop("No zscore_* files found in 'results/' or 'output/'")
}

message("Using zscore files from: ", input_dir)

# =======================================================================
# Loop over each zscore file
# =======================================================================
for (cfin in fin_z) {

  # Corresponding qrate file is in the same directory
  fin_q <- str_replace(cfin, "zscore", "qrate")

  # Check that matching qrate exists
  if (!file.exists(fin_q)) {
    stop(paste("Matching qrate file not found for:", cfin))
  }

  # --- Read zscore ---
  mydata <- read.delim(cfin)
  colnames(mydata) <- c("parameter","zscore")

  # === Extract states & convert to nonzscore ===
  mydf <- mydata %>%
    mutate(
      state_1 = str_sub(parameter, 2, 2),
      state_2 = str_sub(parameter, 3, 3),
      nonzscore = 100 - zscore
    ) %>%
    select(state_1, state_2, nonzscore)

  # === Threshold filtering ===
  mythreshold_lower <- 30
  myedge <- mydf %>% filter(nonzscore > mythreshold_lower)

  if (nrow(myedge) == 0) next

  # --- Build weighted graph ---
  g <- graph_from_data_frame(myedge, directed = TRUE)
  E(g)$weight <- myedge$nonzscore

  # --- Infomap clustering ---
  aaa <- cluster_infomap(g, e.weights = E(g)$weight)

  # --- Save results ---
  mymodule <- tibble(
    label  = names(membership(aaa)),
    module = membership(aaa)
  )

  fout_name <- basename(cfin) %>%
    str_replace("zscore_ch", "modulePathway")

  fout <- file.path(results_dir, fout_name)

  write_tsv(mymodule, fout)

  message("Saved module file: ", fout)
}
