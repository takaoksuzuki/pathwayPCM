library(readr)
library(stringr)
library(dplyr)

# === Directories ===
restrict_dir <- "results"       # output of makerestriction.R
out_dir <- "results"            # output of this script (same folder OK)
dir.create(out_dir, showWarnings = FALSE)

# === List all Restrict files ===
fin <- list.files(restrict_dir, pattern = "^mycommand_", full.names = TRUE)

if (length(fin) == 0){
  stop("No mycommand_*.txt files found in results/. Run makerestriction.R first.")
}

# === Loop through each phylum restriction file ===
for (cf in fin){

  # phylum name
  base <- basename(cf)
  phylum <- str_replace(base, "^mycommand_", "") %>% str_replace(".txt$", "")

  # read restriction lines (Restrict qXX qYY ...)
  restrict_lines <- readLines(cf)

  # output filename (BayesTraits command file)
  fout <- file.path(out_dir, paste0("inputcommands_Res_", phylum, ".txt"))

  # === Write BayesTraits command file ===
  conn <- file(fout, "w")

  writeLines("1", conn)   # index of tree file
  writeLines("2", conn)   # index of data file
  writeLines("RevJumpHP gamma 0 1 0 1", conn)

  # Restrict lines
  writeLines(restrict_lines, conn)

  # RJMCMC settings
  writeLines("Burnin 50000", conn)
  writeLines("Sample 1000", conn)
  writeLines("Iterations 50050000", conn)
  writeLines("Stones 100 10000", conn)

  writeLines("Run", conn)

  close(conn)

  message("Generated: ", fout)
}

