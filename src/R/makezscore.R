library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# === Get the list of target eData files ===
fin <- list.files(getwd(), pattern="^eData_.*\\.tsv\\.tsv$", full.names=FALSE)

# === Get the list of corresponding Restrict files ===
fin_res <- list.files(getwd(), pattern="^mycommand_coding_.*\\.tsv$", full.names=FALSE)

for (cfin in fin) {
  message("Processing: ", cfin)
  
  # Estimate the corresponding mycommand file key
　command_file <- fin_res[grepl("Firmicutes", fin_res)]

  # === Read file ===
  df <- read_tsv(edata_file, show_col_types = FALSE)
  model_matrix <- str_split(df$`Model string`, "\\s+", simplify = TRUE)
  n_rows <- nrow(df)

  # === Retrieve all qXX columns (column names starting with ^q) ===
  all_qxx <- grep("^q", colnames(df), value = TRUE)

  # === Process Restrict statements ===
  group_list <- list()
  group_index <- 1
  used_qxx <- c()
  zero_qxx <- c()

  if (length(command_file) == 1) {
    restrict_lines <- readLines(command_file)
    restrict_lines <- restrict_lines[grepl("^Restrict", restrict_lines)]

    for (line in restrict_lines) {
      tokens <- strsplit(line, "\\s+")[[1]][-1]
      if (tail(tokens, 1) == "0") {
        zero_qxx <- c(zero_qxx, tokens[-length(tokens)])
        next
      }
      used_qxx <- c(used_qxx, tokens)
      group_list[[group_index]] <- tokens
      group_index <- group_index + 1
    }
  }

  used_qxx <- unique(c(used_qxx, zero_qxx))
  free_qxx <- sort(setdiff(all_qxx, used_qxx))

  # Extract the first transition rate from each group
  group_reps <- sapply(group_list, function(g) sort(g)[1])
  group_reps <- unlist(group_reps)
  model_qxx <- sort(c(group_reps, free_qxx))

  # === Calculate Z-scores (in %) ===
  z_counts <- apply(model_matrix, 2, function(col) sum(col == "Z"))
  z_scores <- (z_counts / n_rows) * 100
  zscore_map <- setNames(z_scores, model_qxx)

  # === Record the representative qXX for each group ===
  group_map <- list()
  for (grp in group_list) {
    rep <- sort(grp)[1]
    for (q in grp) group_map[[q]] <- rep
  }
  for (q in free_qxx) group_map[[q]] <- q
  for (q in zero_qxx) group_map[[q]] <- q

  # === Expand Z-scores for all qXX ===
  result <- tibble(
    qXX = sort(all_qxx),
    Z_score = sapply(sort(all_qxx), function(q) {
      if (q %in% zero_qxx) return(100)  # 固定ゼロ
      rep <- group_map[[q]]
      if (is.null(rep) || is.null(zscore_map[[rep]])) return(0)
      round(zscore_map[[rep]], 2)
    })
  )

  # === Generate output file name and save ===
  fout <- paste0("zscore_", str_replace(edata_file, "^eData_", ""))
  write_tsv(result, fout)
}
