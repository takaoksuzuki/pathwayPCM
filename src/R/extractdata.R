library(readr)
library(dplyr)
library(tidyr)
library(stringr)

fin <- list.files(getwd(), pattern = "\\.Log\\.txt$", full.names = FALSE)

for (cfin in fin) {
  message("Processing: ", cfin)
  
  # === Read log file ===
  all_lines <- readLines(cfin)
  
  # Find the line starting with "Iteration"
  start_line <- which(grepl("^Iteration\t", all_lines))[1]
  data_lines <- all_lines[start_line:length(all_lines)]
  
  # Extract column names
  col_names <- strsplit(data_lines[1], "\t")[[1]]
  
  # In the case of fossil models, the "Root" column appears twice - adjust the names
  root_names <- grep("^Root", col_names, value = TRUE)
  n <- length(root_names) / 2
  if (n == floor(n) && identical(root_names[1:n], root_names[(n+1):(2*n)])) {
    second_set <- grep("^Root", col_names)[(n+1):(2*n)]
    col_names[second_set] <- paste0(col_names[second_set], "_F")
  }
  
  # Convert the data into a tibble
  df <- tibble(datacol = data_lines) %>%
    separate(datacol, into = col_names, sep = "\t", extra = "drop") %>%
    slice(-1)   # Remove the header column
  
  # Define output file name
  base_name <- str_remove(cfin, "\\.Log.*$")
  fout <- paste0("eData_", base_name, ".tsv")
  
  # Write output file
  write_tsv(df, fout)
}
