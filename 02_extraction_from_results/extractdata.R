library(readr)
library(dplyr)
library(tidyr)
library(stringr)

fin <- list.files(getwd(), pattern = "\\.Log\\.txt$", full.names = FALSE)

for (cfin in fin) {
  message("Processing: ", cfin)
  
  # === ログファイル読み込み ===
  all_lines <- readLines(cfin)
  
  # Iteration 行を探す
  start_line <- which(grepl("^Iteration\t", all_lines))[1]
  data_lines <- all_lines[start_line:length(all_lines)]
  
  # 列名を抽出
  col_names <- strsplit(data_lines[1], "\t")[[1]]
  
  # Fossil の場合 Root 列が2回出るので調整
  root_names <- grep("^Root", col_names, value = TRUE)
  n <- length(root_names) / 2
  if (n == floor(n) && identical(root_names[1:n], root_names[(n+1):(2*n)])) {
    second_set <- grep("^Root", col_names)[(n+1):(2*n)]
    col_names[second_set] <- paste0(col_names[second_set], "_F")
  }
  
  # データ本体を tibble 化
  df <- tibble(datacol = data_lines) %>%
    separate(datacol, into = col_names, sep = "\t", extra = "drop") %>%
    slice(-1)   # ヘッダー列を除去
  
  # 出力ファイル名
  base_name <- str_remove(cfin, "\\.Log.*$")
  fout <- paste0("eData_", base_name, ".tsv")
  
  # 書き出し
  write_tsv(df, fout)
}
