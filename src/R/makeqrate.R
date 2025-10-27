library(readr)
library(dplyr)
library(stringr)

fin <- list.files(getwd(), pattern="^eData_", full.names=F)

for (cfin in fin){
  message("Processing: ", cfin)

  # === Read file ===
  mydata <- read_tsv(cfin)
  
  # Extract column indices for data section
  aaa <- grep("^q", colnames(mydata))
  mystart <- aaa[1]
  bbb <- grep("^Root", colnames(mydata))
  myend <- bbb[1] - 1
  
  # Extract data columns
  mydata_ex <- mydata[,mystart:myend]
  
  # Calculate posterior means
  mydata_sum <- mydata_ex %>% apply(2,sum) %>% t %>% as.data.frame()
  mydata_ave <- mydata_sum / nrow(mydata)
  
  # Format output data
  # Column 1: transition, Column 2: rates (qAB, etc.)
  mydata_ave <- mydata_ave %>% t %>% as.data.frame
  myrates <- rownames(mydata_ave) %>% as.data.frame
  myresult <- bind_cols(myrates, mydata_ave)
  rownames(myresult) <- NULL
  colnames(myresult) <- c("parameter", "qrates")

  # output file name
  mystart <- regexpr("eData_", cfin) + nchar("eData_")
  myend <- nchar(cfin)
  myfilename <- str_sub(cfin, mystart, myend)
  fout <- str_c("qrate_", myfilename, sep = "")
  
  # Save output
  # Convert the data frame to character type before writing
  df <- apply(myresult, 2, as.character)
  write.table(df, fout, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
}
