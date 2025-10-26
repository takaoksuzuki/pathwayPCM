library(readr)
library(dplyr)
library(stringr)

#setwd("../000_data_IO")
# setwd("../000_data_IO_butterfly/")
fin <- list.files(getwd(), pattern="^eData_", full.names=F)

for (cfin in fin){
  print(cfin)

  #読み込み
  # mydata <- read.delim(cfin)
  mydata <- read_tsv(cfin)
  
  #データ列番号の抽出
  aaa <- grep("^q", colnames(mydata))
  mystart <- aaa[1]
  bbb <- grep("^Root", colnames(mydata))
  myend <- bbb[1] - 1
  
  #データ列の抽出
  mydata_ex <- mydata[,mystart:myend]
  
  #事後確率の計算
  mydata_sum <- mydata_ex %>% apply(2,sum) %>% t %>% as.data.frame()
  mydata_ave <- mydata_sum / nrow(mydata)
  
  #出力へデータ整形
  #１列目：遷移、２列目：rates(qAB,)
  mydata_ave <- mydata_ave %>% t %>% as.data.frame
  myrates <- rownames(mydata_ave) %>% as.data.frame
  myresult <- bind_cols(myrates, mydata_ave)
  rownames(myresult) <- NULL
  colnames(myresult) <- c("parameter", "qrates")

  #ファイル名
  mystart <- regexpr("eData_", cfin) + nchar("eData_")
  myend <- nchar(cfin)
  myfilename <- str_sub(cfin, mystart, myend)
  fout <- str_c("qrate_", myfilename, sep = "")
  
  #保存
  #データフレームを文字型に変換しないと保存できない
  df <- apply(myresult, 2, as.character)
  write.table(df, fout, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
}
